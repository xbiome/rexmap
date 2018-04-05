# For this pipeline to work we will need to provide filenames for forward and
# reverse reads for each sample.

# Load required packages
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(limSolve))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(pso))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(stringr))

# Load HiMAP functions
himap_path = '/Users/igor/cloud/research/microbiome/himap'
source(file.path(himap_path, 'r/mergepairs.R'))
sourceCpp(file.path(himap_path, 'src/hamming.cpp'))
sourceCpp(file.path(himap_path, 'src/fastq_retrieve.cpp'))
# sourceCpp(file.path(himap_path, 'src/fitting_alignment_v2.cpp'))
source(file.path(himap_path, 'r-lib/blast_vs_fasta.R'))
source(file.path(himap_path, 'r-lib/spp_shorten.R'))
source(file.path(himap_path, 'r/read_fastx.R'))
# source(file.path(himap_path, 'r/taxonomy.R'))

#
# s1 r11 ++ r12 -
# s2 r21 =
# s3 r31 x
#
# ++-=x   5 total, s1->3, s2->1, s3->1
#


# Processing parameters

# Minimum BLAST alignment length
blast_aln_len = 200

out_from_input_file = function (new_input, suffix, input) {
  if (new_input == '<<default>>') {
    file.path(ddirname(dirname(input)), suffix)
  } else {
    new_input
  }
}

ie = function(test, yes, no) {
  if (test) {
    yes
  } else {
    no
  }
}


# Check one sequence vs a list of other sequences and return
# a subset of other sequences that are <= N hamming distances away
q = 'ACGTACGTACGT'
ts = c('ACGTACGTACGA', 'ACGTACAAACGT', 'GGTTACGTACGTN', 'ACGTACGGGTGTNN', 'ACGGGTGTACGT', 'TAGGGTGTACGT',
          'TGATGTAAACCC', 'ACGTATGTACGA', 'ACGTACGTGCGT')

find_possible_pcr_errors_query = function (query, targets, hamming_distance=1L) {
  dists = unname(sapply(targets, function (t) hamming(query, t)))
  return(which(dists==hamming_distance))
}

find_possible_pcr_errors = function (sequences, hamming_distance=1L) {
  # Sequences need to be already sorted by abundance, descending.

  for (i in 2:(length(sequences))) {
    targets = sequences[i:length(sequences)][nchar(sequences[i:length(sequences)]) == nchar(sequences[i-1])]
    if (length(targets) == 0) next
    hams = find_possible_pcr_errors_query(sequences[i-1], targets)
    if (length(hams) > 0) {
      cat('Sequence ', i-1, ' has 1-off variants at lower abundance: ')
      # Iterate over hits and return index from the original vector
      for (s in targets[hams]) {
        cat(which(sequences==s), ', ')
      }
      cat('\n')
    }
  }
}

load_copy_number_table = function (ref_cp) {
  # Load copy number information
  cp.dt = fread(ref_cp, select=c(1:3), colClasses=c('character', 'character', 'integer', 'character'))

  # Prepare copy number table columns
  cp.dt[, strain := gsub('_@rrn[0-9]+', '', strain_name)]

  # Generate strain index
  cp.dt[, strain_id := .GRP, by=strain]
  cp.dt[, strain_name := NULL]

  # How many unique 16S rrna variants of each strain?
  cp.dt[, rrn_uniq := length(unique(variant_id)), by=strain]

  # Sort by strain_id, then by
  cp.dt = cp.dt[order(strain_id, variant_id)]

  # Group sequences into OSU
  cp.dt[, spectrum := paste(copy_number, variant_id, sep=':', collapse=','), by=strain]

  # For strains which have only one unique 16S sequence variant, group them into the same
  # same OSU, regardless of the copy number.
  cp.dt[rrn_uniq==1, spectrum := paste0('X:', strsplit(spectrum[1], ':')[[1]][2]), by=strain]

  # OSU is a group of strains all with the same strain spectrum
  cp.dt[, osu_id := .GRP, by=spectrum]

  # Instead of X put mean copy number or 1
  cp.dt[rrn_uniq==1, spectrum := gsub('X', round(mean(copy_number)), spectrum), by=osu_id]
  # ref_cp.dt[rrn_uniq==1, spectrum := gsub('X', '1', spectrum), by=osu_id]

  # Count number of strains binned together
  cp.dt[, no_strains_in_osu := length(unique(strain)), by=osu_id]

  # Generate species names
  cp.dt[, species := gsub('^([^_]+)_([^_]+)_.*', '\\1_\\2', strain)]
  # Count unique number of species per OSU
  cp.dt[, no_species_in_osu := length(unique(species)), by=osu_id]

  cp.dt
}

add_consensus = function (dada_res, derep, fq_tri, fq_fil, truncLen,
                          verbose=T, ncpu=detectCores()-1) {
  if (verbose) cat('Retrieving full-length sequences...\n')
  for (s_id in 1:length(dada_res)) { # For each sample s_id
    if (verbose) cat('Sample ', s_id, '. Load...')
    x = partid_to_fastqid(dada_res[[s_id]]$map-1, derep[[s_id]]$map-1)
    # Load both merged (for retrieval) and filtered sequences for referencing
    # Filtered sequences are enumerated in dada2, so this file is used to
    # extract meta-data information for each read that we need to look up
    # in the merged file.
    fq_fil_data = sfastq_reader(fq_fil[s_id])
    fq_mer_data = sfastq_reader(fq_tri[s_id])
    if (verbose) cat('OK. Consensus...')
    # pid_to_seq = vector("list", length(x))
    pid_to_seq = unlist(mclapply(1:length(x), function (i) {
      metas = fq_fil_data[['meta']][x[i][[1]]+1]
      mer_seqs = fq_mer_data[['seqs']][which(fq_mer_data[['meta']] %in% metas)]
      # Select only the right hanging part
      mer_seqs_right = str_sub(mer_seqs, start=truncLen+1)
      # pid_to_seq[i] = gsub('[-]{1,}.*', '', consensus_sequence(mer_seqs_right))
      mer_seqs_cons = gsub('[-|N]{1,}.*', '', consensus_sequence(mer_seqs_right))
      names(mer_seqs_cons) = as.character(i)
      return(mer_seqs_cons)
    }, mc.cores=ncpu))
    pid_to_seq = pid_to_seq[order(as.integer(names(pid_to_seq)))]
    if (verbose) cat('OK. Update...')
    # for (i in 1:length(x)) { # for each partition i in sample s_id
    #   metas = fq_fil_data[['meta']][x[i][[1]]+1]
    #   mer_seqs = fq_mer_data[['seqs']][which(fq_mer_data[['meta']] %in% metas)]
    #   # Select only the right hanging part
    #   mer_seqs_right = str_sub(mer_seqs, start=truncLen+1)
    #   pid_to_seq[i] = gsub('[-]{1,}.*', '', consensus_sequence(mer_seqs_right))
    #   # if (grepl('-', pid_to_seq[i])) print(i)
    # }
    # sid_to_pid_to_seq[s_id] = pid_to_seq
    # Concatenate left consensus, dada2 middle partition sequence and right consensus
    dada_res[[s_id]]$sequence = paste(dada_res[[s_id]]$sequence,
                                     pid_to_seq,
                                     sep='')
    names(dada_res[[s_id]]$denoised) = dada_res[[s_id]]$sequence
    # cat('.')
    if (verbose) cat('OK.\n')
  }
  return(dada_res)
}



blast_out_to_best_cp = function (
   blast_output, ref_cp, aln_params=c(5L, -4L, -8L, -6L),
   blast_out_fmt='qseqid sseqid qlen length qstart qend sstart send slen qseq sseq',
   verbose=T, ncpu=detectCores()-1, variant_sep='-'
  )
{

  if (verbose) cat('blast out: ')
  blast_out.dt = fread(blast_output)
  names(blast_out.dt) = strsplit(blast_out_fmt, ' ')[[1]]
  # Calculate BLAST scores.
  blast_out.dt[, c('score', 'match', 'mismatch', 'gapopen', 'gapextend') := transpose(unname(mcmapply(function (q,s) {
    aln = compare_alignment(q,s)
    as.integer(c(sum(aln_params*aln), aln[1], aln[2], aln[3], aln[4]))
  }, qseq, sseq, SIMPLIFY=F, mc.cores=ncpu)))]

  # For local alignments that are one-off at each end, reduce the score using
  # a mismatch penalty. (We deliberately use indel score that is very similar.)
  blast_out.dt[qstart >= 2 & sstart >= 2, score := score-aln_params[2]]
  blast_out.dt[qstart >= 2 & sstart >= 2, mismatch := mismatch+1L]
  blast_out.dt[qend <= qlen-1 & send < slen, score := score-aln_params[2]]
  blast_out.dt[qend <= qlen-1 & send < slen, mismatch := mismatch+1L]

  # Filter out alignments with large gaps at beginning or the end
  # blast_out.dt = blast_out.dt[qstart <= 2 & qend >= qlen-1]
  if (verbose) cat('OK. blast best: ')

  # For each dada2 sequence, keep only the best matches by bitscore
  blast_best.dt = blast_best_seq_matches(blast_out.dt, strsplit(blast_out_fmt, ' ')[[1]][1],
                                         strsplit(blast_out_fmt, ' ')[[1]][2])
  blast_best.dt[, variant_id := gsub('^([0-9]+)-.*', '\\1', sseqid)]
  # Remove sequence alignment columns since we already calculated alignment score
  blast_best.dt[, c('qseq', 'sseq') := NULL]
  blast_best.dt[, c('sseqid') := NULL]
  # Calculate percentage similarity
  blast_best.dt[, pctsim := round(100*match/(match+mismatch+gapopen+gapextend), 2)]

  if (verbose) cat('OK. copy number table: ')
  cp.dt = fread(ref_cp, select=c(1:3), colClasses=c('character', 'character', 'integer', 'character'))

  # Prepare copy number table columns
  cp.dt[, strain := gsub('_@rrn[0-9]+', '', strain_name)]

  # Generate strain index
  cp.dt[, strain_id := .GRP, by=strain]
  cp.dt[, strain_name := NULL]

  # How many unique 16S rrna variants of each strain?
  cp.dt[, rrn_uniq := length(unique(variant_id)), by=strain]

  # Sort by strain_id, then by
  cp.dt = cp.dt[order(strain_id, variant_id)]
  if (verbose) cat('OK.\n')
  # Merge with blast hits
  # First do perfect matches

  # Original variant_ids depend how the reduced database is made. That is, they represent unique 16S regions
  # including the overhangs, which do not exactly match the lengths of sequenced regions. Here we just generate
  # new variant_id that pools together all sequences that are different in the extended region but the same
  # in the reduced region that was sequenced (and new variant_id is taken as just the first in the list).
  # blast_best.dt = blast_best.dt[, {
  #   vids = unique(variant_id)
  #   merge(.SD, cp.dt[variant_id %in% vids, .(variant_id, strain, copy_number)], all.x=T, by='variant_id')
  # }, by=qseqid]

  # blast_best.dt[, variant_id := variant_id[1], by=qseqid]

  # Find variant ids from the database that all map to the same shorter sequence in the data and pool together
  # those variant ids into a new one.
  # Find DIFFERENT qseqids with exact same qstart qend
  # blast_best2.dt = blast_best.dt[, , by=.(qseqid, strain)]
  if (verbose) cat('merge: ')
  blast_best.dt = merge(blast_best.dt, cp.dt[, .(variant_id, strain, copy_number)],
                        by='variant_id', allow.cartesian=T)
  if (verbose) cat('OK. Fix overhang differences:')
  hang_diff_qseqids = blast_best.dt[pctsim==100, 1, by=qseqid][, unique(qseqid)]
  for (q in hang_diff_qseqids) {
    vids = blast_best.dt[pctsim==100 & qseqid==q, sort(unique(variant_id))]
    blast_best.dt[pctsim==100 & qseqid==q, variant_id := paste(vids, collapse=variant_sep)]
    cp.dt[variant_id %in% vids, variant_id := paste(vids, collapse=variant_sep)]
  }
  if (verbose) cat('.')
  blast_best.dt[pctsim==100, copy_number := sum(copy_number), by=.(variant_id, strain)]
  blast_best.dt = unique(blast_best.dt)

  if (verbose) cat('.')
  # Add up copy numbers for strain that differ in the undetected overhang
  cp.dt[, copy_number := sum(copy_number), by=.(variant_id, strain_id)]
  cp.dt = unique(cp.dt)
  if (verbose) cat('OK. ')


  # Recalculate rrn_uniq after re-numerating variant_ids (b/c the slightly shorter sequence now has the same
  # or smaller number of variants)
  cp.dt[, rrn_uniq := length(unique(variant_id)), by=strain_id]

  # Group sequences into OSU
  cp.dt[, spectrum := paste(copy_number, variant_id, sep=':', collapse=','), by=strain]

  # For strains which have only one unique 16S sequence variant, group them into the same
  # same OSU, regardless of the copy number.
  cp.dt[rrn_uniq==1, spectrum := paste0('X:', strsplit(spectrum[1], ':')[[1]][2]), by=strain]

  # OSU is a group of strains all with the same strain spectrum
  cp.dt[, osu_id := .GRP, by=spectrum]

  # Instead of X put mean copy number or 1
  cp.dt[rrn_uniq==1, spectrum := gsub('X', round(mean(copy_number)), spectrum), by=osu_id]
  # ref_cp.dt[rrn_uniq==1, spectrum := gsub('X', '1', spectrum), by=osu_id]

  # Count number of strains binned together
  cp.dt[, no_strains_in_osu := length(unique(strain)), by=osu_id]

  # Generate species names
  cp.dt[, species := gsub('^([^_]+)_([^_]+)_.*', '\\1_\\2', strain)]
  # Count unique number of species per OSU
  cp.dt[, no_species_in_osu := length(unique(species)), by=osu_id]
  if (verbose) cat('OK.\n')
  return(list(blast_best.dt, cp.dt))
}



ab_mat_to_dt = function (ab_tab_nochim, fq_prefix_split='_') {
   ab_tab_nochim.dt = as.data.table(unname(ab_tab_nochim))
   ab_tab_nochim.dt[, sample_id := sapply(strsplit(dimnames(ab_tab_nochim)[[1]], fq_prefix_split, fixed=T), `[`, 1)]
   ab_tab_nochim_m.dt = melt(ab_tab_nochim.dt, id.vars = 'sample_id',
                               variable.name = 'dada2_seqid', value.name = 'raw_count')
   ab_tab_nochim_m.dt[, qseqid := as.numeric(gsub('V', '', dada2_seqid))]
   ab_tab_nochim_m.dt[, dada2_seqid := NULL]
   return(ab_tab_nochim_m.dt)
}


pcr_primer_filter = function (fq_in, fq_out, pr_fwd='CCTACGGGNGGCWGCAG',
                              pr_rev='GGATTAGATACCCBDGTAGTCC',
                              pr_fwd_maxoff=10, pr_rev_maxoff=10, return_noprimer=T,
                              multithread=T, ncpu=detectCores()-1, max_mismatch=2) {
  # fq_in = input fastq file
  # fq_out = output fastq file (without primers)
  # pr_fwd = forward primer
  # pr_rev = reverse primer
  #
  # Load Input FASTQ file
  #
  # Count extended DNA symbols
  # Ignore any extended DNA symbols. Shouldn't have too many of them anyway.
  pr_fwd = gsub('[^ACGT-]', 'N', pr_fwd)
  pr_rev = gsub('[^ACGT-]', 'N', pr_rev)
  in_fq = sfastq_reader(fq_in)
  # max_mismatch = 2

  seq = in_fq[['seqs']][1]
  qual = in_fq[['qual']][1]

  fastq_trimmer = function (meta, seq, qual, return_noprimer=return_noprimer) {

    # Search for forward primer
    aln = C_nwalign(pr_fwd, seq, match=1, mismatch=-1, indel=-1)
    pr_fwd_left = max(regexpr('[^-]', aln[1])[1],
                      regexpr('[^-]', aln[2])[1])
    pr_fwd_right = regexpr('[-]{1,}$', aln[1])[1]-1
    pr_fwd_n = lengths(regmatches(pr_fwd, gregexpr('N', pr_fwd)))
    pr_rev_n = lengths(regmatches(pr_rev, gregexpr('N', pr_rev)))
    # Alignment statistics: match, mismatch, gapopen, gapextend
    aln_stat = compare_alignment(str_sub(aln[1], start=pr_fwd_left, end=pr_fwd_right),
                                 str_sub(aln[2], start=pr_fwd_left, end=pr_fwd_right))
    # Forward primer alignment is acceptable if it's near beginning (within first 5 nts)
    # and if it doesn't have more than 2 mismatches for non-N symbols, which includes indels.
    pr_fwd_found = FALSE
    pr_rev_found = FALSE
    if (pr_fwd_left < pr_fwd_maxoff & aln_stat[2]+aln_stat[3]+aln_stat[4]-pr_fwd_n <= max_mismatch) {
      pr_fwd_found = TRUE
    }

    # Search for reverse primer
    aln2 = C_nwalign(pr_rev, seq, match=1, mismatch=-1, indel=-1)
    pr_rev_left = regexpr('[^-]', aln2[1])[1]
    pr_rev_right = min(regexpr('[ACGTN][^ACGTN]*$', aln2[1])[1],
                       regexpr('[ACGTN][^ACGTN]*$', aln2[2])[1])
    aln2_stat = compare_alignment(str_sub(aln2[1], start=pr_rev_left, end=pr_rev_right),
                                  str_sub(aln2[2], start=pr_rev_left, end=pr_rev_right))
    if (pr_rev_left > nchar(seq)-nchar(pr_rev)-pr_rev_maxoff & aln2_stat[2]+aln2_stat[3]+aln2_stat[4]-pr_rev_n <= max_mismatch) {
      # Reverse primer found
      pr_rev_found = TRUE
    }

    # Trim if needed
    start = 1L
    end = -1L
    if (pr_fwd_found) start = pr_fwd_right
    if (pr_rev_found) end = pr_rev_left

    #if ((pr_fwd_found & pr_rev_found) | return_noprimer) {
      return(list('meta' = meta,
                  'seqs' = str_sub(seq, start=start, end=end),
                  'qual' = str_sub(qual, start=start, end=end)))
    #}
  }

  # Apply fastq_trimmer() to each sequence in this file
  out_trimmed = unname(mcmapply(fastq_trimmer, in_fq[['meta']], in_fq[['seqs']], in_fq[['qual']],
                         mc.cores=ncpu, SIMPLIFY=FALSE))

  # Save results in a new file
  fastq_list_writer(out_trimmed, fq_out, ncpu=ncpu)
}

pcr_primer_trimmer = Vectorize(pcr_primer_filter, vectorize.args=c('fq_in', 'fq_out'))





# OSU binning and abundance estimation
blast_cp_to_osu_dt = function (
  blast_best.dt, cp.dt, ab_tab_nochim_m.dt,
  pctsim_min = 100.00, # Min. pct sim for OSU binning
  ncpu = detectCores()-1,
  verbose=T
) {
  if (verbose) cat('OSU table: ')
  # First generate a BLAST best results table with just alignment pctsim
  blast_reduced.dt = unique(blast_best.dt[, .(variant_id, qseqid, pctsim)])[order(qseqid)]
  blast_reduced.dt[, variant_id_new := variant_id[1], by=qseqid]
  # blast_reduced.dt = merge(blast_reduced.dt, ab_tab_nochim_m.dt[, .(qseqid, raw_count, sample_id)],
  #                          by='qseqid', all.x=T)
  blast_reduced.dt = blast_reduced.dt[, .SD[pctsim==max(pctsim)][1], by=variant_id]

  # Generate a separate reference table
  osu.dt = unique(cp.dt[, .(osu_id, spectrum, strain, no_strains_in_osu)])
  osu.dt[, no_variants := sapply(spectrum, function (s) length(strsplit(s, ',', fixed=T)[[1]]))]
  # osu_m.dt = osu.dt[, .(variant_id = sub('.*:([0-9]{5})$', '\\1', strsplit(spectrum, ',')[[1]]),
  #                       spectrum=spectrum, strain, no_variants), by=osu_id]

  # Now generate a smaller OSU table containing only osu id that contain
  # at least one variant ID from the data with > 0 raw count.
  varids_data = blast_reduced.dt[pctsim >= pctsim_min, unique(variant_id)]
  # varids_data = [pctsim > pctsim_min, unique(variant_id)]
  osu_data.dt = osu.dt[, {
    varids = strsplit(gsub('[0-9]+:', '', spectrum), ',')[[1]]
    overlap = intersect(varids_data, varids)
    if (length(overlap) > 0) .SD
  }, by=osu_id]

  if (verbose) cat('OK. Melt table: ')
  sample_ids = ab_tab_nochim_m.dt[, unique(sample_id)]
  osu_data_m.dt = rbindlist(mclapply(sample_ids, function (s) {
    # Melt OSU data spectrum
    osu_data.dt[, { # For each OSU id
      # cat('osu_id: ', osu_id, ', group: ', .GRP, '\n')
      # Extract variant IDs and respective copy numbers from spectrum string
      cp_vid = strsplit(strsplit(spectrum, ',')[[1]], ':')
      cp = as.integer(sapply(cp_vid, `[`, 1))
      vid = sapply(cp_vid, `[`, 2)
      # For each variant id, extract dada2_seqid
      did = unname(sapply(vid, function (v) {
        blast_reduced.dt[variant_id==v & pctsim >= pctsim_min, if (.N==0) NA else qseqid]
      }))
      pctsims = unname(sapply(vid, function (v) {
        blast_reduced.dt[variant_id==v & pctsim >= pctsim_min , if (.N==0) NA else pctsim]
      }))
      raw_c = unname(sapply(did, function (i) {
        if (is.na(i)) {
          0L
        }
        else {
          ab_tab_nochim_m.dt[qseqid==i & sample_id==s, as.integer(raw_count)]
        }
      }))
      .(copy_number=cp, variant_id=vid, raw_count=raw_c, sample_id=s,
        pctsim=as.double(pctsims), no_strains_in_osu=unique(no_strains_in_osu))
    }, by=osu_id]
  }, mc.cores=ncpu))
  if (verbose) cat(' OK.')
  return(osu_data_m.dt)
}

pctsim_range_old = function (p) {
  up = unique(p)
  if (length(up) == 1) return(as.character(up))
  else {
    return(paste0(min(p), ' - ', max(p)))
  }
}

pctsim_range = function (p) {
  return(max(p, na.rm=T))
}


print_strains = function (strains, raw=T) {
  # Input: vector of strains
  # Prints a single string with reduced list of strains
  # such that any non-_bacterium or _sp. strain is shown
  # on a species level.
  if (raw) return(paste(strains, collapse=','))

  # Also add strains for species that occur only once!!
  if (length(strains)==1) return(strains)
  else {
    genuses = gsub('^[ ]*([^_]+)_.*', '\\1', strains)
    species = gsub('^[^_]+_([^_]+)[_]?.*', '\\1', strains)
    na_filter = grepl('^(sp\\.|bacterium)', species)
    sp1 = paste(genuses[!na_filter], species[!na_filter], sep='_')
    sp2 = strains[na_filter]
    ft = as.table(sort(table(c(sp1, sp2)), decreasing=T))
    ft2 = ft[ft>1]
    # Get full strain names for species that occur only once
    sp_once = names(ft[ft==1])
    st_once = c()
    for (sp in sp_once) {
      st_once = c(st_once, grep(sp, strains, fixed=T, value=T))
    }
    out = c()
    if (length(ft2) > 0) {
      out = paste(paste0(paste(names(ft2), ft2, sep='_['), ']'), collapse=',')
    }
    if (length(out) > 0 & length(st_once) > 0) {
      out = paste(out, paste(st_once, collapse=','), sep=',')
    } else if (length(out) == 0 & length(st_once) > 0) {
      out = paste(st_once, collapse=',')
    }
    return(out)
    # return(gsub('_[1]', '', out, fixed=T))
  }
}


osu_cp_to_all_abs = function (osu_data_m.dt, cp.dt, blast_best.dt, ab_tab_nochim_m.dt,
                              pctsim_min=100, ncpu=detectCores()-1, verbose=T, seed=42,
                              osu_offset=1000000L, raw=T) {

  if (verbose) cat('Preparing blast tables...')
  var_count.dt = unique(osu_data_m.dt[, .(variant_id, raw_count, sample_id)])
  common_variant_ids = var_count.dt[, if (all(raw_count > 0)) .SD, by=variant_id][, unique(variant_id)]
  blast_best2.dt = blast_best.dt[
    pctsim < pctsim_min,
    .(pctsim=pctsim_range(pctsim), species=print_strains(strain, raw=raw)), by=qseqid]
  if (verbose) cat('OK.\n')

  # Optimization function
  H = function(x) as.numeric(x>0)
  f = function (x, A, B, x0) {
    # Return an optimization function
    eps = A %*% x - B
    (t(eps) %*% eps) + ((x0^2) %*% H(x))
  }

  osu_data_m_single.dt = unique(osu_data_m.dt[, if (length(unique(variant_id)) == 1) .SD, by=osu_id][raw_count > 0, .(osu_id, variant_id, copy_number)])

  osu_sp.dt = cp.dt[, .(species = print_strains(strain, raw=raw)), by=osu_id]

  sample_ids = ab_tab_nochim_m.dt[, unique(sample_id)]
  all_abs.dt = rbindlist(mclapply(sample_ids, function (s) {

      Ab.dt = dcast(
         osu_data_m.dt[sample_id==s], variant_id ~ osu_id,
         value.var='copy_number', fill=0
      )
      Ab.dt = merge(
         Ab.dt,
         unique(osu_data_m.dt[sample_id==s, .(variant_id, raw_count)]),
         by='variant_id'
      )

      A = as.matrix(Ab.dt[, 2:(ncol(Ab.dt)-1)])
      dimnames(A)[[1]] = Ab.dt[, variant_id]

      # Generate a matrix with B coefficients
      B = as.matrix(Ab.dt[, raw_count])
      dimnames(B)[[1]] = Ab.dt[, variant_id]

      # Solve
      sol = lsei(A, B, fulloutput=T, G=diag(ncol(A)), H=matrix(c(0), nrow=ncol(A), ncol=1), type=2)
      osu_th = 1e-1
      osu_ab = sol$X
      osu_count = as.integer(osu_ab[which(osu_ab > osu_th)])
      osu_ids   = as.integer(names(osu_ab[which(osu_ab > osu_th)]))
      osu_ab.dt = data.table(osu_id=osu_ids, osu_count=osu_count)
      setorder(osu_ab.dt, osu_id)

      # Try optimizing the full table first?
      x0 = osu_ab.dt[, osu_count]
      A_columns = as.character(osu_ab.dt[, osu_id])
      Ar = A[, A_columns, drop=F]
      Ar = Ar[setdiff(names(rowSums(Ar) > 0), names(B[B==0,])), , drop=F]
      Br = B[rownames(Ar), , drop=F]

      # Generate a graph from Ar matrix, then identify all
      # connected clusters.
      g = make_empty_graph()
      A_rows = rownames(Ar)
      # Add vertices like this
      #    1  2  3
      # 4  .  .  .
      # 5  .  .  .
      # 6  .  .  .
      #  nrow=3, ncol=4
      g = add_vertices(g, ncol(Ar))
      g = add_vertices(g, nrow(Ar))
      vxs = c(colnames(Ar), rownames(Ar))
      for (j in 1:ncol(Ar)) {
         g = add_vertices(g, 1)
      }
      for (i in seq(ncol(Ar)+1, ncol(Ar)+1+nrow(Ar))) {
         g = add_vertices(g, 1)
      }
      for (i in 1:nrow(Ar)) {
         for (j in 1:ncol(Ar)) {
            if (Ar[i,j] > 0) {
               # cat('add: ', i+ncol(Ar), '-', j, '\n')
               g = add_edges(g, c(i+ncol(Ar), j))
               g = add_edges(g, c(j, i+ncol(Ar)))
            }
         }
      }
      g = as.undirected(g)
      # Find all connected clusters
      cls = lapply(groups(clusters(g)), function (x) if (length(x[x<=ncol(Ar)]) > 1) x else NA)
      cls = cls[!is.na(cls)]
      osu_ab2.dt = copy(osu_ab.dt)
      # For each cluster i
      if (length(cls) > 0) {
        for (i in 1:length(cls)) {
          # Add back any osu with a single variant_id
          # Sometimes, optimization will omit osu_ids with mapping to single
          # variant_ids so we bring those back here manually.
          cl = cls[i][[1]]
          osu_ids = vxs[cl[cl<=ncol(Ar)]]
          variant_ids = vxs[cl[cl>ncol(Ar)]]
          # osu_ids = c(osu_ids, osu_data_m_single.dt[variant_id %in% variant_ids][, osu_id])
          # variant_ids = c(variant_ids, osu_data_m_single.dt[variant_id %in% variant_ids][, osu_id])
          Ar2 = Ar[variant_ids, osu_ids]
          Br2 = Br[variant_ids]

          x = as.matrix(dcast(osu_data_m_single.dt[variant_id %in% variant_ids], variant_id ~ as.character(osu_id), value.var='copy_number')[, -1])
          rownames(x) = osu_data_m_single.dt[variant_id %in% variant_ids, variant_id]
          x[is.na(x)] = 0
          Ar3.dt = merge(as.data.table(Ar2, keep.rownames=T),
                         as.data.table(x, keep.rownames=T), all.x=T)
          Ar3.dt[is.na(Ar3.dt)] = 0
          Ar3 = as.matrix(Ar3.dt[, -1])
          Ar3 = Ar3[, order(as.integer(colnames(Ar3)))]
          rownames(Ar3) = Ar3.dt[, rn]

          Br3 = as.matrix(B[rownames(Ar3),])
          osu_ab3.dt = merge(osu_ab.dt, data.table(osu_id=as.integer(colnames(Ar3))), by='osu_id', all.y=T)
          osu_ab3.dt[is.na(osu_count), osu_count := 0L]
          # setorder(osu_ab3.dt, osu_id)

          x0 = osu_ab3.dt[, osu_count]

          x0_w = rep(0.3, length(x0))
          ar3 = copy(Ar3)
          for (r in 1:nrow(Ar3)) ar3[r,] = ar3[r,] / Br3[r]
          tmp = lapply(1:1000, function (it) {
            pso = psoptim(
              rep(0, length(x0)),
              f, lower=rep(0, length(x0)), upper=rep(max(Br3), length(x0)),
              A=ar3, B=rep(1, nrow(ar3)),
              x0 = as.numeric(x0_w),
              control=list('maxit'=10,
                           'vectorize'=T
                           #'hybrid'=F, 'rand.order'=F
              ))
            list(round(pso$par, 2),
                 pso$value)
          })
          res = tmp[which.min(sapply(tmp, function (x) x[[2]]))]

          osu_ab3.dt[, osu_count := as.integer(res[[1]][[1]])]
          osu_ab2.dt = merge(osu_ab2.dt[!(osu_id %in% colnames(Ar3))], osu_ab3.dt, all=T, by=names(osu_ab3.dt))
        }
      }

      osu_ab2.dt = osu_ab2.dt[osu_count > 0]

      # Join table to OSU abundances
      osu_ab4.dt = merge(osu_ab2.dt, osu_sp.dt, by='osu_id',
                        all.x=T)
      setorder(osu_ab4.dt, -osu_count)


      ab_tab2.dt = merge(
        ab_tab_nochim_m.dt[sample_id==s],
        blast_best2.dt,
        by='qseqid'
      )
      ab_tab2.dt[, osu_id := osu_offset+qseqid]
      ab_tab2.dt[, osu_count := raw_count]
      ab_tab2.dt[, c('raw_count', 'qseqid') := NULL]
      osu_ab5.dt = merge(osu_ab4.dt,
                        unique(osu_data_m.dt[, .(osu_id, pctsim=pctsim_range(pctsim))]),
                        by='osu_id')

      # Merge OSU analysis with low sim sequences
      all_ab.dt = merge(osu_ab5.dt, ab_tab2.dt, by=intersect(names(osu_ab5.dt), names(ab_tab2.dt)),
                        all=T)
      # Recalculate abundances
      all_ab.dt[, sample_id := s]
      setcolorder(all_ab.dt, c('sample_id', 'osu_id', 'osu_count', 'species',
                               'pctsim'))
      # all_ab.dt[order(-pctsim)]
      all_ab.dt[osu_count>0]
  }, mc.cores=ncpu))
  return(unique(all_abs.dt))
}




collapseNoMismatch2 <- function(seqtab, minOverlap=20, verbose=FALSE) {
  # CollapseNoMismatch but with progress bar
  if (verbose) cat('Start collapseNoMismatchSet\n')
  unqs.srt <- sort(getUniques(seqtab), decreasing=TRUE)
  seqs <- names(unqs.srt) # The input sequences in order of decreasing total abundance
  seqs.out <- character(0) # The output sequences (after collapsing)
  # collapsed will be the output sequence table
  collapsed <- matrix(0, nrow=nrow(seqtab), ncol=ncol(seqtab))
  colnames(collapsed) <- colnames(seqtab) # Keep input ordering for output table
  rownames(collapsed) <- rownames(seqtab)
  i = 0
  for(query in seqs) {
    added=FALSE
    prefix <- substr(query, 1, minOverlap)
    suffix <- substr(query, nchar(query)-minOverlap+1,nchar(query))
    for(ref in seqs.out) { # Loop over the reference sequences already added to output
      # Prescreen to see if costly alignment worthwhile, this all should possibly be C-side
      if(grepl(prefix, ref) || grepl(suffix, ref)) {
        if(nwhamming(query,ref,band=-1) == 0) { # No mismatches/indels, join more abundant sequence
          collapsed[,ref] <- collapsed[,ref] + seqtab[,query]
          added=TRUE
          break
        }
      }
    } # for(ref in seqs.out)
    if(!added) {
      collapsed[,query] <- seqtab[,query]
      seqs.out <- c(seqs.out, query)
    }
    i = i + 1
    if (verbose) cat('\r* processed ', i, ' out of ', length(seqs), ' sequences.')
  } # for(query in seqs)
  cat('\n')
  if(!identical(unname(colSums(collapsed)>0), colnames(collapsed) %in% seqs.out)) {
    stop("Mismatch between output sequences and the collapsed sequence table.")
  }
  collapsed <- collapsed[,colnames(collapsed) %in% seqs.out,drop=FALSE]

  if(verbose) message("Output ", ncol(collapsed), " collapsed sequences out of ", ncol(seqtab), " input sequences.")
  collapsed
}

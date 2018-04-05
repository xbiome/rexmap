# For this pipeline to work we will need to provide filenames for forward and
# reverse reads for each sample.

# Load required packages
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(limSolve))
suppressPackageStartupMessages(library(pso))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(stringr))

# Load HiMAP functions
source('~/cloud/research/microbiome/himap/r/mergepairs.R')
sourceCpp('~/cloud/research/microbiome/himap/src/hamming.cpp')
sourceCpp('~/cloud/research/microbiome/himap/src/fastq_retrieve.cpp')
sourceCpp('~/cloud/research/microbiome/himap/src/fitting_alignment_v2.cpp')
source('~/cloud/research/microbiome/himap/r-lib/blast_vs_fasta.R')
source('~/cloud/research/microbiome/himap/r-lib/spp_shorten.R')
source('~/cloud/research/microbiome/himap/r/read_fastx.R')

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
      mer_seqs_cons = gsub('[-]{1,}.*', '', consensus_sequence(mer_seqs_right))
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
   blast_out_fmt='qseqid sseqid qlen length qstart qend sstart send qseq sseq',
   verbose=T, ncpu=detectCores()-1
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
  blast_out.dt[qstart == 2, score := score-aln_params[2]]
  blast_out.dt[qstart == 2, mismatch := mismatch+1L]
  blast_out.dt[qend == qlen-1, score := score-aln_params[2]]
  blast_out.dt[qend == qlen-1, mismatch := mismatch+1L]

  # Filter out alignments with large gaps at beginning or the end
  blast_out.dt = blast_out.dt[qstart <= 2 & qend >= qlen-1]
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

  # Merge with blast hits
  # First do perfect matches

  # Original variant_ids depend how the reduced database is made. That is, they represent unique 16S regions
  # including the overhangs, which do not exactly match the lengths of sequenced regions. Here we just generate
  # new variant_id that pools together all sequences that are different in the extended region but the same
  # in the reduced region that was sequenced (and new variant_id is taken as just the first in the list).
  blast_best.dt = blast_best.dt[, {
    vids = unique(variant_id)
    merge(.SD, cp.dt[variant_id %in% vids, .(variant_id, strain, copy_number)], all.x=T, by='variant_id')
  }, by=qseqid]
  # blast_best.dt[, variant_id := variant_id[1], by=qseqid]

  # Find variant ids from the database that all map to the same shorter sequence in the data and pool together
  # those variant ids into a new one.
  variant_id_groups = list()
  for (q in blast_best.dt[, unique(qseqid)]) {
    variant_id_groups[[q]] = blast_best.dt[qseqid==q, unique(variant_id)]
    # First, for exact strain hits add the copy numbers
    cp.dt[variant_id %in% blast_best.dt[qseqid==q, unique(variant_id)], variant_id := variant_id[1]]
    blast_best.dt[qseqid==q, variant_id := variant_id[1]]
  }

  # Add up copy numbers for strain that differ in the undetected overhang
  cp.dt[, copy_number := sum(copy_number), by=.(variant_id, strain_id)]
  cp.dt = unique(cp.dt)
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
                              pr_fwd_maxoff=10,
                              pr_rev_maxoff=10,
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

  seq = in_fq[['seq']][1]
  qual = in_fq[['qual']][1]

  fastq_trimmer = function (meta, seq, qual) {

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

    return(list('meta' = meta,
                'seq' = str_sub(seq, start=start, end=end),
                'qual' = str_sub(qual, start=start, end=end)))
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
  blast_reduced.dt = merge(blast_reduced.dt, ab_tab_nochim_m.dt[, .(qseqid, raw_count)],
                           by='qseqid', all.x=T)
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
      cat('osu_id: ', osu_id, ', group: ', .GRP, '\n')
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




osu_cp_to_all_abs = function (osu_data_m.dt, cp.dt, blast_best.dt, ab_tab_nochim_m.dt,
                              pctsim_min=100,
                              ncpu=detectCores()-1, verbose=T, seed=42,
                              osu_offset=1000000L) {
  var_count.dt = unique(osu_data_m.dt[, .(variant_id, raw_count, sample_id)])
  common_variant_ids = var_count.dt[, if (all(raw_count > 0)) .SD, by=variant_id][, unique(variant_id)]
  blast_reduced.dt = unique(blast_best.dt[, .(variant_id, qseqid, pctsim)])[order(qseqid)]
  blast_reduced.dt[, variant_id_new := variant_id[1], by=qseqid]
  blast_nonosu.dt = copy(blast_reduced.dt)
  blast_reduced.dt = merge(blast_reduced.dt, ab_tab_nochim_m.dt[, .(qseqid, raw_count)],
                           by='qseqid', all.x=T)
  blast_reduced.dt = blast_reduced.dt[, .SD[pctsim==max(pctsim) & raw_count==max(raw_count)][1], by=variant_id]

  # blast_nonosu.dt = unique(blast_best.dt[!(qseqid %in% blast_reduced.dt[, qseqid]), .(variant_id, qseqid, pctsim)])
  # blast_nonosu.dt[, variant_id_new := variant_id[1], by=qseqid]

  # osu_data_m2.dt = osu_data_m.dt

  H = function(x) as.numeric(x>0)

  f = function (x, A, B, x0) {
    # Return an optimization function
    eps = A %*% x - B
    (t(eps) %*% eps) + ((x0^2) %*% H(x))
  }

  sample_ids = ab_tab_nochim_m.dt[, unique(sample_id)]
  all_abs.dt = rbindlist(lapply(sample_ids, function (s) {

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
      # osu_ab.dt[, osu_ab := osu_count/sum(osu_count)]
      osu_ab.dt = osu_ab.dt[order(-osu_count)]

      # Try optimizing the full table first?
      # x0 = osu_ab[osu_ab > 0]
      x0 = osu_ab.dt[, osu_count]
      A_columns = as.character(osu_ab.dt[, osu_id])
      Ar = A[, A_columns]
      Ar = Ar[rowSums(Ar) > 0, ]
      Br = B[rownames(Ar), ]


       pso = psoptim(par=x0, fn=f, lower=rep(0, length(x0)),
                     upper=rep(max(x0)*2, length(x0)),
                     A=Ar, B=Br, x0=rep(mean(x0), length(x0)),
                     control=list('maxit'=4000, 'hybrid'=F,
                                  'rand.order'=T))
        osu_ab2.dt = data.table(
           osu_id=as.integer(A_columns[pso$par > 0]),
           osu_count=as.integer(pso$par[pso$par > 0])
         )

      # Generate a list of unique species for each osu id
      # osu_sp.dt = cp.dt[, .(species = paste(unique(species), collapse=',')), by=osu_id]
      osu_sp.dt = cp.dt[, .(species = {
         genuses = gsub('^([^_]+)_.*', '\\1', strain)
         species = gsub('^[^_]+_([^_]+)[_]?.*', '\\1', strain)
         na_species_indexes = species %in% c('sp.', 'bacterium')
         gsp = paste(genuses, species, sep='_')
         gsp[na_species_indexes] = strain[na_species_indexes]
         gsp_t = as.table(sort(table(gsp), decreasing=T))
         gsub('_[1]', '',
              paste(dimnames(gsp_t)[[1]], '_[', as.character(gsp_t), ']', sep='', collapse=','),
              fixed=T)
      }), by=.(osu_id)]

      # Join table to OSU abundances
      osu_ab3.dt = merge(osu_ab2.dt, osu_sp.dt, by='osu_id',
                        all.x=T)
      osu_ab3.dt = osu_ab3.dt[order(-osu_count)]
      osu_ab3.dt[, species := str_sub(species, end=55)]
      osu_ab3.dt[, sample_id := s]

      # Start with the LSEI solution, then numerically minimize cost function f
      # to potentially reduce the number of OSU bins. For Zheng-2015 data this reduces 25
      # OSU to ...

      variants_for_opt.dt = osu_data_m.dt[sample_id==s & osu_id %in% osu_ab.dt[, osu_id]][, if (length(unique(osu_id))>1)
        .(osu_id=osu_id, group=.GRP), by=variant_id]

      if (variants_for_opt.dt[, .N] > 0) {
        var_groups = variants_for_opt.dt[, unique(group)]
        for (g in var_groups) {
          oids = variants_for_opt.dt[group==g, osu_id]
          # oids = osu_data_m.dt[sample_id==s & osu_id %in% oids, if (length(unique(variant_id))==1) osu_id, by=osu_id][, osu_id]

          Ab.dt = dcast(osu_data_m.dt[sample_id==s & osu_id %in% oids], variant_id ~ osu_id, value.var='copy_number', fill=0)
          Ab.dt = merge(Ab.dt, unique(osu_data_m.dt[sample_id==s & osu_id %in% oids, .(variant_id, raw_count)]), by='variant_id')
          A = as.matrix(Ab.dt[, 2:(ncol(Ab.dt)-1)])
          dimnames(A)[[1]] = as.list(Ab.dt[, variant_id])
          B = as.matrix(Ab.dt[, raw_count])
          dimnames(B)[[1]] = as.list(Ab.dt[, variant_id])

          # Now further minimize this with a new cost
          x0 = osu_ab[colnames(A)]
          set.seed(seed)
          pso = psoptim(par=x0, fn=f, lower=rep(0, length(x0)),
                        upper=rep(max(x0)*2, length(x0)),
                        A=A, B=B,
                        x0=x0,
                        control=list('maxit'=10000, 'hybrid'=F,
                                     'rand.order'=F))
          osu_ab[colnames(A)] = pso$par
          no_zeros = length(pso$par[pso$par == 0])
          no_total = length(pso$par)
          if (verbose) cat('Set ', no_zeros, ' zeros out of ', no_total, '\n')
        }

        osu_count = as.integer(osu_ab[which(osu_ab > osu_th)])
        osu_ids = as.integer(names(osu_ab[which(osu_ab > osu_th)]))
        osu_ab2.dt = data.table(osu_id=osu_ids, osu_count=osu_count)
        osu_ab2.dt = osu_ab.dt[osu_count > 0]
      }
     # Add sequences with low pct similarity to known sequences to this table
     # ab_tab2.dt = merge(ab_tab_nochim_m.dt,
     #                    blast_reduced.dt[pctsim < pctsim_min],
     #                    by='qseqid')[sample_id==s & raw_count>0]
     ab_tab2.dt = merge(ab_tab_nochim_m.dt,
                        blast_nonosu.dt[pctsim < pctsim_min],
                        # blast_reduced.dt,
                        by='qseqid')[sample_id==s & raw_count>0]

      # else {
      #   ab_tab2.dt = merge(ab_tab_nochim_m.dt,
      #                      blast_reduced.dt[pctsim < pctsim_min],
      #                      by='qseqid')[sample_id==s & raw_count>0]
      # }
      # ab_tab2.dt is only for sequences that have not gone through optimization
      # due to low pct sim.
      # ab_tab2.dt[, c('adj_count', 'adj_ab', 'variant_id') := NULL]
      cat('Checkpoint 1.')
      ab_tab2.dt[, osu_id := osu_offset+qseqid]
      ab_tab2.dt[, variant_id := variant_id_new]
      # ab_tab2.dt[, c('qseqid') := NULL]
      ab_tab2.dt[, c('qseqid', 'variant_id_new') := NULL]
      cat('Checkpoint 1.1')
      blast_best2.dt = blast_best.dt[variant_id %in% ab_tab2.dt[, variant_id]]
      # Generate species names
      blast_var_sp.dt = blast_best2.dt[, .(species = {
         genuses = gsub('^([^_]+)_.*', '\\1', strain)
         species = gsub('^[^_]+_([^_]+)[_]?.*', '\\1', strain)
         na_species_indexes = species %in% c('sp.', 'bacterium')
         gsp = paste(genuses, species, sep='_')
         gsp[na_species_indexes] = strain[na_species_indexes]
         gsp_t = as.table(sort(table(gsp), decreasing=T))
         gsub('_[1]', '',
              paste(dimnames(gsp_t)[[1]], '_[', as.character(gsp_t), ']', sep='', collapse=','),
              fixed=T)
      }), by=.(qseqid, variant_id)]
      blast_best2.dt = merge(blast_best2.dt, blast_var_sp.dt, by=c('variant_id', 'qseqid'),
                             all.x=T)

      # blast_best2.dt[, species := gsub('^([^_]+)_([^_]+)_.*', '\\1_\\2', strain)]
      ab_tab2.dt = merge(ab_tab2.dt, blast_best2.dt[, .(species = paste(unique(species), collapse=',')), by=variant_id],
                         by='variant_id')
      cat('Checkpoint 1.2')
      ab_tab2.dt[, osu_count := raw_count]
      ab_tab2.dt[, raw_count := NULL]
      ab_tab2.dt[, variant_id := NULL]
      cat('Checkpoint 1.3')
      # osu_ab.dt[, osu_ab := osu_count/sum(osu_count)]
      osu_ab.dt = osu_ab.dt[order(-osu_count)]
      cat('Checkpoint 1.4')
      osu_ab.dt = merge(osu_ab.dt, osu_sp.dt, by='osu_id', all.x=T)
      cat('Checkpoint 1.5')
      osu_ab.dt = merge(osu_ab.dt,
                        # unique(osu_data_m.dt[, .(osu_id, pctsim=unique(pctsim[!is.na(pctsim)]))]),
                        unique(osu_data_m.dt[, .(osu_id, pctsim=round(max(pctsim, na.rm=T), 2))]),
                        by='osu_id')

      # Check if sample_id column is in the osu table
      if (!('sample_id' %in% names(osu_ab.dt))) osu_ab.dt[, sample_id := s]
      cat('Checkpoint 2.\n')
      # Merge OSU analysis with low sim sequences
      all_ab.dt = merge(osu_ab.dt, ab_tab2.dt, by=intersect(names(osu_ab.dt), names(ab_tab2.dt)),
                        all=T)
      cat('Checkpoint 3.\n')
      # cat(str(all_ab.dt))
      # cat('osu_ab.dt, nrow:', nrow(osu_ab.dt), '\n')
      # cat('ab_tab2.dt, nrow:', nrow(ab_tab2.dt), '\n')
      # cat('all_ab.dt, nrow:', nrow(all_ab.dt), '\n')
      # Recalculate abundances
      # all_ab.dt[, osu_ab := osu_count/sum(osu_count)]
      all_ab.dt = all_ab.dt[order(-osu_count)]
      all_ab.dt[, sample_id := s]
      setcolorder(all_ab.dt, c('sample_id', 'osu_id', 'osu_count', 'species',
                               'pctsim'))
      all_ab.dt[order(-pctsim)]
  }))
  return(unique(all_abs.dt))
}

# REMOVE ME: prototyping code (do this for each sample)
# ns = length(dada_res)
# sid_to_pid_to_seq = vector('list', ns)
# for (s_id in 1:length(dada_res)) { # For each sample s_id
#   x = partid_to_fastqid(dada_res[[s_id]]$map-1, derep[[s_id]]$map-1)
#   fq_fil_data = sfastq_reader(fq_fil[s_id])
#   fq_mer_data = sfastq_reader(fq_mer[s_id])
#   pid_to_seq = vector("list", length(x))
#   for (i in 1:length(x)) { # for each partition i in sample s_id
#     metas = fq_fil_data[['meta']][x[i][[1]]+1]
#     mer_seqs = fq_mer_data[['seqs']][which(fq_mer_data[['meta']] %in% metas)]
#     pid_to_seq[i] = gsub('[-]{1,}$', '', consensus_sequence(mer_seqs))
#   }
#   sid_to_pid_to_seq[s_id] = pid_to_seq
# }


# for each partition id, extract all full-length merged sequences
# then summarize them as a consensus sequence
p_id = 95

# /REMOVE ME

himap = function (fq_fwd, fq_rev, # FASTQ input
                  # Output folder and file names
                  out, # Output folder
                  out_prefix = 'himap',
                  fq_mer='<<default>>', fq_fil='<<default>>', # Merged and filtered FASTQ output
                  fq_fil_minlen='<<default>>',
                  fq_prefix_split='_',
                  ab_tab_out='<<default>>', # abundance table
                  # Merging parameters
                  min_pct_sim=0.6, min_aln_len=50,
                  # PCR primer trimming parameters
                  pr_fwd='CCTACGGGNGGCWGCAG', pr_rev='GGATTAGATACCCBDGTAGTCC', # V3-V4 primers from Zheng et al 2015
                  # Filtering parameters
                  minimum_len=20, maxN=0, maxLen=Inf, maxEE=2, truncQ=2, verbose=T, multithread=T,
                  trimLeft=0, truncLen=0,
                  # DADA2 partitioning parameters
                  pvalThreshold=1e-4, # p-value threshold (before FDR correction)
                  # BLAST parameters
                  ref_db='/data1/igor/himap/db/V4_515F-805R_hang1_uniq',
                  ref_cp='/data1/igor/himap/db/V4_515F-805R_hang1_unique_sequences_table.txt',
                  ncpu=39, fa_denoised='<<default>>', pipeline_stats='<<default>>',
                  blast_output='<<default>>', blast_out_best='<<default>>',
                  blast_out_best_sp='<<default>>',
                  max_target_seqs=2000,
                  # Pipeline run parameters, set to false if we already ran this
                  run_merge=TRUE, run_trimmer=TRUE, run_filter=TRUE, run_chim=TRUE,
                  # Do we want the function to directly return output tables (don't need to read files)
                  return_tables=TRUE,
                  # Special parameters
                  mock=FALSE, # Is this a DNA mock community?
                  mock_size_file='/data1/igor/mockrobiota-13/v4_himap_genome_size.txt'
) {

  cat('--- HiMAP pipeline ---\n')
  cat(paste0('Output folder: ', out, '\n'))

  # Automatically generate output filenames if they are not provided
  # fq_prefix = strsplit(sub("(.*)\\..*$", "\\1", basename(fq_fwd)), fq_prefix_split, fixed=T)[[1]][1]
  fq_prefix = unname(
    sapply(
      sub("(.*)\\..*$", "\\1", basename(fq_fwd)),
      function (x) strsplit(x, fq_prefix_split, fixed=T)[[1]][1]
    )
  )
  fq_mer = ie(fq_mer=='<<default>>', file.path(out, 'merged', paste0(fq_prefix, '_merged.fastq')), fq_mer)
  fq_tri = ie(fq_tri=='<<default>>', file.path(out, 'trimmed', paste0(fq_prefix, '_trimmed.fastq')), fq_tri)
  fq_fil = ie(fq_fil=='<<default>>', file.path(out, 'filtered', paste0(fq_prefix, '_filtered.fastq')), fq_fil)
  fq_fil_minlen = ie(fq_fil_minlen=='<<default>>', file.path(out, 'filtered', paste0(fq_prefix, '_minlenfiltered.fastq')),
                     fq_fil_minlen)
  fa_denoised = ie(fa_denoised=='<<default>>', file.path(out, paste0(out_prefix, '_denoised_seqs.fasta')), fa_denoised)
  pipeline_stats = ie(pipeline_stats=='<<default>>', file.path(out, paste0(out_prefix, '_pipeline_reads_each_step.txt')),
                      pipeline_stats)
  blast_output = ie(blast_output=='<<default>>', file.path(out, paste0(out_prefix, '_blast_results.txt')),
                    blast_output)
  blast_out_best = ie(blast_out_best=='<<default>>', file.path(out, paste0(out_prefix, '_blast_best.txt')),
                      blast_out_best)
  blast_out_best_sp = ie(blast_out_best_sp=='<<default>>', file.path(out, paste0(out_prefix, '_blast_best_sp.txt')),
                         blast_out_best_sp)
  ab_tab_out = ie(ab_tab_out=='<<default>>', file.path(out, paste0(out_prefix, '_ab_tab_nochim_m.txt')), ab_tab_out)
  # maxQscore = ie(maxQscore=='<<default>>', )
  cat('Filtered:', fq_fil, '\n')

  # Merge reads
  if (run_merge) {
    cat('Merging reads...')
    if (!dir.exists(dirname(fq_mer[1]))) dir.create(dirname(fq_mer[1]), recursive=T)
    mergestats = merge_pairs(fq_fwd, fq_rev, fq_mer, min_pct_sim=0.6, min_aln_len=50)
    mergestats = as.data.table(t(mergestats))
    cat('OK.\n')
  }

  # PCR primer filter reads
  if (run_trimmer) {
    cat('Trimming PCR primers...')
    if (!dir.exists(dirname(fq_tri[1]))) dir.create(dirname(fq_tri[1]), recursive=T)
    pcr_primer_trimmer(fq_mer, fq_tri, pr_fwd=pr_fwd, pr_rev=pr_rev)
    cat('OK.\n')
  }

  # Filter reads
  if (run_filter) {
    if (!dir.exists(dirname(fq_fil[1]))) dir.create(dirname(fq_fil[1]), recursive=T)
    cat('Filtering reads...')
    filterAndTrim(fq_tri, fq_fil, maxLen=maxLen, maxN=maxN, maxEE=maxEE,
                  truncQ=truncQ, trimLeft=trimLeft,  minLen=minimum_len,
                  truncLen=truncLen, verbose=verbose, multithread=multithread)
    cat('OK.\n')
  }

  # Dereplicate merged reads
  cat('Dereplicating reads...')
  derep = derepFastq(fq_fil)
  if (class(derep)[1] == 'derep') {
    # Only one sample is loaded, so put the derep object into a list
    derep = list(derep)
  }
  cat('OK.\n')

  # Learn errors
  cat('DADA2: Learn errors...')
  dada_errors = learnErrors(fq_fil, multithread=T)
  cat('OK.\n')

  # Plot errors
  # plotErrors(dada_errors, nominalQ=TRUE)

  # Find number of unique sequences
  derep_no_uniq = sapply(derep, function (d) length(d$uniques))
  # derep_no_uniq = length(derep$uniques)

  # Run dada2
  # <= 10^5 unique sequences / sample, max Q score 1e-4, so try OMEGA_A = 1e-9
  cat('Running DADA2...')
  dada_res = dada(derep, err=dada_errors, multithread=multithread, OMEGA_A=pvalThreshold/max(derep_no_uniq))
  cat('OK.\n')

  # Make a copy during debugging of pasting sequences
  # dada_res2 = copy(dada_res)

  # If we ran just 1 sample, then dada_res is a dada-class object, not a list. Make it a list.
  if (class(dada_res)[1] == 'dada') {
    dada_res = list(dada_res)
  }

  # Retrieve original full-length sequences from merged data for further analysis
  # For each partition cluster, we use the dada2 sequence and the consensus of
  # left and right overhangs from actual sequences, concatenated together,
  # if we did any trimming (we need to otherwise dada2 does not work).
  cat('Retrieving full-length sequences...')
  for (s_id in 1:length(dada_res)) { # For each sample s_id
    x = partid_to_fastqid(dada_res[[s_id]]$map-1, derep[[s_id]]$map-1)
    # Load both merged (for retrieval) and filtered sequences for referencing
    # Filtered sequences are enumerated in dada2, so this file is used to
    # extract meta-data information for each read that we need to look up
    # in the merged file.
    fq_fil_data = sfastq_reader(fq_fil[s_id])
    fq_mer_data = sfastq_reader(fq_tri[s_id])
    pid_to_seq = vector("list", length(x))
    for (i in 1:length(x)) { # for each partition i in sample s_id
      metas = fq_fil_data[['meta']][x[i][[1]]+1]
      mer_seqs = fq_mer_data[['seqs']][which(fq_mer_data[['meta']] %in% metas)]
      # Select only the right hanging part
      mer_seqs_right = str_sub(mer_seqs, start=truncLen+1)
      pid_to_seq[i] = gsub('[-]{1,}.*', '', consensus_sequence(mer_seqs_right))
      # if (grepl('-', pid_to_seq[i])) print(i)
    }
    # sid_to_pid_to_seq[s_id] = pid_to_seq
    # Concatenate left consensus, dada2 middle partition sequence and right consensus
    dada_res[[s_id]]$sequence = paste(dada_res[[s_id]]$sequence,
                                      pid_to_seq,
                                      sep='')
    names(dada_res[[s_id]]$denoised) = dada_res[[s_id]]$sequence
    cat('.')
  }
  cat('OK.\n')

  # Generate abundance table
  cat('Generating abundance table...')
  ab_tab = makeSequenceTable(dada_res)
  cat('OK.\n')

  # Remove chimeras
  if (run_chim) {
    cat('Removing chimeras...')
    ab_tab_nochim = removeBimeraDenovo(ab_tab, method='consensus', allowOneOff=T,
                                       multithread=multithread, verbose=verbose)
    cat('OK.\n')
  } else {
    ab_tab_nochim = ab_tab
  }

  # Combine sequences which are just offset shifted
  cat('Adding up counts for offset shifts...')
  ab_tab_nochim = collapseNoMismatch(ab_tab_nochim)
  cat('OK.\n')

  # Track reads through pipeline
  nreads_orig = sapply(fq_fwd, function (fq) length(readFastq(fq)))
  nreads_merg = sapply(fq_mer, function (fq) length(readFastq(fq)))
  nreads_trim = sapply(fq_tri, function (fq) length(readFastq(fq)))
  # nreads_filt = sum(getUniques(dada_res))
  nreads_filt = sapply(dada_res, function (d) sum(getUniques(d)))
  nreads_tabl = rowSums(ab_tab)
  nreads_noch = rowSums(ab_tab_nochim)

  # Save this statistics in a separate file
  # put nreads_* into pipeline_stats.dt then save it to pipeline_stats file
  cat('Saving statistics...')
  nreads.dt = data.table(sample_id=fq_prefix, original=nreads_orig, merged=nreads_merg,
                         trimmed=nreads_trim,
                         filtered=nreads_filt, tabled=nreads_tabl,
                         nonchimeric=nreads_noch)

  write.table(nreads.dt, pipeline_stats, sep='\t', row.names=F, quote=F)
  cat('OK.\n')

  # Write the denoised sequences to FASTA file for Blasting
  write_seqs_to_fasta(dimnames(ab_tab_nochim)[[2]], fa_denoised)

  cat('Running BLAST...')
  blast_out = blast_seqs_to_reference(fa_denoised, ref_db, ncpu=ncpu, output=blast_output,
                                      max_target_seqs=max_target_seqs)
  cat('OK.\n')

  # Now load this output into data table, then find best matches for every dada2
  # denoised sequence.
  cat('Calculating alignment scores...')
  blast_out.dt = fread(blast_output)
  # qseqid sseqid length mismatch gapopen qstart qend sstart send qseq sseq
  names(blast_out.dt) = strsplit(blast_out_fmt, ' ')[[1]]
  #aln = blast_out.dt[, mcmapply(compare_alignment, qseq, sseq, mc.cores=ncpu)]

  # c('score', 'match', 'mismatch', 'gapopen', 'gapextend') :=
  # Calculate BLAST scores.
  blast_out.dt[, c('score', 'match', 'mismatch', 'gapopen', 'gapextend') := transpose(unname(mcmapply(function (q,s) {
    aln = compare_alignment(q,s)
    as.integer(c(sum(aln_params*aln), aln[1], aln[2], aln[3], aln[4]))
  }, qseq, sseq, SIMPLIFY=F, mc.cores=ncpu)))]

  # For local alignments that are one-off at each end, reduce the score using
  # a mismatch penalty. (We deliberately use indel score that is very similar.)
  blast_out.dt[qstart == 2, score := score-aln_params[2]]
  blast_out.dt[qstart == 2, mismatch := mismatch+1L]
  blast_out.dt[qend == qlen-1, score := score-aln_params[2]]
  blast_out.dt[qend == qlen-1, mismatch := mismatch+1L]

  # Filter out alignments with large gaps at beginning or the end
  blast_out.dt = blast_out.dt[qstart <= 2 & qend >= qlen-1]

  cat('OK.\n')

  # Filter out alignments with short alignment length
  cat('Processing BLAST results...')
  # blast_out.dt = blast_out.dt[aln_len >= blast_aln_len]

  # For each dada2 sequence, keep only the best matches by bitscore
  blast_best.dt = blast_best_seq_matches(blast_out.dt, strsplit(blast_out_fmt, ' ')[[1]][1],
                                         strsplit(blast_out_fmt, ' ')[[1]][2])
  blast_best.dt[, variant_id := gsub('^([0-9]+)-.*', '\\1', sseqid)]
  # Remove sequence alignment columns since we already calculated alignment score
  blast_best.dt[, c('qseq', 'sseq') := NULL]
  blast_best.dt[, c('sseqid') := NULL]
  # Calculate percentage similarity
  blast_best.dt[, pctsim := round(100*match/(match+mismatch+gapopen+gapextend), 2)]

  cat('Checkpoint1. Copy number information...')
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

  cat('Checkpoint2.')
  # Merge with blast hits
  # First do perfect matches

  # Original variant_ids depend how the reduced database is made. That is, they represent unique 16S regions
  # including the overhangs, which do not exactly match the lengths of sequenced regions. Here we just generate
  # new variant_id that pools together all sequences that are different in the extended region but the same
  # in the reduced region that was sequenced (and new variant_id is taken as just the first in the list).
  blast_best.dt = blast_best.dt[, {
    vids = unique(variant_id)
    merge(.SD, cp.dt[variant_id %in% vids, .(variant_id, strain, copy_number)], all.x=T, by='variant_id')
  }, by=qseqid]
  # blast_best.dt[, variant_id := variant_id[1], by=qseqid]

  # Find variant ids from the database that all map to the same shorter sequence in the data and pool together
  # those variant ids into a new one.
  variant_id_groups = list()
  for (q in blast_best.dt[, unique(qseqid)]) {
    variant_id_groups[[q]] = blast_best.dt[qseqid==q, unique(variant_id)]
    # First, for exact strain hits add the copy numbers
    cp.dt[variant_id %in% blast_best.dt[qseqid==q, unique(variant_id)], variant_id := variant_id[1]]
    blast_best.dt[qseqid==q, variant_id := variant_id[1]]
  }

  # Add up copy numbers for strain that differ in the undetected overhang
  cp.dt[, copy_number := sum(copy_number), by=.(variant_id, strain_id)]
  cp.dt = unique(cp.dt)
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

  # Regenerate variant_ids for the blast table.


  cat('Checkpoint3.')
  cat('OK.\n')

  # Now convert dada2 results into a data table, each column is a different condition
  # each row
  cat('Converting and saving abundance tables...')
  if (nrow(ab_tab_nochim) == 1) dimnames(ab_tab_nochim)[[1]] = c('1')
  ab_tab_nochim.dt = as.data.table(unname(ab_tab_nochim))
  # ab_tab_nochim.dt[, sample_id := 1]
  ab_tab_nochim.dt[, sample_id := sapply(strsplit(dimnames(ab_tab_nochim)[[1]], fq_prefix_split, fixed=T), `[`, 1)]

  # Now melt the data table
  ab_tab_nochim_m.dt = melt(ab_tab_nochim.dt, id.vars = 'sample_id',
                            variable.name = 'dada2_seqid', value.name = 'raw_count')

  # Cleanup dada2_seqid
  ab_tab_nochim_m.dt[, qseqid := as.numeric(gsub('V', '', dada2_seqid))]
  ab_tab_nochim_m.dt[, dada2_seqid := NULL]

  # Now add raw counts for each sample
  # Find the maximum number of reads across all samples, then set adjusted count such that
  # when we add all fake counts (zeros), they do not add up to more than 1% of the total.
  ab_tab_nochim_m.dt[, adj_count := as.numeric(raw_count)]
  max_sum_raw_counts = ab_tab_nochim_m.dt[, .(N = sum(raw_count)), by=sample_id][, max(N)]
  ab_tab_nochim_m.dt[raw_count==0, adj_count := 0.01*max_sum_raw_counts/.N, by=sample_id]

  # Normalize the adjusted counts for each sample
  ab_tab_nochim_m.dt[, adj_ab := adj_count/sum(adj_count), by=sample_id]

  # Now add raw_count to blast best results, then select strains for each sample id
  # separately.
  sample_ids = ab_tab_nochim_m.dt[, unique(sample_id)]
  ab_blast.dt = rbindlist(lapply(sample_ids, function (s) {
      merge(ab_tab_nochim_m.dt[sample_id==s, .(sample_id, qseqid, raw_count)],
            blast_best.dt,
            by='qseqid', all.x=T)
    }))

  # ab_blast.dt = merge(ab_tab_nochim_m.dt[, .(sample_id, qseqid, raw_count)],
  #                     blast_best.dt,
  #                     by='qseqid', all.x=T)

  # Ignore zero raw_count for each individual sample_id
  ab_blast.dt = ab_blast.dt[, .SD[raw_count > 0], by=sample_id]
  # ab_blast.dt[, rrn := gsub('.*_@rrn([0-9]+)$', '\\1', strain_name)]

  # Save final tables
  write.table(ab_tab_nochim_m.dt, ab_tab_out,
              sep='\t', quote=F, row.names=F)
  # write.table(blast_best_sp.dt, blast_out_best_sp,
  #             sep='\t', quote=F, row.names=F)
  write.table(blast_best.dt, blast_out_best,
              sep='\t', quote=F, row.names=F)
  cat('OK.\n')

  if (mock) {
    size.dt = fread(mock_size_file)
    ab_tab_nochim_m.dt = merge(ab_tab_nochim_m.dt, size.dt, by='dada2_seqid', all.x=T)
    mean_genome_size = size.dt[, mean(genome_sizes_nt)]
    mean_copy_number = size.dt[, mean(copy_number)]
    ab_tab_nochim_m.dt[!is.na(genome_sizes_nt), raw_count_renorm :=
                raw_count*genome_sizes_nt]
    ab_tab_nochim_m.dt[is.na(genome_sizes_nt), raw_count_renorm :=
                raw_count*mean_genome_size]
    ab_tab_nochim_m.dt[, raw_count_renorm := sum(raw_count)*raw_count_renorm/sum(raw_count_renorm)]
    ab_tab_nochim_m.dt[, raw_count := as.integer(raw_count_renorm)]
    # ab_tab.dt[is.na(copy_number), copy_number := 1L]
    #ab_tab.dt[, raw_count := as.integer(raw_count*copy_number)]
    # Remove extra columns
    ab_tab_nochim_m.dt[, c('raw_count_renorm') := NULL]
  }

  ab_tab_nochim_m.dt[, raw_count := as.integer(raw_count)]

  # Abundance estimation part

  #
  blast_reduced.dt = unique(blast_best.dt[, .(variant_id, qseqid, pctsim)])[order(qseqid)]
  blast_reduced.dt[, variant_id_new := variant_id[1], by=qseqid]
  rows_for_updating = blast_reduced.dt[, which(variant_id != variant_id_new)]

  # Combine abundance and BLAST information into one table
  # results.dt = merge(ab_tab.dt[, .(qseqid=dada2_seqid, raw_count)],
  #                    blast_best.dt[, .(variant_id, qseqid, pctsim)])

  cat('Generating reference table...')
  # Add copy number information for sequences with full genome
  # ref_cp.dt = fread(ref_cp, select=1:3, colClasses=c('character', 'character', 'integer',
  #                                                   'character'))
  # ref_cp.dt[, strain := gsub('_@rrn[0-9]+', '', strain_name)]

  # Generate strain index
  #ref_cp.dt[, strain_id := .GRP, by=strain]
  #ref_cp.dt[, strain_name := NULL]

  # Relabel variant ids based on the shorter experimental sequences
  # for (r in rows_for_updating) {
  #   old_varid = blast_reduced.dt[r, variant_id]
  #   new_varid = blast_reduced.dt[r, variant_id_new]
  #   ref_cp.dt[variant_id==old_varid, variant_id := new_varid]
  #   # Update copy number
  #   ref_cp.dt[variant_id==new_varid, copy_number := sum(copy_number), by=strain]
  # }
  # ref_cp.dt = unique(ref_cp.dt)
  #
  # # Update the reduced BLAST reference and cleanup old variant IDs
  # blast_reduced.dt[, variant_id := variant_id_new]
  # blast_reduced.dt[, variant_id_new := NULL]
  # blast_reduced.dt = unique(blast_reduced.dt)
  #
  # # Now we have 1-to-1 mapping between qseqid and variant_id
  #

  cat('OK.\n')


  # Generate a separate reference table
  osu.dt = unique(cp.dt[, .(osu_id, spectrum, strain, no_strains_in_osu)])
  osu.dt[, no_variants := sapply(spectrum, function (s) length(strsplit(s, ',', fixed=T)[[1]]))]

  # Now generate a smaller OSU table containing only osu id that contain
  # at least one variant ID from the data with > 0 raw count.
  pctsim_min = 99.00

  varids_data = blast_reduced.dt[pctsim >= pctsim_min, unique(variant_id)]
  osu_data.dt = osu.dt[, {
    varids = strsplit(gsub('[0-9]+:', '', spectrum), ',')[[1]]
    overlap = intersect(varids_data, varids)
    if (length(overlap) > 0) .SD
  }, by=osu_id]

  cat('Generating OSU table...')

  osu_data_m.dt = rbindlist(lapply(sample_ids, function (s) {
    # Melt OSU data spectrum
    osu_data.dt[, { # For each OSU id
      # Extract variant IDs and respective copy numbers from spectrum string
      cp_vid = strsplit(strsplit(spectrum, ',')[[1]], ':')
      cp = as.integer(sapply(cp_vid, `[`, 1))
      vid = sapply(cp_vid, `[`, 2)
      # For each variant id, extract dada2_seqid
      did = unname(sapply(vid, function (v) {
        blast_reduced.dt[variant_id==v & pctsim==max(pctsim), if (.N==0) NA else qseqid]
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
        no_strains_in_osu=unique(no_strains_in_osu))
    }, by=osu_id]
  }))

  # Some osu_id map to a single variant, just with a different copy number. For these
  # calculate the average copy number.
  # osu_data_m2.dt = osu_data_m.dt[, {
  #
  # }, by=.(sample_id, variant_id)]

  cat('OK.\n')

  # Remove OSU IDs for which we do not have all detected variants
  # osu_data_m.dt = osu_data_m.dt[, if (all(raw_count > 0L)) .SD, by=osu_id]

  # For variant_ids with multiple copy number calculate the average copy number
  # as integer, since there is no way we can tell?

  cat('Estimating OSU abundances...')
  # Now convert this table into 2 matrices, one with copy number coefficients of
  # each variant_id for each OSU id:
  osu_ab_m.dt = rbindlist(lapply(sample_ids, function (s) {
    Ab.dt = dcast(osu_data_m.dt[sample_id==s], variant_id ~ osu_id, value.var='copy_number',
                  fill=0)
    Ab.dt = merge(Ab.dt, unique(osu_data_m.dt[sample_id==s, .(variant_id, raw_count)]), by='variant_id')

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
    osu_ids = as.integer(names(osu_ab[which(osu_ab > osu_th)]))
    osu_ab.dt = data.table(osu_id=osu_ids, osu_count=osu_count)
    osu_ab.dt[, osu_ab := osu_count/sum(osu_count)]
    osu_ab.dt = osu_ab.dt[order(-osu_count)]

    # Generate a list of unique species for each osu id
    osu_sp.dt = cp.dt[, .(species = paste(unique(species), collapse=',')), by=osu_id]

    # Join table to OSU abundances
    osu_ab.dt = merge(osu_ab.dt, osu_sp.dt, by='osu_id',
                      all.x=T)
    osu_ab.dt = osu_ab.dt[order(-osu_count)]
    osu_ab.dt[, sample_id := s]
    # osu_ab.dt

    # Pool together OSUs in which


    # Try manually without constraints
    At = t(A)
    alpha = 0.1
    C = solve(At %*% A + alpha^2 * diag(ncol(A)))
    sol2 = C %*% (At %*% B)
    # Nope. Lots of very negative values.

    # Heaviside theta function
    H = function(x) as.numeric(x>0)

    f = function (x, A, B, lambda) {
      # Return an optimization function
      eps = A %*% x - B
      #t(eps) %*% eps + lambda*sum(x)
      # t(eps) %*% eps + lambda*sum(log(x+1))
      t(eps) %*% eps + lambda^2*sum(H(x))
    }

    # Start with the LSEI solution, then numerically minimize cost function f
    # to potentially reduce the number of OSU bins. For Zheng-2015 data this reduces 25
    # OSU to ...

    # Reduced OSU count vector
    x_red = osu_ab[which(osu_ab > osu_th)]
    # Reduced A matrix (row of A is a sequence, col of A is OSU)
    A_red = A[, which(osu_ab > osu_th)]



    fgrad = function (x, A, B, lambda) {
      # Return a gradient of the optimization function
      # sapply(1:ncol(A), function (j) t(2*(A %*% x - B)) %*% as.matrix(A[,j]) + 2*lambda*x[j])
      sapply(1:ncol(A), function (j) {
        # t(2*(A %*% x - B)) %*% as.matrix(A[,j]) + 2*lambda/(x[j]+1)
        t(2*(A %*% x - B)) %*% as.matrix(A[,j]) + lambda*length(x[x>0])
        #t(2*(A %*% x - B)) %*% as.matrix(A[,j]) + lambda
      })
    }

    # use osu_ab from lsei as good starting values.
    lambda = 1


  }))
  cat('OK.\n')

  # 01-04-2018: 14:39
  # Right now the problem is that the minimization procedure produces many OSUs
  # to fit the exact counts as best as possible. We need to add a cost to adding
  # too many OSUs.



  # Write OSU abundance table together with reference strains table to file

  # Return tables as output if this is set to TRUE
  if (return_tables) {
    return(list('abundance_table'=osu_ab.dt,
                'reference_strains'=ref_cp.dt[osu_id %in% osu_ab.dt[, osu_id],
                                              .(osu_id, no_strains_in_osu, no_species_in_osu,
                                                strain, species, spectrum)],
                'pipeline_stats'=nreads.dt)
    )
  }
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

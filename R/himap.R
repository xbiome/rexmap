# For this pipeline to work we will need to provide filenames for forward and
# reverse reads for each sample.

# Import these functions for everything
#' @importFrom data.table data.table
#' @importFrom data.table setkey
#' @importFrom data.table setcolorder
#' @importFrom data.table key
#' @importFrom dada2 filterAndTrim
#' @importFrom stringr str_sub
#' @useDynLib himap
NULL


.onAttach = function (libname, pkgname) {
  packageStartupMessage('HiMAP v1.0 loaded.')
}

# Default options -------------------------------------------------------------
himap_opts = new.env()
assign('himap_path', '', env=himap_opts)
# FUll current path is obtained by:
# assign()
# BLAST blast() output format
assign('blast_out_fmt',
       'qseqid sseqid qlen length qstart qend sstart send slen qseq sseq',
       env=himap_opts)
# BLAST collapse() format
assign('blast_coll_fmt',
       'qseqid sseqid qlen slen length qstart qend sstart send pident',
       env=himap_opts)
# BLAST paths
assign('path_makeblastdb', 'makeblastdb', env=himap_opts)
assign('path_blastn', 'blastn', env=himap_opts)
# BLAST alignment parameters
assign('aln_params', c(5L, -4L, -8L, -6L), env=himap_opts)
# Autodetect number of available threads for multithreading parts
# Set to 1 if you always  want to use only 1 thread.
assign('ncpu', parallel::detectCores(), env=himap_opts)
# BLAST databases
assign('blast_dbs',
       data.table::fread(system.file('inst', 'database', 'pcr_primers_table.txt',
                   package='himap'), header=T, sep='\t'),
       env=himap_opts)
assign('blast_max_seqs', 2000, env=himap_opts)
assign('blast_word_size', 13, env=himap_opts)
# Read merging
assign('mergepairs_matchqs',
       system.file('inst', 'merge_tables', 'himap_mergepairs_match_qs.txt',
                   package='himap'),
       env=himap_opts)
assign('mergepairs_mismatchqs',
       system.file('inst', 'merge_tables', 'himap_mergepairs_mismatch_qs.txt',
                   package='himap'),
       env=himap_opts)

# Interface to load HiMAP default options
himap_option = function (option_names) {
  if (is.null(option_names)) return(ls(himap_opts))

  if(!all(option_names %in% ls(himap_opts))) {
    warning("Invalid  option: ", option_names[!(option_names %in% ls(himap_opts))])
    option_names = option_names[option_names %in% ls(himap_opts)]
  }
  if (length(option_names) == 0) stop("Invalid options.")
  get(option_names, env=himap_opts)
}

himap_setoption = function (option_name, value) {
  # Simply set option_name to value. Used to change HiMAP defaults. Not finished yet.
  if (option_name == 'ncpu') {
    # Check that it is an integer
    if (!(class(value) == 'integer')) stop('ncpu must be an integer.')
  }
}



ie = function(test, yes, no) {
  if (test) {
    yes
  } else {
    no
  }
}


#' Generate a frequency table of sequence lengths in FASTQ files
#'
#' @param fastq_files A character vector of FASTQ filenames.
#' The resulting "table" object can be visuaized with base plot function.
sequence_length_table = function (fastq_files) {
  return(
    table(unlist(sapply(fastq_files, function (f) nchar(sfastq_reader(f)$seqs))))
  )
}

ftquantile = function (ft, prob) {
  # Return a quantile from a frequency table ft at given probability prob.
  ft_relsum = cumsum(ft)/sum(ft)
  return(as.integer(names(ft_relsum[ft_relsum >= prob])[1]))
}

#' Untrim trimmed sequences in a DADA2 object
#'
#' @param dada_res DADA2 object with sequences to un-trim.
#' @param derep DADA2 derep object that was used to obtain dada_res.
#' @param fq_tri A character vector of FASTQ filenames pre-global-trimming (but after PCR primer trim).
#' @param fq_fil A character vector of FASTQ filenames post-global-trimming and filtering.
#' @param truncLen Length used to trim.
#' @param verbose Boolean specifying whether to display progress bar.
#' @param ncpu Integer specifying number of CPU threads to use. This uses R package "parallel" so works only on macOS and Linux.
#'
add_consensus = function (dada_res, derep, fq_tri, fq_fil, truncLen,
                          verbose=T, ncpu=himap_option('ncpu')) {
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
    pid_to_seq = unlist(parallel::mclapply(1:length(x), function (i) {
      metas = fq_fil_data[['meta']][x[i][[1]]+1]
      mer_seqs = fq_mer_data[['seqs']][which(fq_mer_data[['meta']] %in% metas)]
      # Select only the right hanging part
      mer_seqs_right = str_sub(mer_seqs, start=truncLen+1)
      mer_seqs_cons = gsub('[-|N]{1,}.*', '', consensus_sequence(mer_seqs_right))
      names(mer_seqs_cons) = as.character(i)
      return(mer_seqs_cons)
    }, mc.cores=ncpu))
    pid_to_seq = pid_to_seq[order(as.integer(names(pid_to_seq)))]
    if (verbose) cat('OK. Update...')
    # Concatenate dada2 middle partition sequence and right consensus
    dada_res[[s_id]]$sequence = paste(dada_res[[s_id]]$sequence,
                                     pid_to_seq,
                                     sep='')
    names(dada_res[[s_id]]$denoised) = dada_res[[s_id]]$sequence
    if (verbose) cat('OK.\n')
  }
  return(dada_res)
}

#' BLAST FASTA file against a 16S database
#'
#' Aligns each sequence from the FASTA file against a reference database.
#' Either \code{region} or both \code{ref_db} and \code{ref_cp} arguments must be specificed.
#' If \code{region}
#' is specified, then the \code{ref_db} and \code{ref_cp} arguments are ignored and a
#' reference database is chosen
#' from the pre-computed set matching that hyper-variable region.
#' This function returns a blast-class object (named list; see ?blast-class for info).
#'
#' @param fasta_file FASTA file to BLAST.
#' @param region Hyper-variable region. If this is not NULL, \code{ref_db} is ignored. Possible
#' values: 'V4', 'V3-V4'.
#' @param ref_db Full path to a custom database.
#' @param ref_cp Full path to a copy number table for a custom database.
#' @param max_target_seqs Maximum number of target sequences to return for each query.
#' Note that the temporary output file will be large if this is set > 1000.
#' @param word_size Word size used for BLAST alignment.
#'
#' @export
#'
#'
blast = function (fasta_file, blast_output, region=NULL, ref_db=NULL,
                  ref_cp=NULL, max_target_seqs=himap_option('blast_max_seqs'),
                  word_size=himap_option('blast_word_size')) {
  # Run Blastn
  # Run blast_out_to_best_cp
  # Return himap object
  blast_status = blastn(fa_denoised, blast_output, region=region,
                        max_target_seqs=max_target_seqs,
                        word_size=word_size)
  if (blast_status != 0) stop('blast: error running blastn.')
  # Load BLAST results (can take a while if there are lots of sequences and max_target_seqs
  # is large.
  if (!file.exists(blast_output)) stop('blast: ', blast_output, ' file does not exist.')
  blast_cp = blast_out_to_best_cp(blast_output, region=region)
  names(blast_cp) = c('alignments', 'cp')
  blast_cp$parameters = list(
    'max_target_seqs'=max_target_seqs, 'word_size'=word_size, 'alignment_parameters'=
    paste(c('match', 'mismatch', 'gap_open', 'gap_extend'), himap_option('aln_params'),
          collapse=', ', sep=': '))
  return(as(blast_cp, 'blast'))
}

blast_out_to_best_cp = function (
   blast_output, region=NULL, ref_cp=NULL, aln_params=himap_option('aln_params'),
   blast_out_fmt=himap_option('blast_out_fmt'),
   verbose=T, ncpu=himap_option('ncpu'), variant_sep='-'
  )
{
  # Select region all ref_cp from file
  if (!is.null(region)) { # If region is given, ignore ref_cp
    ref_cp = himap_option('blast_dbs')[Hypervariable_region==region, table]
    ref_cp = system.file('inst', 'database', ref_cp, package='himap')
    if (length(ref_cp) == 0) stop('blast: invalid hyper-variable region.')
  } else { # Region is not given, just check if ref_cp file exist
    if (is.null(ref_cp)) stop('blast: reference copy number table argument missing.')
    if (!file.exists(ref_cp)) stop('blast: reference copy number table file missing.')
  }

  if (verbose) cat('* blast out: ')
  blast_out.dt = data.table::fread(blast_output)
  names(blast_out.dt) = strsplit(blast_out_fmt, ' ')[[1]]
  # Calculate BLAST scores.
  blast_out.dt[, c('score', 'match', 'mismatch', 'gapopen', 'gapextend') := transpose(
    unname(parallel::mcmapply(function (q,s) {
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
  blast_best.dt = blast_best_seq_matches(
    blast_out.dt,
    strsplit(blast_out_fmt, ' ')[[1]][1],
    strsplit(blast_out_fmt, ' ')[[1]][2]
  )
  blast_best.dt[, variant_id := gsub('^([0-9]+)-.*', '\\1', sseqid)]
  # Remove sequence alignment columns since we already calculated alignment score
  blast_best.dt[, c('qseq', 'sseq') := NULL]
  blast_best.dt[, c('sseqid') := NULL]
  # Calculate percentage similarity
  blast_best.dt[, pctsim := round(100*match/(match+mismatch+gapopen+gapextend), 2)]

  if (verbose) cat('OK. copy number table: ')
  cp.dt = data.table::fread(ref_cp, colClasses=c('character', 'character', 'integer'))

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


#' Convert abundance matrix to an abundance data table
#'
#' Each row in abundance matrix is a different sample, each column is a different
#' sequence.
ab_mat_to_dt = function (ab_tab_nochim, fq_prefix_split='_') {
   ab_tab_nochim.dt = as.data.table(unname(ab_tab_nochim))
   ab_tab_nochim.dt[, sample_id := sapply(strsplit(dimnames(ab_tab_nochim)[[1]], fq_prefix_split, fixed=T), `[`, 1)]
   ab_tab_nochim_m.dt = melt(ab_tab_nochim.dt, id.vars = 'sample_id',
                               variable.name = 'dada2_seqid', value.name = 'raw_count')
   ab_tab_nochim_m.dt[, qseqid := as.numeric(gsub('V', '', dada2_seqid))]
   ab_tab_nochim_m.dt[, dada2_seqid := NULL]
   ab_tab_nochim_m.dt = merge(
     ab_tab_nochim_m.dt,
     data.table(qseqid=1:ncol(ab_tab_nochim_coll), sequence=dimnames(ab_tab_nochim_coll)[[2]])
   )
   return(ab_tab_nochim_m.dt[])
}

#' Saves sequences from HiMAP sequence abundance table to FASTA file
#'
#' @param remove_from_table If TRUE, column with sequences (names sequences)
#' is removed from the data table after the FASTA file is written to disk.
sequences_to_fasta = function (abundance_table, fasta_out, remove_from_table=F) {
  if ('sequence' %in% names(abundance_table)) {
    with(unique(abundance_table[, .(qseqid, sequence)])[order(qseqid)],
      fasta_writer(
        1:length(sequence),
        sequence,
        fasta_out
      )
    )
    if (remove_from_table) abundance_table[, sequence := NULL]
  } else {
    warning('Sequences already removed. Nothing to do.')
  }
}


# OSU binning and abundance estimation
blast_cp_to_osu_dt = function (
  blast_best.dt, cp.dt, ab_tab_nochim_m.dt,
  pctsim_min = 100.00, # Min. pct sim for OSU binning
  ncpu = himap_option('ncpu'),
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

  # Optimized osu_data_m2.dt
  osu_data2.dt = copy(osu_data.dt)
  osu_data2.dt = osu_data.dt[,
   transpose(
      strsplit(strsplit(spectrum, ',', fixed=T)[[1]], ':', fixed=T)),
   by=.(osu_id, strain, no_strains_in_osu, no_variants)]
  names(osu_data2.dt)[5:6] = c('copy_number', 'variant_id')

  # Add BLAST alignment information to each variant_id
  osu_data3.dt = merge(
     osu_data2.dt,
     blast_reduced.dt[pctsim >= pctsim_min, .(variant_id, qseqid, pctsim)],
     by='variant_id', all.x=T)

  # Remove strain annotation here
  osu_data4.dt = copy(osu_data3.dt)
  osu_data4.dt[, strain := NULL]
  osu_data4.dt = unique(osu_data4.dt)

  # Fill in zeros for some OSUs where we have some, but not all variants
  sample_ids = ab_tab_nochim_m.dt[, unique(sample_id)]
  ab_tab_nochim_m_fill.dt = rbindlist(list(
     ab_tab_nochim_m.dt,
     data.table(sample_id=sample_ids, raw_count=0, qseqid=NA)))

  # Add abundance information
  osu_data5.dt = merge(
     osu_data4.dt,
     ab_tab_nochim_m_fill.dt,
     by='qseqid', all.x=T
  )

  osu_data5.dt = osu_data5.dt[, .(osu_id, copy_number, variant_id, raw_count, sample_id, pctsim, no_strains_in_osu)]
  osu_data5.dt = osu_data5.dt[order(sample_id, osu_id, variant_id)]

  if (verbose) cat(' OK.')
  return(osu_data5.dt[, .(osu_id, copy_number, variant_id, raw_count, sample_id, pctsim, no_strains_in_osu)])
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

#' Return abundances of each sequence from DADA2 result object
#'
#' @param dada_result dada() result
#' @param remove_bimeras Check and remove bimeric reads? (default: TRUE)
#' @param collapse_sequences Should we check for sequences that differ up to
#' shifts and add up the counts? (default: TRUE)
#' @param remove_bimeras_method Method to remove bimeras. See ?dada2::removeBimeraDenovo
#' for more options.
#'
sequence_abundance = function (dada_result, remove_bimeras=T, collapse_sequences=T,
                               verbose=T, remove_bimeras_method='consensus',
                               remove_bimeras_oneoff=T, fq_prefix_split='.') {
  # Extract sequences and their counts from the dada class result
  ab_tab = dada2::makeSequenceTable(dada_res2)
  if (remove_bimeras) {
    ab_tab_nochim = dada2::removeBimeraDenovo(
      ab_tab, method=remove_bimeras_method, allowOneOff=remove_bimeras_oneoff,
      multithread=ie(himap_option('ncpu') > 1, T, F), verbose=verbose)
  } else {
    ab_tab_nochim = ab_tab
  }
  if (collapse_sequences) ab_tab_nochim_coll = collapse(ab_tab_nochim, verbose=verbose)
  else ab_tab_nochim_coll = ab_tab_nochim
  return(ab_mat_to_dt(ab_tab_nochim_coll, fq_prefix_split=fq_prefix_split))
}

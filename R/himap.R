# For this pipeline to work we will need to provide filenames for forward and
# reverse reads for each sample.

# Import these functions for everything
#' @importFrom data.table data.table
#' @importFrom data.table as.data.table
#' @importFrom data.table setkey
#' @importFrom data.table setcolorder
#' @importFrom data.table setorder
#' @importFrom data.table key
#' @importFrom data.table copy
#' @importFrom dada2 filterAndTrim
#' @importFrom stringr str_sub
#' @useDynLib himap
NULL


.onAttach = function (libname, pkgname) {
  options('datatable.prettyprint.char'=30)
  packageStartupMessage('HiMAP v1.0 loaded.')
}

#' Detect operating system
detect_os = function () {
  switch(
    Sys.info()[['sysname']],
      Windows= {'win'},
      Linux  = {'linux'},
      Darwin = {'osx'}
  )
}

exec_file = function (filename) {
  platform = detect_os()
  if (platform == 'osx') return(filename)
  else if (platform == 'linux') return(paste0(filename, '_linux'))
  else if (platform == 'win') return(paste0(sub('.exe', '', filename),
                                            '_win.exe'))
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
assign('path_makeblastdb', exec_file('makeblastdb'), env=himap_opts)
assign('path_blastn', exec_file('blastn'), env=himap_opts)
# BLAST alignment parameters
assign('aln_params', c(5L, -4L, -8L, -6L), env=himap_opts)
# Autodetect number of available threads for multithreading parts
# Set to 1 if you always  want to use only 1 thread.
assign('ncpu', parallel::detectCores(), env=himap_opts)
# BLAST databases
assign('blast_dbs',
       data.table::fread(
         system.file('extdata', 'pcr_primers_table.txt', package='himap')
         # data('pcr_primers_table.txt', package='himap')
       ),
       env=himap_opts)
assign('blast_max_seqs', 2000, env=himap_opts)
assign('blast_word_size', 13, env=himap_opts)
# Read merging
assign('mergepairs_matchqs',
       system.file('merge_tables', 'himap_mergepairs_match_qs.txt',
                   package='himap'),
       env=himap_opts)
assign('mergepairs_mismatchqs',
       system.file('merge_tables', 'himap_mergepairs_mismatch_qs.txt',
                   package='himap'),
       env=himap_opts)
# Taxonomy
assign('taxonomy_file',
       system.file('database', 'lineages-2017-03-17_bacteria_archaea',
                   package='himap'),
       env=himap_opts)

# Data.table adjustments
assign('string_maxwidth', 50, env=himap_opts)

# assign('maxrows', 12, env=himap_opts)
#' HiMAP options
#'
#' @param option_names A string or a vector of strings with available options. If
#' not given, then the function lists available options.
#'
#' Options:
#'
#'
#' @export
himap_option = function (option_names=NULL) {
  if (is.null(option_names)) {
    cat('HiMAP: available options', fill=T)
    return(ls(himap_opts))
  }
  if(!all(option_names %in% ls(himap_opts))) {
    warning("Invalid  option: ", option_names[!(option_names %in% ls(himap_opts))])
    option_names = option_names[option_names %in% ls(himap_opts)]
  }
  if (length(option_names) == 0) stop("Invalid options.")
  # get(option_names, env=himap_opts)
  if (length(option_names) == 1) get(option_names, env=himap_opts)
  else sapply(option_names, get, env=himap_opts)
}



himap_setoption = function (option_name, value) {
  # Simply set option_name to value. Used to change HiMAP defaults. Not finished yet.
  if (option_name == 'ncpu') {
    # Check that it is an integer
    if (!(class(value) == 'integer')) stop('ncpu must be an integer.')
  }
  if (option_name == 'datatable.prettyprint.char') {
    # Check that value is integer
    options('datatable.prettyprint.char'=value)
  }
  # Else just assign stuff
  assign(option_name, )
}

#' Frequency table of sequence lengths
#'
#' Generate a frequency table of sequence lengths from FASTQ files.
#' The resulting table object can be visuaized with base plot function.
#'
#' @param fastq_files A character vector of FASTQ filenames.
#' @examples
#'
#'
#' @export
sequence_length_table = function (fastq_files) {
  return(
    table(unlist(sapply(fastq_files, function (f) nchar(sfastq_reader(f)$seqs))))
  )
}

#' Quantile from frequency table
#'
#' Outputs a quantile corresponding to a given probability \code{prob}
#' from a frequency table \code{ft}. We use this to find the minimum sequence length
#' above which we have certain fraction of reads.
#'
#' Similar to the base R code \code{\link{quantile}}, but works with
#' frequency table instead of the raw list.
#'
#' @param ft Frequency table of sequence lengths. Output from \code{link{{sequence_length_table}}.
#' @param prob Probability.
#'
#' @examples
#' # Find minimum sequence length, above which we have 99% reads
#' seqlen.ft = sequence_length_table(fq_tri)
#' ftquantile(seqlen.ft, 0.01)
#'
#' @export
ftquantile = function (ft, prob) {
  # Return a quantile from a frequency table ft at given probability prob.
  ft_relsum = cumsum(ft)/sum(ft)
  return(as.integer(names(ft_relsum[ft_relsum >= prob])[1]))
}


#' DADA2 denoising
#'
#' Uses DADA2 algorithm to denoise trimmed FASTQ files \code{fastq_trimmed}.
#' DADA2 partitions sequences based on their quality score profiles, using a
#' \code{pvalue_adjusted} parameter (\code{OMEGA_A} in the paper) as a cutoff
#' for forming new partitions. (See DADA2 paper
#' [Callahan et al. Nature Methods 2016] for details.)
#' If \code{pvalue}
#' is given (default), then the \code{pvalue_threshold} is calculated by dividing
#' the \code{pvalue} by the maximum number of unique sequences from all samples.
#' Otherwise \code{pvalue_threshold} is used directly.
#'
#' After denoising, the untrimmed parts of sequences are retrieved from
#' \code{fastq_untrimmed} files and for each partition a consensus sequence
#' of that part is concatenaded to its end, omitting everything after
#' (and including) any indel and N.
#'
#' @param pvalue P-value threshold before Bonferroni correction for DADA2.
#' @param pvalue_adjusted P-value threshold after Bonferroni correction for DADA2.
#' If left as NULL, adjusted pvalue is calculated by taking \code{pvalue} and
#' dividing it by the maximum number of unique sequences across all samples.
#' @param multithread TRUE/FALSE or the number of CPU threads to use for multithreading.
#' Does not work on Windows due to \code{parallel} package implementation.
#' @param use_intermediate Use saved intermediate files. Used for resuming very
#' long runs.
#' @param verbose TRUE/FALSE: display of status messages.
#'
#' @export
dada_denoise = function (fastq_trimmed, fastq_untrimmed,
                         pvalue=1e-4,
                         pvalue_adjusted=NULL,
                         save_intermediate=TRUE,
                         use_intermediate=FALSE,
                         multithread=himap_option('ncpu'),
                         verbose=T, error_estimation_nsamples=3) {
  # Dereplicate reads into a derep object
  # Check whether intermediate folder exists
  # Load 1 sample at a time, process it, then combine

  # Learn errors
  if (verbose) cat('* learn errors')
  dada_derep = dada2::derepFastq(
    fastq_trimmed[1:min(length(fastq_trimmed), error_estimation_nsamples)])
  if (class(dada_derep) != 'list') dada_derep = list(dada_derep)
  cat('...')
  if (is.null(pvalue_adjusted)) {
    # Find the maximum number of unique sequences across all samples
    max_num_uniques = max(
      sapply(fastq_trimmed, function (f) length(dada2::derepFastq(f)$uniques))
    )
    pvalue_adjusted = pvalue/max_num_uniques
  }
  dada_errors = suppressWarnings(dada2::learnErrors(
    fastq_trimmed, multithread=multithread, OMEGA_A=pvalue_adjusted
  ))
  if (verbose) cat(' OK.', fill=T)
  dada_results = list()
  for (i in 1:length(fastq_trimmed)) {
    fq = fastq_trimmed[i]
    fqt = fastq_untrimmed[i]
    if (verbose) cat('* processing ', fq, fill=T)
    dada_derep = list(dada2::derepFastq(fq))
    if (is.null(pvalue_adjusted)) {
      # Find the maximum number of unique sequences across all samples
      num_uniques = length(dada_derep[[1]]$uniques)
      # Adjust the p-value by an estimate for the total number of comparisons in dada2
      pvalue_adjusted = pvalue/num_uniques^2
    }
    dada_res = list(suppressWarnings(dada2::dada(dada_derep, err=dada_errors, OMEGA_A=pvalue_adjusted, multithread=T)))
    # Extract length to which reads have been truncated
    trunclen = nchar(dada_res[[1]]$sequence[1])
    dada_res = add_consensus(dada_res, dada_derep, fqt, fq, truncLen=trunclen,
                              ncpu=max(1, as.integer(multithread)), verbose=verbose)
    dada_results[[length(dada_results)+1]] = dada_res[[1]]
    names(dada_results)[length(dada_results)] = basename(fq)
  }
  return(dada_results)
}


#' Untrim trimmed sequences in a DADA2 object
#'
#' @param dada_res DADA2 object with sequences to un-trim.
#' @param derep DADA2 derep object that was used to obtain dada_res.
#' @param fq_tri A character vector of FASTQ filenames pre-global-trimming
#' (but after PCR primer trim).
#' @param fq_fil A character vector of FASTQ filenames post-global-trimming and filtering.
#' @param truncLen Length used to trim.
#' @param verbose Boolean specifying whether to display progress bar.
#' @param ncpu Integer specifying number of CPU threads to use. This uses R package "parallel"
#' so works only on macOS and Linux.
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


#' Convert abundance matrix to an abundance data table
#'
#' Each row in abundance matrix is a different sample, each column is a different
#' sequence.
#'
#' @export
ab_mat_to_dt = function (ab_tab_nochim, fq_prefix_split='_') {
   ab_tab_nochim.dt = as.data.table(unname(ab_tab_nochim))
   ab_tab_nochim.dt[, sample_id := sapply(strsplit(dimnames(ab_tab_nochim)[[1]],
                                                   fq_prefix_split, fixed=T), `[`, 1)]
   ab_tab_nochim_m.dt = data.table::melt(ab_tab_nochim.dt, id.vars = 'sample_id',
                               variable.name = 'dada2_seqid', value.name = 'raw_count')
   ab_tab_nochim_m.dt[, qseqid := as.integer(gsub('V', '', dada2_seqid))]
   ab_tab_nochim_m.dt[, dada2_seqid := NULL]
   ab_tab_nochim_m.dt = merge(
     ab_tab_nochim_m.dt,
     data.table(qseqid=1:ncol(ab_tab_nochim), sequence=dimnames(ab_tab_nochim)[[2]])
   )
   setcolorder(ab_tab_nochim_m.dt, c('sample_id', 'qseqid', 'raw_count', 'sequence'))
   return(ab_tab_nochim_m.dt[])
}

#' Saves sequences from HiMAP sequence abundance table to FASTA file
#'
#' @param abundance_table Output from \code{\link{sequence_abundance}} function.
#' @param remove_from_table If TRUE, column with sequences (names sequences)
#' is removed from the data table after the FASTA file is written to disk.
#'
#' @export
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


pctsim_range = function (p) return(max(p, na.rm=T))


#' Combine BLAST object and a sequence abundance table into OSU table
#'
#' @param abundance_table Sequence abundance table
#' @param blast_object BLAST output class, output from \code{\link{blast}} function.
#' @param verbose TRUE/FALSE: display status messages
#' @param raw_strains TRUE/FALSE: whether to report full strain information for each OSU
#' @param ncpu Integer specifying number of CPU threads to use. This uses R package "parallel"
#' (TRUE) or to simplify the output when multiple strains of the same species are in the
#' same OSU (FALSE).
#'
#' @importFrom igraph make_empty_graph
#' @importFrom igraph add_vertices
#' @importFrom igraph add_edges
#' @importFrom igraph as.undirected
#' @importFrom igraph groups
#' @importFrom igraph clusters
#' @importFrom limSolve lsei
#' @importFrom pso psoptim
#' @importFrom data.table dcast
#'
#' @export
abundance = function (abundance_table, blast_object, ncpu=himap_option('ncpu'), verbose=T,
                      raw_strains=FALSE) {
  # Generate OSU data table first
  osu_data_m.dt = blast_cp_to_osu_dt(
    blast_best.dt=blast_object$alignments,
    cp.dt=blast_object$cp,
    ab_tab_nochim_m.dt=abundance_table[, 1:3],
    ncpu=ncpu,
    verbose=verbose
  )
  # Now generate osu abundance table
  osu_ab.dt = osu_cp_to_all_abs(abundance_table[, 1:3],
                                blast_object$alignments,
                                blast_object$cp,
                                osu_data_m.dt,
                                ncpu=ncpu, verbose=verbose,
                                raw=raw_strains)
  setcolorder(osu_ab.dt, c('sample_id', 'osu_id', 'osu_count', 'pctsim', 'species'))
  setorder(osu_ab.dt, sample_id, -osu_count)
  return(osu_ab.dt)
}


osu_cp_to_all_abs = function (ab_tab_nochim_m.dt,
                              blast_best.dt,
                              cp.dt,
                              osu_data_m.dt,
                              ncpu=himap_option('ncpu'), verbose=T, seed=42,
                              raw=T) {

  pctsim_min = 100
  osu_offset = 1000000L

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

  osu_data_m_single.dt = unique(osu_data_m.dt[, if (length(unique(variant_id)) == 1) .SD,
                                              by=osu_id][raw_count > 0, .(osu_id, variant_id, copy_number)])

  osu_sp.dt = cp.dt[, .(species = print_strains(strain, raw=raw)), by=osu_id]

  sample_ids = ab_tab_nochim_m.dt[, unique(sample_id)]
  all_abs.dt = data.table::rbindlist(parallel::mclapply(sample_ids, function (s) {

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
      data.table::setorder(osu_ab.dt, osu_id)

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

          x = as.matrix(dcast(osu_data_m_single.dt[variant_id %in% variant_ids],
                              variant_id ~ as.character(osu_id), value.var='copy_number')[, -1])
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
      data.table::setorder(osu_ab4.dt, -osu_count)


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
      data.table::setcolorder(all_ab.dt, c('sample_id', 'osu_id', 'osu_count', 'species',
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
#' @export
sequence_abundance = function (dada_result, remove_bimeras=T, collapse_sequences=T,
                               verbose=T, remove_bimeras_method='consensus',
                               remove_bimeras_oneoff=T, fq_prefix_split='.') {
  # Extract sequences and their counts from the dada class result
  if (verbose) cat('* generating sequence table...')
  ab_tab = dada2::makeSequenceTable(dada_result)
  if (verbose) cat(' OK.', fill=T)
  if (remove_bimeras) {
    if (verbose) cat('* removing bimeras...')
    ab_tab_nochim = dada2::removeBimeraDenovo(
      ab_tab, method=remove_bimeras_method, allowOneOff=remove_bimeras_oneoff,
      multithread=ie(himap_option('ncpu') > 1, T, F), verbose=verbose)
    if (verbose) cat(' OK.', fill=T)
  } else {
    ab_tab_nochim = ab_tab
  }
  if (collapse_sequences) {
    if (verbose) cat('* adding together sequences that differ in shifts on lengths...')
    ab_tab_nochim_coll = collapse(ab_tab_nochim, verbose=verbose)
    if (verbose) cat(' OK.', fill=T)
  }
  else ab_tab_nochim_coll = ab_tab_nochim
  return(ab_mat_to_dt(ab_tab_nochim_coll, fq_prefix_split=fq_prefix_split))
}


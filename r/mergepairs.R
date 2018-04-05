# suppressPackageStartupMessages(library(ShortRead))

# Merge pairs wrapper


# C_mergepairs = function(read_fwd, read_rev, qual_fwd, qual_rev, match, mismatch, gap_p,
#                         path_posterior) {
#   .Call('C_mergepairs',
#         # PACKAGE = 'dada2',
#         read_fwd, read_rev, qual_fwd, qual_rev, match, mismatch, gap_p,
#         path_posterior)
# }

# HiMAP erge pairs function
himap_merge_pairs = function (fq_fwd, fq_rev, fq_mer, min_pct_sim=0.75, min_aln_len=50,
                              match=5L, mismatch=-5L, gap_p=-7L, rc_reverse=TRUE,
                              threads = 10, verbose=FALSE, timing=FALSE,
                              path_posterior = '/data1/igor/himap/data/') {
  # fq_fwd = vector of filenames (incl. paths) to forward reads
  # fq_rev = vector of filenames (incl. paths) to reverse reads
  # fq_mer = vector of filenames (incl. paths) to write merged reads
  # min_pct_sim = minimum similarity of the overlap aligned region, *100 to get pct (default: 0.75)
  # min_aln_len = minimum overlap aligned region length in nts (default: 50)
  # By default it overwrites existing fq_mer files.
  # Ends-free alignment parameters
  #  match = (int) score for matching ACGT (default: 5)
  #  mismatch = score for mismatching ACGT or alignment vs N (default: -2)
  #  gap_p = gap opening and extension penalty (default: -7)
  #  rc_reverse = (bool) reverse complement reverse reads before alignment? (default: TRUE)
  #  threads = number of cores to use for multi-threading
  #   (tries to detect total number of cores using parallel::detectCores())
  # verbose = (bool) print status messages? (default: FALSE)
  # timing = (bool) print execution time to stderr (default: FALSE)
  # path_precomputed_posterior = path to tab-delimited headerless tables with precomputed
  #  quality scores for posterior probabilities for read merging
  #  output from mergepairs_generate_posterior_probabilities.py
  if (timing) {
    start_time = Sys.time()
  }

  if (is.na(threads) | threads==0) {
    # If we can't figure out number of threads, set it to 2
    threads = 2
  }
  m = function (...) {
    # Print message if verbose is enabled.
    if (verbose) {
      message(..., appendLF=F)
    }
  }

  # Load reads using ShortRead::FastqStreamer
  m('Loading FASTQ reads...')
  f_fwd = FastqStreamer(fq_fwd)
  f_rev = FastqStreamer(fq_rev)
  # Go through each read entry, do alignment
  r_fwd = yield(f_fwd)
  if (length(r_fwd) == 0) break
  if (rc_reverse) {
    r_rev = reverseComplement(yield(f_rev))
  } else {
    r_rev = yield(f_rev)
  }

  # Process chunk
  read_fwd = as.character(sread(r_fwd))
  read_rev = as.character(sread(r_rev))
  qual_fwd = as.character(quality(quality(r_fwd)))
  qual_rev = as.character(quality(quality(r_rev)))
  ids = gsub('^([^ ]+) .*', '\\1', as.character(id(r_fwd)))
  m(' OK.\n')

  # Apply mergepairs
  m('Merging pairs...')
  merged_list = mcmapply(C_mergepairs, read_fwd, read_rev, qual_fwd, qual_rev,
                         match=match, mismatch=mismatch, gap_p=gap_p,
                         min_pct_sim=min_pct_sim, min_aln_len=min_aln_len,
                         posterior_match_file='/data1/igor/himap/data/himap_mergepairs_match_qs.txt',
                         posterior_mismatch_file='/data1/igor/himap/data/himap_mergepairs_mismatch_qs.txt',
                         mc.cores=threads
                       )
  m(' OK.\n')
  merged_aln_filter = as.logical(unname(merged_list[3, ])) &
    as.logical(unname(merged_list[4, ]))

  # Free memory
  close(f_fwd)
  close(f_rev)
  rm(read_fwd, read_rev, qual_fwd, qual_rev, f_fwd, f_rev, r_fwd, r_rev)
  m('Writing output files...')
  merged_sread = ShortReadQ(
    sread = DNAStringSet(unname(merged_list[1, merged_aln_filter])),
    quality = BStringSet(unname(merged_list[2, merged_aln_filter])),
    id = BStringSet(ids[merged_aln_filter])
  )

  # Generate statistics for filtered-out reads
  stats = list(
    'total'=length(merged_aln_filter),
    'low_pct_sim'=length(merged_aln_filter) - sum(as.logical(merged_list[3, ])),
    'low_aln_len'=length(merged_aln_filter) - sum(as.logical(merged_list[4, ]))
  )

  # Free memory and close file connections
  rm(merged_list)

  # If file exists, delete it, then write new one
  if (file.exists(fq_mer)) {
    file_remove_result = file.remove(fq_mer)
  }
  writeFastq(merged_sread, fq_mer, compress = F)
  m(' OK.\n')
  if (timing) {
    end_time = Sys.time()
    m('Finished in ', round(as.numeric(end_time-start_time)/60), ' m ',
      round(as.numeric(end_time-start_time)%%60, 1), ' s.\n')
  }
  return(stats)
}

# Vectorized himap merge pairs
merge_pairs = Vectorize(himap_merge_pairs, c('fq_fwd', 'fq_rev', 'fq_mer'))


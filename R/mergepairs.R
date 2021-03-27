#' Merge paired-end reads
#'
#' Merges reads using overlap version of the Needleman-Wunsch alignment algorithm and calculate
#' merged quality scores based on the posterior probabilities (cite: Edgar, Fyelberg)
#'
#' @param fq_fwd A character vector for input forward FASTQ files (can be gzipped).
#' @param fq_rev A character vector for input reverse FASTQ files (can be gzipped).
#' @param fq_mer A character vector for output merged FASTQ files.
#' @param min_sim Minimum similarity to accept an alignment. A floating point number between 0 and 1. Default: 0.75, corresponding to minimum 75\% similarity in the aligned overlap region.
#' @param min_aln_len Minimum alignment length. Ignore alignments below this alignment length.
#' @param match Score for character match in the alignment (default: 5.
#' @param mismatch Score for character mismatch in the alignment (default: -5)
#' @param gap_p Score for gap penalty, either gap opening on gap extension (default: -7).
#' @param rc_forward TRUE/FALSE Reverse complement forward reads? (default: FALSE)
#' @param rc_reverse TRUE/FALSE Reverse complement reverse reads? (default: TRUE)
#' @param ncpu Number of CPU threads to use for multithreading.
#' @param force Do not stop merging on fatal errors. (default: FALSE)
#' @param verbose TRUE/FALSE Display of status messages.
#' @param timing TRUE/FALSE Time merging.
#'
#' @export
merge_pairs = function (
  fq_fwd, fq_rev, fq_mer, min_sim=0.75, min_aln_len=50,
  match=5L, mismatch=-5L, gap_p=-7L, rc_forward=FALSE, rc_reverse=TRUE,
  ncpu=rexmap_option('ncpu'), ncpu_samples=1,
  force=FALSE, verbose=TRUE, timing=FALSE
  ) {
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

  m('------ Merge reads -------', verbose=verbose)
  m('Running with parameters:', verbose=verbose)
  m('Min:', 100*min_sim, '% sim |', 'Min aln len:', min_aln_len, 'nt | RC reverse:', rc_reverse, verbose=verbose)
  m('Run', ncpu, 'files in parallel, ', ncpu_samples, 'threads/file.', verbose=verbose)

  # m2 = function (..., fill=TRUE, time_stamp=TRUE, verbose=TRUE) {
  #   if (ncpu > 1) {
  #
  #   } else {
  #     m(.., fill=fill, time_stamp=time_stamp, verbose=verbose)
  #   }
  # }
  empty_result = data.table(
    fq_fwd=NA, fq_rev=NA,
    total=NA, merged=NA, low_pct_sim=NA, low_aln_len=NA)

  out_per_sample = parallel::mcmapply(
    function (fq_fwd_i, fq_rev_i, fq_mer_i) {
      start_time = Sys.time()
      # Does output folder exist?
      output_folders = unique(dirname(fq_mer_i))
      for (out_f in output_folders) {
        if (!dir.exists(out_f)) dir.create(out_f, recursive=T)
      }

      # Load reads using ShortRead::FastqStreamer
      # m(paste0('Loading FASTQ reads: ', basename(fq_fwd_i), ', ',
      #          basename(fq_rev_i), ' ...'))
      m_buffer = ''
      if (ncpu == 1) {
        m(paste0('* ', basename(fq_fwd_i), '+', basename(fq_rev_i), ':'), fill=F, time_stamp=F)
      } else {
        m_buffer = paste0(m_buffer, '* ', basename(fq_fwd_i), '+', basename(fq_rev_i), ':')
      }

      f_fwd = tryCatch(
        ShortRead::FastqStreamer(fq_fwd_i),
        error = function (x) NA
      )
      f_rev = tryCatch(
        ShortRead::FastqStreamer(fq_rev_i),
        error = function (x) NA
      )
      if (class(f_fwd) != 'FastqStreamer') {
        if (is.na(f_fwd)) {
          m('Warning: corrupted FWD file', basename(fq_fwd_i), ': Skipping.', verbose=verbose)
          return(empty_result)
        }
      }
      if (class(f_rev) != 'FastqStreamer') {
        if (is.na(f_rev)) {
          m('Warning: corrupted FWD file', basename(fq_fwd_i), ': Skipping.', verbose=verbose)
          return(empty_result)
        }
      }
      # Go through each read entry, do alignment
      if (rc_forward) {
        r_fwd = tryCatch(
          ShortRead::reverseComplement(ShortRead::yield(f_fwd)),
          error = function (x) NA
        )
      } else {
        r_fwd = tryCatch(
          ShortRead::yield(f_fwd),
          error = function (x) NA
        )
      }
      if (length(r_fwd) == 0) {
        m('Warning: No reads in the forward file ', basename(fq_fwd_i), 'found. Skipping.', verbose=verbose)
        return(empty_result)
      }
      if (rc_reverse) {
        r_rev = tryCatch(
          ShortRead::reverseComplement(ShortRead::yield(f_rev)),
          error = function (x) NA
        )
      } else {
        r_rev = tryCatch(
          ShortRead::yield(f_rev),
          error = function (x) NA
        )
      }
      if (class(r_fwd) != 'ShortReadQ') {
        if (is.na(r_fwd)) {
          m('Warning: corrupted FWD file', fq_fwd_i, ': Skipping.', verbose=verbose)
          return(empty_result)
        }
      }
      if (class(r_rev) != 'ShortReadQ') {
        if (is.na(r_rev)) {
          m('Warning: corrupted REV file', fq_rev_i, ': Skipping.', verbose=verbose)
          return(empty_result)
        }
      }

      # Process chunk
      read_fwd = as.character(ShortRead::sread(r_fwd))
      read_rev = as.character(ShortRead::sread(r_rev))
      qual_fwd = as.character(Biostrings::quality(Biostrings::quality(r_fwd)))
      qual_rev = as.character(Biostrings::quality(Biostrings::quality(r_rev)))
      ids = gsub('^([^ ]+) .*', '\\1', as.character(ShortRead::id(r_fwd)))
      # m('  Loaded.', fill=F, time_stamp=F)
      if (ncpu == 1) {
        m(' Loaded.', fill=F, time_stamp=F, verbose=verbose)
      } else {
        m_buffer = paste0(m_buffer, ' Loaded.')
      }
      # Filter out useless/invalid reads at this point before sending it to C++
      # function.
      if (length(read_fwd) != length(read_rev)) {
        m(' * Error: Unequal number of reads in forward and reverse files.', verbose=verbose)
        m('   Forward file: ', length(read_fwd), ' reads | Reverse file: ',
          length(read_rev), ' reads.', verbose=verbose)
        if (!force) {
          m('   Stop.\n', verbose=verbose)
          return()
        } else {
          m('   Proceeding by ignoring extra reads.', verbose=verbose)
          read_cap = min(length(read_fwd), length(read_rev))
          read_fwd = read_fwd[1:read_cap]
          read_rev = read_rev[1:read_cap]
          qual_fwd = qual_fwd[1:read_cap]
          qual_rev = qual_rev[1:read_cap]
          ids = ids[1:read_cap]
        }

      }


      # m('Removing invalid reads...')
      # NNNNN reads
      n_mask = grepl('^[N]+$', read_fwd) & grepl('^[N]+$', read_rev)
      if (sum(n_mask) > 0) {
        read_fwd = read_fwd[!n_mask]
        read_rev = read_rev[!n_mask]
        qual_fwd = qual_fwd[!n_mask]
        qual_rev = qual_rev[!n_mask]
        ids = ids[!n_mask]
      }

      # Non-ACGTN reads
      x_mask = grepl('[^ACGTN]+', read_fwd) & grepl('[^ACGTN]+', read_rev)
      if (sum(x_mask) > 0) {
        read_fwd = read_fwd[!x_mask]
        read_rev = read_rev[!x_mask]
        qual_fwd = qual_fwd[!x_mask]
        qual_rev = qual_rev[!x_mask]
        ids = ids[!x_mask]
      }
      if (ncpu == 1) {
        m(' QC.', fill=F, time_stamp=F, verbose=verbose)
      } else {
        m_buffer = paste0(m_buffer, ' QC.')
      }

      # print(head(read_fwd))
      # print(head(read_rev))
      # print(head(qual_fwd))
      # print(head(qual_rev))


      # Apply mergepairs
      # m('Merging pairs...')
      merged_list = parallel::mcmapply(C_mergepairs, read_fwd, read_rev, qual_fwd, qual_rev,
                                       match=match, mismatch=mismatch, gap_p=gap_p,
                                       min_pct_sim=min_sim, min_aln_len=min_aln_len,
                                       posterior_match_file=rexmap_option('mergepairs_matchqs'),
                                       posterior_mismatch_file=rexmap_option('mergepairs_mismatchqs'),
                                       mc.cores=ncpu_samples
      )
      if (ncpu == 1) {
        m(' Merged ', fill=F, time_stamp=F, verbose=verbose)
      } else {
        m_buffer = paste0(m_buffer, ' Merged ')
      }

      # Filter out low % sim and low aln length alignments
      merged_aln_filter = as.logical(unname(merged_list[3, ])) &
        as.logical(unname(merged_list[4, ]))

      pct_merged = 100*sum(merged_aln_filter)/length(merged_aln_filter)
      if (ncpu == 1) {
        m(round(pct_merged, 1), '%.', fill=F, time_stamp=F, verbose=verbose)
      } else {
        m_buffer = paste0(m_buffer, round(pct_merged, 1), '%.')
      }

      # Free memory
      close(f_fwd)
      close(f_rev)
      rm(read_fwd, read_rev, qual_fwd, qual_rev, f_fwd, f_rev, r_fwd, r_rev)
      # m('Writing output files')
      # sequences
      final_seqs = unname(merged_list[1, merged_aln_filter])
      # m('.')
      final_qual = unname(merged_list[2, merged_aln_filter])
      # m('.')
      final_ids  = ids[merged_aln_filter]
      # m('.')

      final_filter = !is.null(final_seqs) & !is.null(final_qual)
      # if (!all(final_filter)) m('.')
      final_seqs = final_seqs[final_filter]
      final_qual = final_qual[final_filter]
      final_ids  = final_ids[final_filter]


      # Generate statistics for filtered-out reads
      stats = data.table(
        fq_fwd=fq_fwd_i, fq_rev=fq_rev_i,
        total=length(merged_aln_filter),
        merged=sum(merged_aln_filter),
        low_pct_sim=length(merged_aln_filter) - sum(as.logical(merged_list[3, ])),
        low_aln_len=length(merged_aln_filter) - sum(as.logical(merged_list[4, ]))
      )

      # Sometimes there are no sequences left, in which case just return
      # stats

      if (length(final_seqs) > 1) {
        merged_sread = ShortRead::ShortReadQ(
          sread = Biostrings::DNAStringSet(final_seqs),
          quality = Biostrings::BStringSet(final_qual),
          id = Biostrings::BStringSet(final_ids)
        )
      }


      # Free memory and close file connections
      rm(merged_list)

      # If file exists, delete it, then write new one
      if (file.exists(fq_mer_i)) {
        file_remove_result = file.remove(fq_mer_i)
      }
      if (length(final_seqs) > 1) {
        ShortRead::writeFastq(merged_sread, fq_mer_i, compress = F)
      } else {
        # m(' (No reads merged!) ', fill=F, time_stamp=F)
        if (ncpu == 1) {
          m(' (No reads merged) ', fill=F, time_stamp=F, verbose=verbose)
        } else {
          m_buffer = paste0(m_buffer, ' (No reads merged!)')
        }

      }
      # m(' OK.\n')
      end_time = Sys.time()
      dt = end_time - start_time


      if (ncpu == 1) {
        m(' [', round(dt, 1), ' ', attr(dt, 'units'), ']', fill=T, time_stamp=F, verbose=verbose)
      } else {
        m_buffer = paste0(m_buffer, ' [', round(dt, 1), ' ', attr(dt, 'units'), ']')
        m(m_buffer)
      }


      # if (timing) {
      #   end_time = Sys.time()
      #   diff_time = difftime(end_time, start_time, units='secs')
      #   m('Finished in ', round(as.numeric(diff_time)/60), ' m ',
      #     round(as.numeric(diff_time)%%60, 1), ' s.\n')
      # }
      return(stats)
    }, fq_fwd, fq_rev, fq_mer, SIMPLIFY=F, USE.NAMES=F, mc.cores=ncpu
  )
  out.dt = data.table::rbindlist(out_per_sample)
  # fq_fwd_processed = names(fq_fwd)
  # processed_filter = sapply(fq_fwd, function (x) x %in% fq_fwd_processed)
  # out.dt[, fq_fwd := basename(fq_fwd[processed_filter])]
  # out.dt[, fq_rev := basename(fq_rev[processed_filter])]
  # setcolorder(out.dt, c('fq_fwd', 'fq_rev', 'total', 'merged', 'low_pct_sim', 'low_aln_len'))
  return(out.dt)
}

#' Convert mergestats table to a normal data.table
#'
#' @export
mergeout_to_table = function (mergestats) {
  mergestats.dt = data.table::data.table(
    matrix(
      unlist(t(mergestats)),
      ncol=nrow(mergestats),
      dimnames=list(sampleids_from_filenames(colnames(mergestats), '_'), rownames(mergestats))
    ),
    keep.rownames='sample_id'
  )
  return(mergestats.dt)
}

#' Auto-detect overlap length in paired-end reads
#'
#' Tests different minimum overlap lengths for `merge_pairs()`
#'
#' @param min_sim Minimum % similarity in the overlap region to
#' accept the overlap alignment as valid
#' @param ncpu Number of threads
#'
#' @export
detect_overlap_length = function (
  fq_fwd, fq_rev,
  min_sim=0.75,
  match=5L, mismatch=-5L, gap_p=-7L,
  rc_forward=FALSE, rc_reverse=TRUE,
  ncpu = rexmap_option('ncpu'),
  nseqs=1000, nseqs_seed=42,
  minalnlen_min=10, minalnlen_max=70, minalnlen_step=5,
  minalnlen_drop_pct=5,
  force=T,
  verbose=T
) {
  # Try different min_aln_len and rc_reverse TRUE/FALSE if we cant merge much
  # fq_fwd and fq_rev are characters vectors with forward and reverse reads
  # sample nseqs each run to speed it up
  if (length(fq_fwd) != length(fq_rev)) {
    stop('fq_fwd and fq_rev have different lengths.')
  }
  fq_fwd_exist = file.exists(fq_fwd)
  fq_rev_exist = file.exists(fq_rev)
  if (any(!fq_fwd_exist)) {
    stop('Some forward read files do not exist: ', paste(fq_fwd[!fq_fwd_exist], collapse=','))
  }
  if (any(!fq_rev_exist)) {
    stop('Some reverse read files do not exist: ', paste(fq_rev[!fq_rev_exist], collapse=','))
  }

  m = function (..., fill=TRUE, time_stamp=TRUE) {
    if (!time_stamp) {
      cat(..., append=TRUE, fill=fill)
    } else {
      cat(as.character(as.POSIXlt(Sys.time())), ' | ', ...,
          append=TRUE, fill=fill)
    }
  }
  m('---- Detect forward/reverse overlap lengths ----')
  min_aln_lens = seq(minalnlen_min, minalnlen_max, minalnlen_step)
  empty_result = NA
  out = parallel::mcmapply(function(fqf, fqr) {
    # Load files fqf and fqr
    m(paste0('Loading FASTQ reads: ', basename(fqf), ', ', basename(fqr), ' ...'),
      fill=F)

    f_fwd = tryCatch(
      ShortRead::FastqStreamer(fqf),
      error = function (x) NA
    )
    f_rev = tryCatch(
      ShortRead::FastqStreamer(fqr),
      error = function (x) NA
    )
    if (class(f_fwd) != 'FastqStreamer') {
      if (is.na(f_fwd)) {
        m('Warning: corrupted FWD file', basename(fqf), ': Skipping.')
        return(empty_result)
      }
    }
    if (class(f_rev) != 'FastqStreamer') {
      if (is.na(f_rev)) {
        m('Warning: corrupted FWD file', basename(fqf), ': Skipping.')
        return(empty_result)
      }
    }
    # Go through each read entry, do alignment
    if (rc_forward) {
      r_fwd = tryCatch(
        ShortRead::reverseComplement(ShortRead::yield(f_fwd)),
        error = function (x) NA
      )
    } else {
      r_fwd = tryCatch(
        ShortRead::yield(f_fwd),
        error = function (x) NA
      )
    }
    if (length(r_fwd) == 0) {
      m('Warning: No reads in the forward file ', basename(fqf), 'found. Skipping.')
      return(empty_result)
    }
    if (rc_reverse) {
      r_rev = tryCatch(
        ShortRead::reverseComplement(ShortRead::yield(f_rev)),
        error = function (x) NA
      )
    } else {
      r_rev = tryCatch(
        ShortRead::yield(f_rev),
        error = function (x) NA
      )
    }
    if (class(r_fwd) != 'ShortReadQ') {
      if (is.na(r_fwd)) {
        m('Warning: corrupted FWD file', fqf, ': Skipping.')
        return(empty_result)
      }
    }
    if (class(r_rev) != 'ShortReadQ') {
      if (is.na(r_rev)) {
        m('Warning: corrupted REV file', fqr, ': Skipping.')
        return(empty_result)
      }
    }

    f_fwd = ShortRead::FastqStreamer(fqf)
    f_rev = ShortRead::FastqStreamer(fqr)
    # Go through each read entry, do alignment
    r_fwd = ShortRead::yield(f_fwd)
    if (length(r_fwd) == 0) {
      m('Warning: corrupted files. Skipping.')
      return(NA)
    }
    # if (rc_reverse) {
    #   r_rev = ShortRead::reverseComplement(ShortRead::yield(f_rev))
    # } else {
    #   r_rev = ShortRead::yield(f_rev)
    # }

    # Process chunk
    read_fwd = as.character(ShortRead::sread(r_fwd))
    read_rev = as.character(ShortRead::sread(r_rev))
    qual_fwd = as.character(Biostrings::quality(Biostrings::quality(r_fwd)))
    qual_rev = as.character(Biostrings::quality(Biostrings::quality(r_rev)))
    ids = gsub('^([^ ]+) .*', '\\1', as.character(ShortRead::id(r_fwd)))
    m(' OK.', time_stamp=F)

    # Filter out useless/invalid reads at this point before sending it to C++
    # function.
    if (length(read_fwd) != length(read_rev)) {
      m(' * Error: Unequal number of reads in forward and reverse files.')
      m('   Forward file: ', length(read_fwd), ' reads | Reverse file: ',
        length(read_rev), ' reads.')
      if (!force) {
        m('   Stop.')
        return(list('total'=NA, 'low_pct_sim'=NA, 'low_aln_len'=NA))
      } else {
        m('   Proceeding by ignoring extra reads.')
        read_cap = min(length(read_fwd), length(read_rev))
        read_fwd = read_fwd[1:read_cap]
        read_rev = read_rev[1:read_cap]
        qual_fwd = qual_fwd[1:read_cap]
        qual_rev = qual_rev[1:read_cap]
        ids = ids[1:read_cap]
      }

    }


    m('Removing invalid reads...', fill=F)
    # NNNNN reads
    n_mask = grepl('^[N]+$', read_fwd) & grepl('^[N]+$', read_rev)
    if (sum(n_mask) > 0) {
      read_fwd = read_fwd[!n_mask]
      read_rev = read_rev[!n_mask]
      qual_fwd = qual_fwd[!n_mask]
      qual_rev = qual_rev[!n_mask]
      ids = ids[!n_mask]
    }

    # Non-ACGTN reads
    x_mask = grepl('[^ACGTN]+', read_fwd) & grepl('[^ACGTN]+', read_rev)
    if (sum(x_mask) > 0) {
      read_fwd = read_fwd[!x_mask]
      read_rev = read_rev[!x_mask]
      qual_fwd = qual_fwd[!x_mask]
      qual_rev = qual_rev[!x_mask]
      ids = ids[!x_mask]
    }
    m(' OK.', time_stamp=F)

    # Subsample reads to nseqs or all reads if there are less than nseqs
    nseqs_max = min(nseqs, length(read_fwd))
    set.seed(nseqs_seed)
    read_sampler = sample(1:length(read_fwd), nseqs_max)
    # Sample reads
    read_fwd = read_fwd[read_sampler]
    read_rev = read_rev[read_sampler]
    qual_fwd = qual_fwd[read_sampler]
    qual_rev = qual_rev[read_sampler]

    # For each min_aln_len run this
    # Apply mergepairs
    pct_merged_list = lapply(min_aln_lens, function (min_aln_len) {
      m('* Minimum overlap length: ', min_aln_len, fill=F)
      # m('  Merging pairs...', fill=F)
      merged_list = parallel::mcmapply(
        C_mergepairs, read_fwd, read_rev, qual_fwd, qual_rev,
        match=match, mismatch=mismatch, gap_p=gap_p,
        min_pct_sim=min_sim, min_aln_len=min_aln_len,
        posterior_match_file=rexmap_option('mergepairs_matchqs'),
        posterior_mismatch_file=rexmap_option('mergepairs_mismatchqs'),
        mc.cores=ncpu
      )
      # m(' OK.', time_stamp=F)
      # Filter out low % sim and low aln length alignments
      merged_aln_filter = as.logical(unname(merged_list[3, ])) &
        as.logical(unname(merged_list[4, ]))

      final_seqs = unname(merged_list[1, merged_aln_filter])
      final_qual = unname(merged_list[2, merged_aln_filter])
      final_filter = !is.null(final_seqs) & !is.null(final_qual)
      final_seqs = final_seqs[final_filter]
      final_qual = final_qual[final_filter]
      pct_merged = 100*length(final_seqs)/nseqs_max
      m(' | ', round(pct_merged, 1), '% reads merged.', time_stamp=F)
      return(pct_merged)
    })

    pct_merged.dt = data.table(min_aln_len=min_aln_lens,
                               pct_merged=as.numeric(pct_merged_list))

    # Free memory
    close(f_fwd)
    close(f_rev)
    rm(read_fwd, read_rev, qual_fwd, qual_rev, f_fwd, f_rev, r_fwd, r_rev)
    return(pct_merged.dt)

  }, fq_fwd, fq_rev, mc.cores=1, SIMPLIFY=F, USE.NAMES=F)
  if (class(out) == 'logical') {
    if (is.na(out)) return(empty_result)
  }
  min_aln_lens_best = sapply(out, function (pm.dt) {
    max_pct_merged = pm.dt[, max(pct_merged)]
    max2_pct_merged = (1-minalnlen_drop_pct/100)*max_pct_merged
    best_min_aln_len = pm.dt[
      abs(pct_merged-max2_pct_merged)==min(
        abs(pct_merged-max2_pct_merged))
    ][, max(min_aln_len)]
    return(best_min_aln_len)
  })
  min_aln_len_best = round(mean(min_aln_lens_best, na.rm=T))
  m(' * min_aln_len with largest overlap: ', min_aln_len_best)

  return(min_aln_len_best)
}

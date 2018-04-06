
pcr_primer_filter = function (fq_in, fq_out, pr_fwd='CCTACGGGNGGCWGCAG',
                              pr_rev='GGATTAGATACCCBDGTAGTCC',
                              pr_fwd_maxoff=10, pr_rev_maxoff=10, return_noprimer=T,
                              multithread=T, ncpu=detectCores()-1, max_mismatch=2,
                              timing=F) {
  # fq_in = input fastq file
  # fq_out = output fastq file (without primers)
  # pr_fwd = forward primer 5'->3'
  # pr_rev = reverse primer 5'->3'
  #
  # Load Input FASTQ file
  if (timing) start_time = Sys.time()
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
    aln_stat = compare_alignment(
      stringr::str_sub(aln[1], start=pr_fwd_left, end=pr_fwd_right),
      stringr::str_sub(aln[2], start=pr_fwd_left, end=pr_fwd_right)
    )
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
    aln2_stat = compare_alignment(
      stringr::str_sub(aln2[1], start=pr_rev_left, end=pr_rev_right),
      stringr::str_sub(aln2[2], start=pr_rev_left, end=pr_rev_right)
    )
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
                'seqs' = stringr::str_sub(seq, start=start, end=end),
                'qual' = stringr::str_sub(qual, start=start, end=end),
                'trim_fwd' = pr_fwd_found,
                'trim_rev' = pr_rev_found))
  }

  # Apply fastq_trimmer() to each sequence in this file
  out_trimmed = unname(
    parallel::mcmapply(
      fastq_trimmer, in_fq[['meta']], in_fq[['seqs']], in_fq[['qual']],
      mc.cores=ncpu, SIMPLIFY=FALSE
    )
  )

  fwd_trimmed = sum(sapply(out_trimmed, function (x) x$trim_fwd))
  rev_trimmed = sum(sapply(out_trimmed, function (x) x$trim_rev))

  # Save results in a new file
  fastq_list_writer(out_trimmed, fq_out, ncpu=ncpu)
  if (timing) {
    end_time = Sys.time()
    diff_time = difftime(end_time, start_time, units='secs')
    cat('Finished in ', round(as.numeric(diff_time)/60), ' m ',
      round(as.numeric(diff_time)%%60, 1), ' s.\n')
  }
  return(c('fwd_trim'=fwd_trimmed, 'rev_trim'=rev_trimmed))
}

pcr_primer_trimmer = Vectorize(pcr_primer_filter, vectorize.args=c('fq_in', 'fq_out'))

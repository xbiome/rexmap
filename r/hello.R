# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

# this.dir = dirname(parent.frame(2)$ofile)

hello = function() {
  print("Hello, world!")
  #system2(file.path(this.dir, 'vsearch'))
}



# Generate merged reads for merged-first check
vsearch_mergepairs = function (fwd_file, rev_file, out_file, max_dif=20, minov=100, minlen=240,
                               vs_path='/home/isegota/bin/vsearch') {
  # Check if vsearch exists
  if (!file.exists(vs_path)) {
    cat(paste0('vsearch not found in ', vs_path))
    return(FALSE)
  }
  # Check if output folder exists (if not make it)
  if (!dir.exists(dirname(out_file))) {
    dir.create(dirname(out_file), recursive = TRUE)
  }
  # Run vsearch
  x = system2(vs_path, args = c('--fastq_mergepairs', fwd_file, '--reverse', rev_file,
                                '--fastq_maxdiffs', as.character(max_dif), '--fastq_minovlen',
                                as.character(minov), '--fastq_allowmergestagger',
                                '--fastq_minmergelen', as.character(minlen), '--fastqout', out_file,
                                '--quiet'), stdout=T, stderr=T)
  return(x)
}

vsearch_mergestats = function (vs_out) {
  # Take the result vector of vsearch_mergepairs() and return a data frame with
  # merging statistics. Use apply to apply this to all results.
  list('total'=as.numeric(str_extract(vs_out[1], '[0-9]+')),
       'merged'=as.numeric(str_extract(vs_out[2], '[0-9]+')),
       'not_merged'=as.numeric(str_extract(vs_out[3], '[0-9]+')),
       'merge_pct'=as.numeric(str_extract(vs_out[2], '[0-9]+\\.[0-9]+')))
}

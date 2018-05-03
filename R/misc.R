
#' Shortcut IF/ELSE function
#'
#' if \code{test} evaluates to TRUE, return \code{yes}. Otherwise return \code{no}.
#'
#' @export
#'
#' @examples
#' ie(1==2, '1 equals 2', '1 does not equal 2')
#' # Output: '1 does not equal 2'
ie = function(test, yes, no) {
  if (test) yes
  else no
}


#' Extract sample identifiers from FASTQ filenames
#'
#' @param filenames A character vector of filenames.
#' @param separator A single character used delimiting sample id from the
#' rest of FASTQ filename.
#'
#' @export
sampleids_from_filenames = function (filenames, separator='_') {
  # Get sample ids from FASTQ filenames
  sapply(strsplit(basename(filenames), separator, fixed=T), `[`, 1)
}

#' Read files from a folder satisfying a pattern.
#'
#' Shortcut for \code{sort(dir(path, pattern, full.names=T))}
#'
#' @param path Full path to the folder to list files in.
#' @param pattern Pattern for pattern matching. If '' then
#' list all files.
#'
#' @examples
#' read_files('~/data/diabimmune/fastq_tutorial', 'R1') # Get forward reads
#' read_files('~/data/diabimmune/fastq_tutorial', 'R2') # Get reverse reads
#'
#'
#' @export
read_files = function (path, pattern='') {
  sort(dir(path.expand(path), pattern, full.names=T))
}

#' Strain vector to shorter strain string
#'
#' Condense a character vector of strain names to a more managable
#' shorter string.
#'
#' @param strains (Required) A character vector of strain names, with underscore _ instead of space
#' to separate genus, species and strain designations.
#' @param raw FALSE/TRUE If TRUE, then the output is just a list of concatenated strain
#' names (can be very large!). Default: FALSE.
#' @param deduplicate FALSE/TRUE If strain name occurs multiple times, count only unique names.
#' Default: TRUE.
#' @param nmax Integer. If the number of strains is below \code{nmax}, just list all the strain names.
#'
#' @export
print_strains = function (strains, raw=T, deduplicate=T,
                          nmax=himap_option('print_strains_nmax')) {
  # Input: vector of strains
  # Prints a single string with reduced list of strains
  # such that any non-_bacterium or _sp. strain is shown
  # on a species level.

  if (deduplicate) uniq_strains = unique(strains)
  else uniq_strains = strains

  if (raw | length(uniq_strains) <= nmax) return(paste(uniq_strains, collapse=','))

  # Also add strains for species that occur only once!!
  if (length(strains)==1) return(strains)
  else {
    genuses = gsub('^[ ]*([^_]+)_.*', '\\1', uniq_strains)
    species = gsub('^[^_]+_([^_]+)[_]?.*', '\\1', uniq_strains)
    na_filter = grepl('^(sp\\.|bacterium)', species)
    # sp1 is the list of species that does not end with sp. or bacterium
    sp1 = paste(genuses[!na_filter], species[!na_filter], sep='_')
    # sp2 is the list of strains which do not have a species name (sp. or bacterium)
    sp2 = uniq_strains[na_filter]
    ft = as.table(sort(table(c(sp1, sp2)), decreasing=T))
    ft2 = ft[ft>1]
    # Get full strain names for species that occur only once
    sp_once = names(ft[ft==1])
    st_once = c()
    for (sp in sp_once) {
      st_once = c(st_once, grep(sp, uniq_strains, fixed=T, value=T))
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


print_fixed_length_string = function (x, len=30) {
   if (nchar(x) >= len+3) {
      ndots=3
      left = floor((len-ndots)/2)
      right = ceiling((len-ndots)/2)
      return(paste0(substr(x, 1, left), '...', substr(x, nchar(x)-right+1, nchar(x))))
   } else if (nchar(x) == len+2) {
      ndots=2
      left = floor((len-ndots)/2)
      right = ceiling((len-ndots)/2)
      return(paste0(substr(x, 1, left), '..', substr(x, nchar(x)-right+1, nchar(x))))
   } else if (nchar(x) == len+1) {
      ndots=1
      left = floor((len-ndots)/2)
      right = ceiling((len-ndots)/2)
      return(paste0(substr(x, 1, left), '.', substr(x, nchar(x)-right+1, nchar(x))))
   } else {
      return(x)
   }
}

string_fixed_len = function (x, len=30) {
  if (is.na(x)) return(NA)
  if (nchar(x) > len) return(paste0(substr(x, 1, len-3), '...'))
  else return(x)
}

print.data.table2 = function (dt, width=himap_option('string_maxwidth'),
                             topn=himap_option('maxrows')) {
   dt2 = head(dt, topn)
   extra_rows = nrow(dt) - topn
   for (j in names(dt2)) {
      if (dt2[, class(get(j))] == 'character') {
         values = dt2[[j]]
         set(dt2, j = j, value = sapply(values, string_fixed_len, len=width))
      }
   }
   print.data.frame(dt2)
   if (extra_rows > 0) cat('...', paste0('(', nrow(dt), ' rows)'), fill=T)
}

#' Reverse complement sequence string
#'
#' Take a character string \code{seq_string} and reverse complement it using
#' Biostrings package function \code{reverseComplement}. Supports extended
#' nucleotide code, so good for reverse complementing PCR primer sequences.
#'
#' Shortcut for the monstrostity of:
#' \code{as.character(Biostrings::reverseComplement(Biostrings::DNAString(seq_string)))}
#'
#' @export
reverse_complement = function (seq_string) {
   as.character(Biostrings::reverseComplement(Biostrings::DNAString(seq_string)))
}

all_exist = function (files) all(file.exists(files))

#' Generate random sequence of given length
#'
#' @export
random_sequences = function(len=50, n=10, dictionary=c('A', 'C', 'G', 'T')) {
  do.call(paste0, replicate(len, sample(dictionary, n, TRUE), FALSE))
}


#' Write data table to file
#'
#' Like base function \code{\link{write.table}} but with normal defaults.
#'
#' @export
write_table = function (table, output, sep='\t', verbose=F) {
  if (verbose) cat('* writing ', output, '...')
  write.table(table, output, sep=sep, quote=F, row.names=F)
  if (verbose) cat('OK.', fill=T)
}

#' Filter and trim reads. Wrapper for dada2::filterAndTrim
#'
#' @export
filter_and_trim = function (fq_in, fq_out, maxEE=2,
                            multithread=himap_option('ncpu'),
                            verbose=himap_option('verbose'),
                            compress=FALSE,
                            ...) {
  dada2::filterAndTrim(fq_in, fq_out, maxEE=maxEE, multithread=multithread,
                       verbose=verbose, compress=compress, ...)
}



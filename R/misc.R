
#' Shortcut IF/ELSE function
ie = function(test, yes, no) {
  if (test) yes
  else no
}


#' Extract sample identifiers from FASTQ filenames
#'
#' @param filenames A character vector of filenames.
#' @param separator A single character used delimiting sample id from the
#' rest of FASTQ filename.
sampleids_from_filenames = function (filenames, separator='_') {
  # Get sample ids from FASTQ filenames
  sapply(strsplit(basename(filenames), separator, fixed=T), `[`, 1)
}


read_files = function (path, pattern='') {
  sort(dir(path, pattern, full.names=T))
}

#' Strain vector to shorter strain string
#'
#' Condense a character vector of strain names to a more managable
#' shorter string.
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
reverse_complement = function (seq_string) {
   as.character(Biostrings::reverseComplement(Biostrings::DNAString(seq_string)))
}

all_exist = function (files) all(file.exists(files))



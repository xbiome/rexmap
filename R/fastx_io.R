# Read FASTA and FASTQ files regardless of number of characters in a row
# (a limitation of ShortRead::readFastq()).

file_reader = function (file_name) {
  # read entire file at once into a list. used for reading our FASTA database
  # Basically same runtime as C++ code.
  size = file.info(file_name)$size
  strsplit(readChar(file_name, size, useBytes=T), "\n", fixed=T, useBytes=T)[[1]]
}

fasta_reader = function (file_name) {
  x = file_reader(file_name)
  list('meta'=sub('>', '', x[c(T,F)], fixed=T, useBytes=T), 'seqs'=x[c(F,T)])
}

fastq_reader = function (file_name) {
  x = file_reader(file_name)
  list('meta'=x[c(T,F,F,F)], 'seqs'=x[c(F,T,F,F)], 'qual'=x[c(F,F,F,T)])
}

sfastq_reader = function (file_name) {
  x = ShortRead::readFastq(file_name)
  xx = ShortRead::sread(x)
  list('meta'= as.character(ShortRead::id(x)),
       'seqs'=as.character(xx),
       'qual'=as.character(Biostrings::quality(Biostrings::quality(x))))
}

fastq_list_writer = function (fastq_list, output, ncpu=detectCores()-1) {
  # Write FASTQ file at once. Input is a list where each element is a list
  # with names 'meta', 'seq' and 'qual'. I use this to save the result from (mc)mapply.
  list_paste = function (l) {
    paste(paste0('@', l$meta), l$seq, '+', paste0(l$qual, '\n'), sep='\n')
  }
  writeChar(paste(parallel::mclapply(fastq_list, list_paste, mc.cores=ncpu),
                  sep='\n'), output, eos=NULL)
}

fasta_writer = function (meta, seqs, output) {
  # Write FASTA sequences from meta data and sequences
  writeChar(paste(paste0('>', meta), seqs, sep='\n', collapse='\n'),
            output, eos=NULL)
}


# Old functions that might not be used?
strains_from_meta_data = function (meta_data_v, sep='$') {
  out = c()
  for (meta_data in meta_data_v) {
    meta_data = gsub('[\\$]{2,}', '\\$', meta_data)
    if (sep == '$') sep = '\\$'
    out = c(out, gsub('_@rrn[0-9]+', '', unlist(strsplit(meta_data, sep))))
  }
  unique(out)
}

get_blast_path = function () {
  # Find blastn path (assume its in system PATH)
  # First look for the package binary
  bin = system.file('exec', 'blastn', package='himap')
  # Binary file from the package not found. Fall back to system blastn
  if (bin == '') bin = suppressWarnings(system2('which', 'blastn', stdout=T))
  # If this is not found either, then just stop.
  if (attr(bin, 'status') == 1) stop('Error: makeblastdb command not found.')
  else return(bin)
}

blast_seqs_to_reference = function (seqs_fa, ref_db, blast_path=get_blast_path(),
                                    match=aln_params[1], mismatch=aln_params[2],
                                    gapopen=-aln_params[3], gapextend=-aln_params[4],
                                    word_size=13, ncpu=4, max_target_seqs=5000,
                                    query_coverage_pct=100, outfmt=paste0('6 ', blast_out_fmt),
                                    output='/data1/igor/zr_16s2/blast_test.txt',
                                    output_err=F) {
   x = system2(blast_path, args = c(
      '-dust', 'no', '-word_size', word_size,
      '-reward', match, '-penalty', mismatch, '-gapopen', gapopen, '-gapextend', gapextend,
      '-no_greedy', '-outfmt', shQuote(outfmt),
      '-query', seqs_fa, '-db', ref_db, '-num_threads', ncpu,
      '-max_target_seqs', max_target_seqs
      # '-qcov_hsp_perc', query_coverage_pct
   ), stdout = output, stderr = output_err)
   return(x)
}

# Use fasta_writer instead. This is very slow in comparison.
# write_seqs_to_fasta = function (seqs, fa_out) {
#    # Write a vector of sequences seqs to fa_out fasta file and index them
#    # sequentally.
#    for (i in seq_along(seqs)) {
#       if (i == 1) app = F
#       else app = T
#       cat('>', i, '\n', seqs[i], '\n', file=fa_out, sep='', append=app)
#    }
# }

out_from_input_file = function (new_input, suffix, input) {
  if (new_input == '<<default>>') {
    file.path(dirname(dirname(input)), suffix)
  } else {
    new_input
  }
}


species_from_metadata = function (meta_data, sep='$', space_char='_', collapse=F,
                                  collapse_str=', ') {
  # Process FASTA meta-data entry and extract all (unique) species names
  # into a character vector.
  # First replace any double $$ with single
  meta_data = gsub('[\\$]{2,}', '\\$', meta_data)
  if (sep == '$') sep = '\\$'
  strains = unlist(strsplit(meta_data, sep))
  species = sapply(strains, function (s) {
    x = unlist(strsplit(s, space_char))
    paste0(x[1], ' ', x[2])
  })
  if (collapse) {
    return(paste(sort(unique(species)), collapse=collapse_str))
  } else {
    return(sort(unique(species)))
  }
}



filter_rrn_from_best_matches = function (dt, id_col='dada_seqid',
                                         strain_col='strain') {
   # For each dada_seqid check if all entries are just different rRNA genes
   # do this by extracting strain names, then find number of uniques
   dt[, {
      .SD[]
   }, by=c(id_col)]
}

strain_name_from_match_strain = function (match_strain) {
   x = gsub('\\$$', '', gsub('_@rrn[0-9]', '', match_strain))
   gsub('_', ' ', strsplit(x, '\\$')[[1]][1])
}

species_from_strain_name = function (strain_name, sep='_') {
   # Parse out species name from strain name
   items = strsplit(strain_name, '_')

}

pctsim_range_old = function (p) {
  up = unique(p)
  if (length(up) == 1) return(as.character(up))
  else {
    return(paste0(min(p), ' - ', max(p)))
  }
}

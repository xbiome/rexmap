# Run BLAST with denoised sequences from dada2 against a reference FASTA file
# Record all matches.

# We will use this to check how many exact matches, one off etc we have against
# sequences in the reference file.


blast_best_seq_matches = function (dt, id_col='dada_seqid', match_col='match_strains',
                                   exclude_match=FALSE) {
  # For each strain sequence, find the best match (by bitscore) for
  # fwd and rev primers.
  # exclude_match=T will not show the match_col in the output and will just show
  # statistics for the best matched hit.
  dt2 = dt[, {
    if (.N > 0) {
      best_match = .SD[score == max(score)]
      if (exclude_match) {
        unique(best_match[, setdiff(names(best_match), match_col), with=F])
      } else {
        best_match
      }
    }
  }, by=c(id_col)]
  return(dt2)
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

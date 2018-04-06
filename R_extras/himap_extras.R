# Check one sequence vs a list of other sequences and return
# a subset of other sequences that are <= N hamming distances away
q = 'ACGTACGTACGT'
ts = c('ACGTACGTACGA', 'ACGTACAAACGT', 'GGTTACGTACGTN', 'ACGTACGGGTGTNN', 'ACGGGTGTACGT', 'TAGGGTGTACGT',
          'TGATGTAAACCC', 'ACGTATGTACGA', 'ACGTACGTGCGT')

find_possible_pcr_errors_query = function (query, targets, hamming_distance=1L) {
  dists = unname(sapply(targets, function (t) hamming(query, t)))
  return(which(dists==hamming_distance))
}

find_possible_pcr_errors = function (sequences, hamming_distance=1L) {
  # Sequences need to be already sorted by abundance, descending.

  for (i in 2:(length(sequences))) {
    targets = sequences[i:length(sequences)][nchar(sequences[i:length(sequences)]) == nchar(sequences[i-1])]
    if (length(targets) == 0) next
    hams = find_possible_pcr_errors_query(sequences[i-1], targets)
    if (length(hams) > 0) {
      cat('Sequence ', i-1, ' has 1-off variants at lower abundance: ')
      # Iterate over hits and return index from the original vector
      for (s in targets[hams]) {
        cat(which(sequences==s), ', ')
      }
      cat('\n')
    }
  }
}

load_copy_number_table = function (ref_cp) {
  # Load copy number information
  cp.dt = fread(ref_cp, select=c(1:3), colClasses=c('character', 'character', 'integer', 'character'))

  # Prepare copy number table columns
  cp.dt[, strain := gsub('_@rrn[0-9]+', '', strain_name)]

  # Generate strain index
  cp.dt[, strain_id := .GRP, by=strain]
  cp.dt[, strain_name := NULL]

  # How many unique 16S rrna variants of each strain?
  cp.dt[, rrn_uniq := length(unique(variant_id)), by=strain]

  # Sort by strain_id, then by
  cp.dt = cp.dt[order(strain_id, variant_id)]

  # Group sequences into OSU
  cp.dt[, spectrum := paste(copy_number, variant_id, sep=':', collapse=','), by=strain]

  # For strains which have only one unique 16S sequence variant, group them into the same
  # same OSU, regardless of the copy number.
  cp.dt[rrn_uniq==1, spectrum := paste0('X:', strsplit(spectrum[1], ':')[[1]][2]), by=strain]

  # OSU is a group of strains all with the same strain spectrum
  cp.dt[, osu_id := .GRP, by=spectrum]

  # Instead of X put mean copy number or 1
  cp.dt[rrn_uniq==1, spectrum := gsub('X', round(mean(copy_number)), spectrum), by=osu_id]
  # ref_cp.dt[rrn_uniq==1, spectrum := gsub('X', '1', spectrum), by=osu_id]

  # Count number of strains binned together
  cp.dt[, no_strains_in_osu := length(unique(strain)), by=osu_id]

  # Generate species names
  cp.dt[, species := gsub('^([^_]+)_([^_]+)_.*', '\\1_\\2', strain)]
  # Count unique number of species per OSU
  cp.dt[, no_species_in_osu := length(unique(species)), by=osu_id]

  cp.dt
}

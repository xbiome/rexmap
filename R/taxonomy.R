# Add taxonomy annotation to the final

# test

# all_ab.dt = readRDS('~/cloud/research/microbiome/diabimmune/all_ab.dt_top4')
# blast_best.dt = readRDS('~/cloud/research/microbiome/diabimmune/blast_best.dt')
#
# var.dt = fread('~/cloud/research/microbiome/genomes/data/vregions_db/V3-V4_337F-805R_hang21_wrefseq_table.txt',
#                select=c())


#' Read lineages*.csv file and convert it to compressed R object
#'
#' This will be used for taxonomy matching.
lineagecsv_to_robj = function (lincsv, out) {
   taxonomy.dt = data.table::fread(lincsv, select=1:8)
   taxonomy.dt = taxonomy.dt[genus != '']
   taxonomy.dt[, species := gsub('^Candidatus ', '', species)]
   taxonomy.dt = unique(taxonomy.dt)
   saveRDS(taxonomy.dt, file=out, compress=T)
}

load_taxonomy = function (tax_file=himap_option('taxonomy_file')) {
   readRDS(tax_file)
}

#' Parse genus names and strain counts from each OSU species string
#'
osuab_genuses = function (osuab) {
  osuab[, {
    species_list = strsplit(species[1], ',', fixed=T)[[1]]
    genus_counts = lapply(species_list, function (s) {
      genus = gsub('^([^_]+)_.*', '\\1', s)
      strain_count = gsub('.*_\\[([0-9]+)\\]$', '\\1', s)
      if (strain_count == s) strain_count = '1'
      return(c(genus, strain_count))
    })
    genus_counts_unique = list()
    for (i in 1:length(genus_counts)) {
      genus_i = genus_counts[[i]][1]
      counts_i = as.integer(genus_counts[[i]][2])
      if (genus_i %in% names(genus_counts_unique)) {
        genus_counts_unique[[genus_i]] = genus_counts_unique[[genus_i]] + counts_i
      } else {
        genus_counts_unique[[genus_i]] = counts_i
      }
    }
    # Now pool together same genuses and add up the counts
    list('genus'=names(genus_counts_unique), 'genus_strain_count'=unname(genus_counts_unique))
  }, by=osu_id]
}

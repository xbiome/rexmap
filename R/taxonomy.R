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

#' Load HiMAP taxonomy reference file
#'
#' This is a NCBI Taxonomy database, exported as a tab-delimited file.
#'
#' @export
load_taxonomy = function (tax_file=himap_option('taxonomy_file')) {
   readRDS(tax_file)
}

#' Parse genus names and strain counts from each OSU species string
#'
osuab_genuses = function (osuab) {
  dt = osuab[, {
    species_list = strsplit(species[1], ',', fixed=T)[[1]]
    genus_counts = lapply(species_list, function (s) {
      genus = gsub('^([^_]+)_.*', '\\1', s)
      if (genus == 'Candidatus') genus = gsub('^[^_]+_([^_]+)_.*', '\\1', s)
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
    list('genus'=names(genus_counts_unique),
         # 'order_from_genus'=as.character(NA),
         'strain_count'=as.integer(unname(genus_counts_unique)),
         'pctsim'=unique(pctsim))
  }, by=.(osu_id)]
  # dt[genus %like% 'les$', c('order_from_genus', 'genus') := list(genus, NA)]
  return(dt)
}

#' Generate taxonomy for each OSU
#'
#' @param osu_abundance_table Data table with OSU abundances. Output from abundance().sssss
#'
#' @export
taxonomy = function (osu_abundance_table) {
  taxonomy.dt = load_taxonomy()
  taxonomy.dt[, c('tax_id', 'species') := NULL]
  taxonomy.dt = unique(taxonomy.dt)
  tax_family.dt = unique(taxonomy.dt[, .(superkingdom, phylum, class, order, family)])
  tax_order.dt = unique(taxonomy.dt[, .(superkingdom, phylum, class, order, family=NA)])
  tax_class.dt = unique(taxonomy.dt[, .(superkingdom, phylum, class, order=NA, family=NA)])
  tax_phylum.dt = unique(taxonomy.dt[, .(superkingdom, phylum, class=NA, order=NA, family=NA)])

  osu_ab_g.dt = osuab_genuses(osu_abundance_table)
  osu_ab_g2.dt = merge(osu_ab_g.dt, taxonomy.dt, by='genus', all.x=T) # Genus matches
  col_order = c('osu_id', 'strain_count', 'pctsim', 'superkingdom', 'phylum',
                'class', 'order', 'family', 'genus')
  setcolorder(osu_ab_g2.dt, col_order)
  osu_ab_g2.dt[is.na(superkingdom) & genus %in% taxonomy.dt[, family],
               c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus') := {
                 dt = merge(.SD[, .(osu_id, genus)],
                            tax_family.dt,
                            by.x='genus', by.y='family')
                 list(dt[, superkingdom], dt[, phylum], dt[, class], dt[, order],
                      dt[, genus], NA)
               }]
  osu_ab_g2.dt[is.na(superkingdom) & genus %in% taxonomy.dt[, order],
               c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus') := {
                 dt = merge(.SD[, .(osu_id, genus)],
                            tax_order.dt,
                            by.x='genus', by.y='order')
                 list(dt[, superkingdom], dt[, phylum], dt[, class], dt[, genus],
                      NA, NA)
               }]
  osu_ab_g2.dt[is.na(superkingdom) & genus %in% taxonomy.dt[, class],
               c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus') := {
                 dt = merge(.SD[, .(osu_id, genus)],
                            tax_class.dt,
                            by.x='genus', by.y='class')
                 list(dt[, superkingdom], dt[, phylum], dt[, genus],
                      NA, NA, NA)
               }]
  osu_ab_g2.dt[is.na(superkingdom) & genus %in% taxonomy.dt[, phylum],
               c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus') := {
                 dt = merge(.SD[, .(osu_id, genus)],
                            tax_phylum.dt,
                            by.x='genus', by.y='phylum')
                 list(dt[, superkingdom], dt[, genus],
                      NA, NA, NA, NA)
               }]
  # Now count unique taxonomic ranks for each osu_id
  osu_ab_ranks.dt = osu_ab_g2.dt[, {
    phylum_num = .SD[!is.na(phylum), .(phylum_num = sum(strain_count)),
                     by=phylum][order(-phylum_num)][, paste(phylum,
                        phylum_num, sep='_[', collapse='],')]
    class_num = .SD[!is.na(class), .(class_num = sum(strain_count)),
                     by=class][order(-class_num)][, paste(class,
                        class_num, sep='_[', collapse='],')]
    order_num = .SD[!is.na(order), .(order_num = sum(strain_count)),
                     by=order][order(-order_num)][, paste(order,
                        order_num, sep='_[', collapse='],')]
    family_num = .SD[!is.na(family), .(family_num = sum(strain_count)),
                     by=family][order(-family_num)][, paste(family,
                        family_num, sep='_[', collapse='],')]
    genus_num = .SD[!is.na(genus) & !grepl('^[a-z]', genus),
                    .(genus_num = sum(strain_count)),
                     by=genus][order(-genus_num)][, paste(genus,
                        genus_num, sep='_[', collapse='],')]
    list('phylum'=paste0(phylum_num, ']'), 'class'=paste0(class_num, ']'),
         'order'=paste0(order_num, ']'), 'family'=paste0(family_num, ']'),
         'genus'=paste0(genus_num, ']'))
  }, by=.(osu_id, pctsim)]
  osu_ab_ranks.dt[!(grepl('^[A-Z]', phylum)), phylum := NA]
  osu_ab_ranks.dt[!(grepl('^[A-Z]', class)), class := NA]
  osu_ab_ranks.dt[!(grepl('^[A-Z]', order)), order := NA]
  osu_ab_ranks.dt[!(grepl('^[A-Z]', family)), family := NA]
  osu_ab_ranks.dt[!(grepl('^[A-Z]', genus)), genus := NA]
  return(osu_ab_ranks.dt[])
}


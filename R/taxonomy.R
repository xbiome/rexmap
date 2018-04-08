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
   dt = data.table::fread(lincsv, select=1:8)
   dt = dt[genus != '']
   dt[, species := gsub('^Candidatus ', '', species)]
   dt = unique(dt)
   save(dt, file='inst/database/lineages-2017-03-17_bacteria_archaea', compress=T)
}

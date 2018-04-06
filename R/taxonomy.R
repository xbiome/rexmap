# Add taxonomy annotation to the final

# test

all_ab.dt = readRDS('~/cloud/research/microbiome/diabimmune/all_ab.dt_top4')
blast_best.dt = readRDS('~/cloud/research/microbiome/diabimmune/blast_best.dt')

var.dt = fread('~/cloud/research/microbiome/genomes/data/vregions_db/V3-V4_337F-805R_hang21_wrefseq_table.txt',
               select=c())

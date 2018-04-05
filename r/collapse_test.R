# To test our collapse function first take some data from
# diabimmune samples where collapseNoMismatch finds stuff to be
# collapsed.

ab_tab = readRDS('~/cloud/research/microbiome/diabimmune/ab_tab')

# Select some sample_id (rows), then remove all sequences that have
# zero counts in those sample_ids
rows = 1:2
ab_tab2 = ab_tab[rows,]
cols_zero = colSums(ab_tab2)
cols_zero_filter = cols_zero==0
ab_tab3 = ab_tab2[, !cols_zero_filter]

# Do we have anything to collapse?
ab_tab3_coll = collapseNoMismatch2(ab_tab3, verbose=T)
# No.

# Try QIIME1 analysis of zheng 2015
# Load FASTA and the table and generate abundance matrix
q1_fa = fasta_reader('~/cloud/research/microbiome/zheng-2015/qiime1_analysis/V3V4Rep1_denovo_otus/rep_set/V3V4Rep1_nochim_rep_set.fasta')
q1_tb = fread('~/cloud/research/microbiome/zheng-2015/qiime1_analysis/V3V4Rep1_denovo_otus/otu_counts.txt')
names(q1_tb) = c('id', 'count')
q1_fa$meta = unname(sapply(q1_fa$meta, function (x) strsplit(x, ' ', fixed=T)[[1]][1]))
q1_fa.dt = data.table(id=q1_fa$meta, seq=q1_fa$seqs)
q1.dt = merge(q1_fa.dt, q1_tb, by='id')
q1.dt[, sample_id := 'V3V4Rep1']
q1.m = matrix(q1.dt[, count], nrow=1)
rownames(q1.m) = 'V3V4Rep1'
colnames(q1.m) = q1.dt[, seq]

# Sample 2
q1_fa = fasta_reader('~/cloud/research/microbiome/zheng-2015/qiime1_analysis/V3V4Rep2_denovo_otus/rep_set/V3V4Rep2_nochim_rep_set.fasta')
q1_tb = fread('~/cloud/research/microbiome/zheng-2015/qiime1_analysis/V3V4Rep2_denovo_otus/otu_counts.txt')
names(q1_tb) = c('id', 'count')
q1_fa$meta = unname(sapply(q1_fa$meta, function (x) strsplit(x, ' ', fixed=T)[[1]][1]))
q1_fa.dt = data.table(id=q1_fa$meta, seq=q1_fa$seqs)
q1.dt2 = merge(q1_fa.dt, q1_tb, by='id')
q1.dt2[, sample_id := 'V3V4Rep2']
q1.m2 = matrix(q1.dt2[, count], nrow=1)
rownames(q1.m2) = 'V3V4Rep2'
colnames(q1.m2) = q1.dt2[, seq]

# Combine into a single matrix
q1_all.dt = merge(q1.dt, q1.dt2, all=T, by=names(q1.dt))
q1_all.dt[, id := NULL]
q1_all_d.dt = dcast(q1_all.dt, sample_id ~ seq, 
                    value.var='count', fill=0)
q1_all.m = as.matrix(q1_all_d.dt[, 2:ncol(q1_all_d.dt)])
rownames(q1_all.m) = q1_all_d.dt[, sample_id]

q1.m_coll = collapseNoMismatch2(q1_all.m, verbose=T)

# ok check which one

# Load test data from diabimmune
ab_tab_nochim = readRDS('~/cloud/research/microbiome/diabimmune/ab_tab_nochim')
ab_tab_nochim_colsums = colSums(ab_tab_nochim)
ab_tab_nochim_maxab = apply(ab_tab_nochim, 2, max)
# ab_tab_nochim_red = ab_tab_nochim[, which(ab_tab_nochim_colsums > 10)]
ab_tab_nochim_red = ab_tab_nochim[, which(ab_tab_nochim_maxab > 20)]



ab_tab_nochim_coll = readRDS('~/cloud/research/microbiome/diabimmune/ab_tab_nochim_coll')


write_seqs_to_fasta(dimnames(ab_tab_nochim_coll)[[2]], '~/cloud/research/')


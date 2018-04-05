source('~/cloud/research/microbiome/himap/r-lib/himap_pipeline_v4.R')

tax_file = file.path(himap_path, 'data/lineages-2017-03-17_bacteria_archaea.csv')
gen_a_file = file.path(himap_path, 'data/archaea_assembly_summary_filter_2018-02-23.txt')
gen_b_file = file.path(himap_path, 'data/bacteria_assembly_summary_filter_2018-02-23.txt')
vregion_tab = file.path(himap_path, 'data/V4-V5-2_515F-926R_hang22_wrefseq_table.txt')
tax_out = file.path(himap_path, 'data/taxonomy_2018-04-05.txt')

tax.dt = fread(tax_file, select=c(1:8))
names(tax.dt) = c('taxid', 'kingdom', 'phylum', 'class', 'order',
                  'family', 'genus', 'species')
tax.dt[, species := gsub(' ', '_', species, fixed=T)]
tax.dt[, species := gsub('sp\\..*', 'sp.', species)]
tax.dt[, species := gsub('^([^_]+_[^_]+)[_]?.*$', '\\1', species)]
# tax.dt[, ncbi_tax_id := NULL]
# tax.dt = unique(tax.dt)

# tax_genus.dt = unique(tax.dt[, .(kingdom, phylum, class, order, family, genus)])

# Now get species and tax id for full genome sequences
gen.dt = rbindlist(lapply(c(gen_a_file, gen_b_file), fread))

# Add taxonomy information
gen_tax.dt = merge(gen.dt, tax.dt, by='taxid', all.x=T)
# Basically we only need assembly_accesion, taxid and taxonomic ranks
gen_tax.dt = gen_tax.dt[, .(ncbi_id=get('# assembly_accession'), strain_name, taxid, kingdom, phylum,
               class, order, family, genus, species)]
gen_tax.dt[, strain_name := gsub(' ', '_', strain_name)]

# Now load NCBI RefSeq 16S search results
rs.dt = fread(vregion_tab, select=1:2)
rs.dt = rs.dt[!(assembly_id %like% '^GCF')]
rs.dt[, species := sub('^([^_]+_[^_]+)[_]?.*', '\\1', strain_name), by=assembly_id]
rs_tax.dt = merge(rs.dt, tax.dt, by='species', all.x=T)
rs_tax2.dt = rs_tax.dt[, .(strain_name=strain_name[1],
                           taxid=if (length(unique(taxid)) == 1) taxid[1] else as.integer(NA),
                           kingdom=paste(unique(kingdom[kingdom!='']), collapse=','),
                           phylum=paste(unique(phylum[phylum!='']), collapse=','),
                           class=paste(unique(class[class!='']), collapse=','),
                           order=paste(unique(order[order!='']), collapse=','),
                           family=paste(unique(family[family!='']), collapse=','),
                           genus=paste(unique(genus[genus!='']), collapse=','),
                           species=species[1]), by=assembly_id]
names(rs_tax2.dt)[1] = 'ncbi_id'
# This should give unique assignments for every species (like it is at a time I am testing),
# but in case it doesn't multiple assignments will be delimited by ','

# Combine these into a single table
all_tax.dt = rbindlist(list(gen_tax.dt, rs_tax2.dt))

# Write output
write.table(all_tax.dt, tax_out, sep='\t', quote=F, row.names=F)

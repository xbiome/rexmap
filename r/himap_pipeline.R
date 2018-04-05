# For this pipeline to work we will need to provide filenames for forward and
# reverse reads for each sample.

# Load required packages
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ShortRead))

# Load HiMAP functions
source('/data1/igor/himap/r/mergepairs.R')
sourceCpp('/data1/igor/himap/src/hamming.cpp')
sourceCpp('/data1/igor/himap/src/fastq_retrieve.cpp')
source('/data1/igor/himap/r-lib/blast_vs_fasta.R')
source('/data1/igor/himap/r-lib/spp_shorten.R')
source('/data1/igor/himap/r/read_fastx.R')

#
# s1 r11 ++ r12 -
# s2 r21 =
# s3 r31 x
# 
# ++-=x   5 total, s1->3, s2->1, s3->1
# 


# Processing parameters

# Minimum BLAST alignment length
blast_aln_len = 200

out_from_input_file = function (new_input, suffix, input) {
  if (new_input == '<<default>>') {
    file.path(ddirname(dirname(input)), suffix)
  } else {
    new_input
  }
}

ie = function(test, yes, no) {
  if (test) {
    yes
  } else {
    no
  }
}


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

# REMOVE ME: prototyping code (do this for each sample)
# ns = length(dada_res)
# sid_to_pid_to_seq = vector('list', ns)
# for (s_id in 1:length(dada_res)) { # For each sample s_id 
#   x = partid_to_fastqid(dada_res[[s_id]]$map-1, derep[[s_id]]$map-1)
#   fq_fil_data = sfastq_reader(fq_fil[s_id])
#   fq_mer_data = sfastq_reader(fq_mer[s_id])
#   pid_to_seq = vector("list", length(x))
#   for (i in 1:length(x)) { # for each partition i in sample s_id
#     metas = fq_fil_data[['meta']][x[i][[1]]+1]
#     mer_seqs = fq_mer_data[['seqs']][which(fq_mer_data[['meta']] %in% metas)]
#     pid_to_seq[i] = gsub('[-]{1,}$', '', consensus_sequence(mer_seqs))
#   }
#   sid_to_pid_to_seq[s_id] = pid_to_seq
# }


# for each partition id, extract all full-length merged sequences
# then summarize them as a consensus sequence
p_id = 95

# /REMOVE ME

himap = function (fq_fwd, fq_rev, # FASTQ input 
                  # Output folder and file names
                  out, # Output folder
                  out_prefix = 'himap',
                  fq_mer='<<default>>', fq_fil='<<default>>', # Merged and filtered FASTQ output
                  fq_fil_minlen='<<default>>',
                  fq_prefix_split='_',
                  ab_tab_out='<<default>>', # abundance table
                  # Merging parameters
                  min_pct_sim=0.6, min_aln_len=50,
                  # Filtering parameters
                  minimum_len=20, maxN=0, maxLen=Inf, maxEE=2, truncQ=2, verbose=T, multithread=T,
                  trimLeft=0, truncLen=0,
                  # DADA2 partitioning parameters
                  pvalThreshold=1e-4, # p-value threshold (before FDR correction)
                  # BLAST parameters
                  ref_db='/data1/igor/zr_16s2/16s_uniq_refseq_fullgen_bact_arch_2017-08-24',
                  ref_cp='/data1/igor/himap/db/V4_515F-805R_hang1_unique_sequences_table.txt',
                  ncpu=39, fa_denoised='<<default>>', pipeline_stats='<<default>>',
                  blast_output='<<default>>', blast_out_best='<<default>>',
                  blast_out_best_sp='<<default>>', 
                  max_target_seqs=2000,
                  # Pipeline run parameters, set to false if we already ran this
                  run_merge=TRUE, run_filter=TRUE, run_chim=TRUE,
                  # Do we want the function to directly return output tables (don't need to read files)
                  return_tables=TRUE
) {
  
  cat('--- HiMAP pipeline ---\n')
  cat(paste0('Output folder: ', out, '\n'))
  
  # Automatically generate output filenames if they are not provided
  # fq_prefix = strsplit(sub("(.*)\\..*$", "\\1", basename(fq_fwd)), fq_prefix_split, fixed=T)[[1]][1]
  fq_prefix = unname(
    sapply(
      sub("(.*)\\..*$", "\\1", basename(fq_fwd)), 
      function (x) strsplit(x, fq_prefix_split, fixed=T)[[1]][1]
    )
  )
  fq_mer = ie(fq_mer=='<<default>>', file.path(out, 'merged', paste0(fq_prefix, '_merged.fastq')), fq_mer)
  fq_fil = ie(fq_fil=='<<default>>', file.path(out, 'filtered', paste0(fq_prefix, '_filtered.fastq')), fq_fil)
  fq_fil_minlen = ie(fq_fil_minlen=='<<default>>', file.path(out, 'filtered', paste0(fq_prefix, '_minlenfiltered.fastq')), 
                     fq_fil_minlen)
  fa_denoised = ie(fa_denoised=='<<default>>', file.path(out, paste0(out_prefix, '_denoised_seqs.fasta')), fa_denoised)
  pipeline_stats = ie(pipeline_stats=='<<default>>', file.path(out, paste0(out_prefix, '_pipeline_reads_each_step.txt')),
                      pipeline_stats)
  blast_output = ie(blast_output=='<<default>>', file.path(out, paste0(out_prefix, '_blast_results.txt')), 
                    blast_output)
  blast_out_best = ie(blast_out_best=='<<default>>', file.path(out, paste0(out_prefix, '_blast_best.txt')), 
                      blast_out_best)
  blast_out_best_sp = ie(blast_out_best_sp=='<<default>>', file.path(out, paste0(out_prefix, '_blast_best_sp.txt')), 
                         blast_out_best_sp)
  ab_tab_out = ie(ab_tab_out=='<<default>>', file.path(out, paste0(out_prefix, '_ab_tab_nochim_m.txt')), ab_tab_out)
  # maxQscore = ie(maxQscore=='<<default>>', )
  cat('Filtered:', fq_fil, '\n')
  
  # Merge reads
  if (run_merge) {
    cat('Merging reads...')
    if (!dir.exists(dirname(fq_mer[1]))) dir.create(dirname(fq_mer[1]), recursive=T)
    mergestats = merge_pairs(fq_fwd, fq_rev, fq_mer, min_pct_sim=0.6, min_aln_len=50)
    mergestats = as.data.table(t(mergestats))
    cat('OK.\n')
  }
  
  # Filter reads
  if (run_filter) {
    if (!dir.exists(dirname(fq_fil[1]))) dir.create(dirname(fq_fil[1]), recursive=T)
    cat('Filtering reads...')
    # Do we need minimum length filtering first?
    # if (minimum_len < Inf) {
    #   filterAndTrim(fq_mer, fq_fil_minlen, minLen=minimum_len, multithread=multithread)
    #   filterAndTrim(fq_fil_minlen, fq_fil, maxLen=maxLen, maxN=maxN, maxEE=maxEE, 
    #                 truncQ=truncQ, trimLeft=trimLeft, 
    #                 truncLen=truncLen, verbose=verbose, multithread=multithread)
    #   
    # } else {
    #   
      filterAndTrim(fq_mer, fq_fil, maxLen=maxLen, maxN=maxN, maxEE=maxEE, 
                    truncQ=truncQ, trimLeft=trimLeft,  minLen=minimum_len,
                    truncLen=truncLen, verbose=verbose, multithread=multithread)
    # }
    cat('OK.\n')
  }
  
  # Dereplicate merged reads
  cat('Dereplicating reads...')
  derep = derepFastq(fq_fil)
  if (class(derep)[1] == 'derep') {
    # Only one sample is loaded, so put the derep object into a list
    derep = list(derep)
  }
  cat('OK.\n')
  
  # Learn errors
  cat('DADA2: Learn errors...')
  dada_errors = learnErrors(fq_fil, multithread=T)
  cat('OK.\n')
  
  # Plot errors
  # plotErrors(dada_errors, nominalQ=TRUE)
  
  # Find number of unique sequences
  derep_no_uniq = sapply(derep, function (d) length(d$uniques))
  # derep_no_uniq = length(derep$uniques)
  
  # Run dada2
  # <= 10^5 unique sequences / sample, max Q score 1e-4, so try OMEGA_A = 1e-9
  cat('Running DADA2...')
  dada_res = dada(derep, err=dada_errors, multithread=multithread, OMEGA_A=pvalThreshold/max(derep_no_uniq))
  cat('OK.\n')
  
  # If we ran just 1 sample, then dada_res is a dada-class object, not a list. Make it a list.
  if (class(dada_res)[1] == 'dada') {
    dada_res = list(dada_res)
  }
  
  # Retrieve original full-length sequences from merged data for further analysis
  # For each partition cluster, we use the dada2 sequence and the consensus of
  # left and right overhangs from actual sequences, concatenated together,
  # if we did any trimming (we need to otherwise dada2 does not work).
  cat('Retrieving full-length sequences...')
  for (s_id in 1:length(dada_res)) { # For each sample s_id 
    x = partid_to_fastqid(dada_res[[s_id]]$map-1, derep[[s_id]]$map-1)
    fq_fil_data = sfastq_reader(fq_fil[s_id])
    fq_mer_data = sfastq_reader(fq_mer[s_id])
    pid_to_seq = vector("list", length(x))
    for (i in 1:length(x)) { # for each partition i in sample s_id
      metas = fq_fil_data[['meta']][x[i][[1]]+1]
      mer_seqs = fq_mer_data[['seqs']][which(fq_mer_data[['meta']] %in% metas)]
      pid_to_seq[i] = gsub('[-]{1,}$', '', consensus_sequence(mer_seqs))
    }
    # sid_to_pid_to_seq[s_id] = pid_to_seq
    # Concatenate left consensus, dada2 middle partition sequence and right consensus
    dada_res[[s_id]]$sequence = paste(str_sub(pid_to_seq, end=trimLeft), 
                                      dada_res[[s_id]]$sequence, 
                                      str_sub(pid_to_seq, start=truncLen+1), 
                                      sep='')
    names(dada_res[[s_id]]$denoised) = dada_res[[s_id]]$sequence
    cat('.')
  }
  cat('OK.\n')
  
  # Generate abundance table
  cat('Generating abundance table...')
  ab_tab = makeSequenceTable(dada_res)
  cat('OK.\n')
  
  # Remove chimeras
  if (run_chim) {
    cat('Removing chimeras...')
    ab_tab_nochim = removeBimeraDenovo(ab_tab, method='consensus', allowOneOff=T,
                                       multithread=multithread, verbose=verbose)
    cat('OK.\n')
  } else {
    ab_tab_nochim = ab_tab
  }
  
  # Combine sequences which are just offset shifted
  cat('Adding up counts for offset shifts...')
  ab_tab_nochim = collapseNoMismatch(ab_tab_nochim)
  cat('OK.\n')
  
  # Track reads through pipeline
  nreads_orig = sapply(fq_fwd, function (fq) length(readFastq(fq)))
  nreads_merg = sapply(fq_mer, function (fq) length(readFastq(fq)))
  # nreads_filt = sum(getUniques(dada_res))
  nreads_filt = sapply(dada_res, function (d) sum(getUniques(d)))
  nreads_tabl = rowSums(ab_tab)
  nreads_noch = rowSums(ab_tab_nochim)
  
  # Save this statistics in a separate file
  # put nreads_* into pipeline_stats.dt then save it to pipeline_stats file
  cat('Saving statistics...')
  nreads.dt = data.table(sample_id=fq_prefix, original=nreads_orig, merged=nreads_merg,
                         filtered=nreads_filt, tabled=nreads_tabl,
                         nonchimeric=nreads_noch)
  
  write.table(nreads.dt, pipeline_stats, sep='\t', row.names=F, quote=F)
  cat('OK.\n')
  
  # Write the denoised sequences to FASTA file for Blasting
  write_seqs_to_fasta(dimnames(ab_tab_nochim)[[2]], fa_denoised)
  
  cat('Running BLAST...')
  blast_out = blast_seqs_to_reference(fa_denoised, ref_db, ncpu=ncpu, output=blast_output,
                                      max_target_seqs=max_target_seqs)
  cat('OK.\n')
  
  # Now load this output into data table, then find best matches for every dada2
  # denoised sequence.
  cat('Calculating alignment scores...')
  blast_out.dt = fread(blast_output)
  # qseqid sseqid length mismatch gapopen qstart qend sstart send qseq sseq
  names(blast_out.dt) = strsplit(blast_out_fmt, ' ')[[1]]
  #aln = blast_out.dt[, mcmapply(compare_alignment, qseq, sseq, mc.cores=ncpu)]
  
  # c('score', 'match', 'mismatch', 'gapopen', 'gapextend') := 
  # Calculate BLAST scores.
  blast_out.dt[, c('score', 'match', 'mismatch', 'gapopen', 'gapextend') := transpose(unname(mcmapply(function (q,s) {
    aln = compare_alignment(q,s)
    as.integer(c(sum(aln_params*aln), aln[1], aln[2], aln[3], aln[4]))
  }, qseq, sseq, SIMPLIFY=F, mc.cores=ncpu)))]
  
  # For local alignments that are one-off at each end, reduce the score using
  # a mismatch penalty. (We deliberately use indel score that is very similar.)
  blast_out.dt[qstart == 2, score := score-aln_params[2]]
  blast_out.dt[qstart == 2, mismatch := mismatch+1L]
  blast_out.dt[qend == qlen-1, score := score-aln_params[2]]
  blast_out.dt[qend == qlen-1, mismatch := mismatch+1L]
  
  # Filter out alignments with large gaps at beginning or the end
  blast_out.dt = blast_out.dt[qstart <= 2 & qend >= qlen-1]
  
  cat('OK.\n')
  
  # Filter out alignments with short alignment length
  cat('Processing BLAST results...')
  # blast_out.dt = blast_out.dt[aln_len >= blast_aln_len]
  
  # For each dada2 sequence, keep only the best matches by bitscore
  blast_best.dt = blast_best_seq_matches(blast_out.dt, strsplit(blast_out_fmt, ' ')[[1]][1], 
                                         strsplit(blast_out_fmt, ' ')[[1]][2])
  blast_best.dt[, variant_id := gsub('^([0-9]+)-.*', '\\1', sseqid)]

  cat('Checkpoint1.')
  # Load copy number information
  cp.dt = fread(ref_cp, select=c(1:3), colClasses=c('character', 'character', 'integer', 'character'))
  
  cat('Checkpoint2.')
  # Merge with blast hits
  blast_best.dt = merge(blast_best.dt, cp.dt[, .(variant_id, strain_name, copy_number)], 
                        by='variant_id', all.x=T)
  
  cat('Checkpoint3.')
  # Remove sequence alignment columns since we already calculated alignment score
  blast_best.dt[, c('qseq', 'sseq') := NULL]
  blast_best.dt[, c('sseqid') := NULL]
  blast_best.dt[, strain := gsub('^(.*)_@rrn[0-9]+$', '\\1', strain_name)]
  
  # For each strain_name (unique 16S variant) record its best hits. Sometimes, the same
  # strain_name will hit multiple qseqid, some with a lower score. 
  # best_hit column record whether 
  cat('Checkpoint4.')
  # blast_best2.dt = merge(blast_best.dt, 
  #                        blast_best.dt[, .(qseqid=.SD[score==max(score), qseqid],
  #                                    best_hit=TRUE), by=strain_name],
  #                        all.x=T,
  #                        by=c('strain_name', 'qseqid'))
  # blast_best2.dt[is.na(best_hit), best_hit := FALSE]
  
  
  # Now from the best matches filter out multiple entries that differ only
  # in @rrnX. We just select the first entry since it doesn't matter -- they have
  # equally good BLAST alignment.
  
  # blast_best.dt = blast_best.dt[, {
  #   if (length(unique(gsub('_@rrn[0-9]+', '', sseqid))) == 1 & .N > 1) {
  #     .SD[1]
  #   } else {
  #     .SD
  #   }
  # }, by=c(strsplit(blast_out_fmt, ' ')[[1]][1])]
  
  # Now we still have multiple hits for each sequence. Keep that, but also summarize
  # the best we can. If there are multiple species 
  # First generate species names.
  # blast_best_sp.dt = blast_best.dt[, {
  #   .(species_str = spp_shorten(species_from_metadata(sseqid, collapse=T, collapse_str=',')),
  #     score = unique(score),
  #     match = unique(match)[1],
  #     mismatch = unique(mismatch)[1],
  #     gapopen = unique(gapopen)[1],
  #     gapextend = unique(gapextend)[1])
  # }, by=c(strsplit(blast_out_fmt, ' ')[[1]][1])]
  cat('OK.\n')
  
  # Now convert dada2 results into a data table, each column is a different condition
  # each row 
  cat('Converting and saving abundance tables...')
  if (nrow(ab_tab_nochim) == 1) dimnames(ab_tab_nochim)[[1]] = c('1')
  ab_tab_nochim.dt = as.data.table(unname(ab_tab_nochim))
  # ab_tab_nochim.dt[, sample_id := 1]
  ab_tab_nochim.dt[, sample_id := sapply(strsplit(dimnames(ab_tab_nochim)[[1]], ".", fixed=T), `[`, 1)]
  
  # Now melt the data table
  ab_tab_nochim_m.dt = melt(ab_tab_nochim.dt, id.vars = 'sample_id', 
                            variable.name = 'dada2_seqid', value.name = 'raw_count')
  
  # Cleanup dada2_seqid
  ab_tab_nochim_m.dt[, dada2_seqid := as.numeric(gsub('V', '', dada2_seqid))]
  
  # Now add raw counts for each sample 
  # Find the maximum number of reads across all samples, then set adjusted count such that
  # when we add all fake counts (zeros), they do not add up to more than 1% of the total.
  ab_tab_nochim_m.dt[, adj_count := as.numeric(raw_count)]
  max_sum_raw_counts = ab_tab_nochim_m.dt[, .(N = sum(raw_count)), by=sample_id][, max(N)]
  ab_tab_nochim_m.dt[raw_count==0, adj_count := 0.01*max_sum_raw_counts/.N, by=sample_id]
  
  # Normalize the adjusted counts for each sample
  ab_tab_nochim_m.dt[, adj_ab := adj_count/sum(adj_count), by=sample_id]
  
  # Now add raw_count to blast best results, then select strains for each sample id
  # separately.
  ab_blast.dt = merge(ab_tab_nochim_m.dt[, .(sample_id, qseqid=dada2_seqid, raw_count)],
                      blast_best.dt,
                      by='qseqid', all.x=T)
  
  # Ignore zero raw_count for each individual sample_id
  ab_blast.dt = ab_blast.dt[, .SD[raw_count > 0], by=sample_id]
  ab_blast.dt[, rrn := gsub('.*_@rrn([0-9]+)$', '\\1', strain_name)]
  
  # Save final tables
  write.table(ab_tab_nochim_m.dt, ab_tab_out,
              sep='\t', quote=F, row.names=F)
  # write.table(blast_best_sp.dt, blast_out_best_sp,
  #             sep='\t', quote=F, row.names=F)
  write.table(blast_best.dt, blast_out_best,
              sep='\t', quote=F, row.names=F)
  cat('OK.\n')
  
  # Return tables as output if this is set to TRUE
  if (return_tables) {
    return(list('abundance_table'=ab_tab_nochim_m.dt,
                'blast_best_hits'=blast_best.dt,
                'blast_best_hits_sp'=ab_blast.dt,
                'pipeline_stats'=nreads.dt)
    )
  }
}




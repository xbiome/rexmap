rm(list=ls(all=T))
library(tcR)
# suppressPackageStartupMessages(library(Rcpp))
# suppressPackageStartupMessages(library(Biostrings))
# suppressPackageStartupMessages(library(parallel))
sourceCpp('/data1/igor/himap/src/nw_fitting_align.cpp')
source('/data1/igor/himap/r/read_fastx.R')

# Test
query_fasta_filename = '/data1/igor/mockrobiota-13/himap_denoised_seqs.fasta'
subject_fasta_filename = '/data1/igor/himap/16s_uniq_refseq_fullgen_bact_arch_2017-08-24.fasta'
q1 = 'ACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGCGCGCGTAGGCGGTTTTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGAAAACTTGAGTGCAGAAGAGGAAAGTGGAATTCCATGTGTAGCGGTGAAATGCGCAGAGATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACGCTGATGTGCGAAAGCGTGGGGATCAAAC'

# Load FASTA files
# query_list = list('meta'=c('query1'), 'seq'=c(q1))
query_list = fasta_reader(query_fasta_filename)
subjects_list = fasta_reader(subject_fasta_filename)
subjects_list2 = list('meta'=subjects_list[['meta']][1],
                      'seqs'=subjects_list[['seqs']][1])
timing = TRUE
# Auxiliary message function
m = function (...) {
  message(..., appendLF=FALSE)
}

# Default parameters
MATCH = 1
MISMATCH = -1
GAP_OPEN = -2
GAP_EXTEND = -1
INDEL=-1
timing=TRUE
sample_size = 93
subjects_list = fasta_reader(subject_fasta_filename)
tests = grepl('^[ACGT]+$', subjects_list[['seqs']])
subjects_list2 = list('meta'=subjects_list[['meta']][1:sample_size],
                      'seqs'=subjects_list[['seqs']][1:sample_size])

subjects = subjects_list[['seqs']][tests][1:1000]
metas = subjects_list[['meta']][tests][1:1000]


if (timing) start_time = Sys.time()
i = 0
# We can try BLAST first to pre-scan potentially good targets?
# seq 94 has Y. change every non-ACGT letter into N, during database assembly in Python.
query_hits = mapply(function (meta, seq) {
  # aln = fitting_alignment_fast(query_list[['seq']], seq, match=MATCH, mismatch=MISMATCH, indel=INDEL)
  # data.table('query' = query_list$meta, 'subject' = meta,
  #            # 'score' = aln$score,
  #            'qseq' = aln[1], 'sseq' = aln[2])
  aln = C_nwalign(query_list[['seq']], seq, match=MATCH, mismatch=-MISMATCH, 
                  indel=INDEL)
  data.table('query' = query_list$meta, 'subject' = meta,
             'score' = 1,
             'qseq' = aln[1], 'sseq' = aln[2])

}, metas, subjects, SIMPLIFY=FALSE)
# query_hits.dt = rbindlist(query_hits)

if (timing) {
  end_time = Sys.time()
  m('Execution time: ', round(as.numeric(end_time-start_time)/60), ' m ', 
    round(as.numeric(end_time-start_time)%%60, 2), ' s.\n')
}



























fitting_alignment = function (query_list, subjects_list, match, mismatch, gap_open, gap_extend,
                              threads=detectCores()-1, timing=FALSE) {
  # subjects_list[['seqs']] are sequences.
  # query_list[['seq']] = sequence, query_list[['meta']] is the ID of the query sequence.
  subjects_list2 = list('meta'=subjects_list[['meta']][1:100],
                       'seqs'=subjects_list[['seqs']][1:100])
  if (timing) start_time = Sys.time()
  i = 0
  # We can try BLAST first to pre-scan potentially good targets?
  query_hits = mapply(function (meta, seq) {
    aln = fit_align_fast(query_list[['seq']], seq, match=MATCH, mismatch=MISMATCH, indel=INDEL)
    # aln = fit_align(query_list[['seq']], seq, match=MATCH, mismatch=MISMATCH, gap_open=GAP_OPEN, gap_extend=GAP_EXTEND)
    # list('query' = query_list[['meta']], 'subject' = meta,
    #      'no_diff' = nchar(aln[['query']]) - aln[['right']] - aln[['left']] + 1 -
    #        aln[['no_matches']],
    #      'q_gapopen' = aln[['q_no_gapopen']], 'q_gapextend' = aln[['q_no_gapext']],
    #      's_gapopen' = aln[['s_no_gapopen']], 's_gapextend' = aln[['s_no_gapext']],
    #      'query_seq' = substring(aln[['query']], aln[['left']], nchar(aln[['query']]) - 
    #                              aln[['right']]),
    #      'subject_seq' = substring(aln[['subject']], aln[['left']], nchar(aln[['subject']]) - 
    #                                aln[['right']])
    # )
    data.table('query' = query_list$meta, 'subject' = meta,
               #'no_diff' = nchar(aln$query) - aln$right -aln$left + 1 - aln$no_matches
               'score' = aln$score
               )
  }, subjects_list2[['meta']], subjects_list2[['seqs']], SIMPLIFY=FALSE)
  # query_hits.dt = rbindlist(query_hits)
  
  if (timing) {
    end_time = Sys.time()
    m('Execution time: ', round(as.numeric(end_time-start_time)/60), ' m ', 
      round(as.numeric(end_time-start_time)%%60, 1), ' s.\n')
  }
  
  return list('query_id'=query_list[['meta']], 'query_')
}



# tempdir() to get temporary dir path

get_makeblastdb_path = function () {
  # Find makeblastdb location
  out = suppressWarnings(system2('which', 'makeblastdb', stdout=T))
  # If makeblastdb is not found return an error
  if (is.null(attr(out, 'status'))) return(out)
  if (attr(out, 'status') == 1) stop('Error: makeblastdb command not found.')
}

makeblastdb = function (fasta_in, db_out, verbose=T) {
  # Generate a BLAST database from fasta_in
  # db_out contains a full location and prefix for BLASTDB files
  # This function does not return anything, but can optionally
  # (if verbose == T) write the status to console.
  makeblastdb_path = get_makeblastdb_path()
  if (verbose) output = ''
  else output = FALSE
  system2(makeblastdb_path, c('-dbtype', 'nucl', '-in', fasta_in,
                              '-out', db_out), stdout=output)
}

cleanup_blastdb = function (db_out) {
  # Delete blastdb files
  db_files = paste(db_out, c('.nin', '.nhr', '.nsq'), sep='')
  for (f in db_files) {
    if (file.exists(f)) file.remove(f)
  }
}

blast_coll_fmt = 'qseqid sseqid qlen slen length qstart qend sstart send pident'

blast = function (
  seqs_fa, ref_db,
  blast_path=blast_path_def,
  match=aln_params[1], mismatch=aln_params[2], gapopen=-aln_params[3],
  gapextend=-aln_params[4], word_size=13, ncpu=4, max_target_seqs=5000,
  perc_identity=75, outfmt=paste0('6 ', blast_out_fmt), dust='20 64 1',
  output='blast_output.txt',
  output_err=F) {
   x = system2(blast_path, args = c(
      '-dust', shQuote(dust), '-word_size', word_size,
      '-reward', match, '-penalty', mismatch, '-gapopen', gapopen,
      '-gapextend', gapextend, '-outfmt', shQuote(outfmt),
      '-query', seqs_fa, '-db', ref_db, '-num_threads', ncpu,
      '-max_target_seqs', max_target_seqs, '-perc_identity', perc_identity
      # '-qcov_hsp_perc', query_coverage_pct
   ), stdout = output, stderr = output_err)
   return(x)
}

# Normally formatted timestamp (to be used in file names)
ts = function () return(sub(' ', '_', gsub(':', '-', Sys.time())))

# Find minimal sequence length among the first N sequences from the fasta file
min_seq_len = function (fasta_files, n=50) {
   return(min(sapply(fasta_files, function (f) {
     return(nchar(fasta_reader(f)$seqs[1:n]))
   })))
}

collapse = function (ab_in, verbose=T) {
  # Collapse sequences that are exact matches, up to shifts and/or length
  # Provide ab_tab_nochim as an input argument

  # Generate temp files, this automatically generates file names
  if (verbose) cat('collapse:')
  if (verbose) cat('- generating temporary files...')
  ab = copy(ab_in)
  out_files = ab_to_files(ab)
  fasta = out_files[1]
  db = out_files[2]
  if (verbose) cat('OK.\n')

  # Generate temp output file
  if (verbose) cat('- blast word size: ')
  blast_out = paste0(db, '_blast_output.txt')
  ws = min_seq_len(fasta) - 50
  if (verbose) cat(ws, '\n')

  # BLAST fasta_in vs db
  if (verbose) cat('- running blast...')
  blast_status = blast(fasta, db, blast_path=get_blast_path(),
                       output=blast_out, max_target_seqs=50,
                       outfmt=paste0('6 ', blast_coll_fmt),
                       perc_identity=100, word_size=ws)
  if (blast_status != 0) {
     stop('\nError in BLAST alignment step.')
  }
  if (verbose) cat('OK.\n')
  # Load blast output and remove the temp blast output file
  if (verbose) cat('- selecting ends-free alignments...')
  blast.dt = fread(blast_out)

  names(blast.dt) = strsplit(blast_coll_fmt, ' ', fixed=T)[[1]]
  blast2.dt = blast.dt[(qseqid != sseqid) & (pident==100)]
  # Select only ends-free alignments
  blast3.dt = blast2.dt[(qstart==1 & qend==qlen) | (sstart==1 & send==slen) |
                        (qstart==1 & send==slen) | (sstart==1 & qend==qlen)]
  blast3.dt[, pair := paste0(min(.BY[[1]], .BY[[2]]), ',',
                             max(.BY[[1]], .BY[[2]])), by=.(qseqid, sseqid)]
  # Select only 1 alignment from each pair
  blast3.dt = blast3.dt[, .SD[1], by=pair]
  if (verbose) cat('OK.\n')

  # Addup counts to the sequence with more total reads
  if (verbose) cat('- adding up abundances...')
  qs_names = names(blast.dt)[1:2]
  filtered_columns = c()
  for (i in 1:nrow(blast3.dt)) {
    num_reads = c(sum(ab[, blast3.dt[i, qseqid]]), sum(ab[, blast3.dt[i, sseqid]]))
    seq1 = qs_names[which.max(num_reads)]
    seq2 = qs_names[which.min(num_reads)]
    if (seq1 == seq2) { # In the unlikely case of tie, doesn't matter which to use
       seq1 = qs_names[1]
       seq2 = qs_names[2]
    }
    ab[, blast3.dt[i, get(seq1)]] = ab[, blast3.dt[i, get(seq1)]] + ab[, blast3.dt[i, get(seq1)]]
    filtered_columns = c(filtered_columns, blast3.dt[i, get(seq2)])
  }
  ab = ab[, setdiff(1:ncol(ab), sort(filtered_columns))]
  ab_colsums = unname(colSums(ab))
  if (verbose) cat('OK.\n')


  # How to collapse
  # Add the first sequence to the queue (or list)


  # Cleanup the temp database
  if (verbose) cat('- cleaning up temporary files...')
  file.remove(blast_out)
  file.remove(fasta)
  cleanup_blastdb(db)
  if (verbose) cat('OK.\n')

  return(ab[, order(-ab_colsums), drop=F])
}


ab_to_files = function (ab) {
   # Convert the abundance table matrix into a FASTA file, together
   # with the read counts.
   meta = colnames(ab)
   colnames(ab) = 1:ncol(ab)

   rand_id = sample(LETTERS, 10)
   # Generate a temp fasta file and blast database prefix
   fasta_out = file.path(
      tempdir(),
      paste(c('collapse_', rand_id, '.fasta'), collapse='')
   )
   db_prefix = file.path(dirname(fasta_out), paste(c('collapse_', rand_id), collapse=''))

   # Generate a temp abundance table.
   tab_out = file.path(
      tempdir(),
      paste(c('collapse_', rand_id, '.txt'), collapse='')
   )

   tab.dt = data.table(ab)
   tab.dt[, sample_id := rownames(ab)]
   tab_m.dt = melt(tab.dt, variable.name='seq_id', value.name='count',
                   id.vars='sample_id')

   # Write a temp fasta
   fasta_writer(colnames(ab), meta, fasta_out)

   # Generate a BLAST database from this fasta file
   makeblastdb(fasta_out, db_prefix)

   # Return file names
   return(c(fasta_out, db_prefix))
}


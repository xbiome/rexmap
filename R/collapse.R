# tempdir() to get temporary dir path

get_makeblastdb_path = function () {
  # Find makeblastdb location
  out = suppressWarnings(system2('which', 'makeblastdb', stdout=T))
  # If makeblastdb is not found return an error
  if (is.null(attr(out, 'status'))) return(out)
  if (attr(out, 'status') == 1) stop('Error: makeblastdb command not found.')
}

makeblastdb = function (fasta_in, db_out, verbose=T,
                        makeblastdb_path=himap_option('path_makeblastdb')) {
  # Generate a BLAST database from fasta_in
  # db_out contains a full location and prefix for BLASTDB files
  # This function does not return anything, but can optionally
  # (if verbose == T) write the status to console.
  # makeblastdb_path = get_makeblastdb_path()
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

blastn = function (
  seqs_fa, output, region=NULL, ref_db=NULL,
  blast_path=himap_option('path_blastn'),
  match=himap_option('aln_params')[1], mismatch=himap_option('aln_params')[2],
  gapopen=-himap_option('aln_params')[3], gapextend=-himap_option('aln_params')[4],
  word_size=himap_option('blast_word_size'), ncpu=himap_option('ncpu'),
  max_target_seqs=himap_option('blast_max_seqs'),
  perc_identity=75, outfmt=paste0('6 ', himap_option('blast_out_fmt')),
  dust='20 64 1',
  output_err=F) {

  dbs = himap_option('blast_dbs')
  if (!is.null(region)) {
    # Region is specificed, ignore ref_db
    if (!(region %in% dbs[, Hypervariable_region])) {
      stop('blast: wrong hypervariable region specified.')
    } else {
      if (dbs[Hypervariable_region==region, DB] == '') {
        stop('blast: no database listed for this hypervariable region.')
      } else {
        ref_db = system.file('database',
                             paste0(dbs[Hypervariable_region==region, DB], '.nhr'),
                             package='himap')
        if (ref_db == '') stop('blast: missing database.')
        else ref_db = sub('.nhr$', '', ref_db)
      }
    }
  } else {
    # ref_db is specified so just use that and check if the files are there.
    if (!all(file.exists(paste0(ref_db, c('.nin', '.nhr', '.nsq'))))) {
      stop('blast: missing database file(s).')
    }
  }
  blast_args = c(
      '-dust', shQuote(dust), '-word_size', word_size,
      '-reward', match, '-penalty', mismatch, '-gapopen', gapopen,
      '-gapextend', gapextend, '-outfmt', shQuote(outfmt),
      '-query', seqs_fa, '-db', ref_db, '-num_threads', ncpu,
      '-max_target_seqs', max_target_seqs, '-perc_identity', perc_identity
      # '-qcov_hsp_perc', query_coverage_pct
  )
  if (output_err != F) cat(paste(blast_args, sep=' '), fill=T)
  x = system2(blast_path, args=blast_args, stdout=output, stderr=output_err)
  return(x)
}

# Normally formatted timestamp (to be used in file names)
ts = function () return(sub(' ', '_', gsub(':', '-', Sys.time())))

# Find minimal sequence length among the first N sequences from the fasta file
min_seq_len = function (fasta_files, n=50) {
   return(min(sapply(fasta_files, function (f) {
     return(nchar(fasta_reader(f)$seqs[1:n]))
   }), na.rm=T))
}

#' Collapse sequences that are exact matches, up to shifts and/or length
#'
#' @importFrom igraph make_graph
#' @importFrom igraph groups
#' @importFrom igraph clusters
#' @export
collapse = function (ab_in, verbose=himap_option('verbose')) {
  # Provide ab_tab_nochim as an input argument

  # To do: Check that collapse pulls together more than 2 sequences.

  # Generate temp files, this automatically generates file names
  if (verbose) cat('collapse:', fill=T)
  if (verbose) cat('* generating temporary files...')
  ab = copy(ab_in)
  out_files = ab_to_files(ab, verbose=verbose)
  fasta = out_files[1]
  db = out_files[2]
  if (verbose) cat('OK.\n')

  # Generate temp output file
  if (verbose) cat('* blast word size: ')
  blast_out = paste0(db, '_blast_output.txt')
  ws = max(round(min_seq_len(fasta) * 0.8), 7)
  if (verbose) cat(ws, '\n')

  # BLAST fasta_in vs db
  if (verbose) cat('* running blast...')
  blast_status = blastn(fasta, ref_db=db,
                        output=blast_out, max_target_seqs=50,
                        outfmt=paste0('6 ', himap_option('blast_coll_fmt')),
                        perc_identity=100, word_size=ws)
  cat('blast status: ', blast_status, fill=T)
  if (blast_status != 0) {
     stop('\nError in BLAST alignment step.')
  }
  if (verbose) cat('OK.\n')
  # Load blast output and remove the temp blast output file
  if (verbose) cat('* selecting ends-free alignments...')
  blast.dt = data.table::fread(blast_out)

  names(blast.dt) = strsplit(himap_option('blast_coll_fmt'), ' ', fixed=T)[[1]]
  blast2.dt = blast.dt[(qseqid != sseqid) & (pident==100)]
  # Select only ends-free alignments
  blast3.dt = blast2.dt[(qstart==1 & qend==qlen) | (sstart==1 & send==slen) |
                        (qstart==1 & send==slen) | (sstart==1 & qend==qlen)]
  if (nrow(blast3.dt) == 0) {
    # This means there aren't any sequences that need collapsing.
    # Return the input back then.
    if (verbose) cat('OK.\n')
    if (verbose) cat('* no sequences need collapsing.\n')
    if (verbose) cat('* cleaning up temporary files...')
    file.remove(blast_out)
    file.remove(fasta)
    cleanup_blastdb(db)
    if (verbose) cat('OK.\n')
    if (verbose) cat('* returning input.\n')
    return(ab_in)
  }

  # Find all connected clusters of sequence IDs
  g = make_graph(as.vector(t(blast3.dt[, .(qseqid, sseqid)])), directed=F)
  # Find all connected clusters
  cls = groups(clusters(g))
  filtered_columns = c()
  for (cl in cls) {
    if (length(cl) > 1) {
      column_sums = colSums(ab[, cl])
      # max_column is the sequence with max total number of reads
      max_column = cl[which(column_sums==max(column_sums))]
      # add all other columns to that one
      for (id in cl[cl != max_column]) {
        ab[, max_column] = ab[, max_column] + ab[, id]
        filtered_columns = c(filtered_columns, id)
      }
    }
  }
  ab = ab[, setdiff(1:ncol(ab), sort(filtered_columns))]

  ab_colsums = unname(colSums(ab))
  if (verbose) cat('OK.\n')

  # Cleanup the temp database
  if (verbose) cat('* cleaning up temporary files...')
  file.remove(blast_out)
  file.remove(fasta)
  cleanup_blastdb(db)
  if (verbose) cat('OK.\n')

  # Output abundance table sorted by total abundance (add up samples) in
  # descending order.
  return(ab[, order(-ab_colsums), drop=F])
}


ab_to_files = function (ab, verbose=T) {
   # Convert the abundance table matrix into a FASTA file, together
   # with the read counts. Internal function used by collapse().
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
   tab_m.dt = data.table::melt(tab.dt, variable.name='seq_id', value.name='count',
                   id.vars='sample_id')

   # Write a temp fasta
   fasta_writer(colnames(ab), meta, fasta_out)

   # Generate a BLAST database from this fasta file
   makeblastdb(fasta_out, db_prefix, verbose=verbose)

   # Return file names
   return(c(fasta_out, db_prefix))
}




# Pool two experiments (final HiMAP outputs) and generate a new HiMAP output
# Input consists of two lists:
#   osu_ab.list = list of osu_ab.dt data tables (output from abundance())
#   osu_seq.list = list of osu_seq.dt data tables (output from osu_sequences())
# Optionally
#   dataset_names = a list or a vector of dataset names. If it is not provided
#      it will be automatically generated as 'dataset001', 'dataset002', etc.
# Output
# verbose = verbose output
# temp_dir = directory to store temporary files
# osu_offset = Integer offset for non-exact matching OSU IDs (default: 1000000)

# Input arguments
# verbose=T
# temp_dir = tempdir()
#
# # Example FINRISK will be first
# # osu_ab_1.dt = fread('~/data/finrisk/data/himap-analysis/osu_ab.dt.txt')
# # osu_seq_1.dt = readRDS('~/data/finrisk/data/himap-analysis/osu_seq.dt')
# osu_ab_1.dt = readRDS('~/data/skorea/Rdata/osu_ab.dt')
# osu_seq_1.dt = readRDS('~/data/skorea/Rdata/osu_seq.dt')
#
# # TwinsUK will be second
# osu_ab_2.dt = readRDS('~/data/twinsuk_2014/Rdata/osu_ab.dt')
# osu_seq_2.dt = readRDS('~/data/twinsuk_2014/Rdata/osu_seq.dt')
#
# osu_offset = himap_option('osu_offset')
# dataset_names = c('SKOREA', 'TwinsUK_2014')
# Dataset names need to be added as columns to each data table before we begin
# pairwise merging

pool_himap_results = function (result_list, dataset_names = NULL) {
   # Add dataset name columns
   if (!is.null(dataset_names)){
      if (length(result_list) != length(dataset_names)) {
         stop('Error: The number of input data tables (', length(result_list),
              ') doesn\'t match the number of dataset names (',
              length(dataset_names), ').')
      } else {
         if (length(unique(dataset_names)) < 2) {
            stop('Error: Dataset names must be unique.')
         }
         for (i in 1:length(dataset_names)) {
            result_list[[i]][[1]]$dataset = dataset_names[i]
         }
      }
   }
   pooled_list = Reduce(pool_two_himap_outputs, result_list)
   return(pooled_list)
}

fremove = function (filename) {
   if (file.exists(filename)) {
      file.remove(filename)
   }
}

vcat = function (..., verbose=himap_option('verbose')) {
   if (verbose) {
      cat(...)
   }
}

lu = function (x) length(unique(x))

pool_two_himap_outputs = function (osu_1, osu_2,
                                   osu_offset=himap_option('osu_offset'),
                                   verbose=T, temp_dir=tempdir()) {

   # SET _1 = SUBJECT
   # SET _2 = QUERY
   pool_blast = function (qseq_1.dt, qseq_2.dt, temp_dir=tempdir(), verbose=T) {

      # Generate temp FASTA filenames
      rand_id = sample(LETTERS, 10)
      fasta_out_1 = file.path(
         temp_dir, paste(c('pool_', rand_id, '_1.fasta'), collapse='')
      )
      fasta_out_2 = file.path(
         temp_dir, paste(c('pool_', rand_id, '_2.fasta'), collapse='')
      )
      db_prefix_1 = file.path(dirname(fasta_out_1),
                              paste(c('pool_', rand_id, '_1'), collapse=''))

      # Generate temp BLAST DB
      fasta_writer(qseq_1.dt[, qseqid], qseq_1.dt[, sequence], fasta_out_1)
      fasta_writer(qseq_2.dt[, qseqid], qseq_2.dt[, sequence], fasta_out_2)

      # Generate a BLAST database from this fasta file
      makeblastdb(fasta_out_1, db_prefix_1, verbose=verbose)

      # Automatically generate a good BLAST word size for set _1
      # Find the dataset that has shorter average sequence length
      nch_1 = qseq_1.dt[, mean(nchar(sequence))]
      nch_2 = qseq_2.dt[, mean(nchar(sequence))]
      shorter_set = 0
      if (nch_1 <= nch_2) {
         ws = max(round(min_seq_len(fasta_out_1) * 0.5), 7)
         shorter_set = 1
      } else {
         ws = max(round(min_seq_len(fasta_out_2) * 0.5), 7)
         shorter_set = 2
      }

      # Run BLAST and generate output table
      blast_out = file.path(
         temp_dir, paste(c('pool_', rand_id, '_blast_output.txt'), collapse='')
      )
      # QUERY IS SET _2 (qseqid = qseqid_2)
      # SUBJECT IS SET _1 (sseqid = qseqid_1)
      blast_status = blastn(fasta_out_2, ref_db=db_prefix_1,
                            output=blast_out, max_target_seqs=NULL,
                            match=1, mismatch=-2, gapopen=0, gapextend=2,
                            outfmt=paste0('6 ', himap_option('blast_coll_fmt')),
                            perc_identity=100, word_size=ws, verbose=T)
      vcat('blast status: ', blast_status, fill=T)

      # Cleanup temporary filenames
      fremove(fasta_out_1)
      fremove(fasta_out_2)
      fremove(paste0(db_prefix_1, '.nhr'))
      fremove(paste0(db_prefix_1, '.nin'))
      fremove(paste0(db_prefix_1, '.nseq'))

      # Load and parse BLAST results
      if (blast_status != 0) {
         stop('\nError in BLAST alignment step.')
      }
      vcat('OK.\n')

      # Load blast output and remove the temp blast output file
      vcat('* selecting ends-free alignments')
      blast.dt = data.table::fread(blast_out)
      vcat('.')
      names(blast.dt) = strsplit(himap_option('blast_coll_fmt'), ' ', fixed=T)[[1]]

      # Remove last temporary file (blast output)
      fremove(blast_out)

      # First pool out qseqids that do not show up in the alignment from both sets
      # unaligned_qseqids_1 = setdiff(qseq_1.dt[, qseqid], blast.dt[, qseqid])
      # unaligned_qseqids_2 = setdiff(qseq_2.dt[, qseqid], blast.dt[, sseqid])

      # For the pairs of qseqids that do match, keep only the longest alignments
      #### blast2.dt = blast.dt[, .SD[length==max(length)], by=.(qseqid)]
      blast2.dt = copy(blast.dt)
      vcat('.')
      # Filter out partial alignments !!!
      #
      #   Keep:
      #
      #   ------------------------------------------------------ 1 (q)
      #   ------------------------------------------------------ 2 (s)
      #
      #
      #      -----------------------   qstart=1, qend=qlen (fully aligned 1st)
      #     -------------------------
      #
      #     -------------------------
      #        --------------------     sstart=1, send=slen (fully aligned 2nd)
      #
      #       --------------------------- qstart=1, send=slen
      #    -------------------------
      #
      #   -----------------------------
      #      -----------------------------  qend=qlen, sstart=1
      #
      #  ------------------------------
      #  ---------------------------    sstart=1, send=slen
      #
      # Select only ends-free alignments
      blast3.dt = blast2.dt[
         (qstart==1 & qend==qlen) | (sstart==1 & send==slen) |
            (qstart==1 & send==slen) | (sstart==1 & qend==qlen)
         ]
      vcat('.')
      if (shorter_set == 2) {
         # Set _2 is shorter (qseqids in this table)
         blast3.dt[, qseqid_new := .GRP, by=qseqid]
      } else {
         # Set _1 is shorter (sseqids) so collapse to that one
         blast3.dt[, qseqid_new := .GRP, by=sseqid]
      }

      vcat(' OK.\n')


      # if (shorter_set == 2) {
      #    # _1 = sseqid (in debug this is SKOREA, longer)
      #    # _2 = qseqid (in debug example this is TwinsUK, shorter)
      #    # If longer sequence maps to multiple shorter, then we have a situation
      #    # like this:
      #    #                          -------------  (shorter)
      #    #  --------------------------------       (longer)
      #    # In these cases, let's not combine the reads, because then everything
      #    # else will have to be shortened to:
      #    #                          --------
      #    blast4.dt = blast3.dt[, if (length(unique(qseqid)) == 1) .SD, by=sseqid]
      #    # Now sort out the OSUs
      # } else {
      #    # Shorter_seq == 1
      #    # this is just like prev if block, but swap qseqid <--> sseqid
      #    blast4.dt = blast3.dt[, if (length(unique(sseqid)) == 1) .SD, by=qseqid]
      # }
      blast4.dt = copy(blast3.dt)
      return(list(blast4.dt, shorter_set))
   }


   # osu_1 = list(osu_ab_1.dt, osu_seq_1.dt)
   vcat('\n')
   vcat('POOLING: \n')
   vcat('- Dataset 1: ', paste(unique(osu_1[[1]]$dataset), collapse=','), '\n')
   vcat('- Dataset 2: ', paste(unique(osu_2[[1]]$dataset), collapse=','), '\n')
   vcat('\n')

   osu_ab_1.dt = osu_1[[1]]
   osu_seq_1.dt = osu_1[[2]]
   osu_ab_2.dt = osu_2[[1]]
   osu_seq_2.dt = osu_2[[2]]


   osu_seq_1.dt[, sequence := as.character(sequence)]
   # osu_seq_1.dt[, nqs := lu(qseqid), by=osu_id]
   osu_seq_2.dt[, sequence := as.character(sequence)]
   # osu_seq_2.dt[, nqs := lu(qseqid), by=osu_id]

   osu_ab_1_colorder = names(osu_ab_1.dt)

   # To update the counts of all species, we will have to expand things like
   # Escherichia_coli_[3] to Escherichia_coli,Escherichia_coli,Escherichia_coli

   species_expand = function (x, sep=',', coll=T) {
      # Internal function for expanding collapsed strains from print_strains
      # Escherichia_coli_[3] to Escherichia_coli,Escherichia_coli,Escherichia_coli
      x_split = strsplit(x, sep, fixed=T)[[1]]
      x_split_expanded = lapply(x_split, function (x_i) {
         # Extract the number
         regex_count = sub('.*_\\[([0-9]+)\\]$', '\\1', x_i)
         if (regex_count == x_i) {
            # There is no count for this species/strain so just return as is
            return(x_i)
         } else {
            x_i_name = sub('^(.*)_\\[[0-9]+\\]$', '\\1', x_i)
            x_i_new = rep(x_i_name, as.integer(regex_count))
            x_i_new = paste(
               x_i_new,
               1:length(x_i_new),
               sep='_'
            )
            if (coll==TRUE) {
               return(paste(x_i_new, collapse=sep))
            } else {
               return(x_i_new)
            }

         }
      })
      if (coll==TRUE) {
         return(paste(x_split_expanded, collapse=sep))
      } else {
         return(unname(x_split_expanded))
      }
   }

   vector_split = function (x, sep=',') {
      strsplit(x, ',')[[1]]
   }

   species_relabel = function (x, sep=',') {
      print_strains(vector_split(species_expand(x, sep=sep)))
   }


   # Extract sequences
   qseq_1.dt = na.omit(unique(osu_seq_1.dt[, .(qseqid, sequence)])[order(qseqid)])
   qseq_2.dt = na.omit(unique(osu_seq_2.dt[, .(qseqid, sequence)])[order(qseqid)])

   pool_blast_result = pool_blast(qseq_1.dt, qseq_2.dt, temp_dir=temp_dir,
                                  verbose=verbose)
   blast4.dt = pool_blast_result[[1]]
   shorter_set = pool_blast_result[[2]]
   rm(pool_blast_result)

   # Generate BLAST DB from the first set, FASTA File from the second set then run
   # BLAST.


   vcat('* preparing blast tables..')
   osu_seq_map_1.dt = unique(osu_seq_1.dt[
      , .(osu_id, qseqid, copy_number, sequence)])

   osu_seq_map_2.dt = unique(osu_seq_2.dt[
      , .(osu_id, qseqid, copy_number, sequence)])


   # Merge Set 1 with BLAST results, find shared OSUs
   osu_seq_map_blast_1.dt = merge(
      osu_seq_map_1.dt,
      blast4.dt[, .(sseqid, qseqid2=qseqid, qseqid_new, qlen, slen,
                    qstart, qend, sstart, send, length)],
      # blast4.dt,
      by.x='qseqid', by.y='sseqid', all.x=T, allow.cartesian=T
   )
   # osu_seq_map_blast_1.dt[, sequence := substring(sequence, sstart[1], send[1]),
   #                        by=qseqid]
   osu_seq_map_blast_1.dt = osu_seq_map_blast_1.dt[
      , .SD[qlen==max(qlen)], by=.(qseqid, osu_id)]
   osu_seq_map_blast_1.dt[
      , sequence := substring(sequence, max(sstart), min(send)),
                          by=qseqid]
   ## osu_seq_map_blast_1.dt = osu_seq_map_blast_1.dt[, if(all(!is.na(sequence))) .SD, by=osu_id]
   osu_seq_map_blast_1a.dt = copy(na.omit(osu_seq_map_blast_1.dt))

   vcat('.')
   osu_seq_map_blast_2.dt = merge(
      osu_seq_map_2.dt,
      blast4.dt[, .(qseqid, qseqid1=sseqid, qseqid_new, slen, qlen,
                    qstart, qend, sstart, send, length)],
      by.x='qseqid', by.y='qseqid', all.x=T, allow.cartesian=T
   )
   osu_seq_map_blast_2.dt = osu_seq_map_blast_2.dt[
      , .SD[slen==max(slen)], by=.(qseqid, osu_id)]

   osu_seq_map_blast_2.dt[
      , sequence := substring(sequence,  max(qstart), min(qend)),
      by=qseqid]

   osu_seq_map_blast_2a.dt = copy(na.omit(osu_seq_map_blast_2.dt))
   vcat('OK.\n')

   vcat('* collapsing OSUs in set _1 ')
   seqs_with_diff_qseqids = osu_seq_map_blast_1a.dt[
      , if (lu(qseqid) > 1) .SD, by=sequence][, unique(sequence)]

   #-------------------- COLLAPSE SEQ FOR SET _1 -=============================
   seqs_with_diff_qseqids = osu_seq_map_blast_2a.dt[
      , if (lu(qseqid) > 1) .SD, by=sequence][, unique(sequence)]


   # osu_seq_map_blast_2.dt[, osu_n_qseqids := lu(qseqid), by=osu_id]
   # osu_seq_map_blast_2.dt[, osu_n_qseqid_news := lu(qseqid_new), by=osu_id]
   # osu_seq_map_blast_2.dt[, copy_number := sum(copy_number), by=.(osu_id, qseqid_new)]
   # # osu_seq_map_blast_2.dt[osu_n_qseqid_news==1, osu_id , by=sequence]
   #
   # qseqid_new_collapse = osu_seq_map_blast_2.dt[
   #    osu_n_qseqid_news==1 & osu_n_qseqids != osu_n_qseqid_news,
   #    unique(qseqid_new)
   #    ]
   #
   # for (qn in qseqid_new_collapse) {
   #    qn_collapse_osus = osu_seq_map_blast_2.dt[
   #       osu_n_qseqid_news==1 & qseqid_new==qn,
   #       sort(unique(osu_id))
   #       ]
   #    qn_collapse_osus_best = qn_collapse_osus[qn_collapse_osus < osu_offset]
   #    if (length(qn_collapse_osus_best) == 0) {
   #       qn_collapse_osus_best = qn_collapse_osus
   #    }
   #    qn_new_osu = min(qn_collapse_osus_best)
   #    # Update the sequence map
   #    osu_seq_map_blast_2.dt[
   #       osu_n_qseqid_news==1 & qseqid_new==qn,
   #       osu_id := qn_new_osu
   #       ]
   #    # Update the sequence descriptions for that OSU
   #    qn_new_pctsim = osu_seq_2.dt[osu_id %in% qn_collapse_osus_best,
   #                                 max(pctsim, na.rm=T)]
   #    old_coll_species = osu_seq_2.dt[osu_id %in% qn_collapse_osus_best,
   #                                    paste(species, collapse=',')]
   #    old_species = strsplit(old_coll_species, ',')[[1]]
   #    qn_new_species = print_strains(old_species)
   #    osu_seq_2.dt[
   #       osu_id==qn_new_osu,
   #       c('species', 'pctsim') := list(qn_new_species, qn_new_pctsim)]
   # }

   for (seq in seqs_with_diff_qseqids) {
      # vcat('.')
      min_id = osu_seq_map_blast_1a.dt[sequence==seq, which(osu_id==min(osu_id))[1]]
      # min_osuid = osu_seq_map_blast_1a.dt[sequence==seq, min(osu_id)]
      all_osuids = osu_seq_map_blast_1a.dt[sequence==seq, osu_id]
      min_cp = osu_seq_map_blast_1a.dt[sequence==seq, .SD[min_id, copy_number]]
      # min_cp = max(
      #    osu_seq_map_blast_1a.dt[sequence==seq, .SD[min_id, copy_number]],
      #    osu_seq_map_blast_2a.dt[
      #       sequence==seq &
      #          qseqid==osu_seq_map_blast_1a.dt[sequence==seq, .SD[min_id, qseqid2]],
      #       copy_number
      #       ]
      # )

      # min_pctsim = osu_seq_map_blast_2.dt[sequence==seq, .SD[min_id, pctsim]]
      # min_species = osu_seq_map_blast_2.dt[sequence==seq, .SD[min_id, species]]
      min_osuid = min(all_osuids)
      osu_seq_map_blast_1.dt[
         sequence==seq, c('osu_id', 'copy_number', 'qseqid_new') := list(
            .SD[min_id, osu_id],
            .SD[min_id, copy_number],
            # min_cp,
            .SD[min_id, qseqid_new]
         )
         ]
      # Update abundance table by pooling OSU counts
      # If we already have single OSU with multiple sequences (spectra)
      # then we collapse to the first sequence
      # Update sequence reference table
      # osu_seq_1.dt[
      #    osu_id %in% all_osuids,
      #    c('osu_id', 'copy_number') := list(min_osuid, min_cp)
      # ]
   }
   osu_seq_map_blast_1.dt = unique(osu_seq_map_blast_1.dt)
   vcat(' OK.\n')

   osu_spectra_1.dt = osu_seq_map_blast_1.dt[
      , .(
         spectrum = paste(
            sort(unique(
               paste(sequence, copy_number, sep=':')
            )),
            collapse=','
          ),
          len_1=slen[1]),
      by=osu_id
   ]

   #------------------ Collapse sequences for Set _2 ---------------------------
   vcat('* collapsing OSUs in set _2 ')

   seqs_with_diff_qseqids = osu_seq_map_blast_2a.dt[
      , if (lu(qseqid) > 1) .SD, by=sequence][, unique(sequence)]

#
#    osu_seq_map_blast_2.dt[, osu_n_qseqids := lu(qseqid), by=osu_id]
#    osu_seq_map_blast_2.dt[, osu_n_qseqid_news := lu(qseqid_new), by=osu_id]
#    osu_seq_map_blast_2.dt[, copy_number := sum(copy_number), by=.(osu_id, qseqid_new)]
#    # osu_seq_map_blast_2.dt[osu_n_qseqid_news==1, osu_id , by=sequence]
#
#    qseqid_new_collapse = osu_seq_map_blast_2.dt[
#       osu_n_qseqid_news==1 & osu_n_qseqids != osu_n_qseqid_news,
#       unique(qseqid_new)
#    ]
#
#    for (qn in qseqid_new_collapse) {
#       qn_collapse_osus = osu_seq_map_blast_2.dt[
#          osu_n_qseqid_news==1 & qseqid_new==qn,
#          sort(unique(osu_id))
#       ]
#       qn_collapse_osus_best = qn_collapse_osus[qn_collapse_osus < osu_offset]
#       if (length(qn_collapse_osus_best) == 0) {
#          qn_collapse_osus_best = qn_collapse_osus
#       }
#       qn_new_osu = min(qn_collapse_osus_best)
#       # Update the sequence map
#       osu_seq_map_blast_2.dt[
#          osu_n_qseqid_news==1 & qseqid_new==qn,
#          osu_id := qn_new_osu
#       ]
#       # Update the sequence descriptions for that OSU
#       qn_new_pctsim = osu_seq_2.dt[osu_id %in% qn_collapse_osus_best,
#                                    max(pctsim, na.rm=T)]
#       old_coll_species = osu_seq_2.dt[osu_id %in% qn_collapse_osus_best,
#                                       paste(species, collapse=',')]
#       old_species = strsplit(old_coll_species, ',')[[1]]
#       qn_new_species = print_strains(old_species)
#       osu_seq_2.dt[
#          osu_id==qn_new_osu,
#          c('species', 'pctsim') := list(qn_new_species, qn_new_pctsim)]
#    }

   # osu_seq_map_blast_2.dt[osu_n_qseqid_news]

      for (seq in seqs_with_diff_qseqids) {
         # vcat('.', which(seqs_with_diff_qseqids == seq), '.')
         min_id = osu_seq_map_blast_2.dt[sequence==seq, which(osu_id==min(osu_id))[1]]
         all_osuids = osu_seq_map_blast_2.dt[sequence==seq, unique(osu_id)]
         min_cp = osu_seq_map_blast_2.dt[sequence==seq, .SD[min_id, copy_number]]
         # Always choose the larger copy number between the two matched sequences
         # min_cp = max(
         #    osu_seq_map_blast_2.dt[sequence==seq, .SD[min_id, copy_number]],
         #    osu_seq_map_blast_1.dt[
         #       sequence==seq &
         #       qseqid==osu_seq_map_blast_2.dt[sequence==seq, .SD[min_id, qseqid1]],
         #       copy_number
         #    ]
         # )
         # min_pctsim = osu_seq_map_blast_2.dt[sequence==seq, .SD[min_id, pctsim]]
         # min_species = osu_seq_map_blast_2.dt[sequence==seq, .SD[min_id, species]]
         min_osuid = min(all_osuids)
         osu_seq_map_blast_2.dt[
            sequence==seq, c('osu_id', 'copy_number', 'qseqid_new') := list(
               .SD[min_id, osu_id],
               .SD[min_id, copy_number],
               # min_cp,
               .SD[min_id, qseqid_new]
            )
         ]
         # Update abundance table by pooling OSU counts
         # If we already have single OSU with multiple sequences (spectra)
         # then we collapse to the first sequence
         # osu_ab_2.dt[
         #    osu_id %in% all_osuids,
         #    c('osu_id', 'osu_count', 'pctsim', 'species') := list(
         #       min_osuid, sum(.SD[, osu_count[1], by=osu_id]),
         #       .SD[osu_id==min_osuid, pctsim],
         #       .SD[osu_id==min_osuid, species]), by=.(sample_id, dataset)]
         # Update sequence reference table
         # osu_seq_2.dt[
         #    osu_id %in% all_osuids,
         #    c('osu_id', 'copy_number') := list(min_osuid, min_cp)
         # ]
      }
   osu_seq_map_blast_2.dt = unique(osu_seq_map_blast_2.dt)
   vcat(' OK.\n')

   ## osu_seq_map_blast_2.dt = osu_seq_map_blast_2.dt[, if(all(!is.na(sequence))) .SD, by=osu_id]
   osu_spectra_2.dt = osu_seq_map_blast_2.dt[
      , .(
         spectrum = paste(
            sort(unique(
               paste(sequence, copy_number, sep=':')
            )),
            collapse=','
         ),
         len_2=qlen[1]
      ),by=osu_id
   ]

   osu_matching_spectra.dt = merge(
      osu_spectra_1.dt[, .(osu_id_1=osu_id, spectrum, len_1)],
      osu_spectra_2.dt[, .(osu_id_2=osu_id, spectrum, len_2)],
      by='spectrum'
   )[order(osu_id_1)]

   # Renumerate matching OSUs
   # First exact matches
   if (shorter_set == 1) {
      # osu_matching_spectra.dt[osu_id_1 < osu_offset,
      #                         osu_id_new := .GRP, by=spectrum]
      # osu_matching_spectra.dt[osu_id_1 >= osu_offset,
      #                         osu_id_new := .GRP + as.integer(osu_offset),
      #                         by=spectrum]
      osu_matching_spectra.dt[
         , osu_id_new := ifelse(
            osu_id_1 < osu_offset,
            .GRP,
            .GRP + as.integer(osu_offset)
         )
         , by=spectrum]
   } else {  # shorter_set == 2
      # osu_matching_spectra.dt[osu_id_2 < osu_offset,
      #                         osu_id_new := .GRP, by=spectrum]
      # osu_matching_spectra.dt[osu_id_2 >= osu_offset,
      #                         osu_id_new := .GRP + as.integer(osu_offset),
      #                         by=spectrum]
      osu_matching_spectra.dt[
         , osu_id_new := ifelse(
            osu_id_2 < osu_offset,
            .GRP,
            .GRP + as.integer(osu_offset)
         )
         , by=spectrum]
   }

   # After sequence trimming sometimes < 100% match can collapse to 100% match
   # so the mapping will go to the higher pctsim OSU
   osu_matching_spectra.dt[, osu_id_new := min(osu_id_new), by=osu_id_1]
   osu_matching_spectra.dt[, osu_id_new := min(osu_id_new), by=osu_id_2]

   # Renumerate OSUs from _1 that are not matching
   vcat('* renumerating OSUs _1 ...')
   osu_renum_1.dt = unique(osu_seq_1.dt[, .(osu_id_1 = osu_id)])[
      !(osu_id_1 %in% osu_matching_spectra.dt[, osu_id_1])]
   osu_100_last = osu_matching_spectra.dt[osu_id_new < as.integer(osu_offset), max(osu_id_new)]
   osu_non100_last = osu_matching_spectra.dt[osu_id_new >= as.integer(osu_offset), max(osu_id_new)]
   osu_renum_1.dt[osu_id_1 < osu_offset,
                  osu_id_new := .GRP + osu_100_last, by=osu_id_1]
   osu_renum_1.dt[osu_id_1 >= osu_offset,
                  osu_id_new := .GRP + osu_non100_last, by=osu_id_1]
   osu_renum_all_1.dt = rbindlist(list(
      osu_renum_1.dt, unique(osu_matching_spectra.dt[, .(osu_id_1, osu_id_new)])
   ))[order(osu_id_new)]
   vcat(' OK.\n')
   # Renumerate OSUs from _2 that are not matching
   vcat('* renumerating OSUs _2 ...')
   osu_renum_2.dt = unique(osu_seq_2.dt[, .(osu_id_2 = osu_id)])[
      !(osu_id_2 %in% osu_matching_spectra.dt[, osu_id_2])]
   osu_100_last = osu_renum_1.dt[osu_id_new < osu_offset, max(osu_id_new)]
   osu_non100_last = osu_renum_1.dt[osu_id_new >= osu_offset, max(osu_id_new)]
   osu_renum_2.dt[osu_id_2 < as.integer(osu_offset),
                  osu_id_new := .GRP + osu_100_last, by=osu_id_2]
   osu_renum_2.dt[osu_id_2 >= as.integer(osu_offset),
                  osu_id_new := .GRP + osu_non100_last, by=osu_id_2]
   osu_renum_all_2.dt = rbindlist(list(
      osu_renum_2.dt, unique(osu_matching_spectra.dt[, .(osu_id_2, osu_id_new)])
   ))[order(osu_id_new)]
   vcat(' OK.\n')


   # Add new OSU IDs to abundance tables
   vcat('* adding new OSUs .')
   osu_ab_1_new.dt = merge(osu_ab_1.dt, osu_renum_all_1.dt, by.x='osu_id',
                           by.y='osu_id_1', all.x=T, all.y=T)

   osu_seq_1_new.dt = merge(osu_seq_1.dt, osu_renum_all_1.dt, by.x='osu_id',
                           by.y='osu_id_1', all.x=T, all.y=T)

   # Update non-matching OSUs with trimmed sequences, e.g.
   # we have a matching spectra 5:s1, 1:s2, 1:s3
   # but _2 contains 1:s1 spectra that is not matched in _1
   osu_seq_1_new.dt = merge(
      osu_seq_1_new.dt,
      unique(osu_seq_map_blast_1.dt[, .(qseqid, seq_new=sequence)]),
      by='qseqid', all.x=T
   )
   osu_seq_1_new.dt[!is.na(seq_new) & sequence != seq_new,
                    sequence := seq_new]
   osu_seq_1_new.dt[, seq_new := NULL]
   setcolorder(osu_seq_1_new.dt, c('osu_id', 'species', 'qseqid', 'copy_number',
                                   'pctsim', 'sequence', 'osu_id_new'))


   # If _2 set is shorter then update the sequences here in _1 with the shorter
   # ones from _2, for matching OSUs.

   if (osu_ab_1_new.dt[is.na(osu_id), .N] > 0 |
       osu_ab_1_new.dt[is.na(osu_id_new), .N] > 0) {
      stop('Error: Missing or extra OSU IDs in _1 after renumeration.')
   }
   osu_ab_1_new.dt[, osu_id := osu_id_new]
   osu_ab_1_new.dt[, osu_id_new := NULL]
   osu_ab_1_new.dt[, osu_count_new := sum(osu_count), by=.(osu_id, sample_id)]
   # Pool osu_counts for OSUs that have been collapsed because they exactly
   # matched shorter sequence from the other dataset.
   osu_ab_1_new.dt[, osu_count := osu_count_new]
   osu_ab_1_new.dt[, c('osu_count_new', 'pctsim', 'species') := NULL]
   osu_ab_1_new.dt = unique(osu_ab_1_new.dt)

   vcat('.')

   # DO THE SAME WITH OSU_SEQs
   osu_ab_2_new.dt = merge(osu_ab_2.dt, osu_renum_all_2.dt, by.x='osu_id',
                           by.y='osu_id_2', all.x=T, all.y=T)
   osu_seq_2_new.dt = merge(osu_seq_2.dt, osu_renum_all_2.dt, by.x='osu_id',
                            by.y='osu_id_2', all.x=T, all.y=T)

   # Update non-matching OSUs with trimmed sequences, e.g.
   # we have a matching spectra 5:s1, 1:s2, 1:s3
   # but _2 contains 1:s1 spectra that is not matched in _1
   osu_seq_2_new.dt = merge(
      osu_seq_2_new.dt,
      unique(osu_seq_map_blast_2.dt[, .(qseqid, seq_new=sequence)]),
      by='qseqid', all.x=T
   )
   osu_seq_2_new.dt[!is.na(seq_new) & sequence != seq_new,
                    sequence := seq_new]
   osu_seq_2_new.dt[, seq_new := NULL]
   setcolorder(osu_seq_2_new.dt, c('osu_id', 'species', 'qseqid', 'copy_number',
                                   'pctsim', 'sequence', 'osu_id_new'))

   if (osu_ab_2_new.dt[is.na(osu_id), .N] > 0 |
       osu_ab_2_new.dt[is.na(osu_id_new), .N] > 0) {
      stop('Error: Missing or extra OSU IDs in _2 after renumeration.')
   }
   osu_ab_2_new.dt[, osu_id := osu_id_new]
   osu_ab_2_new.dt[, osu_id_new := NULL]
   osu_ab_2_new.dt[, osu_count_new := sum(osu_count), by=.(osu_id, sample_id)]
   # Pool osu_counts for OSUs that have been collapsed because they exactly
   # matched shorter sequence from the other dataset.
   osu_ab_2_new.dt[, osu_count := osu_count_new]
   osu_ab_2_new.dt[, c('osu_count_new', 'pctsim', 'species') := NULL]
   osu_ab_2_new.dt = unique(osu_ab_2_new.dt)

   vcat('.')

   # Now sort out OSU_SEQ tables
   # From the matching table
   # If _1 is longer than _2 then add shorter seqs of _2 to _1

   osu_seq_1_update.dt = merge(
      osu_matching_spectra.dt[len_1 >= len_2],
      osu_seq_2_new.dt,
      by.x=c('osu_id_2', 'osu_id_new'),
      by.y=c('osu_id', 'osu_id_new')
   )

   osu_seq_1_update.dt = osu_seq_1_update.dt[
      , .(osu_id=osu_id_1, species, qseqid, copy_number, pctsim, sequence,
          osu_id_new)]

   osu_seq_1_new_all.dt = data.table::rbindlist(list(
      osu_seq_1_new.dt[!(osu_id %in% osu_seq_1_update.dt[, unique(osu_id)])],
      osu_seq_1_update.dt
   ))

   vcat('.')


   # 1 is shorter, so update 2 with seqs from 1
   osu_seq_2_update.dt = merge(
      osu_matching_spectra.dt[len_1 < len_2],
      osu_seq_1_new.dt,
      by.x=c('osu_id_1', 'osu_id_new'),
      by.y=c('osu_id', 'osu_id_new')
   )

   osu_seq_2_update.dt = osu_seq_2_update.dt[
      , .(osu_id=osu_id_2, species, qseqid, copy_number, pctsim, sequence,
          osu_id_new)]

   osu_seq_2_new_all.dt = data.table::rbindlist(list(
      osu_seq_2_new.dt[!(osu_id %in% osu_seq_2_update.dt[, unique(osu_id)])],
      osu_seq_2_update.dt
   ))
   vcat('.')
   osu_seq_1_new_all_f.dt = copy(osu_seq_1_new_all.dt)
   osu_seq_1_new_all_f.dt[, osu_id := osu_id_new]
   osu_seq_1_new_all_f.dt[, osu_id_new := NULL]
   osu_seq_2_new_all_f.dt = copy(osu_seq_2_new_all.dt)
   osu_seq_2_new_all_f.dt[, osu_id := osu_id_new]
   osu_seq_2_new_all_f.dt[, osu_id_new := NULL]

   #------------------ Finalize OSU_SEQ TABLE ----------------------------------
   osu_seq_new.dt = unique(data.table::rbindlist(list(
      osu_seq_1_new_all_f.dt,
      osu_seq_2_new_all_f.dt
   )))
   # Finally, renumerate all qseqid by sequence
   osu_seq_new.dt[, qseqid := NULL]
   osu_seq_new.dt[, qseqid := .GRP, by=sequence]

   osu_seq_1_ref.dt = unique(osu_seq_1_new_all_f.dt[, .(osu_id, pctsim, species)])
   osu_seq_1_ref.dt = osu_seq_1_ref.dt[
      , .(species = print_strains(species),
          pctsim = unique(pctsim[!is.na(pctsim)])[1]), by=osu_id]
   osu_seq_2_ref.dt = unique(osu_seq_2_new_all_f.dt[, .(osu_id, pctsim, species)])
   osu_seq_2_ref.dt = osu_seq_2_ref.dt[
      , .(species = print_strains(species),
          pctsim = unique(pctsim[!is.na(pctsim)])[1]), by=osu_id]
   # Add new OSU annotations to OSU AB table
   osu_ab_1_new.dt = merge(
      osu_ab_1_new.dt,
      osu_seq_1_ref.dt,   # <-------------- change this
      by='osu_id' #, allow.cartesian=T
   )
   setcolorder(osu_ab_1_new.dt, osu_ab_1_colorder)
   osu_ab_2_new.dt = merge(
      osu_ab_2_new.dt,
      osu_seq_2_ref.dt,
      #unique(osu_seq_2_new_all_f.dt[, .(osu_id, pctsim, species)]),
      by='osu_id'# , allow.cartesian=T
   )
   setcolorder(osu_ab_2_new.dt, osu_ab_1_colorder)

   vcat(' OK.\n')

   vcat('* merging output...')
   osu_ab_new.dt = data.table::rbindlist(list(osu_ab_1_new.dt, osu_ab_2_new.dt))
   # DO THE SAME WITH OSU_SEQ TABLES
   vcat(' OK.\n')

   vcat('* finalizing output...')

   osu_seq_new.dt[, nqs := lu(qseqid), by=osu_id]
   seqs_collapse_final = osu_seq_new.dt[
      nqs==1, if (lu(osu_id) > 1) 1, by=sequence][, unique(sequence)]
   for (s in seqs_collapse_final) {
      vcat('\n process seq ', which(seqs_collapse_final==s), '...')
      osu_seq_new_sub.dt = osu_seq_new.dt[sequence==s & nqs==1]
      s_osuids = osu_seq_new_sub.dt[, sort(unique(osu_id))]
      s_final_pctsim = osu_seq_new_sub.dt[, max(pctsim, na.rm=T)]
      s_final_osuid = osu_seq_new_sub.dt[pctsim==s_final_pctsim,
         min(osu_id)
      ]
      s_final_species = osu_seq_new_sub.dt[pctsim==s_final_pctsim,
         species_relabel(species)
      ]
      s_final_cp = osu_seq_new_sub.dt[
         pctsim==s_final_pctsim & osu_id==s_final_osuid,
         max(copy_number)
      ]
      osu_seq_new.dt[
         osu_id %in% s_osuids & nqs==1 & sequence==s,
         c('osu_id', 'species', 'pctsim', 'copy_number') := list(
            s_final_osuid, s_final_species, s_final_pctsim, s_final_cp
         )
      ]
      osu_ab_new.dt[
         osu_id %in% s_osuids,
         c('osu_id', 'osu_count', 'pctsim', 'species') := list(
            s_final_osuid, sum(osu_count), s_final_pctsim, s_final_species
         ),
         by=.(dataset, sample_id)
      ]
   }
   osu_seq_new.dt = unique(osu_seq_new.dt)
   osu_seq_new.dt[, nqs := NULL]
   osu_ab_new.dt = unique(osu_ab_new.dt)

   vcat(' OK.\n')

   return(list(osu_ab_new.dt, osu_seq_new.dt))
}

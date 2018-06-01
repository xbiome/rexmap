#' Download blastn and makeblastdb executables
#'
#' Downloads NCBI blastn and makeblastdb executables for the Operating System
#' running R and copies it to the HiMAP package folder. This is executed automatically
#' if blastn and makeblastdb executables are missing upon attaching HiMAP package.
#'
#' @param verbose (default: TRUE) Whether or not to print progress during download.
#'
#' @export
download_blast = function (verbose=T) {
   # First find the platform
   if (verbose) cat('HiMAP: Download blastn+makeblastdb command line utilities\n')
   sys_os = detect_os()
   if (sys_os == 'linux') {
      blast_url = 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+//2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz'
   } else if (sys_os == 'macos') {
      blast_url = 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+//2.7.1/ncbi-blast-2.7.1+-x64-macosx.tar.gz'
   } else {
      stop('Error: Unsupported OS.')
   }

   # Check internet connection
   if (!internet_connection()) {
      stop('Error: Internet connection unavailable.')
   }

   # Download into temp folder
   if (verbose) cat('* Downloading executables...')
   blast_file_gz = file.path(tempdir(), basename(blast_url))
   download.file(blast_url, blast_file_gz, quiet=T)
   if (verbose) cat('OK.\n')

   # Unpack
   if (verbose) cat('* Unpacking...')
   blast_files = untar(blast_file_gz, verbose=F, list=T, exdir = tempdir())
   untar(blast_file_gz, verbose=F, exdir = tempdir())
   if (verbose) cat(' OK.\n')

   # Copy to HiMAP folder
   if (verbose) cat('* Copying...')
   blast_exec = file.path(dirname(blast_file_gz),
                          grep('/bin/blastn$', blast_files, value=T))
   makedb_exec = file.path(dirname(blast_file_gz),
                          grep('/bin/makeblastdb$', blast_files, value=T))
   himap_path = find.package('himap')
   himap_blastn_path = file.path(himap_path, 'exec', paste0('blastn_', sys_os))
   himap_makedb_path = file.path(himap_path, 'exec', paste0('makeblastdb_', sys_os))
   copy_success = file.copy(blast_exec, himap_blastn_path, overwrite=T)
   if (!copy_success) stop('Error. Failed: ', blast_exec, ' --[copy]--> ', himap_blastn_path)
   copy_success = file.copy(makedb_exec, himap_makedb_path, overwrite=T)
   if (!copy_success) stop('Error. Failed: ', makedb_exec, ' --[copy]--> ', himap_makedb_path)
   if (verbose) cat(' OK.\n')

   # Clean up files in the temp folder
   if (verbose) cat('* Cleaning temporary files...')
   temp_dirs = grep('/$', blast_files, value=T)
   temp_files = blast_files[!(blast_files %in% temp_dirs)]
   file.remove(file.path(tempdir(), temp_files))
   unlink(file.path(tempdir(), temp_dirs), recursive=T)
   file.remove(blast_file_gz)
   if (verbose) cat(' OK.\nDone.\n')
}

#' Check if there's an active internet connection
#'
#' Return TRUE or FALSE
internet_connection = function(url='https://cran.r-project.org') {
    out = tryCatch(
       readChar(con=url, nchars=1),
       error = function(e) return(FALSE),
       warning = function(w) return(FALSE)
    )
    if (class(out) == 'logical') return(FALSE)
    else return(TRUE)
}

#' curl path for system2 call
curl_path = function () {
   # Try system path first
   curl_sys_path = suppressWarnings(system2(c('which', 'curl'), stdout=T, stderr=F))
   if (length(curl_sys_path) > 0) return(curl_sys_path)
   # Curl not found in system path use the one bundled in HiMAP
   return(file.path(find.package('himap'), 'exec', paste0('curl_', detect_os())))
}


#' Update the HiMAP database
#'
#' Updates the database from HiMAP GitHub repository.
#'
#' @export
update_database = function (verbose=T) {
   if (verbose) cat('HiMAP database update\n')
   if (verbose) cat('* check internet connection... ')
   if (!internet_connection()) stop('\nError: Internet connection unavailable.')
   if (verbose) cat('OK.\n')
   if (verbose) cat('* downloading database files:\n')
   himap_database_path = file.path(find.package('himap'), 'database')
   # Implements this method:
   # https://medium.com/@caludio/how-to-download-large-files-from-github-4863a2dbba3b
   json_out = paste(system2(curl_path(), c(
      '-H', '"Authorization: token 99f22e14f4ed6ec6899bebe79dbf6fd7fbf9bac6"',
      '-L', 'https://api.github.com/repos/taolonglab/himap/contents/inst/database/'
   ), stdout=T, stderr=F), collapse='')

   db.dt = as.data.table(as.list(jsonlite::fromJSON(json_out)))[, c(1,4,11)]
   # Select only database files for hypervariable regions
   db.dt = db.dt[name %like% 'V[0-9][-]?(V[0-9])?']
   for (i in 1:nrow(db.dt)) {
      f = db.dt[i, `_links.git`]
      if (verbose) cat('* -', db.dt[i, name], '\n')
      # download.file(f, file.path(himap_database_path, basename(f)), extra=c(
      system2(curl_path(), c(
         '-s',
         '-H', '"Accept: application/vnd.github.v3.raw"',
         '-H', '"Authorization: token 99f22e14f4ed6ec6899bebe79dbf6fd7fbf9bac6"',
         '-L', f,
         '-o', file.path(himap_database_path, db.dt[i, name])
      ), stdout=F, stderr=F)
   }
   if (verbose) cat('* OK.\n')
   if (verbose) cat('* downloading reference table...')
   download.file(
      url = 'https://api.github.com/repos/taolonglab/himap/contents/inst/extdata/pcr_primers_table.txt',
      destfile = file.path(find.package('himap'), 'inst', 'extdata', 'pcr_primers_table.txt'),
      method = 'curl',
      extra = c(
         '-s',
         '-H', '"Accept: application/vnd.github.v3.raw"',
         '-H', '"Authorization: token 99f22e14f4ed6ec6899bebe79dbf6fd7fbf9bac6"'
      )
   )
   if (verbose) cat('OK.\n')
   if (verbose) cat('Done.\n')

}


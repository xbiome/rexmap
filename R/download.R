#' Download blastn and makeblastdb executables
#'
#' Downloads NCBI blastn and makeblastdb executables for the Operating System
#' running R and copies it to the HiMAP package folder.
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








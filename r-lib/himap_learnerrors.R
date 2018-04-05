is.list.of <- function(x, ctype) {
  if(!is.list(x)) return(FALSE)
  else return(all(sapply(x, is, ctype)))
}

C_isACGT <- function(seqs) {
  .Call('_dada2_C_isACGT', PACKAGE = 'dada2', seqs)
}

learnErrors <- function(fls, nreads=1e6, errorEstimationFunction = loessErrfun, multithread=FALSE, 
                        randomize=FALSE, MAX_CONSIST=10, verbose=FALSE, ...) {
  NREADS <- 0
  if(is(fls, "derep")) { fls <- list(fls) } # A single derep=class object
  drps <- vector("list", length(fls))
  if(randomize) { fls <- sample(fls) }
  for(i in seq_along(fls)) {
    if (is.list.of(fls, "derep")){
      drps[[i]] <- fls[[i]]
    } else {
      drps[[i]] <- derepFastq(fls[[i]])
    }
    NREADS <- NREADS + sum(drps[[i]]$uniques)
    if(NREADS > nreads) { break }
  }
  drps <- drps[1:i]
  # Run dada in self-consist mode on those samples
  dds <- dada(drps, err=NULL, errorEstimationFunction=errorEstimationFunction, selfConsist=TRUE, 
              multithread=multithread, verbose=verbose, MAX_CONSIST=MAX_CONSIST, ...)
  if(is.logical(verbose) || verbose > 0) cat("Total reads used: ", NREADS, "\n")
  return(getErrors(dds, detailed=TRUE))
}
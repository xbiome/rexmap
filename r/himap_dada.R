himap_dada <- function(derep,
                 err,
                 errorEstimationFunction = loessErrfun,
                 selfConsist = FALSE, 
                 pool = FALSE,
                 multithread = FALSE, 
                 verbose=FALSE, test=FALSE,...) {
  
  if (test) message('my data active.')
  call <- sys.call(1)
  # Read in default opts and then replace with any that were passed in to the function
  opts <- getDadaOpt()
  args <- list(...)
  # Catch the deprecated VERBOSE option
  if("VERBOSE" %in% names(args)) {
    stop("DEPRECATED: The VERBOSE option has been replaced by the verbose argument. Please update your code.")
    verbose <- args[["VERBOSE"]]
    args[["VERBOSE"]] <- NULL
  }
  for(opnm in names(args)) {
    if(opnm %in% names(opts)) {
      opts[[opnm]] <- args[[opnm]]
    } else {
      warning(opnm, " is not a valid DADA option.")
    }
  }
  
  # Parse verbose
  if(is.logical(verbose)) {
    if(verbose == FALSE) { verbose <- 1 }
    else { verbose <- 2 }
  }
  
  # If a single derep object, make into a length 1 list
  if(class(derep) == "derep") { derep <- list(derep) }
  if(!is.list.of(derep, "derep")) { stop("The derep argument must be a derep-class object or list of derep-class objects.") }
  if(opts$USE_QUALS && any(is.null(lapply(derep, function(x) x$quals)))) { stop("The input derep-class object(s) must include quals if USE_QUALS is TRUE.") }
  
  # Validate derep object(s)
  for(i in seq_along(derep)) {
    if(!(is.integer(derep[[i]]$uniques))) {
      stop("Invalid derep$uniques vector. Must be integer valued.")
    }
    if(!(all(C_isACGT(names(derep[[i]]$uniques))))) {
      stop("Invalid derep$uniques vector. Names must be sequences made up only of A/C/G/T.")
    }
  }
  
  # Validate quals matrix(es)
  qmax <- 0
  for(i in seq_along(derep)) {
    if(nrow(derep[[i]]$quals) != length(derep[[i]]$uniques)) {
      stop("derep$quals matrices must have one row for each derep$unique sequence.")
    }
    if(any(sapply(names(derep[[i]]$uniques), nchar) > ncol(derep[[i]]$quals))) { ###ITS
      stop("derep$quals matrices must have as many columns as the length of the derep$unique sequences.")
    }
    if(any(sapply(seq(nrow(derep[[i]]$quals)), 
                  function(row) any(is.na(derep[[i]]$quals[row,1:nchar(names(derep[[i]]$uniques)[[row]])]))))) { ###ITS
      stop("NAs in derep$quals matrix. Check that all input sequences had valid associated qualities assigned.")
    }
    if(min(derep[[i]]$quals, na.rm=TRUE) < 0) {
      stop("Invalid derep$quals matrix. Quality values must be positive integers.")
    }
    qmax <- max(qmax, max(derep[[i]]$quals, na.rm=TRUE))
  }
  
  qmax <- ceiling(qmax) # Only getting averages from derep$quals
  # if(qmax > 45) {
  #   if(qmax > 62) {
  #     stop("derep$quals matrix has an invalid maximum Phred Quality Scores of ", qmax) 
  #   }
  #   warning("derep$quals matrix has Phred Quality Scores >45. For Illumina 1.8 or earlier, this is unexpected.")
  # }
  # 
  # Pool the derep objects if so indicated
  if(length(derep) <= 1) { pool <- FALSE }
  if(pool) { # Make derep a length 1 list of pooled derep object
    derep.in <- derep
    derep <- list(combineDereps2(derep))
  }
  
  # Validate err matrix
  initializeErr <- FALSE
  if(selfConsist && (missing(err) || is.null(err))) {
    err <- NULL
    initializeErr <- TRUE
  } else {
    err <- getErrors(err, enforce=TRUE)
    if(ncol(err) < qmax+1 && verbose) { # qmax = 0 if USE_QUALS = FALSE
      message("The supplied error matrix does not extend to maximum observed Quality Scores in derep (", qmax, ").
              Extending error rates by repeating the last column of the Error Matrix (column ", ncol(err), ").
              In selfConsist mode this should converge to the proper error rates, otherwise this may not be what you want.")
      for (q in seq(ncol(err), qmax)) { 
        err <- cbind(err, err[1:16, q])
        colnames(err)[q+1] <- q
      }
    }
}
  
  # Might want to check for summed transitions from NT < 1 also.
  
  # Validate errorEstimationFunction
  if(!opts$USE_QUALS) {
    if(!missing(errorEstimationFunction) && verbose) message("The errorEstimationFunction argument is ignored when USE_QUALS is FALSE.")
    errorEstimationFunction <- noqualErrfun  # NULL error function has different meaning depending on USE_QUALS
  } else {
    if(!is.function(errorEstimationFunction)) stop("Must provide a function for errorEstimationFunction.")
  }
  
  # Validate alignment parameters
  if(opts$GAP_PENALTY>0) opts$GAP_PENALTY = -opts$GAP_PENALTY
  if(is.null(opts$HOMOPOLYMER_GAP_PENALTY)) { # Set gap penalties equal
    opts$HOMOPOLYMER_GAP_PENALTY <- opts$GAP_PENALTY
  }
  if(opts$HOMOPOLYMER_GAP_PENALTY > 0) opts$HOMOPOLYMER_GAP_PENALTY = -opts$HOMOPOLYMER_GAP_PENALTY
  
  if(opts$HOMOPOLYMER_GAP_PENALTY != opts$GAP_PENALTY) { # Use homopolymer gapping
    opts$VECTORIZED_ALIGNMENT <- FALSE # No homopolymer gapping in vectorized aligner
  }
  if(opts$VECTORIZED_ALIGNMENT) {
    if(length(unique(diag(opts$SCORE)))!=1 || 
       length(unique(opts$SCORE[upper.tri(opts$SCORE) | lower.tri(opts$SCORE)]))!=1) {
      if(verbose) message("The vectorized aligner requires that the score matrix reduces to match/mismatch. Turning off vectorization.")
      opts$VECTORIZED_ALIGNMENT=FALSE
    }
    if(opts$BAND_SIZE > 0 && opts$BAND_SIZE<8) {
      if(verbose) message("The vectorized aligner is slower for very small band sizes.")
    }
    if(opts$BAND_SIZE == 0) opts$VECTORIZED_ALIGNMENT=FALSE
  }
  
  # Parse multithreading argument
  if(is.logical(multithread)) {
    if(multithread==TRUE) { RcppParallel::setThreadOptions(numThreads = "auto") }
  } else if(is.numeric(multithread)) {
    RcppParallel::setThreadOptions(numThreads = multithread)
    multithread <- TRUE
  } else {
    if(verbose) message("Invalid multithread parameter. Running as a single thread.")
    multithread <- FALSE
  }
  
  # Initialize
  cur <- NULL
  if(initializeErr) { nconsist <- 0 } else { nconsist <- 1 }
  errs <- list()
  # The main loop, run once, or repeat until err repeats if selfConsist=T
  
  repeat{
    clustering <- list()
    clusterquals <- list()
    birth_subs <- list()
    trans <- list()
    map <- list()
    #    exp <- list()
    prev <- cur
    if(nconsist > 0) errs[[nconsist]] <- err
    
    for(i in seq_along(derep)) {
      qi <- unname(t(derep[[i]]$quals)) # Need transpose so that sequences are columns
      
      if(nconsist == 1 && verbose) {
        if(pool) {
          cat(length(derep.in), "samples were pooled:", sum(derep[[i]]$uniques), "reads in", 
              length(derep[[i]]$uniques), "unique sequences.\n")
        } else {
          cat("Sample", i, "-", sum(derep[[i]]$uniques), "reads in", 
              length(derep[[i]]$uniques), "unique sequences.\n")
        }
      } else if(i==1 && verbose) {
        if(nconsist == 0) {
          cat("Initializing error rates to maximum possible estimate.\n")
        } else {
          cat("   selfConsist step", nconsist, "\n")
        }
      }
      # Initialize error matrix if necessary
      if(initializeErr) {
        err <- matrix(1, nrow=16, ncol=max(41,qmax+1))
      }
      res <- dada_uniques(names(derep[[i]]$uniques), unname(derep[[i]]$uniques), 
                          err,
                          qi, 
                          opts[["SCORE_MATRIX"]], opts[["GAP_PENALTY"]],
                          opts[["USE_KMERS"]], opts[["KDIST_CUTOFF"]],
                          opts[["BAND_SIZE"]],
                          opts[["OMEGA_A"]], 
                          if(initializeErr) { 1 } else { opts[["MAX_CLUST"]] }, ###!
                          opts[["MIN_FOLD"]], opts[["MIN_HAMMING"]], opts[["MIN_ABUNDANCE"]],
                          TRUE, #opts[["USE_QUALS"]],
                          FALSE,
                          opts[["VECTORIZED_ALIGNMENT"]],
                          opts[["HOMOPOLYMER_GAP_PENALTY"]],
                          multithread,
                          (verbose>=2),
                          opts[["SSE"]])
      
      # Augment the returns
      res$clustering$sequence <- as.character(res$clustering$sequence)
      
      # List the returns
      clustering[[i]] <- res$clustering
      clusterquals[[i]] <- t(res$clusterquals) # make sequences rows and positions columns
      birth_subs[[i]] <- res$birth_subs
      trans[[i]] <- res$subqual
      map[[i]] <- res$map
      #      exp[[i]] <- res$exp
      rownames(trans[[i]]) <- c("A2A", "A2C", "A2G", "A2T", "C2A", "C2C", "C2G", "C2T", "G2A", "G2C", "G2G", "G2T", "T2A", "T2C", "T2G", "T2T")
      colnames(trans[[i]]) <- seq(0, ncol(trans[[i]])-1)  # Assumes C sides is returning one col for each integer starting at 0
    }
    # Accumulate the sub matrix
    cur <- Reduce("+", trans) # The only thing that changes is err(trans), so this is sufficient
    
    # Estimate the new error model (if applicable)
    if(is.null(errorEstimationFunction)) {
      err <- NULL
    } else {
      err <- tryCatch(suppressWarnings(errorEstimationFunction(cur)),
                      error = function(cond) {
                        if(verbose) message("Error rates could not be estimated.")
                        return(NULL)
                      })
    }
    if(initializeErr) {
      initializeErr <- FALSE
      err[c(1,6,11,16),] <- 1.0 # Set self-transitions (A2A, C2C, G2G, T2T) to max of 1
    }
    
    if(selfConsist) { # Validate err matrix
      if(!is.numeric(err)) stop("Error matrix returned by errorEstimationFunction not numeric.")
      if(!(nrow(err)==16)) stop("Error matrix returned by errorEstimationFunction does not have 16 rows.")
      if(!all(err>=0)) stop("Error matrix returned by errorEstimationFunction has entries <0.")
      if(!all(err<=1)) stop("Error matrix returned by errorEstimationFunction has entries >1.")
      if(any(err==0)) warning("Error matrix returned by errorEstimationFunction has 0s in some entries.")      
    }
    
    # Termination condition for selfConsist loop
    if((!selfConsist) || any(sapply(errs, identical, err)) || (nconsist >= opts$MAX_CONSIST)) {
      break
    } 
    nconsist <- nconsist+1
  } # repeat
  
  if(selfConsist && verbose) {
    if(nconsist >= opts$MAX_CONSIST) {
      message("Self-consistency loop terminated before convergence.")
    } else {
      cat("Convergence after ", nconsist, " rounds.\n")
    }
  }
  
  # Construct return object
  # A single dada-class object if one derep object provided.
  # A list of dada-class objects if multiple derep objects provided.
  rval2 = replicate(length(derep), list(denoised=NULL, clustering=NULL, sequence=NULL, quality=NULL, birth_subs=NULL, trans=NULL, map=NULL,
                                        err_in=NULL, err_out=NULL, opts=NULL, call=NULL), simplify=FALSE)
  for(i in seq_along(derep)) {
    rval2[[i]]$denoised <- getUniques(clustering[[i]])
    rval2[[i]]$clustering <- clustering[[i]]
    rval2[[i]]$sequence <- names(rval2[[i]]$denoised)
    rval2[[i]]$quality <- clusterquals[[i]]
    rval2[[i]]$birth_subs <- birth_subs[[i]]
    rval2[[i]]$trans <- trans[[i]]
    rval2[[i]]$map <- map[[i]]
    #    rval2[[i]]$exp <- exp[[i]]
    # Return the error rate(s) used as well as the final estimated error matrix
    if(selfConsist) { # Did a self-consist loop
      rval2[[i]]$err_in <- errs
    } else {
      rval2[[i]]$err_in <- errs[[1]]
    }
    rval2[[i]]$err_out <- err
    
    # Store the call and the options that were used in the return object
    rval2[[i]]$opts <- opts
    rval2[[i]]$call <- call
  }
  
  # If pool=TRUE, expand the rval and prune the individual return objects
  if(pool) {
    # Expand rval into a list of the proper length
    rval1 <- rval2[[1]]
    rval2 = replicate(length(derep.in), list(denoised=NULL, clustering=NULL, sequence=NULL, quality=NULL, birth_subs=NULL, trans=NULL, map=NULL,
                                             err_in=NULL, err_out=NULL, opts=NULL, call=NULL), simplify=FALSE)
    # Make map named by the pooled unique sequence
    map <- map[[1]]
    names(map) <- names(derep[[1]]$uniques)
    for(i in seq_along(derep.in)) {
      rval2[[i]] <- rval1
      # Identify which output clusters to keep
      keep <- unique(map[names(derep[[1]]$uniques) %in% names(derep.in[[i]]$uniques)])
      keep <- seq(length(rval2[[i]]$denoised)) %in% keep # -> logical
      newBi <- cumsum(keep) # maps pooled cluster index to individual index
      # Prune $denoised, $clustering, $sequence, $quality
      rval2[[i]]$denoised <- rval2[[i]]$denoised[keep]
      rval2[[i]]$clustering <- rval2[[i]]$clustering[keep,] # Leaves old (char of integer) rownames!
      rownames(rval2[[i]]$clustering) <- as.character(newBi[as.integer(rownames(rval2[[i]]$clustering))])
      rval2[[i]]$sequence <- rval2[[i]]$sequence[keep]
      rval2[[i]]$quality <- rval2[[i]]$quality[keep,,drop=FALSE] # Not the qualities for this sample alone!
      # Prune birth_subs and remap its $clust column
      rval2[[i]]$birth_subs <- rval2[[i]]$birth_subs[keep[rval2[[i]]$birth_subs$clust],,drop=FALSE]
      rval2[[i]]$birth_subs$clust <- newBi[rval2[[i]]$birth_subs$clust]      
      # Remap $map
      rval2[[i]]$map <- newBi[map[names(derep.in[[i]]$uniques)]]
      # Recalculate abundances (both $denoised and $clustering$abundance)
      rval2[[i]]$denoised[] <- tapply(derep.in[[i]]$uniques, rval2[[i]]$map, sum)
      rval2[[i]]$clustering$abundance <- rval2[[i]]$denoised
    }
    derep <- derep.in
    rm(derep.in)
  }
  
  names(rval2) <- names(derep)
  if(length(rval2) == 1) {  # Unlist if just a single derep object provided
    rval2 <- rval2[[1]]
    rval2 <- as(rval2, "dada")
  } else {
    for(i in seq_along(rval2)) {
      rval2[[i]] <- as(rval2[[i]], "dada")
    }
  }
  
  return(rval2)
  }

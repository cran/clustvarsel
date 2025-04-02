#############################################################################
## Sequential & parallel backward greedy search
#############################################################################
   
clvarselgrbkw <- function(X, G = 1:9,
                          emModels1 = c("E","V"), 
                          emModels2 = mclust.options("emModelNames"),
                          samp = FALSE, sampsize = 2000, 
                          hcModel = "VVV", allow.EEE = TRUE, forcetwo = TRUE, 
                          BIC.diff = 0, itermax = 100,
                          parallel = FALSE, 
                          verbose = interactive())
{

  X <- as.matrix(X)
  n <- nrow(X) # number of rows = number of observations
  d <- ncol(X) # number of columns = number of variables
  if(is.null(colnames(X))) 
    colnames(X) <- paste("X", 1:d, sep = "")
  G <- setdiff(G, 1)
  
  # If needed, sample the subset of observations for hierarchical clustering
  if(samp) { sub <- sample(1:n, min(sampsize,n), replace = FALSE) }
  else     { sub <- seq.int(1,n) }

  hcModel1 <- if(any(grep("V", hcModel))) "V" else "E"

  # info records the proposed variable, BIC for clustering using the S matrix
  # and BICdiff for BIC difference of clustering versus no clustering on S
  info <- as.data.frame(matrix(as.double(NA), nrow = 0, ncol = 7))
  #info <- data.frame(Var = NULL, BIC = NULL, BICdiff = NULL, 
  #                   Step = NULL, Decision = NULL,
  #                   stringsAsFactors = FALSE)
  S <- X
  NS <- matrix(as.double(NA), n, 0)

  mod <- NULL
  try(mod <- Mclust(S, G = G, modelNames = emModels2,
                    initialization = list(hcPairs = hc(hcModel, data = S[sub,],
                                                       use = mclust.options("hcUse")), 
                                          subset = sub),
                    verbose = FALSE),
      silent = TRUE)
  # If we get all NA's from above starting hierarchical values use "EEE"
  if((allow.EEE) & (sum(is.finite(mod$BIC))==0))
    { try(mod <- Mclust(S, G = G, modelNames = emModels2,
                        initialization = list(hcPairs = hc("EEE", data = S[sub,],
                                                           use = mclust.options("hcUse")),
                                              subset = sub),
                        verbose = FALSE),
          silent = TRUE)
    }
  # BICS is the BIC for the clustering model with all variables in S
  BICS <- if(sum(is.finite(mod$BIC))==0) NA 
          else max(mod$BIC[is.finite(mod$BIC)])

  # Start "parallel backend" if needed
  if(is.logical(parallel))
    { if(parallel) 
        { parallel <- startParallel(parallel)
          stopCluster <- TRUE }
      else
      { parallel <- stopCluster <- FALSE } 
    }
  else
    { stopCluster <- if(inherits(parallel, "cluster")) FALSE else TRUE
      parallel <- startParallel(parallel) 
    }
  on.exit(if(parallel & stopCluster)
          parallel::stopCluster(attr(parallel, "cluster")) )

  # define operator to use depending on parallel being TRUE or FALSE
  `%DO%` <- if(parallel) `%dopar%` else `%do%`
  i <- NULL # dummy to trick R CMD check 
  
  criterion <- 1
  iter <- 0
  
  while((criterion == 1) & (iter < itermax) & (ncol(S) > 1))
  {
   iter <- iter+1
   check1 <- colnames(S)
   
   if(verbose) cat(paste("iter", iter, "\n"))
   # removing step
   if(verbose) cat("- removing step\n")
   out <- foreach(i = 1:ncol(S)) %DO% 
   {
     # Calculate the BIC for the regression of the proposed variable on 
     # the other variable(s) in S
     BICreg <- BICreg(y = S[,i], x = S[,-i,drop=FALSE])        
     # Fit the cluster model on the S/{i} variables for 2 to G groups 
     mod <- NULL
     hcPairs <- hc(modelName = if(ncol(S) > 2) hcModel else hcModel1,
                   data = S[,-i,drop=FALSE][sub,],
                   use = ifelse(ncol(S[,-i,drop=FALSE]) > 1,
                                mclust.options("hcUse"), "VARS"))
     try(mod <- Mclust(S[,-i,drop=FALSE], G = G, 
                       modelNames = if(ncol(S) > 2) emModels2
                                               else emModels1,
                       initialization = list(hcPairs = hcPairs, 
                                             subset = sub),
                       verbose = FALSE),
         silent = TRUE)
     # If we get all NA's from above starting hierarchical values use "EEE"
     if((allow.EEE) & (sum(is.finite(mod$BIC))==0))
       { hcPairs <- hc(if(ncol(S) > 2) "EEE" else "E",
                       data = S[,-i,drop=FALSE][sub,])
         try(mod <- Mclust(S[,-i,drop=FALSE], G = G, 
                           modelNames = if(ncol(S) > 2) emModels2
                                                   else emModels1,
                           initialization = list(hcPairs = hcPairs, 
                                                 subset = sub),
                           verbose = FALSE),
             silent = TRUE) 
       }
     # BIC for the no clustering model on S/{j} 
     BICnotclust <- if(sum(is.finite(mod$BIC))==0) NA
                    else (max(mod$BIC[is.finite(mod$BIC)]) + BICreg)
     #
     return(list(BICreg, BICnotclust, mod$modelName, mod$G))
   }
   BICreg <- sapply(out, "[[", 1)
   BICnotclust <- sapply(out, "[[", 2)
   # cdiff is the difference between BIC for clustering on S versus {j} 
   # being conditionally independent of the clustering
   cdiff <- BICS - BICnotclust
   # Choose the variable with the smallest difference 
   m <- min(cdiff[is.finite(cdiff)])
   arg <- which(cdiff == m, arr.ind=TRUE)[1]
   if(cdiff[arg] < BIC.diff)
     { # if this difference is smaller than cut-off remove this variable 
       # from S and update the clustering model BICS
       BICS <- BICS - BICreg[arg] - cdiff[arg]
       nks <- c(colnames(NS), colnames(S)[arg])
       k <- colnames(S)[-arg]
       if(nrow(info) > 0)
         info <- rbind(info,
                       c(colnames(S)[arg], BICS, cdiff[arg], 
                       "Remove", "Accepted", out[[arg]][3:4])) 
       else
         info <- data.frame(Var = colnames(S)[arg], 
                            BIC = BICS, BICdiff = cdiff[arg], 
                            Step = "Remove", Decision = "Accepted",
                            Model = out[[arg]][[3]],
                            G = out[[arg]][[4]],
                            stringsAsFactors = FALSE)
       NS <- cbind(NS, S[,arg])
       S <- S[,-arg,drop=FALSE]
       colnames(NS) <- nks
       colnames(S) <- k
     } 
   else
     { if(nrow(info) > 0)
         info <- rbind(info, 
                       c(colnames(S)[arg], 
                         (BICnotclust-BICreg)[arg], cdiff[arg],
                         "Remove", "Rejected", out[[arg]][3:4]))
       else
         info <- data.frame(Var = colnames(S)[arg], 
                            BIC = BICS, BICdiff = cdiff[arg], 
                            Step = "Remove", Decision = "Rejected",
                            Model = out[[arg]][[3]],
                            G = out[[arg]][[4]],
                            stringsAsFactors = FALSE)
     }

   if(ncol(NS) > 2)
     { # adding step 
       if(verbose) cat("+ adding step\n")
       out <- foreach(i = 1:ncol(NS)) %DO% 
       {
         # Fit clustering model with proposed variable
         hcPairs <- hc(modelName = if(ncol(S) > 0) hcModel else hcModel1,
                       data = cbind(S,NS[,i])[sub,],
                       use = mclust.options("hcUse"))
         try(mod <- Mclust(cbind(S,NS[,i]), G = G, 
                           modelNames = if(ncol(S) > 0) emModels2 
                                        else            emModels1,
                           initialization = list(hcPairs = hcPairs,
                                                 subset = sub),
                           verbose = FALSE),
             silent = TRUE)
         # If we get all NA's from above starting hierarchical values use "EEE"
         if((allow.EEE) & (sum(is.finite(mod$BIC))==0))
           { 
             hcPairs <- hc(if(ncol(S) > 0) "EEE" else "E",
                           data = cbind(S,NS[,i])[sub,])
             try(mod <- Mclust(cbind(S,NS[,i]), G = G, 
                               modelNames = if(ncol(S) > 0) "EEE" else "E",
                               initialization = list(hcPairs = hcPairs,
                                                     subset = sub),
                               verbose = FALSE),
                silent = TRUE)
           }
         # BICS is the BIC for the clustering model with all variables in S
         BICNS <- if(sum(is.finite(mod$BIC))==0) NA 
                  else max(mod$BIC[is.finite(mod$BIC)])
         # Calculate the BIC for the regression of the proposed variable 
         # on the other variable(s) in S
         BICreg <- BICreg(y = NS[,i], x = S)
         #         
         return(list(BICreg, BICNS, mod$modelName, mod$G))
       }
       BICreg <- sapply(out, "[[", 1)
       BICNS <- sapply(out, "[[", 2)
       # cdiff is the difference between BIC for clustering on S u NS versus
       # being conditionally independent of the clustering
       cdiff <- BICNS - (BICS + BICreg)
       # Choose the variable with the largest difference 
       m <- max(cdiff[is.finite(cdiff)])
       arg <- which(cdiff==m,arr.ind=TRUE)[1]
       if(cdiff[arg] > BIC.diff)
         { # if this difference is larger than cut-off add this variable 
           # to S and update the clustering model BICS
           BICS <- BICNS[arg]
           k <- c(colnames(S), colnames(NS)[arg])
           nks <- colnames(NS)[-arg]
           info <- rbind(info, c(colnames(NS)[arg], BICS, cdiff[arg], 
                                 "Add", "Accepted", out[[arg]][3:4]))
           S <- cbind(S,NS[,arg])
           NS <- NS[,-arg,drop=FALSE]
           colnames(S) <- k
           colnames(NS) <- nks
         } 
       else
         { info <- rbind(info, c(colnames(NS)[arg], BICNS[arg], cdiff[arg],
                                 "Add", "Rejected", out[[arg]][3:4]))
         }
   }
   
   info$BIC <- as.numeric(info$BIC)
   info$BICdiff <- as.numeric(info$BICdiff)

   if(verbose) 
     { if(iter > 2) print(info[seq(nrow(info)-1,nrow(info)),c(1,3:5)])
       else         print(info[nrow(info),c(1,3:5),drop=FALSE]) }
   
   # Check if the variables in S have changed or not
   check2 <- colnames(S)
   if(is.null(check2))
     # all variables have been removed
     { criterion <- 0 }
   else
     # if they have changed (either added one or removed one or changed one)
     # then continue the algorithm (criterion is 1) otherwise stop 
     # (criterion is 0)
     { if(length(check2) != length(check1))
         { criterion <- 1 }
       else
         { criterion <- if(sum(check1==check2) != length(check1)) 1 else 0 }
     }
        
  }

  if(iter >= itermax) 
    warning("Algorithm stopped because maximum number of iterations was reached")
  
  # List the selected variables and the matrix of steps' information
  info$BIC <- as.numeric(info$BIC)
  info$BICdiff <- as.numeric(info$BICdiff)
  # reorder steps.info
  info <- info[,c(1,4,2,6,7,3,5),drop=FALSE]
  colnames(info) <- c("Variable proposed", "Type of step",
                      "BICclust", "Model", "G", "BICdiff", "Decision")
  varnames <- colnames(X)
  subset <- if(is.null(S)) NULL else
               sapply(colnames(S), function(x) which(x == varnames))

  out <- list(variables = varnames,
              subset = subset,
              steps.info = info,
              search = "greedy",
              direction = "backward")
  class(out) <- "clustvarsel"  
  return(out)
}

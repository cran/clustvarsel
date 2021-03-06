clustvarsel <- function(data, G = 1:9, 
                        search = c("greedy", "headlong"),
                        direction = c("forward", "backward"),
                        emModels1 = c("E", "V"), 
                        emModels2 = mclust.options("emModelNames"),
                        samp = FALSE, 
                        sampsize = round(nrow(data)/2), 
                        hcModel = "VVV", 
                        allow.EEE = TRUE, 
                        forcetwo = TRUE,
                        BIC.diff = 0, 
                        BIC.upper = 0, 
                        BIC.lower = -10, 
                        itermax = 100, 
                        parallel = FALSE,
                        fit = TRUE,
                        verbose = interactive())
{

  call <- match.call()
  mc <- match.call(expand.dots = TRUE)
  
  # supress warning temporarily
  warn <- getOption("warn")
  on.exit(options("warn" = warn))
  options("warn" = -1)

  X <- data.matrix(data)
  # Check whether there are variable names to identify selected variables
  if(is.null(colnames(X))) 
    colnames(X) <- paste("X", 1:ncol(X), sep = "")
  mc$data <- NULL
  
  # search    <- match.arg(search, choices = eval(formals(clustvarsel)$search))
  # direction <- match.arg(direction, choices = eval(formals(clustvarsel)$direction))
  search    <- match.arg(search, choices = c("greedy", "headlong"))
  direction <- match.arg(direction, choices = c("forward", "backward"))
  
  if(search == "greedy" & direction == "forward")
    { mc[[1]] <- as.name("clvarselgrfwd")
      mc$X <- X
      mc[c("search", "direction", "BIC.upper", "BIC.lower", "fit")] <- NULL
      out <- eval(mc, parent.frame())
    }
  else if(search == "greedy" & direction == "backward")
    { mc[[1]] <- as.name("clvarselgrbkw")
      mc$X <- X
      mc[c("search", "direction", "forcetwo", "BIC.upper", "BIC.lower", "fit")] <- NULL
      out <- eval(mc, parent.frame())
    }
  else if(search == "headlong" & direction == "forward")
    { mc[[1]] <- as.name("clvarselhlfwd")
      mc$X <- X
      mc[c("search", "direction", "BIC.diff", "parallel", "fit")] <- NULL
      out <- eval(mc, parent.frame())
    }   
  else if(search == "headlong" & direction == "backward")
    { mc[[1]] <- as.name("clvarselhlbkw")
      mc$X <- X
            mc[c("search", "direction",  "forcetwo", 
                 "BIC.diff", "parallel", "fit")]  <- NULL
      out <- eval(mc, parent.frame())
    }
  else stop("selected search and/or direction not available")
  
  if(!is.null(out) & fit)
  { 
    if(verbose)
      cat(paste("final iter\n* fitting model on selected subset\n"))
    out$model <- Mclust(X[,out$subset,drop=FALSE], 
                        G = G, 
                        modelNames = if(length(out$subset) > 1) emModels2 
                                     else emModels1,
                        verbose = FALSE)
  }
  class(out) <- "clustvarsel"
  return(out)
}

print.clustvarsel <- function(x, digits = getOption("digits"), ...) 
{
  # cat("'", class(x)[1], "' model object:\n\n", sep = "")
  header1 <- c("Variable selection for Gaussian model-based clustering")
  dir <- switch(x$direction,
                "forward" = "(forward/backward)",
                "backward" = "(backward/forward)",
                NULL)
  header2 <- switch(x$search,
                    "greedy"   = paste("Stepwise", dir, "greedy search"),
                    "headlong" = paste("Headlong", dir, "search"),
                    NULL)
  sep <- paste(rep("-", max(nchar(header1),nchar(header2))), collapse ="") 
  cat(sep, "\n")
  cat(header1, sep = "\n")
  cat(header2, sep = "\n")
  cat(sep, "\n\n")
  #
  print(x$steps.info, na.print = "", row.names = FALSE, digits = digits, ...)
  cat("\n")
  # cat("\nSelected subset:",
  #    paste(names(x$subset), collapse = ", "), fill = TRUE)
  footer <- strwrap(paste("Selected subset:", 
                          paste(names(x$subset), collapse = ", ")),
                    width = getOption("width"), simplify = TRUE, 
                    indent = 0, exdent = 2)
  cat(footer, sep = "\n")
  invisible()
}

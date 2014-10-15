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
                        parallel = FALSE)
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
  search <- match.arg(search)
  mc$data <- NULL
  
  if(search == "greedy" & direction == "forward")
    { mc[[1]] <- as.name("clvarselgrfwd")
      mc$X <- X
      mc$search <- mc$direction <- mc$BIC.upper <- mc$BIC.lower <- NULL
      out <- eval(mc, parent.frame())
    }
  else if(search == "greedy" & direction == "backward")
    { mc[[1]] <- as.name("clvarselgrbkw")
      mc$X <- X
      mc$search <- mc$direction <- mc$forcetwo <- mc$BIC.upper <- mc$BIC.lower <- NULL
      out <- eval(mc, parent.frame())
    }
  else if(search == "headlong" & direction == "forward")
    { mc[[1]] <- as.name("clvarselhlfwd")
      mc$X <- X
      mc$search <- mc$direction <- mc$BIC.diff <- mc$parallel <- NULL
      out <- eval(mc, parent.frame())
    }   
  else if(search == "headlong" & direction == "backward")
    { mc[[1]] <- as.name("clvarselhlbkw")
      mc$X <- X
      mc$search <- mc$direction <- mc$forcetwo <- mc$BIC.diff <- mc$parallel <- NULL
      out <- eval(mc, parent.frame())
    }
  else stop("selected search and/or direction not available")
  
  if(!is.null(out)) 
    class(out) <- "clustvarsel"  
  return(out)
}

print.clustvarsel <- function(x, digits = getOption("digits"), ...) 
{
  cat("'", class(x)[1], "' model object:\n\n", sep = "")
  switch(x$search,
         "greedy"   = cat("Stepwise (", x$direction, ") greedy search:\n", sep=""),
         "headlong" = cat("Headlong (", x$direction, ") search:\n", sep="")
         )
  print(x$steps.info, na.print = "", print.gap = 2, digits = digits, ...)
  cat("\nSelected subset:",
      paste(names(x$subset), collapse = ", "), fill = TRUE)
  invisible()
}

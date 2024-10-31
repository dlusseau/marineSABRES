

networkfromDf <- function (file, bodySeparator = ",", lowercaseGenes = FALSE, 
          symbolic = FALSE) 
{
  
  # require(glue) For trim function
  
  func <- readLines(file, -1)
  func <- gsub("#.*", "", glue::trim(func))
  func <- func[nchar(func) > 0]
  if (length(func) == 0) 
    stop("Header expected!")
  header <- func[1]
  header <- tolower(trim(strsplit(header, bodySeparator)[[1]]))
  if (length(header) < 2 || header[1] != "targets" || !(header[2] %in% 
                                                        c("functions", "factors")) || (length(header) == 3 && 
                                                                                       header[3] != "probabilities")) 
    stop(paste("Invalid header:", func[1]))
  func <- func[-1]
  if (lowercaseGenes) 
    func <- tolower(func)
  func <- gsub("[^\\[\\]a-zA-Z0-9_\\|\\&!\\(\\) \t\\-+=.,]+", 
               "_", func, perl = TRUE)
  
  
  
  
  tmp <- unname(lapply(func, function(x) {
    bracketCount <- 0
    lastIdx <- 1
    chars <- strsplit(x, split = "")[[1]]
    res <- c()
    if (length(chars) > 0) {
      for (i in seq_along(chars)) {
        if (chars[i] == "(") 
          bracketCount <- bracketCount + 1
        else if (chars[i] == ")") 
          bracketCount <- bracketCount - 1
        else if (chars[i] == bodySeparator && bracketCount == 
                 0) {
          res <- c(res, glue::trim(paste(chars[lastIdx:(i - 
                                                    1)], collapse = "")))
          lastIdx <- i + 1
        }
      }
      res <- c(res, glue::trim(paste(chars[lastIdx:length(chars)], 
                               collapse = "")))
    }
    return(res)
  }))
  
  
  targets <- sapply(tmp, function(rule) glue::trim(rule[1]))
  targets1 <- as.character(boolean.df$targets)
  

  factors <- sapply(tmp, function(rule) trim(rule[2]))
  
  factors1 <- as.character(boolean.df$factors)
  

  probabilities <- sapply(tmp, function(rule) {
    if (length(rule) >= 3) 
      as.numeric(trim(rule[3]))
    else 1
  })
  
  
  # factors.tmp <- lapply(factors1, matchNames)
  # genes <- unique(c(targets, unname(unlist(factors.tmp))))

    # This should be just a list of all nodes, which we here assume equals all targets  
  genes <- targets1
  
  isProbabilistic <- (length(unique(targets1)) < length(targets1))
  
  # It is not symbolic, so this if statement can go
  
  # if (symbolic) {
  #   if (isProbabilistic) 
  #     stop("Probabilistic networks cannot be loaded with symbolic=TRUE!")
  #   interactions <- lapply(factors, function(rule) parseBooleanFunction(rule, 
  #                                                                       genes))
  #   onlyInputs <- setdiff(genes, targets)
  #   if (length(onlyInputs) > 0) {
  #     for (gene in onlyInputs) {
  #       warning(paste("There is no transition function for gene \"", 
  #                     gene, "\"! Assuming an input!", sep = ""))
  #       interactions[[gene]] = parseBooleanFunction(gene, 
  #                                                   genes)
  #     }
  #   }
  #   delays <- apply(sapply(interactions, maxTimeDelay, genes = genes), 
  #                   1, max)
  #   names(interactions) <- genes
  #   fixed <- as.integer(sapply(interactions, function(int) {
  #     if (int$type == "constant") int$value else -1L
  #   }))
  #   names(fixed) <- genes
  #   res <- list(genes = genes, interactions = interactions, 
  #               fixed = fixed)
  #   res$internalStructs <- .Call("constructNetworkTrees_R", 
  #                                res)
  #   res$timeDelays <- delays
  #   class(res) <- "SymbolicBooleanNetwork"
  #   return(res)
  # }
  
  # Generate an interaction by parsing <expressionString>
  # and building the corresponding truth table.
  # Here, <genes> is a list of all genes in the network.
  # Returns an interaction as used in the BooleanNetwork class.
  generateInteraction <- function(expressionString, genes)
  {
    res <- .Call("getTruthTable_R", parseBooleanFunction(expressionString, genes), length(genes))
    names(res) <- c("input","func")
    res$expression <- expressionString
    return(res)
  }
  
  else {
    fixed <- rep(-1, length(genes))
    names(fixed) <- genes
    interactions <- list()
    for (i in seq_along(targets)) {
      target <- targets[i]
      interaction <- generateInteraction(factors[i], genes)
      if (isProbabilistic) {
        interaction$probability <- probabilities[i]
        interactions[[target]][[length(interactions[[target]]) + 
                                  1]] <- interaction
      }
      else {
        if (length(interaction$func) == 1) {
          fixed[target] <- interaction$func
        }
        interactions[[target]] <- interaction
      }
    }
    onlyInputs <- setdiff(genes, targets)
    if (length(onlyInputs) > 0) {
      for (gene in onlyInputs) {
        warning(paste("There is no transition function for gene \"", 
                      gene, "\"! Assuming an input!", sep = ""))
        if (isProbabilistic) 
          interactions[[gene]][[1]] = list(input = length(interactions) + 
                                             1, func = c(0, 1), expression = gene, probability = 1)
        else interactions[[gene]] = list(input = length(interactions) + 
                                           1, func = c(0, 1), expression = gene)
      }
    }
    if (isProbabilistic) {
      wrongProb <- sapply(interactions, function(interaction) {
        abs(1 - sum(sapply(interaction, function(func) func$probability))) > 
          1e-04
      })
      if (any(wrongProb)) 
        stop(paste("The probabilities of gene(s) ", 
                   paste(genes[wrongProb], collapse = ", "), 
                   " do not sum up to 1!", sep = ""))
    }
    res <- list(interactions = interactions, genes = genes, 
                fixed = fixed)
    if (isProbabilistic) 
      class(res) <- "ProbabilisticBooleanNetwork"
    else class(res) <- "BooleanNetwork"
    return(res)
  }
}

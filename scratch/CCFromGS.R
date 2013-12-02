library(flowWorkspace)

setwd("/shared/silo_researcher/Gottardo_R/mike_working/HVTN/086/merged")
gs <- load_gslist("1287")
parent <- "4+"

COMPASSContainerFromGatingSet <- function(gs, parent, markers, meta, individual_id, sample_id, fixed=TRUE) {
  
  if (!require(flowWorkspace)) {
    stop("This function depends on the suggested package \"flowWorkspace\".")
  }
  
  if (!(inherits(gs, "GatingSet") || inherits(gs, "GatingSetList"))) {
    stop("'gs' must be an object inheritting from 'GatingSet' or 'GatingSetList'")
  }
  
  ## some guessing of the parent path
  paths <- getNodes( gs[[1]], isPath=TRUE )
  if (!(parent %in% paths)) {
    pr <- parent ## parent_regex
    n <- nchar(parent)
    if (substring(pr, 1, 1) != "/")
      pr <- paste0("/", pr)
    if (substring(pr, n, n) != "$")
      pr <- paste0(pr, "$")
    pr <- gsub("+", "\\+", pr, fixed=TRUE)
    path <- grep(pr, paths, value=TRUE)
    if (length(path) != 1) {
      stop("Could not uniquely match 'parent' to a path in the gating set.")
    }
  } else {
    path <- parent
  }
  
  message("Taking node path as '", path, "'.")
  
  ## Initialize memory to store the intensities and counts
  n <- length(gs)
  counts <- vector("list", n)
  sn <- sampleNames(gs)
  intensities <- vector("list", n)
  names(counts) <- sn
  names(intensities) <- sn
  
  for (i in seq_along(gs)) {
    
    ff <- getData(gs[[i]], path)
    
    ## Get the cell counts for the parent
    ind <- flowWorkspace:::.getNodeInd(gs[[i]], path)
    stats <- flowWorkspace:::.getPopStat(gs[[i]], ind)
    counts[[i]] <- stats$flowCore["count"]
    
    children <- getChildren( gs[[i]], path, isPath=TRUE )
    if (missing(markers)) {
      markers <- gsub(".*/|\\+|-", "", children)
    }
    
    ## get the channel name associated with each marker
    pm <- parameters(ff)@data
    map <- sapply(markers, function(marker) {
      
      m <- matrix( nrow=2, byrow=TRUE, c(
        adist(marker, pm$name),
        adist(marker, pm$desc)
      ))
      rownames(m) <- c("name", "desc")
      
      return (pm$desc[which.min(apply(m, 2, min))])
    })
    names(map) <- gsub(".*/", "", children)
    
    ## construct a boolean expression
    ## should be a function of the channel names (?)
    expr <- as.name( paste(map, collapse="|") )
    
    map <- as.list(map)
    
    res <- getData( gs[[i]], parent=path, expr, pop_marker_list=map )
    ## the parent-level data
    res <- exprs(getData( gs[[i]], path ))
    
    res_child <- lapply(children, function(child) {
      exprs(getData( gs[[i]], child ))
    })
    names(res_child) <- children
    
    ## the idea: we have the cell intensities for the parent population,
    ## and we have the cell intensities for cells that made it through a particular
    ## 'marker' gate. all the cells in the parent population that
    ## did not make it through the 'marker' gate should get 0 for that marker
    
    for (j in seq_along(children)) {
      res[, bgi[j]][ !(rownames(res) %in% rownames(res_child[[j]])) ] <- 0
    }
    
    ## keep only the elements where at least one is positive
    res <- res[, bgi, drop=FALSE]
    res <- res[ apply(res, 1, sum) != 0, , drop=FALSE]
    intensities[[i]] <- res
    
  }
  
  ## try to get the metadata, if available
  if (missing(meta)) {
    tryCatch(meta <- pData(gs),
      error=function(e) {
        stop("No metadata available!")
      })
  }
  
  CC <- COMPASSContainer(
    data=intensities,
    counts=counts,
    meta=meta,
    individual_id=individual_id,
    sample_id=sample_id
  )
  
  return(CC)
  
}
# list files that fit a pattern,
# store the contente within a list,
# handle EMPTY FILES !
allFilters <- list.files(paste0(PATH_ROOT, "/03_qc"), pattern = DATASET, full.names = TRUE)
removedCells <- lapply(allFilters, function(file) {
    x <- try(read.table(file, head = FALSE), silent = TRUE)
    if(inherits(x, "try-error"))
        return(NULL)
    else
        return(x)
})
names(removedCells) <- basename(allFilters)

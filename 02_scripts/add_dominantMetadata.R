add_dominantMetadata <- function(SO, metadata = NULL){
    metaSO <- grep(paste0("^", metadata,"$"), names(SO@meta.data), value = TRUE)
    cluster_mtx <- table(SO$clusters, SO@meta.data[[metaSO]])
    
    pct_mtx <- prop.table(cluster_mtx, 1)
    max_val <- apply(pct_mtx, 1, max)
    idents_idx <- apply(pct_mtx, 1, which.max)
    max_idents <- colnames(cluster_mtx)[idents_idx]
    
    cluster_idents <- data.frame(rownames(cluster_mtx), max_idents, round(max_val*100, digits = 2))
    colnames(cluster_idents) <- c("cluster", paste0("most.expressed.", metadata), paste0("pct.", metadata))
    
    return(cluster_idents)
}
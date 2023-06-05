extract_topn <- function(markers, topn=20){

    # Rearrange the columns to be more intuitive
    # Order values to be more intuitive
    if ("cluster" %in% names(markers)){
        markers <- markers %>%
            dplyr::arrange(cluster, p_val_adj)
        # Extract top markers per cluster
        topnPos <- markers %>%
            group_by(cluster) %>%
            top_n(n = topn,
                  wt = avg_log2FC)
        
        topnMarkers <-  topnPos %>%
            dplyr::arrange(cluster, p_val_adj)

    } else {
        
        markers <- markers %>%
            dplyr::arrange(p_val_adj)
        # Extract top markers per cluster
        topnPos <- markers %>%
            top_n(n = topn,
                  wt = avg_log2FC)
        
        topnMarkers <-  topnPos %>%
            dplyr::arrange(p_val_adj)
    }
    return(topnMarkers)
}

# ##################################################
# Global declarations and libraries for the analysis
# ##################################################


######## R Libraries

# library( digest)        # digest (hashing)
# library( DT)            # datatable
# library( forcats)       # fct_inorder (factors)
# library( fs)            # path_sanitize
# library( future)        # plan (multicore)
# library( ggplot2)
# library( ggpubr)
# library( ggrepel)
# #library( grid)
# #library( gridExtra)
# library( htmltools)     # browsable
# library( htmlwidgets)   # JS (for datatable)
# library( iheatmapr)     # iheatmap
# library( kableExtra)
# library( knitr)
# #library( magrittr)
# library( pander)        # pander
# library( patchwork)     # +/ (ggplot layout)
# library( pheatmap)      # pheatmap
# library( plotly)
# #library( plyr)
# #library( dplyr)
# library( rmarkdown)
# library( scales)        # hue_pal
# 
# # Single-cell technology
# library( Seurat)
# 
# # Functional Enrichment analysis
# library( biomaRt)
# library( clusterProfiler)
# library( org.Mm.eg.db)
# library( org.Hs.eg.db)
# library( rrvgo)
# library( rstatix)

######## DEFINITION OF FUNCTIONS #########################

# ########################################################
# Convert a list of Human genes to their homolog in Mouse
# ########################################################

convertHumanGeneList <- function(x){
  
  human = useEnsembl( biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse = useEnsembl( biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  return( genesV2)
}

######## Violin/Jitter plot using plotly

# Function that creates a violin plot with jitter points
# To be used for showing #UMIs, #Genes, and %Mito
plotViolinJitter = function(plotData, xAxisFormula, yAxisFormula, colorFormula = NULL, pointSize = 3, pointsColorFormula = "#444444", pointsOpacity = 0.5, xAxisTitle = "", yAxisTitle = "", hoverText = NULL, traceName = "", xTicklabels = c(""), panelWidth = 300, panelHeight = 600, thresholdHigh = NULL, thresholdLow = NULL)
{
  # Generate a jitter that will not be recomputed at each redraw
  constantJitter = jitter(rep(0, nrow(plotData)), amount = 0.15);
  # Add it to data so it can be accessed with updated formula (markers)
  plotData = cbind(plotData, jitterCoords = constantJitter + 0.17);
  
  p = plot_ly( plotData,
               x = xAxisFormula,
               y = yAxisFormula,
               name = traceName,
               fillcolor = colorFormula, # 'color' argument does not work with custom values (?)
               height = panelHeight,
               width = panelWidth) %>%
    add_trace( type='violin',
               showlegend=FALSE,
               points = FALSE,
               line = list( width = 1.5, 
                            color = "#44444444"), # impossible to specify individual colors (?)
               meanline = list(visible = TRUE),
               scalemode = 'width',
               width = 0.6,
               side = 'negative',
               hoveron = 'kde',
               hoverinfo = 'skip') %>%
    add_trace( type = "scatter",
               mode = "markers",
               x = update(xAxisFormula, ~ . + jitterCoords), # Update formula adding jitter values
               marker = list(size = pointSize,
                             color = pointsColorFormula,
                             opacity = pointsOpacity),
               text = hoverText,
               hoverinfo = 'text+name') %>%
    layout( yaxis = list( title = yAxisTitle,
                          rangemode = "nonnegative"),
            xaxis = list( title = xAxisTitle,
                          ticktext = xTicklabels,
                          tickvals = ((1:length( xTicklabels)) + 0.01)),
            margin = list( l = 0, r = 0, t = 0, b = 0, pad = 5));
  
  
  # Function creating a horizontal line as a shape to be added to layout (use 'paper' coords ref for x axis)
  hline = function(y, color = "blue") 
  {
    if (is.null( y)) return( NULL);
    list( list( type = "line", 
                x0 = 0, 
                x1 = 1, 
                xref = "paper",
                y0 = y, 
                y1 = y, 
                line = list( color = color)))
  }
  
  # Shapes to be added to layout (for eventual threshold lines)
  shapes = c( hline( thresholdHigh, "red"), hline( thresholdLow, "blue"));
  
  if(length( shapes)) return(p %>% layout( shapes = shapes));
  return(p);
}




######## Highlight cells from a specific cluster on a dimreduc plot
highlightClusterPlot = function(clusterName, seuratObject, reduction = "umap")
{
  clusterCells = which( Idents( seuratObject) == clusterName);
  # Create a dimreduc plot with current cluster highlighted (useful for large number of clusters)
  print( DimPlot( seuratObject, 
                  reduction = reduction, 
                  cols="#44444422", 
                  cells.highlight = clusterCells, 
                  cols.highlight = "#FF000088", 
                  sizes.highlight = Seurat:::AutoPointSize(seuratObject)*1.5, 
                  order = clusterCells,  # Plot highlighted cells last
                  group.by=NULL) + 
           ggtitle( paste( "Cluster", clusterName)) +
           theme( axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  plot.title = element_text( face = "bold",
                                             size = rel( 16/14), 
                                             hjust = 0.5,
                                             vjust = 1, 
                                             margin = margin( b = 7)),
                  legend.position = "none"));

}




######## Create a violinplot showing feature (gene or metadata) values for each cluster (+ number of 'zero' and 'not zero' cells)
violinFeatureByCluster = function(currentFeature, seuratObject, slot = "counts", clustersColor = NULL, yLabel = "Counts", addStats = TRUE, addTitle = TRUE, trimTitle = FALSE) 
{
  if(length( currentFeature)>1) warning( "Several features in 'currentFeature' argument, corresponding values will be averaged...");

  # Extract feature values from Seurat object (average if several features given)
  featureValues = apply( FetchData( object = seuratObject, vars = currentFeature, slot = slot), 1, mean);

  # Plot expression values of current annotation for each cluster as violin + jitter
  ggFigure = ggplot( data.frame( Counts = featureValues, Cluster = Idents( seuratObject)), 
                                 aes( x = Cluster, y = Counts)) +
               geom_jitter( width = 0.2, height = 0.3, size = 0.4, colour = "#44444444") +
               geom_violin( aes( col = Cluster, fill = Cluster), scale = "width", alpha = 0.4, draw_quantiles = 0.5) +
               theme_minimal() +
               ylab( label = yLabel) +
               theme( legend.position = "None");

  # Change ggplot default clusters colors ('clustersColor' should be a named vector of colors)
  if(!is.null(clustersColor))
  {
    ggFigure = ggFigure + 
      scale_color_manual( values = clustersColor) +
      scale_fill_manual( values = clustersColor);
  }

  if(addStats)
  {
    # Get the number of cells in which score is zero (or not) for this annotation
    statsZero = tapply( featureValues, 
                        Idents( seuratObject), 
                        function(x) { paste( table( factor( as.logical( x), levels = c( TRUE, FALSE))), collapse="\n") });

    # Create the text layer with stats to be added
    ggFigure = ggFigure + geom_text( data = data.frame( Stat = statsZero, Cluster = names( statsZero)), 
                                     aes( x = Cluster, y = max( featureValues)*1.1, label = Stat), 
                                     size = 2.5) +
                          geom_text( data = data.frame( text = "Pos.\nNull"), 
                                     aes( label = text), 
                                     x = -0.4, y = max( featureValues)*1.1, hjust = 0, size = 2.5) +
                          coord_cartesian( clip = 'off', xlim = c( 0, length(statsZero))); # Allow drawing outside plot region
  } 

  if(addTitle)
  {
    # Remove 'trimTitle' number of ending character(s) (module names have a number suffix added automatically)
    title = if(trimTitle) substr( currentFeature, 1, nchar( currentFeature)-trimTitle) else currentFeature;

    # Create the title layer to be added (copy 'FeaturePlot' title style from cowplot)
    ggFigure = ggFigure + ggtitle( label = title) + 
                          theme( plot.title =  element_text( face = "bold",
                                                             size = rel( 16/14), 
                                                             hjust = 0.5,
                                                             vjust = 1, 
                                                             margin = margin( b = 7)));
  }

  # Render figure
  #print( ggFigure)
  suppressWarnings(print( ggFigure));
}




######## Datatable default buttons
# Configure a set of buttons for exporting data from datatables
# (NOTE: an orthogonal named 'export' must be defined in the concerned columns for proper formatting)
exportButtonsListDT = list( 
  list( # Collection containing export buttons with altered 'exportOptions' to use non-formatted/non-ordered data and save full table
    extend = "collection",
    text = "<center>Export<br>Full</center>",
    buttons = list( list( extend  = "excel",
                          exportOptions = list( modifier = list(order    = 'index', # give result in the same order as original table
                                                                page     = 'all',   # do not restrict to currently displayed page
                                                                search   = 'none'), # do not restrict to current filtered rows
                                                # do not define 'selected' here as 'undefined' gives expected behavior (selection ignored, not sure how to set to null or undefined using R)
                                                orthogonal = 'export')), # use 'export' orthogonal source (no formatting, defined in 'columnsDef')
                    list( extend  = "csv",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'all',
                                                                search   = 'none'),
                                                orthogonal = 'export')),
                    list( extend  = "copy",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'all',
                                                                search   = 'none'),
                                                orthogonal = 'export')))),
  list( # Collection containing export buttons with altered 'exportOptions' to use non-formatted/non-ordered data and save filtered table
    extend = "collection",
    text = "<center>Export<br>Filtered</center>",
    buttons = list( list( extend  = "excel",
                          exportOptions = list( modifier = list(order    = 'index',      # give result in the same order as original table
                                                                page     = 'all',        # do not restrict to currently displayed page
                                                                search   = 'applied'),   # restrict to current filtered rows
                                                # do not define 'selected' here as 'undefined' gives expected behavior (selection ignored, not sure how to set to null or undefined using R)
                                                orthogonal = 'export')),
                    list( extend  = "csv",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'all',
                                                                search   = 'applied'),
                                                orthogonal = 'export')),
                    list( extend  = "copy",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'all',
                                                                search   = 'applied'),
                                                orthogonal = 'export')))),
  list( # Collection containing export buttons with altered 'exportOptions' to use non-formatted/non-ordered data and save rows from current page only
    extend = "collection",
    text = "<center>Export<br>Page</center>",
    buttons = list( list( extend  = "excel",
                          exportOptions = list( modifier = list(order    = 'index',      # give result in the same order as original table
                                                                page     = 'current',    # restrict to currently displayed page
                                                                search   = 'none'),      # do not restrict to current filtered rows
                                                # do not define 'selected' here as 'undefined' gives expected behavior (selection ignored, not sure how to set to null or undefined using R)
                                                orthogonal = 'export')),
                    list( extend  = "csv",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'current',
                                                                search   = 'none'),
                                                orthogonal = 'export')),
                    list( extend  = "copy",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'current',
                                                                search   = 'none'),
                                                orthogonal = 'export')))),
  list( # Collection containing export buttons with altered 'exportOptions' to use non-formatted/non-ordered data and save selected rows only (independently of the filtering)
    extend = "collection",
    text = "<center>Export<br>Selected</center>",
    buttons = list( list( extend  = "excel",
                          exportOptions = list( modifier = list(order    = 'index',      # give result in the same order as original table
                                                                page     = 'all',        # do not restrict to currently displayed page
                                                                search   = 'none',       # do not restrict to current filtered rows
                                                                selected = TRUE),        # extension Select
                                                orthogonal = 'export')),
                    list( extend  = "csv",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'all',
                                                                search   = 'none',
                                                                selected = TRUE),
                                                orthogonal = 'export')),
                    list( extend  = "copy",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'all',
                                                                search   = 'none',
                                                                selected = TRUE),
                                                orthogonal = 'export')))),
  # # -- Print buttons
  list( # Collection containing print buttons with altered 'exportOptions' to render visible columns only in displayed rows order, and select desired rows (see comments in export buttons, here it uses default orthogonal that gives formatted data)
    extend = "collection",
    text = "<center>Print<br>(as seen)</center>",
    buttons = list( list( extend  = "print",
                          text = "Full",
                          exportOptions = list( modifier = list(order    = 'current',
                                                                page     = 'all',
                                                                search   = 'none'),
                                                columns = ':visible')),
                    list( extend  = "print",
                          text = "Filtered",
                          exportOptions = list( modifier = list(order    = 'current',
                                                                page     = 'all',
                                                                search   = 'applied'),
                                                columns = ':visible')),
                    list( extend  = "print",
                          text = "Page",
                          exportOptions = list( modifier = list(order    = 'current',
                                                                page     = 'current',
                                                                search   = 'none'),
                                                columns = ':visible')),
                    list( extend  = "print",
                          text = "Selected",
                          exportOptions = list( modifier = list(order    = 'current',
                                                                page     = 'all',
                                                                search   = 'none',
                                                                selected = TRUE),
                                                columns = ':visible')))),
  list( extend = "selectNone",
        text = "<center>Deselect<br>all</center>"));




######## Generate a datatable summarizing values 
# For environments (parameters), all values/variables are shown

showSimpleDT = function( dataToShow, rownames = TRUE, colnames = "Value")
{
  valuesDF = NULL;
  
  if(is.environment( dataToShow))
  {
    # Extract values from environment as character strings
    envValues = sapply(lapply(dataToShow, function(x) {ifelse(is.null(x), "NULL", x)}), paste, collapse = ", ");
    # Sort them by name and convert to data.frame
    valuesDF = data.frame("Value" = envValues[order(names(envValues))]);
  } else
  {
    valuesDF = dataToShow;
  }
  
  # Create a datatable with collected information
  # Create datatable
  datatable( as.data.frame(valuesDF), 
             class = "compact",
             rownames = rownames,
             colnames = colnames,
             options = list(dom = "<'row'rt>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                            autoWidth = FALSE,
                            columnDefs = list( # Center all columns
                              list( targets = "_all",
                                    className = 'dt-center')),
                            orderClasses = FALSE, # Disable flag for CSS to highlight columns used for ordering (for performance)
                            ordering = FALSE,
                            paging = FALSE, # Disable pagination (show all)
                            processing = TRUE, 
                            scrollCollapse = TRUE,
                            scroller = TRUE,  # Only load visible data
                            scrollX = TRUE,
                            scrollY = "525px",
                            stateSave = TRUE));
}


######## 

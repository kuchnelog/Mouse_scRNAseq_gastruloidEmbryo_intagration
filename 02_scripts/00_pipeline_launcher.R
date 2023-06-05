#!/usr/bin/env Rscript

library(optparse)

source("load_parameters.R")

args <- commandArgs(trailingOnly = TRUE)

main_help <- "
Usage:
Rscript pipeline_launcher.R <pipelineName> [--flag <flag_arg>]

Required Argument:
    <pipelineName>      This argument is required and must be one of the
                        following: qc, process, deg

Flag Argument:
    --flag <flag_arg>   There are optional arguments according to the required
                        argument

Example:
# For more details to each pipeline
Rscript pipeline_launcher.R <pipelineName> --help
"

# make the call to the main help working:
# Rscript pipeline_launcher.R -h         # -h takes the place of args[1] -_-'
if (!is.na(args[1]) && grepl("^-", args[1])) {
  args[1] <- NA
}

# handle the multiple pipelineName input call main help:
# Rscript pipeline_launcher.R qc process
if (length(args) > 1 && !grepl("^-", args[1]) && !grepl("^-", args[2])) {
  cat(main_help)
  quit(status = 1)
}

# call main help message when user call unknown pipelineName:
if (!args[1] %in% c("NA", "qc", "process", "filters", "deg")) {
  cat(main_help)
  quit(status = 1)
}

# defines the parameters to show in the displayed help menu
switch(args[1],
  "NA" = {
    cat(main_help)
    quit(status = 1)
  },
  "qc" = {
    option_list1 <- list(
        make_option(
            c("-i", "--input_dataset"),
            action = "store",
            default = NA,
            type = "character",
            help = ""
        ),
        make_option(
        c("-m", "--mito_high"),
        action = "store",
        default = 10,
        type = "numeric",
        help = "", metavar = "[1:100]"
      ),
      make_option(
        c("-q", "--ribo_low"),
        action = "store",
        default = 25,
        type = "numeric",
        help = ""
      ),
      make_option(
        c("-f", "--min_feat"),
        action = "store",
        default = 200,
        type = "integer",
        help = "minimum number of features detected in every cells"
      ),
      make_option(
        c("-c", "--min_cells"),
        action = "store",
        default = 3,
        type = "integer",
        help = "minimum number in which a feature must be detected"
      ),
      make_option(
        c("--min_counts"),
        action = "store",
        default = 1500,
        type = "integer",
        help = ""
      ),
      make_option(
        c("--max_counts"),
        action = "store",
        default = 150000,
        type = "integer",
        help = ""
      )
    )
    parsed <- OptionParser(
      usage = "Usage: \n\t%prog qc [--flag <flag_arg>]",
      option_list = option_list1,
      add_help_option = TRUE,
      description = "\nMay add description here",
      epilogue = "Add some details, examples, ...",
      formatter = IndentedHelpFormatter # TitleHelpFormatter
    )
  },
  "process" = {
    option_list2 <- list(
        make_option(
            c("-i", "--input_dataset"),
            action = "store",
            default = NA,
            type = "character",
            help = ""
        ),
        make_option(
            c("-f", "--filter"),
            action = "store",
            default = "filtered",
            type = "character",
            help = "indicate whether to use complete or filtered dataset.
                Default value is 'filtered'. Set to 'complete' to study the whole
                dataset."
        ),
      make_option(
        c("-n", "--norm_method"),
        action = "store",
        default = "LogNormalize",
        type = "character",
        help = "normalization method. Imply to not use the same Seurat function in the pipeline.
                The choice is between:
                \t- NormalizeData, with the parameter normalization.method = \"LogNormalize\"
                \t- SCTransform"
      ),
      make_option(
        c("-v", "--hvg_method"),
        action = "store",
        default = "mvp",
        type = "character",
        help = "FindVariableFeatures method. The choice is given between:
                \t- vst, what requires to set the number of features to select --hvg_number
                \t- mvp, where the --hvg_number is useless
                \t- disp, that also requires --hvg_number to be set"
      ),
      make_option(
        c("-w", "--hvg_number"),
        action = "store",
        default = FALSE,
        type = "integer",
        help = "number of highly variable genes to use for the downstream analysis.
                Useful when --hvg_method is set to \"vst\" or \"disp\""
      ),
      make_option(
        c("-d", "--do_scaling"),
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "whether to compute scaling or not. If not mentionned in the command line,
                it will not be done"
      ),
      make_option(
        c("-p", "--pca_npcs"),
        action = "store",
        default = 50,
        type = "integer",
        help = "number of PCs to compute for the Seurat function RunPCA"
      ),
      make_option(
        c("--pca_print"),
        action = "store",
        default = 10,
        type = "integer",
        help = "number of features to print from the top of each PC"
      ),
      make_option(
        c("-t", "--top_pcs"),
        action = "store",
        default = 30,
        type = "integer",
        help = "number of PCs to select for the downstream analysis"
      ),
      make_option(
        c("-o", "--doublet_rate"),
        action = "store",
        default = 8,
        type = "integer",
        help = "fraction of doublets estimated by 10XGenomics,
                according to the technical specifiactions"
      ),
      make_option(
        c("-s", "--selected_resolution"),
        action = "store",
        default = 1,
        type = "double",
        help = "value of the resolution for the Seurat function FindClusters"
      ),
      make_option(
        c("-a", "--algo_clustering"),
        action = "store",
        default = 4,
        type = "integer",
        help = "select the algorithm to run the Seurat function FindClusters
                1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement;
                3 = SLM algorithm; 4 = Leiden algorithm (recommended)"
      )
    )
    parsed <- OptionParser(
      usage = "Usage: \n\t%prog process [--flag <flag_arg>]",
      option_list = option_list2,
      add_help_option = TRUE,
      description = "\nMay add description here",
      epilogue = "Add some details, examples, ...",
      formatter = IndentedHelpFormatter # TitleHelpFormatter
    )
  },
  "filters" = {
      option_list3 <- list(
          make_option(
              c("-i", "--input_dataset"),
              action = "store",
              default = NA,
              type = "character",
              help = ""
          ),
          make_option(
              c("--manual"),
              action = "store",
              default = FALSE,
              type = "logical",
              help = ""
          )
      )
      parsed <- OptionParser(
          usage = "Usage: \n\t%prog filter? [--flag <flag_arg>]",
          option_list = option_list3,
          add_help_option = TRUE,
          description = "\nMay add description here",
          epilogue = "Add some details, examples, ...",
          formatter = IndentedHelpFormatter # TitleHelpFormatter
      )
  },
  "deg" = {
    option_list3 <- list(
        make_option(
            c("-i", "--input_dataset"),
            action = "store",
            default = NA,
            type = "character",
            help = ""
        ),
      make_option(
        c("-k", "--markers_number"),
        action = "store",
        default = 20,
        type = "integer",
        help = "reduce the list of markers to the top [-k] for each cluster at the
				[--selected_resolution]"
      )
    )
    parsed <- OptionParser(
      usage = "Usage: \n\t%prog deg [--flag <flag_arg>]",
      option_list = option_list3,
      add_help_option = TRUE,
      description = "\nMay add description here",
      epilogue = "Add some details, examples, ...",
      formatter = IndentedHelpFormatter # TitleHelpFormatter
    )
  }
)

opt <- parse_args(parsed, positional_arguments = TRUE)

load_parameters(paste0("../01_requirements/globalParameters_", opt$options$input_dataset,".param"))

general_seed = as.numeric(general_seed)
# qc pipeline variables ----------------------
mito_high = if (exists("MITO_HIGH")) as.numeric(MITO_HIGH) else opt$options$mito_high
ribo_low = if (exists("RIBO_LOW")) as.numeric(RIBO_LOW) else opt$options$ribo_low
min_feat = if (exists("MIN_FEAT")) as.numeric(MIN_FEAT) else opt$options$min_feat
min_cells = if (exists("MIN_CELLS")) as.numeric(MIN_CELLS) else opt$options$min_cells
min_counts = if (exists("MIN_COUNTS")) as.numeric(MIN_COUNTS) else opt$options$min_counts
max_counts = if (exists("MAX_COUNTS")) as.numeric(MAX_COUNTS) else opt$options$max_counts
# process pipeline variables -----------------
FILTER = if (exists("FILTER")) FILTER else opt$options$filter
norm_meth = if (exists("NORM_METH")) NORM_METH else opt$options$norm_method
hvg_meth = if (exists("HVG_METH")) HVG_METH else opt$options$hvg_method
hvg_num = if (exists("HVG_NUM")) as.numeric(HVG_NUM) else opt$options$hvg_number
do_scale = if (exists("DO_SCALE")) as.logical(DO_SCALE) else opt$options$do_scaling
pca_npcs = if (exists("PCA_NPCS")) as.numeric(PCA_NPCS) else opt$options$pca_npcs
pca_print = if (exists("PCA_PRINT")) as.numeric(PCA_PRINT) else opt$options$pca_print
top_pcs = if (exists("TOP_PCS")) as.numeric(TOP_PCS) else opt$options$top_pcs
#dblt_rate = if (exists("MITO_HIGH")) as.numeric(MITO_HIGH) else opt$options$doublet_rate
clust_res = if (exists("CLUST_RES")) as.numeric(CLUST_RES) else opt$options$selected_resolution
clust_meth = if (exists("CLUST_METH")) as.numeric(CLUST_METH) else opt$options$algo_clustering
# deg pipeline variables ---------------------
top_markers = opt$options$markers_number
# unclassified variables ---------------------
do_integ = opt$options$do_integ
do_merge = opt$options$do_merge
combine_meth = opt$options$combine_meth
manual = opt$options$manual


switch(args[1],
       "qc" = {
           rmarkdown::render(
               "01_qc_pipeline.Rmd",
               #params = param_list,
               output_file = paste0(PATH_OUT_HTML, "01_qc_", DATASET, ".html")
           )
       },
       "process" = {
           rmarkdown::render(
               "02_process_pipeline.Rmd",
               #params = param_list,
               output_file = paste0(PATH_OUT_HTML, "02_process_", DATASET, "_", FILTER, ".html")
           )
       },
       "filters" = {
           if (opt$options$manual){
               rmarkdown::render(
                   "03_1_manualControls_pipeline.Rmd",
                   #params = param_list,
                   output_file = paste0(PATH_OUT_HTML, "03_1_manualControls_", DATASET, ".html")
               )
           } else {
               rmarkdown::render(
                   "03_filtersControls_pipeline.Rmd",
                   #params = param_list,
                   output_file = paste0(PATH_OUT_HTML, "03_filtersControls_", DATASET, ".html")
               )   
           }
       },
       "deg" = {
           rmarkdown::render(
               "03_deg_pipeline.Rmd",
               #params = param_list,
               output_file = paste0(PATH_OUT_HTML, "03_deg_", DATASET, ".html")
           )
       }
)


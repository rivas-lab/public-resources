#########
# tsvd_cmdarg_parser.R
#   Yosuke Tanigawa (2017/10/15)
#
#  This script has command line arguement parser
#########


tsvd_cmdarg_parser <- function(bim_filter_names, use_default = FALSE){
    args <- commandArgs(trailingOnly = TRUE)
    
    if(use_default){
        param.suffix             <- "yt_debug"
        param.bim_filter_type    <- "codingNonMHC"
        param.z_or_lor           <- "z"
        param.center_tsvd_str    <- "nonCenter"
        param.pvalue_str         <- "0.001"
        param.num_components_str <- "100"        
        param.algorithm          <- "tsvd"
        param.value_filter_mode  <- "default"
        
    }else if(length(args) < 6){
        print("[error] wrong number of args")
        print("  usage: Rscript tsvd.R <suffix (dataset name))> <bim_filter_name> <z | lor> <center | non-center> <p-val> <num_components>")
        print("The following filters are available:")
        print(names(bim_filters))
        quit("no", 1)

    }else{                
        param.suffix             <- args[1]
        param.bim_filter_type    <- args[2]
        param.z_or_lor           <- args[3]
        param.center_tsvd_str    <- args[4]
        param.pvalue_str         <- args[5]
        param.num_components_str <- args[6]      
        if(length(args) > 6){
            param.algorithm <- args[7]
        }else{
            param.algorithm <- "tsvd"
        }        
        if(length(args) > 7){
            param.value_filter_mode <- args[8]
        }else{
            param.value_filter_mode <- "default"
        }
    }


    if(! param.bim_filter_type %in% bim_filter_names){
      print(paste0("[error] bim_filter '", param.bim_filter_type, "' is not available."))
      print("The following filters are available:")
      print(bim_filter_names)
      quit("no", 1)
    }
    
    if(! param.z_or_lor %in% c("z", "lor")){
        print(sprintf("[error] can't recognize %s. Please pass either z or lor", param.z_or_lor))
        quit("no", 1)        
    }

    if(! param.center_tsvd_str %in% c("center", "nonCenter")){
        print(sprintf("[error] can't recognize %s. Please pass either center or nonCenter", param.center_tsvd_str))
        quit("no", 1)        
    }
    
    if(! param.algorithm %in% c("tsvd", "ssvd")){
        print(sprintf("[error] can't recognize %s. Please pass either tsvd or ssvd", param.algorithm))
        quit("no", 1)        
    }        
        
    dummy_pvalue         <- as.numeric(param.pvalue_str)
    dummy_num_components <- as.integer(param.num_components_str)

    if ( param.z_or_lor == "z" ) {
      param.z_or_lor_se <- "lor/se"
    } else if ( param.z_or_lor == "lor" ) {
      param.z_or_lor_se <- "lor"
    }

    if(! param.value_filter_mode %in% c("default", "pvalue", "se")){
        print(sprintf("[error] can't recognize %s. Please specify one of the followings: default, pvalue, or se.", param.value_filter_mode))
        quit("no", 1)
    }

    params <- c(
        param.suffix, 
        param.bim_filter_type,
        param.z_or_lor,
        param.z_or_lor_se,
        param.center_tsvd_str,
        param.pvalue_str,
        param.num_components_str,
        param.algorithm,
        param.value_filter_mode
    )
    
    names(params) <- c(
        "suffix", 
        "bim_filter_type",
        "z_or_lor",
        "z_or_lor_se",
        "center_tsvd_str",
        "pvalue_str",
        "num_components_str",
        "algorithm",
        "value_filter_mode"
    )
    
    print("the following arguments are used for this computation:")
    print(params)
    
    return(params)    
}


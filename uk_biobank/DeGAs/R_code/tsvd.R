suppressMessages(library(scidb))
suppressMessages(library(irlba))
DB = scidbconnect()

#Exclude the tail of NULL ICDs
#523 originally
# 7020 as of 2018/5/15 (for the production)
end_icd_idx=7551


use_default_command_line_arguments <- TRUE
use_default_command_line_arguments <- FALSE

source('tsvd_filters.R')
source('tsvd_cmdarg_parser.R')

##########################################################################################
# parse command line arguments. 
#   one can use 'debug = TRUE' to set the default params
##########################################################################################
params <- tsvd_cmdarg_parser(get_bim_filter_names(), use_default = use_default_command_line_arguments)
pvalue <- as.numeric(params[["pvalue_str"]])
if (params[["center_tsvd_str"]] == "center") {
  center_tsvd <- TRUE
} else if (params[["center_tsvd_str"]] == "nonCenter") {
  center_tsvd <- FALSE
}

##########################################################################################
array_suffix=paste0("_", params[["suffix"]])

icd_idx_array       = paste0('ap_icd_idx',          array_suffix)
icd_prefilter_array = paste0('ap_icd_pf',           array_suffix)
variant_idx_array   = paste0('ap_variant_idx',      array_suffix)
var_matrix_array    = paste0('ap_icd_var_matrix',   array_suffix)
var_matrix_t_array  = paste0('ap_icd_var_matrix_t', array_suffix)
initial_array       = paste0('ap_icd_initial',      array_suffix)
left_array          = paste0('ap_icd_left',         array_suffix)
right_array         = paste0('ap_icd_right',        array_suffix)
tsvd_array          = paste0('ap_icd_var_tsvd',     array_suffix)
cm_array            = paste0('ap_icd_var_centered', array_suffix)
cor1_array          = paste0('ap_icd_cor1',         array_suffix)
cor2_array          = paste0('ap_icd_cor2',         array_suffix)

##########################################################################################
# set filters
##########################################################################################
bim_filters <- set_bim_filters()
icd_filter  <- set_icd_filter()
value_filter <- set_value_filters(pvalue = pvalue, mode = params[["value_filter_mode"]])
bim_filter <- bim_filters[[params[["bim_filter_type"]]]]

print(c("bim_filter:", bim_filter))
print(c("icd_filter:", icd_filter))
print(c("value_filter:", value_filter))

##########################################################################################

remove_versions = function(array)
{
  mv = max(iquery(DB, sprintf("versions(%s)", array), return=TRUE)$version_id)
  iquery(DB, sprintf("remove_versions(%s, %i)", array, mv))
}

build_icd_matrix = function()
{
  tryCatch({ iquery(DB,  paste0("remove(", icd_idx_array,         ")" )) }, error =invisible)
  tryCatch({ iquery(DB,  paste0("remove(", icd_prefilter_array,   ")" )) }, error =invisible)
  tryCatch({ iquery(DB,  paste0("remove(", variant_idx_array,     ")" )) }, error =invisible)
  tryCatch({ iquery(DB,  paste0("remove(", var_matrix_array,      ")" )) }, error =invisible)
  tryCatch({ iquery(DB,  paste0("remove(", var_matrix_t_array,    ")" )) }, error =invisible)
  print("Making pre-filter array")
  iquery(DB, sprintf("
   store(
    project(
     apply(
      cross_join(
       cross_join(
        filter(icd,      %s) as A,
        filter(join(icd_info, aggregate(icd, count(*) as num_cells, icd_idx)), (%s) and icd_idx<=%i) as B,
        A.icd_idx, B.icd_idx
       ) as C,
       project(
        join(
         filter(variant, filter='PASS'),
         filter(bim, %s)
        ),
        ld
       ) as D,
       C.chrom, D.chrom,
       C.pos, D.pos
      ),
      z, %s
     ),
     z
    ), 
    %s
   )", value_filter, icd_filter, end_icd_idx, bim_filter, params[["z_or_lor_se"]], icd_prefilter_array))
  print("Making variant index")
  iquery(DB, sprintf("
  store(
   cast(
    sort(
     project(
      grouped_aggregate(
       apply(
        %s, 
        signature, string(chrom) + '_' + string(pos)
       ),
       count(*), signature
      ),
      signature
     ),
     signature, 50000
    ),
    <signature:string> [variant_idx]
   ),
   %s
  )",
  icd_prefilter_array, variant_idx_array))
  nvar = iquery(DB, sprintf("op_count(%s)", variant_idx_array), return=T)$count
  print(paste("got", nvar, "variants"))
  print("Making ICD index")
  iquery(DB, sprintf("
   store(
    project(
     unpack(
      join(
       icd_info,
       aggregate(
        %s,
        count(*) as num_cells,
        icd_idx
       )
      ),
      new_icd_idx
     ),
     icd_idx, icd, Case, Name, num_cells
    ),
    %s
   )",
   icd_prefilter_array, icd_idx_array))
  nicd = iquery(DB, sprintf("op_count(%s)", icd_idx_array), return=T)$count
  print(paste("got", nicd, "icd codes"))
  iquery(DB, sprintf("create array %s
                    <z:double not null>
                    [new_icd_idx=0:%i,      1, 0, 
                     variant_idx = 0:%i, %i, 0]", 
                     var_matrix_array, nicd-1, nvar-1, nvar))
  iquery(DB, sprintf("create array %s
                    <z:double not null>
                    [variant_idx = 0:%i, 1024, 0,
                     new_icd_idx     = 0:%i,%i,   0]",
                     var_matrix_t_array, nvar-1, nicd-1, nicd))
  #TSVD does not need a zero pad
  #print("Performing zero-pad")
  #iquery(DB, "store(build(ap_icd_var_matrix, 0), ap_icd_var_matrix)")
  #iquery(DB, "store(build(ap_icd_var_matrix_t, 0), ap_icd_var_matrix_t)")
  icd_chunk =80
  icd_idx = 0
  while(icd_idx< end_icd_idx)
  {
    segment_icd_idx= min(icd_idx + icd_chunk -1, end_icd_idx)
    print(paste("Inserting ICD", icd_idx, "to", segment_icd_idx))
    t1=proc.time()
    redim_query = sprintf("
         substitute(
          index_lookup(
           index_lookup(
            project(
             apply(
              between(
               %s, 
               %i, null, null, null, null,
               %i, null, null, null, null
              ),
              signature,string(chrom) + '_' + string(pos),
              old_icd_idx, icd_idx
             ),
             old_icd_idx, z, signature
            ) as A,
            %s,
            A.signature,
            variant_idx
           ),
           project(%s, icd_idx),
           A.old_icd_idx,
           new_icd_idx
          ),
          build(<val:double>[i=0:0,1,0],0),
          z
         )",icd_prefilter_array, icd_idx, segment_icd_idx,  
            variant_idx_array, icd_idx_array)
    iquery(DB,sprintf(
      "insert(
        redimension(
         %s,
         %s,
         false
        ),
        %s
      )", redim_query, var_matrix_array, var_matrix_array
     )
    ) 
    remove_versions(var_matrix_array)
    print("Transposing") 
    iquery(DB,sprintf(
      "insert(
        redimension(
         %s,
         %s,
         false
        ),
        %s
      )", redim_query, var_matrix_t_array, var_matrix_t_array
     )
    )    
    remove_versions(var_matrix_t_array)   
    icd_idx = icd_idx + icd_chunk
    print(proc.time()-t1)
  }
}

compute_tsvd = function()
{
  nvar = iquery(DB, sprintf("op_count(%s)", variant_idx_array), return=T)$count
  nicd = iquery(DB, sprintf("op_count(%s)", icd_idx_array), return=T)$count
  end_icd_idx=nicd-1
  print(paste("Running TSVD on",nicd,"ICDs by",nvar,"variants")) 
  tryCatch({ iquery(DB,  paste0("remove(", initial_array,     ")" )) }, error =invisible)
  tryCatch({ iquery(DB,  paste0("remove(", left_array,        ")" )) }, error =invisible)
  tryCatch({ iquery(DB,  paste0("remove(", right_array,       ")" )) }, error =invisible)
  tryCatch({ iquery(DB,  paste0("remove(", tsvd_array,        ")" )) }, error =invisible)
  t1=proc.time()
  print("Running....")
  iquery(DB, sprintf("store(redimension(filter(tsvd_random_seed, variant_idx<=%i),<val:double not null> [variant_idx=0:%i:0:%i]), %s)", 
                      nvar-1, nvar-1, nvar, initial_array))
  iquery(DB, sprintf("store(build(<v:double not null> [i=0:%i,%i,0], 1), %s)",end_icd_idx, nicd, left_array))
  iquery(DB, sprintf("store(
       project(
         substitute(
           apply(
             aggregate(
               %s,
               sum(z) as colsum,
               variant_idx
             ),
             colmean,
             colsum * 1.0 / %i
           ),
           build(<v:double NOT NULL>[i=0:0], 0.0),
           colmean
         ),
         colmean
       ),
       %s
     )", var_matrix_t_array, nicd, right_array))
  if( center_tsvd ) {
    #CENTERED: NOTE the number of singular vectors
    iquery(DB, sprintf("store(tsvd(%s, %s, %s, 0.0001, 2000, %s, %s, %s), %s)",
                        var_matrix_array, var_matrix_t_array, params[["num_components_str"]], initial_array, left_array, right_array, tsvd_array))
  } else {
    #NOT CENTERED:
    iquery(DB, sprintf("store(tsvd(%s, %s, %s, 0.0001, 2000, %s), %s)", 
                        var_matrix_array, var_matrix_t_array, params[["num_components_str"]], initial_array, tsvd_array))
  }
  print(paste("Computed ", tsvd_array))
  print(proc.time()-t1)
}

compute_cor = function()
{
  nvar = iquery(DB, sprintf("op_count(%s)", variant_idx_array), return=T)$count
  nicd = iquery(DB, sprintf("op_count(%s)", icd_idx_array), return=T)$count
  end_icd_idx=nicd-1
  print(paste("Running post-TSVD cor on",nicd,"ICDs by",nvar,"variants")) 
  tryCatch({ iquery(DB,  paste0("remove(", cm_array,     ")" )) }, error =invisible)
  tryCatch({ iquery(DB,  paste0("remove(", cor1_array,        ")" )) }, error =invisible)
  tryCatch({ iquery(DB,  paste0("remove(", cor2_array,       ")" )) }, error =invisible)
  mat_schema = paste0("<z:double not null>[new_icd_idx=0:",nicd-1,":0:512; variant_idx=0:",nvar-1,":0:512]")
  print("Creating Centered Matrix")
  iquery(DB,  paste0("create array ", cm_array, " ", mat_schema))
  iquery(DB,  sprintf("store(build( %s, 0), %s)", cm_array, cm_array))
  icd_row=0
  while(icd_row<nicd)
  {
    icd_row_end =icd_row+511
    icd_row_end = min(icd_row_end, nicd-1)
    print(paste("Inserting idx",icd_row, "to", icd_row_end))
    iquery(DB, sprintf("insert(redimension(between(%s, %i, null, %i, null), %s), %s)",
                       var_matrix_array, icd_row, icd_row_end, cm_array, cm_array))
    remove_versions(cm_array)
    icd_row = icd_row_end +1
  }
  print("Centering")
  iquery(DB, sprintf("
   store(
    substitute(
     project(
      apply(
       cross_join(
        %s as A, 
        aggregate(%s, avg(z) as mean, variant_idx) as B,
        A.variant_idx, B.variant_idx
       ),
       z_centered, z - mean
      ),
      z_centered
     ),
     build(<val:double>[i=0:0,1,0],0)
    ),
    %s
   )", cm_array, cm_array, cm_array
  ))
  remove_versions(cm_array)
  print("Computing cor1")
  iquery(DB, sprintf(
   "store(
     dmetric(
      %s,
      transpose(redimension(apply(slice(%s, matrix, 2), z, value), %s)),
      'metric=pearson',
      'rightReplicate=1'
     ),
     %s
    )", cm_array, tsvd_array, cm_array, cor1_array))
  print(paste("Computed ", cor1_array))
  print("Computing cor2")
  iquery(DB, sprintf(
   "store(
     dmetric(
      transpose(%s),
      redimension(apply(slice(%s, matrix, 0), z, value), %s),
      'metric=pearson',
      'rightReplicate=1'
     ),
     %s
    )", cm_array, tsvd_array, cm_array, cor2_array))
  print(paste("Computed ", cor2_array))
}


get_matrix_from_scidb <- function(db, matrix_name) {
  matrix_data <- as.R(scidb(db, matrix_name))
  matrix <- Matrix::sparseMatrix(
    matrix_data$variant_idx+1, 
    matrix_data$new_icd_idx+1, 
    x=matrix_data$z
  )
  return(matrix)
}


compute_ssvd <- function(db, matrix_name, num_components, rds_file){
    set.seed(1)
    print(paste("Computing SSVD"))    
    ssvd_results <- ssvd(
      get_matrix_from_scidb(db, matrix_name),
      k=as.numeric(num_components),
      n=50
    )
    saveRDS(ssvd_results, rds_file)
}

## run the computation.

build_icd_matrix()

if( params[["algorithm"]] == "tsvd" ){
    compute_tsvd()
    #compute_cor()    
}else if( params[["algorithm"]] == "ssvd" ){
    compute_ssvd(
        db = DB,
        matrix_name = var_matrix_array,
        num_components = params[["num_components_str"]], 
        rds_file = paste('/home/scidb/R_code/results', params[["suffix"]], 'ssvd.rds', sep='/')
    )    
}


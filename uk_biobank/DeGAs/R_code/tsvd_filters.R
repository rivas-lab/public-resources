#########
# tsvd_filters.R
#   Yosuke Tanigawa (2017/10/15, updated on 2018/01/22)
#
#  This script has functions to define bim_filters and icd_filter
#########
set_icd_filter <- function(){
    icd_filter <- "(Case > 100) and (num_cells > 10000) and ((icd_idx<>6585) and (icd_idx<>6819) and (icd_idx<>6820) and (icd_idx<>4794) and (icd_idx<>4795) and (icd_idx<>4796) and (icd_idx<>4797) and (icd_idx<>4790) and (icd_idx<>4791) and (icd_idx<>4792) and (icd_idx<>5210) and (icd_idx<>4788) and (icd_idx<>4789) and (icd_idx<>4947) and (icd_idx<>4946) and (icd_idx<>4945) and (icd_idx<>4944) and (icd_idx<>4951) and (icd_idx<>4950) and (icd_idx<>4949) and (icd_idx<>4948) and (icd_idx<>6181) and (icd_idx<>6178) and (icd_idx<>4389) and (icd_idx<>4390) and (icd_idx<>4387) and (icd_idx<>4388) and (icd_idx<>4393) and (icd_idx<>4394) and (icd_idx<>4391) and (icd_idx<>4392) and (icd_idx<>5557) and (icd_idx<>4382) and (icd_idx<>4580) and (icd_idx<>4579) and (icd_idx<>4582) and (icd_idx<>4581) and (icd_idx<>4584) and (icd_idx<>4583) and (icd_idx<>4586) and (icd_idx<>4585) and (icd_idx<>4572) and (icd_idx<>4571) and (icd_idx<>5427) and (icd_idx<>5428) and (icd_idx<>5429) and (icd_idx<>5430) and (icd_idx<>6963) and (icd_idx<>6962) and (icd_idx<>1988) and (icd_idx<>1453) and (icd_idx<>992) and (icd_idx<>2044) and (icd_idx<>2485) and (icd_idx<>2486) and (icd_idx<>2837) and (icd_idx<>1486) and (icd_idx<>725) and (icd_idx<>6927) and (icd_idx<>6945) and (icd_idx<>6942) and (icd_idx<>6896) and (icd_idx<>1469) and (icd_idx<>2526) and (icd_idx<>1993) and (icd_idx<>2969) and (icd_idx<>3003) and (icd_idx<>2967) and (icd_idx<>2966) and (icd_idx<>2972) and (icd_idx<>2270) and (icd_idx<>3424) and (icd_idx<>2802) and (icd_idx<>2976) and (icd_idx<>2975) and (icd_idx<>2488) and (icd_idx<>2489) and (icd_idx<>2490) and (icd_idx<>977) and (icd_idx<>2492) and (icd_idx<>2493) and (icd_idx<>983) and (icd_idx<>4109) and (icd_idx<>4110) and (icd_idx<>4121) and (icd_idx<>4108) and (icd_idx<>4118) and (icd_idx<>4061) and (icd_idx<>3368) and (icd_idx<>2888) and (icd_idx<>493) and (icd_idx<>497) and (icd_idx<>1426) and (icd_idx<>1423) and (icd_idx<>1419) and (icd_idx<>1420) and (icd_idx<>2924) and (icd_idx<>2923) and (icd_idx<>2452) and (icd_idx<>935) and (icd_idx<>936) and (icd_idx<>934) and (icd_idx<>939) and (icd_idx<>937) and (icd_idx<>2237) and (icd_idx<>2242) and (icd_idx<>2241) and (icd_idx<>2230) and (icd_idx<>504) and (icd_idx<>471) and (icd_idx<>714) and (icd_idx<>764) and (icd_idx<>762) and (icd_idx<>759) and (icd_idx<>758) and (icd_idx<>1206) and (icd_idx<>1203) and (icd_idx<>2828) and (icd_idx<>2827) and (icd_idx<>2832) and (icd_idx<>1970) and (icd_idx<>3174) and (icd_idx<>2218) and (icd_idx<>718) and (icd_idx<>1771) and (icd_idx<>1769) and (icd_idx<>1768) and (icd_idx<>710) and (icd_idx<>2908) and (icd_idx<>2911) and (icd_idx<>2909) and (icd_idx<>2916) and (icd_idx<>709) and (icd_idx<>2460) and (icd_idx<>2457) and (icd_idx<>2456) and (icd_idx<>2454))"
    return(icd_filter)
}

get_bim_filter_names <- function(){
    return(c(
      "PTVs", "coding", "nonCoding", "all",
      "PTVsNonMHC","codingNonMHC", "nonCodingNonMHC", "allNonMHC",
      "PTVsNoLDpruning","codingNoLDpruning", "nonCodingNoLDpruning", "allNoLDpruning",
      "PTVsNoLDpruningNonMHC","codingNoLDpruningNonMHC", "nonCodingNoLDpruningNonMHC", "allNoLDpruningNonMHC",
      "allMafOnly", "allMafOnlyNonMHC"
    ))
}

set_bim_filters <- function(){
    #Exclude variants that don't match this filter

    filter_base      = " (allfilter=0 and ld=True and maf >= .0001) "

    filter_noLDPruning      = " (allfilter=0 and maf >= .0001) "

    filter_mafOnly      = " (maf >= .0001) "

    filter_nonMHC    = " ((chrom = 6 and pos < 25477797 or pos > 36448354) or (chrom <> 6)) "

    filter_PTVs      = " (consequence='stop_gained' or consequence='frameshift_variant' or consequence='splice_acceptor_variant' or consequence='splice_donor_variant') "

    filter_coding    = " (consequence='missense_variant' or consequence='stop_gained' or consequence='frameshift_variant' or consequence='splice_acceptor_variant' or consequence='splice_donor_variant' or consequence='splice_region_variant' or consequence='start_lost' or consequence='stop_lost') "

    filter_noncoding = " (consequence<>'missense_variant' and consequence<>'stop_gained' and consequence<>'frameshift_variant' and consequence<>'splice_acceptor_variant' and consequence<>'splice_donor_variant' and consequence<>'splice_region_variant' and consequence<>'start_lost' and consequence<>'stop_lost') "

    bim_filters <- c(
      paste(filter_base, filter_PTVs, sep=" and "),
      paste(filter_base, filter_coding, sep=" and "),
      paste(filter_base, filter_noncoding, sep=" and "),  
      paste(filter_base, sep=" and "),  
      paste(filter_base, filter_nonMHC, filter_PTVs, sep=" and "),
      paste(filter_base, filter_nonMHC, filter_coding, sep=" and "),
      paste(filter_base, filter_nonMHC, filter_noncoding, sep=" and "),  
      paste(filter_base, filter_nonMHC, sep=" and "),
      paste(filter_noLDPruning, filter_PTVs, sep=" and "),
      paste(filter_noLDPruning, filter_coding, sep=" and "),
      paste(filter_noLDPruning, filter_noncoding, sep=" and "),  
      paste(filter_noLDPruning, sep=" and "),  
      paste(filter_noLDPruning, filter_nonMHC, filter_PTVs, sep=" and "),
      paste(filter_noLDPruning, filter_nonMHC, filter_coding, sep=" and "),
      paste(filter_noLDPruning, filter_nonMHC, filter_noncoding, sep=" and "),  
      paste(filter_noLDPruning, filter_nonMHC, sep=" and "),
      paste(filter_mafOnly, sep=" and "),
      paste(filter_mafOnly, filter_nonMHC, sep=" and ")
    )

    names(bim_filters) <- get_bim_filter_names()
    
    return(bim_filters)
}

set_value_filters <- function(pvalue = 1, mode = 'default'){
    pvalue_filter  = sprintf('(pvalue<%g)', pvalue)
    se_filter      = '((se < .08 and lor = or_val) or (se < .2 and lor <> or_val))'
    default_filter = sprintf('(%s and %s)', pvalue_filter, se_filter)

    if(mode == "pvalue"){
        return(pvalue_filter)
    }else if(mode == "se"){
        return(se_filter)
    }else{
        return(default_filter)
    }
}


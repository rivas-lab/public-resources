library(tidyverse)

## usage of the script
## 
## $ Rscript circular_bar_great.R ../private_data/results/dev_codingNonMHC_lor_center_p001_100PCs_20180123.npz/phenotypes/asthma MGIPhenotype BFold
## 

source('circular_bar.R')

args <- commandArgs(TRUE)
phenotype_dir <- args[1]
ontology <- args[2]
score <- args[3]

circos_data  <- file.path(phenotype_dir, 'great', ontology, 'circos-data.csv')
cosine_score <- file.path(phenotype_dir, 'squared_cosine_scores.tsv') 
figure_out <-   file.path(phenotype_dir, 'great', ontology, 'circos.pdf')  

## read data
print('readling data')

if(score == 'BFold'){
    data <- read_csv(circos_data)  %>% 
    mutate(Score = BFold, Alpha = -log10(BPval)) %>% 
    select(-PC_rank, -BFold, -BPval) %>%
    rename(
        Group_id = PC,
        Label = Term
    )    
}else if(score == 'BPval'){
    data <- read_csv(circos_data)  %>% 
    mutate(Alpha = BFold, Score = -log10(BPval)) %>% 
    select(-PC_rank, -BFold, -BPval) %>%
    rename(
        Group_id = PC,
        Label = Term
    )    
    
}

print('readling cosine scores')    
groups <- read.table(cosine_score) %>% 
rename(
    Group_order = V1,
    PC_zero_based = V2,
    Group_fraction = V3
) %>% mutate(
    Group_id = paste0('PC', PC_zero_based + 1)
) %>% select(-PC_zero_based)
    
## plot
    
circular_bar_plot(data, groups, loc_margin=0.05, quantile_thr=0.05, alpha_min=0.3)    

ggsave(figure_out)

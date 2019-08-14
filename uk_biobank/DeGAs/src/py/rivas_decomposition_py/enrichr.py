import os, logging, collections
import numpy as np
import pandas as pd

from logging.config import dictConfig
from logging import getLogger

dictConfig(dict(
    version = 1,
    formatters = {'f': {'format': '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'}},
    handlers = {
        'h': {'class': 'logging.StreamHandler','formatter': 'f',
              'level': logging.DEBUG}},
    root = {'handlers': ['h'], 'level': logging.DEBUG,},
))

logger = getLogger('enrichr')
import enrichr_py


def generate_enrichr_query(d, out_dir, topk=None,
    contribution_thr=None, contribution_cumsum_thr=None):
    logger = logging.getLogger('generate_enrichr_query')
    if(not os.path.exists(os.path.join(out_dir, 'enrichr'))):
        os.mkdir(os.path.join(out_dir, 'enrichr'))    
    if(not os.path.exists(os.path.join(out_dir, 'enrichr', 'query'))):
        os.mkdir(os.path.join(out_dir, 'enrichr', 'query'))
    
    logger.debug('\t'.join(['PC', 'num_genes', 'sum(contribution_score)']))
    nums = np.zeros(d.d['n_PCs'], dtype=np.int)
    scores = np.zeros(d.d['n_PCs'], dtype=np.float)
    
    logger.info(d.d['n_PCs'])
    for i in range(d.d['n_PCs']):    
        try:
            txt, num, score=d.text_data_contribution_gene(
                pc_index=i, topk=topk, 
                contribution_thr=contribution_thr, 
                contribution_cumsum_thr=contribution_cumsum_thr)
            logger.info('PC{:02d}\t{}\t{:.4f}'.format(i+1, num, score))        
            nums[i]=num
            scores[i]=score
            with open(            
                os.path.join(out_dir, 'enrichr', 'query', 'PC{}.txt'.format(i)), 
                'w') as f:
                f.write(txt)
        except Exception as e:
            logger.warning('Error at PC{:02d}'.format(i))
            logger.warning(e)
    pd.DataFrame({
        '#PC_(zero_based)': np.arange(d.d['n_PCs']),
        'num_genes': nums,
        'sum_contribution_scores': scores
    })[['#PC_(zero_based)', 'num_genes', 'sum_contribution_scores']].to_csv(
        os.path.join(out_dir, 'enrichr', 'query', 'Enrichr_query_meta.tsv'), 
        index=False, sep='\t'
    )

    
def run_enrichr(d, out_dir, gene_set_library_list):
    logger = logging.getLogger('run_enrichr')
    if(not os.path.exists(os.path.join(out_dir, 'enrichr', 'results'))):
        os.mkdir(os.path.join(out_dir, 'enrichr', 'results'))
    logger.info('we are testing the following gene sets:')
    for gene_set_library in gene_set_library_list:
        logger.info('  {}'.format(gene_set_library))
        if(not os.path.exists(os.path.join(out_dir, 'enrichr', 'results', gene_set_library))):
            os.mkdir(os.path.join(out_dir, 'enrichr', 'results', gene_set_library))
            
    for pc_index in range(d.d['n_PCs']):  
        query_str_f = os.path.join(out_dir, 'enrichr', 'query', 'PC{}.txt'.format(pc_index))
        if(os.path.exists(query_str_f) and 
           (not os.path.exists(os.path.join(
               out_dir, 'enrichr', 'results', gene_set_library, 'PC{}.txt'.format(pc_index))))):
            logger.info(pc_index)
            with open(query_str_f) as f:
                query_str = f.read()
            try:
                e = enrichr_py.Enrichr(query_str)
            except Exception as E:
                print(E)

            for gene_set_library in gene_set_library_list:
                try:
                    e.enrich(
                        gene_set_library, 
                        outFileName = os.path.join(
                            out_dir, 'enrichr', 'results', 
                            gene_set_library, 'PC{}.txt'.format(pc_index)))
                except Exception as E:
                    print(E)
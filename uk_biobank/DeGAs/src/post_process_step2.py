from __future__ import print_function

_README_ = '''
-------------------------------------------------------------------------
Post_process_step2.py

After running GREAT and place them on a proper directory, 
apply this script to generate circular bar plots

Usage: python post_process_step2.py -d ../private_data/results/dev_PTVsNonMHC_z_nonCenter_p0001_100PCs.npz -p ../public_data/phenotype_of_interest.lst

Author: Yosuke Tanigawa (ytanigaw@stanford.edu)
Date: 2018/01/22
-------------------------------------------------------------------------
'''

import numpy as np
import pandas as pd
import itertools as it
import os
import logging
import subprocess
import argparse


from logging.config import dictConfig


import great_io

######

dictConfig(dict(
    version = 1,
    formatters = {
        'f': {'format':
              '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'}
        },
    handlers = {
        'h': {'class': 'logging.StreamHandler',
              'formatter': 'f',
              'level': logging.DEBUG}
        },
    root = {
        'handlers': ['h'],
        'level': logging.INFO,
        },
))
    
def read_great_results(great_res_dir, pc_index_list, ontology, BFold = None, BPval = None, sort_by='BFold'):
    def read_great_results_sub(great_res_dir, pc_index, ontology, BFold = None, BPval = None, sort_by='BFold'):
        df=great_io.read_great_res(
            filename=os.path.join(great_res_dir, 'PC{}'.format(pc_index), '{}.tsv'.format(ontology)),
            BFold=BFold, BPval=BPval, sort_by=sort_by
        )[['Desc', 'BFold', 'BPval']]
        df['pc_index'] = pc_index
        return df
       
    assert(len(pc_index_list) > 0)
    if(len(pc_index_list) == 1):
        return read_great_results_sub(great_res_dir, pc_index_list[0], ontology, BFold, BPval, sort_by=sort_by)
    else:
        df=read_great_results_sub(great_res_dir, pc_index_list[0], ontology, BFold, BPval, sort_by=sort_by)
        for i in range(len(pc_index_list) - 1):
            df = df.append(
                read_great_results_sub(great_res_dir, pc_index_list[i + 1], ontology, BFold, BPval, sort_by=sort_by)
            )
        return df

    

def prepare_circos_csv(
    out_dir, phe, BFold = None, BPval = None, 
    sort_by='BFold', squared_cos_cumsum_threshold = 0.80,
    ontology='MGIPhenotype'
):
    logger = logging.getLogger('prepare_circos_csv')
    logger.info(phe)
    
    great_res_dir = os.path.join(out_dir, 'great', 'results')
    if(not os.path.exists(great_res_dir)):
        try:
            subprocess.run(
                ('tar', '-xzvf', 'results.tar.gz'), check=True, cwd=os.path.join(out_dir, 'great')
            )
        except Exception as e:
            logger.warning(out_dir)
            logger.warning(e)
            raise e
    

    phe_safe = ''.join([
        c if (c.isalnum() or c in ['_', '.'] ) else '_' for c in phe])
    phe_dir = os.path.join(out_dir, 'phenotypes', phe_safe)


    squared_cos_all = pd.read_table(os.path.join(phe_dir, 'squared_cosine_scores.tsv'))
    squared_cos_all['PC'] = squared_cos_all['PC_(zero_based)'].map(
        lambda x: 'PC{}'.format(x + 1)
    )

    squared_cos = squared_cos_all.head(
        np.sum(squared_cos_all['squared_cosine_score'].cumsum() <= squared_cos_cumsum_threshold) + 1
    )
    
    if(not os.path.exists(os.path.join(phe_dir, 'great'))):
        os.mkdir(os.path.join(phe_dir, 'great'))    
    if(not os.path.exists(os.path.join(phe_dir, 'great', ontology))):
        os.mkdir(os.path.join(phe_dir, 'great', ontology))    
        

    circular_plot_data = read_great_results(
        great_res_dir, squared_cos['PC_(zero_based)'], ontology,
        BFold=BFold, BPval=BPval, sort_by=sort_by
    ).merge(
        squared_cos, left_on='pc_index', right_on='PC_(zero_based)'
    ).rename(columns={
        '#Rank': 'PC_rank', 'Desc': 'Term',
    })[['PC_rank', 'PC' , 'Term', 'BFold', 'BPval']]
    
    circular_plot_data.to_csv(
        os.path.join(phe_dir, 'great', ontology, 'circos-data.csv'), index=False
    )
    
    if(len(circular_plot_data) > 0):
        subprocess.run(
            ('Rscript', os.path.join(os.path.dirname(__file__), 'circular_bar_great.R'), phe_dir, ontology, sort_by), 
            check=True)
    


def post_process_step2(out_dir, phe_list, BFold, BPval, sort_by, squared_cos_cumsum_threshold, ontology='MGIPhenotype'):
    logger = logging.getLogger('post_process_step2')
    for phe in phe_list:
        try:
            prepare_circos_csv(out_dir, phe, BFold, BPval, sort_by, squared_cos_cumsum_threshold, ontology)
        except ValueError as e:
            logger.warning(e)            
            logger.warning("phenotype {} is not found in the data. Skipping".format(phe))


def main():
    logger_main = logging.getLogger('main')    
    
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=_README_
    )
    
    parser.add_argument('-d', metavar='d', required=True,
                        help='out_dir')
    parser.add_argument('-p', metavar='p', required=True,
                        help='phenotype list file')
    parser.add_argument('--onto', metavar='o', default='MGIPhenotype',
                        help='ontology (default: MGIPhenotype)')

    args = parser.parse_args()

    out_dir=os.path.abspath(args.d)
    phe_list_f=os.path.abspath(args.p)
    ontology=args.onto
    
    logger_main.info('  out_dir   : {}'.format(out_dir))
    logger_main.info('  phenotype : {}'.format(phe_list_f))
    logger_main.info('  ontology  : {}'.format(ontology))
           
    with open(phe_list_f) as f:
        phe_list=[x for x in f.read().splitlines() if len(x) > 0]        
        
    post_process_step2(
        out_dir, phe_list, BFold=2, BPval=float('5e-7'), sort_by='BFold',
        squared_cos_cumsum_threshold=0.80, 
	ontology=ontology
	)    

    
if __name__ == "__main__":
    main()     
    

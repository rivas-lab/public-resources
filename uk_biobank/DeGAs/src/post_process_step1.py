from __future__ import print_function

_README_ = '''
-------------------------------------------------------------------------
Post_process_step1.py

Given a npz file, generate plots and GREAT query bed file

npz_file, out_dir_root, phe_list

Usage: python post_process_step1.py -d ../private_data/results/dev_PTVsNonMHC_z_nonCenter_p0001_100PCs.npz -p ../public_data/phenotype_of_interest.lst

Author: Yosuke Tanigawa (ytanigaw@stanford.edu)
Date: 2018/01/22 (update on 2018/5/18)
-------------------------------------------------------------------------
'''

import os, logging, collections, itertools
import numpy as np
import pandas as pd
import argparse
import matplotlib

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

matplotlib.rc('font',**{'size':16, 'family':'sans-serif','sans-serif':['HelveticaNeue', 'Helvetica']})

logger = getLogger('post_process')


# https://github.com/yk-tanigawa/yt_misc_py
import yt_misc_py as yt_misc

# https://github.com/rivas-lab/decomposition/tree/master/src/py
import rivas_decomposition_py as decomposition

#npz_file='../private_data/npz/dev_PTVsNonMHC_z_nonCenter_p0001_100PCs.npz'
#phe_list_f='../public_data/phenotype_of_interest.lst'
#phe_list=['asthma']
#with open(phe_list_f) as f:
#    phe_list=[x for x in f.read().splitlines() if len(x) > 0]
#out_dir_root='../private_data/results'

######

                  
def analyze_phenotype(d, out_dir, phe, num_top_components=5):
    def get_safe_phe_name(phe):
        return ''.join([c if (c.isalnum() or c in ['_', '.'] ) else '_' for c in phe])    
    def cos2score_tbl(d, phe):
        pcs, score = d.get_topk_pcs_for_phe_with_scores_by_label(phe)
        return pd.DataFrame(collections.OrderedDict([
            ('#Rank', np.arange(d.d['n_PCs']) + 1),
            ('PC_(zero_based)', pcs),
            ('squared_cosine_score', score)
        ]))
    def cont_plot_filename(phe_dir, rank, pc, target):
        return os.path.join(
            phe_dir, 'contribution', 
            '{:02d}_PC{:02d}_{}'.format(rank, pc, target)
        )    
    
    logger.info(phe)
    phe_dir=os.path.join(out_dir, 'phenotypes', get_safe_phe_name(phe))
    if(not os.path.exists(os.path.join(phe_dir, 'contribution'))):
        os.makedirs(os.path.join(phe_dir, 'contribution'))
        
    cos2score_tbl=cos2score_tbl(d, phe)
    cos2score_tbl.to_csv(
        os.path.join(phe_dir, 'squared_cosine_scores.tsv'),
        index=False, sep='\t'
    )
    yt_misc.plot_scatter(
        d.plot_data_cos_phe_by_label(phe),
        save=os.path.join(phe_dir, 'squared_cosine_scores.pdf')
    )       
    
    for rank, pc in enumerate(d.get_topk_pcs_for_phe_by_label(phe, topk=num_top_components)):
        try:
            decomposition.plot_contribution_legend_gene(
                d, pc, topk=20, save_exts=['pdf'],
                save=cont_plot_filename(phe_dir, rank, pc, 'gene'),                
            )
        except Exception as e:  
            logger.warning('Gene contribution score plot failed! Rank = {}, PC = {}'.format(rank, pc))
            logger.warning(e)
        try:            
            decomposition.plot_contribution_legend_phe(
                d, pc, topk=20, save_exts=['pdf'],
                save=cont_plot_filename(phe_dir, rank, pc, 'phe')
            )
        except Exception as e:  
            logger.warning('Phenotype contribution score plot failed! Rank = {}, PC = {}'.format(rank, pc))
            logger.warning(e)
        matplotlib.pyplot.close('all')
    
    
    
    

    
            
def post_process_step1(
    npz_file, out_dir_root, phe_list, 
    great_topk=5000, 
    enrichr_contribution_thr=0.005, 
    gene_set_library_list=['JENSEN_Tissues']
):
    def pca_plot_filename(out_dir, pc1, pc2, ext=None):
        if(ext is None):
            ext_str=""
        else:
            ext_str='.{}'.format(ext)
        return os.path.join(out_dir, 'pca', 'pca_PC{:02d}x{:02d}{}'.format(pc1, pc2, ext_str))
    def cont_plot_filename(out_dir, pc, target, ext=None):
        if(ext is None):
            ext_str=""
        else:
            ext_str='.{}'.format(ext)        
        return os.path.join(
            out_dir, 'contribution', 
            'PC{:02d}_{}{}'.format(pc, target, ext_str)
        )    
    
    
    logger = logging.getLogger('post_process_step1')
        
    d = decomposition.decomposition(npz_file)    
    out_dir = os.path.join(
        os.path.abspath(out_dir_root), 
        os.path.basename(npz_file))
    for subdir in ['pca', 'contribution']:
        if(not os.path.exists(os.path.join(out_dir, subdir))):
            os.makedirs(os.path.join(out_dir, subdir))                        
    if False:
        yt_misc.plot(
            d.plot_data_scree(), 
            save=os.path.join(out_dir, 'scree.pdf'))                            
    
    if False:
        decomposition.generate_great_file(d, out_dir, topk=great_topk, )
    
    if True:
        for pc1 in range(d.d['n_PCs']):
            if(pc1 % 10 == 0):
                logger.info('Generating PCA plot & contribution score plots for PC{}'.format(pc1))
            if(pc1 < d.d['n_PCs'] - 1):                                            
                if(not os.path.exists(pca_plot_filename(out_dir, pc1, pc1 + 1, ext='pdf'))):
                    decomposition.plot_pca(
                        d, pc1, pc1 + 1, save_exts=['pdf'],
                        save=pca_plot_filename(out_dir, pc1, pc1 + 1),                        
                    )

            if(not os.path.exists(cont_plot_filename(out_dir, pc1, 'phe', ext='pdf'))):
                decomposition.plot_contribution_legend_phe(
                    d, pc1, topk=20, save_exts=['pdf'],
                    save=cont_plot_filename(out_dir, pc1, 'phe')
                )

            if(not os.path.exists(cont_plot_filename(out_dir, pc1, 'gene', ext='pdf'))):
                decomposition.plot_contribution_legend_gene(
                    d, pc1, topk=20, save_exts=['pdf'],
                    save=cont_plot_filename(out_dir, pc1, 'gene')
                )
                
            matplotlib.pyplot.close('all')
        logger.info('PCA plot  & contribution score plots done!')
            
    if False:
        logger.info(out_dir)
        try:
            decomposition.generate_enrichr_query(
                d, out_dir, contribution_thr=enrichr_contribution_thr)        
        except Exception as e:
            logger.warning(e)
        if False:
            decomposition.run_enrichr(d, out_dir, gene_set_library_list)        

    if True:
        for phe in phe_list:
            try:
                analyze_phenotype(d, out_dir, phe)
            except ValueError as e:
                logger.warning(e)            
                logger.warning("phenotype {} is not found in the data. Skipping".format(phe))

def main():
    logger_main = logging.getLogger('main')    
    
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=_README_
    )

    parser.add_argument('-i', metavar='i', required=True,
                        help='npz file')    
    parser.add_argument('-o', metavar='o', required=True,
                        help='out_dir_root')
    parser.add_argument('-p', metavar='p', required=True,
                        help='phenotype list file')

    args = parser.parse_args()

    npz_file=os.path.abspath(args.i)    
    out_dir_root=os.path.abspath(args.o)
    phe_list_f=os.path.abspath(args.p)
    
    logger_main.info('  npz file     : {}'.format(npz_file))
    logger_main.info('  out_dir_root : {}'.format(out_dir_root))
    logger_main.info('  phenotype    : {}'.format(phe_list_f))
    
    with open(phe_list_f) as f:
        phe_list=[x for x in f.read().splitlines() if len(x) > 0]
    
    post_process_step1(npz_file, out_dir_root, phe_list)    
    
    logger_main.info(' '.join([
        'python', 'post_process_step2.py', 
        '-p', phe_list_f,
        '-d', os.path.join(out_dir_root, os.path.basename(npz_file)),
    ]))
    
if __name__ == "__main__":
    main()     
    

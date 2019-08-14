from __future__ import print_function

_README_ = '''
-------------------------------------------------------------------------
Generate JSON files for GBE decomposition page.
-p option outputs python numpy npz file (compressed format) for python

Author: Yosuke Tanigawa (ytanigaw@stanford.edu)
Date: 2017/12/01
-------------------------------------------------------------------------
'''

import pandas as pd
import numpy as np
import os, sys, json, re, gzip, argparse, logging, collections
from datetime import datetime
from functools import reduce
from scipy.sparse import dok_matrix
from logging.config import dictConfig
import rpy2.robjects as robjects

logging_config = dict(
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
        #'level': logging.DEBUG,
        },
)
dictConfig(logging_config)


def parse_label_phe(label_phe_f):
    label_phe_df = pd.read_csv(label_phe_f, sep='\t', compression='gzip')
    label_phe_code = label_phe_df['icd'].as_matrix()   
    label_phe      = label_phe_df['Name'].map(lambda x: re.sub('_', ' ', re.sub('_/_', '/', x))).as_matrix()
    return label_phe, label_phe_code


def parse_label_var(label_var_f):
    label_var_df = pd.read_csv(label_var_f, sep='\t', compression='gzip')
    return label_var_df['signature'].map(lambda x: re.sub('_', '-', x)).as_matrix()


def read_eigen_values(tsvd_f):
    eigen_v_dict = dict([])
    with gzip.open(tsvd_f) as f:
        for line in f:
            l = line.split('\t')
            if(l[0] == '1'):
                eigen_v_dict[int(l[2])] = float(l[3])
    return np.array([eigen_v_dict[x] for x in sorted(eigen_v_dict.keys())])


def read_eigen_vectors(tsvd_f, n_PCs, n_phes, n_vars):
    eigen_phe_dok = dok_matrix((n_phes, n_PCs), dtype = np.float)
    eigen_var_dok = dok_matrix((n_vars, n_PCs), dtype = np.float)    
    with gzip.open(tsvd_f) as f:
        for line in f:
            l = line.split('\t')
            if(  l[0] == '0' and int(l[1]) < n_phes and int(l[2]) < n_PCs):
                eigen_phe_dok[int(l[1]), int(l[2])] = float(l[3])
            elif(l[0] == '2' and int(l[2]) < n_vars and int(l[1]) < n_PCs):
                eigen_var_dok[int(l[2]), int(l[1])] = float(l[3]) 
    return np.array(eigen_phe_dok.todense()), np.array(eigen_var_dok.todense())


def dok_from_tsv(tsv_f, dtype=np.float):
    logger = logging.getLogger('dok_from_tsv')
    logger.info('reading {}'.format(tsv_f))
    df = pd.read_csv(tsv_f, sep='\t', compression='gzip')
    logger.info('constructing a dok matrix of size {} x {}'.format(len(set(df.ix[:, 0])), len(set(df.ix[:, 1]))))
    dok_mat = dok_matrix(
        (len(set(df.ix[:, 0])), len(set(df.ix[:, 1]))),
        dtype = dtype
    )
    dok_mat.update(
        dict(
            zip(
                zip(
                    df.ix[:, 0].tolist(), 
                    df.ix[:, 1].tolist()
                ),
                df.ix[:, 2].tolist()
            )
        )
    )
    return dok_mat


def read_ssvd_rds(rds_file):
    '''
    Read RDS file that contains the results of ssvd
    '''
    r_funcs = collections.OrderedDict()
    for func in ['readRDS', 'as.matrix']:
        r_funcs[func] = robjects.r[func]

    res = r_funcs['readRDS'](rds_file)
    return dict(zip(
        res.names,
        [x[0,0] if x.shape == (1,1) else x 
         for x in 
         [np.array(r_funcs['as.matrix'](x)) 
          for x in list(res)]]
    ))    


def compute_factor(eigen_vec, eigen_values):
    return np.dot(eigen_vec, np.diag(eigen_values))


def compute_contribution(factor):
    return (factor ** 2) / (np.sum(factor ** 2, axis = 0).reshape((1, factor.shape[1])))


def compute_cos(factor):
    return (factor ** 2) / (np.sum(factor ** 2, axis = 1).reshape((factor.shape[0], 1)))


def compute_contribution_gene(
    var2gene_dict, label_var, contribution_var
):
    contribution_var_df = pd.DataFrame(contribution_var)    
    contribution_var_df['gene'] = [var2gene_dict[x] for x in label_var]
    contribution_gene_df = contribution_var_df.groupby('gene').sum()
    
    return contribution_gene_df.as_matrix(), np.array(contribution_gene_df.index)


def generate_data_mat_for_stacked_bar(contribution_scores, label, threshold):    
    def generate_mask_for_contribution_scores(contribution_scores, threshold):
        return np.apply_along_axis(
            lambda l: reduce(lambda x, y: x or y, l), 1, 
            np.vectorize(lambda z: z > threshold)(np.array(contribution_scores))
        )    
       
    mask = generate_mask_for_contribution_scores(contribution_scores, threshold)
    stacked_bar_label = np.hstack([np.array(label[mask]), ['others']])
    truncated_data = contribution_scores[mask, :]    
    #truncated_data[truncated_data < threshold] = 0    
    stacked_bar_data  = np.vstack([truncated_data, 1 - truncated_data.sum(axis = 0)])

    return stacked_bar_data, stacked_bar_label


def stacked_bar_per_pc(stacked_bar_data, stacked_bar_label, pc):
    sort_order = (-stacked_bar_data[:-1, pc]).argsort() 
    data  = stacked_bar_data[sort_order, pc][:50].tolist()
    label = stacked_bar_label[sort_order][:50].tolist()
    return data, label


def sparsify_contributoin_scores(contribution_mat, label, pci, threshold=0.0001):
#    mask = contribution_mat[:, pci] > (0.1 / contribution_mat.shape[0])
    mask = contribution_mat[:, pci] > threshold
    xs = np.arange(contribution_mat.shape[0])[mask]
    ys = contribution_mat[mask, pci]
    ls = label[mask]
    return xs, ys, ls, mask


def get_label_var_to_label_gene_dict(tsv_file):
    df = pd.read_csv(tsv_file, sep='\t')
    return dict(zip(df['label_var'], df['label_gene']))


def write_json_misc(
    out_dir, dataset, metadata, n_PCs, total_inertia, eigen_v, 
    label_phe, label_var, label_phe_code, label_gene,
    label_phe_stackedbar, label_gene_stackedbar,
    stackedbar_phe, stackedbar_gene):
    
#    eigen_relative = eigen_v ** 2 / np.sum(eigen_v ** 2)
    eigen_relative = eigen_v ** 2 / total_inertia

    if not os.path.exists(os.path.join(out_dir, dataset)):
        os.makedirs(os.path.join(out_dir, dataset))
        
    stackedbar_phe_json = [
        {
            'x':['PC{}'.format(pc + 1) for pc in range(n_PCs)],
            'y':stackedbar_phe[i].tolist(),
            'name': label_phe_stackedbar[i],
            'type': 'bar',
            'hoverinfo': 'none'
        } for i in range(stackedbar_phe.shape[0])
    ]
    
    stackedbar_gene_json = [
        {
            'x':['PC{}'.format(pc + 1) for pc in range(n_PCs)],
            'y':stackedbar_gene[i].tolist(),
            'name': label_gene_stackedbar[i],
            'type': 'bar',
            'hoverinfo': 'none'
        } for i in range(stackedbar_gene.shape[0])
    ]                    
        
    with open(os.path.join(out_dir, dataset, '{}_misc.json'.format(dataset)), 'w') as f:
            json.dump({
                'metadata' : metadata,
                'total_inertia' : total_inertia,
                'eigen_v'  : eigen_v.tolist(),
                'eigen_r'  : eigen_relative.tolist(),
                'label_phe': label_phe,
                'label_var': label_var,
                'label_phe_code' : label_phe_code,
                'label_phe_code_idx' : dict(zip(label_phe_code, range(len(label_phe_code)))),
                'label_gene' : label_gene,
                'label_pc':  ['PC{}'.format(pci + 1) for pci in range(n_PCs)],
                'label_pc_idx':  dict(zip(['PC{}'.format(pci + 1) for pci in range(n_PCs)], range(n_PCs))),
                'label_phe_stackedbar':  label_phe_stackedbar,
                'label_gene_stackedbar': label_gene_stackedbar,
                'stackedbar_phe':  stackedbar_phe_json,
                'stackedbar_gene': stackedbar_gene_json
                }, f) 
            

def write_json_data(
    out_dir, dataset, n_PCs, n_phes, n_vars, 
    label_phe, label_var, label_phe_code, gene2Ensembl_dict,
    factor_phe, factor_var, 
    contribution_phe, contribution_var, 
    cos_phe, cos_var,
    label_phe_stackedbar, label_gene_stackedbar,
    stackedbar_phe, stackedbar_gene,
    loading_phe, loading_var
):

    loading_sq_phe = np.array(loading_phe) ** 2
    loading_sq_var = np.array(loading_var) ** 2

    if not os.path.exists(os.path.join(out_dir, dataset)):
        os.makedirs(os.path.join(out_dir, dataset))    
    for pci in range(n_PCs):
        with open(os.path.join(out_dir, dataset, '{}_factor_phe_{}.json'.format(dataset, pci)), 'w') as f:
            json.dump(factor_phe[:, pci].tolist(), f)
        with open(os.path.join(out_dir, dataset, '{}_factor_var_{}.json'.format(dataset, pci)), 'w') as f:
            json.dump(factor_var[:, pci].tolist(), f)    

        contribution_phe_x, contribution_phe_y, contribution_phe_l, _ = sparsify_contributoin_scores(contribution_phe, label_phe, pci, 0.0001)
        contribution_var_x, contribution_var_y, contribution_var_l, _ = sparsify_contributoin_scores(contribution_var, label_var, pci, 0.001)

        with open(os.path.join(out_dir, dataset, '{}_contribution_phe_x_{}.json'.format(dataset, pci)), 'w') as f:
            json.dump(contribution_phe_x.tolist(), f)
        with open(os.path.join(out_dir, dataset, '{}_contribution_phe_y_{}.json'.format(dataset, pci)), 'w') as f:
            json.dump(contribution_phe_y.tolist(), f)
        with open(os.path.join(out_dir, dataset, '{}_contribution_phe_l_{}.json'.format(dataset, pci)), 'w') as f:
            json.dump(contribution_phe_l.tolist(), f)
        with open(os.path.join(out_dir, dataset, '{}_contribution_var_x_{}.json'.format(dataset, pci)), 'w') as f:
            json.dump(contribution_var_x.tolist(), f)
        with open(os.path.join(out_dir, dataset, '{}_contribution_var_y_{}.json'.format(dataset, pci)), 'w') as f:
            json.dump(contribution_var_y.tolist(), f)
        with open(os.path.join(out_dir, dataset, '{}_contribution_var_l_{}.json'.format(dataset, pci)), 'w') as f:
            json.dump(contribution_var_l.tolist(), f)

        loading_phe_x, loading_phe_y, loading_phe_l, _ = sparsify_contributoin_scores(loading_sq_phe, label_phe, pci, 0.0001)
        loading_var_x, loading_var_y, loading_var_l, _ = sparsify_contributoin_scores(loading_sq_var, label_var, pci, 0.001)

        with open(os.path.join(out_dir, dataset, '{}_loading_phe_x_{}.json'.format(dataset, pci)), 'w') as f:
            json.dump(loading_phe_x.tolist(), f)
        with open(os.path.join(out_dir, dataset, '{}_loading_phe_y_{}.json'.format(dataset, pci)), 'w') as f:
            json.dump(loading_phe_y.tolist(), f)
        with open(os.path.join(out_dir, dataset, '{}_loading_phe_l_{}.json'.format(dataset, pci)), 'w') as f:
            json.dump(loading_phe_l.tolist(), f)
        with open(os.path.join(out_dir, dataset, '{}_loading_var_x_{}.json'.format(dataset, pci)), 'w') as f:
            json.dump(loading_var_x.tolist(), f)
        with open(os.path.join(out_dir, dataset, '{}_loading_var_y_{}.json'.format(dataset, pci)), 'w') as f:
            json.dump(loading_var_y.tolist(), f)
        with open(os.path.join(out_dir, dataset, '{}_loading_var_l_{}.json'.format(dataset, pci)), 'w') as f:
            json.dump(loading_var_l.tolist(), f)

        bar_phe,  bar_phe_label  = stacked_bar_per_pc(stackedbar_phe,  label_phe_stackedbar,  pci)
        bar_gene, bar_gene_label = stacked_bar_per_pc(stackedbar_gene, label_gene_stackedbar, pci)
        bar_phe_code = [dict(zip(label_phe, label_phe_code))[x] for x in bar_phe_label]
        bar_gene_code = [gene2Ensembl_dict[x] for x in bar_gene_label]

        with open(os.path.join(out_dir, dataset, '{}_bar_phe_{}.json'.format(dataset, pci)), 'w') as f:
            json.dump(bar_phe, f)
        with open(os.path.join(out_dir, dataset, '{}_bar_phe_label_{}.json'.format(dataset, pci)), 'w') as f:
            json.dump(bar_phe_label, f)
        with open(os.path.join(out_dir, dataset, '{}_bar_phe_code_{}.json'.format(dataset, pci)), 'w') as f:
            json.dump(bar_phe_code, f)
        with open(os.path.join(out_dir, dataset, '{}_bar_gene_{}.json'.format(dataset, pci)), 'w') as f:
            json.dump(bar_gene, f)
        with open(os.path.join(out_dir, dataset, '{}_bar_gene_label_{}.json'.format(dataset, pci)), 'w') as f:
            json.dump(bar_gene_label, f)
        with open(os.path.join(out_dir, dataset, '{}_bar_gene_code_{}.json'.format(dataset, pci)), 'w') as f:
            json.dump(bar_gene_code, f)
                    
    for phe in range(n_phes):
        with open(os.path.join(out_dir, dataset, '{}_cos_phe_{}.json'.format(dataset, phe)), 'w') as f:
            json.dump(cos_phe[phe, :].tolist(), f)    
        with open(os.path.join(out_dir, dataset, '{}_loading_phe_{}.json'.format(dataset, phe)), 'w') as f:
            json.dump(loading_sq_phe[phe, :].tolist(), f)

    for var in range(n_vars):
        with open(os.path.join(out_dir, dataset, '{}_cos_var_{}.json'.format(dataset, var)), 'w') as f:
            json.dump(cos_var[var, :].tolist(), f)       
        with open(os.path.join(out_dir, dataset, '{}_loading_var_{}.json'.format(dataset, var)), 'w') as f:
            json.dump(loading_sq_var[var, :].tolist(), f)
            

def data_prep_for_gbe_main(
    in_data_dir, out_dir, dataset, variant_tsv, 
    stackedbar_threshold = 0.1, json_out = False, python_out = False, loading=False
):
    
    logger_main = logging.getLogger('data_prep_for_gbe_main') 
    
    # filenames
    label_phe_f     = os.path.join(in_data_dir, dataset, 'ap_icd_idx.tsv.gz')
    label_var_f     = os.path.join(in_data_dir, dataset, 'ap_variant_idx.tsv.gz')
    tsvd_f          = os.path.join(in_data_dir, dataset, 'ap_icd_var_tsvd.tsv.gz')
    loading_phe_f   = os.path.join(in_data_dir, dataset, 'ap_icd_svd_cor1.tsv.gz')
    loading_var_f   = os.path.join(in_data_dir, dataset, 'ap_icd_svd_cor2.tsv.gz')
    ssvd_f          = os.path.join(in_data_dir, dataset, 'ssvd.rds')    
    
    meta_f          = os.path.join(in_data_dir, dataset, 'metadata.txt')
    total_inertia_f = os.path.join(in_data_dir, dataset, 'total_inertia.txt')

    # meta data
    with open(meta_f) as f:
        metadata = f.read().splitlines()
    metadata.append('data conversion script has started on {}'.format(str(datetime.now())))
    
    # total inertia
    with open(total_inertia_f) as f:
        total_inertia = float(f.read().splitlines()[0])

    # read dict to convert label_var to label_gene 
    variant_df = pd.read_csv(variant_tsv, sep='\t', compression='gzip')
    var2gene_dict = dict(zip(variant_df['label_var'], variant_df['label_gene']))    
    gene2Ensembl_dict = dict(zip(variant_df['label_gene'], variant_df['Gene']))    
        
    # read the data (1) labels and eigen values
    logger_main.info('reading labels and eigen values ...')    
    label_phe, label_phe_code = parse_label_phe(label_phe_f)
    label_var_unsorted = parse_label_var(label_var_f)    
    
    # sort variant labels
    label_var_argsort = np.argsort(
        [int(x.split('-')[1]) + 1000000000 * int(x.split('-')[0]) for x in label_var_unsorted]
    )    
    label_var = label_var_unsorted[label_var_argsort]

    if(os.path.isfile(tsvd_f)):
        # results from TSVD        
        is_ssvd = False
        eigen_v = read_eigen_values(tsvd_f)
        # get the number of PCs, variants, and phenotyps
        n_phes = len(label_phe)
        n_vars = len(label_var_unsorted)
        n_PCs = len(eigen_v)

        # read the data (2) eigen vectors
        logger_main.info('reading eigen vectors ...')        
        eigen_phe, eigen_var_unsorted = read_eigen_vectors(tsvd_f, n_PCs, n_phes, n_vars)
        eigen_var = eigen_var_unsorted[label_var_argsort, :]

        # read the data (3) loadings (correlation of phenotype/variant vector and PCs)
        if(loading and os.path.isfile(loading_phe_f) and os.path.isfile(loading_var_f)):
            logger_main.info('reading phenotype loading (correlation) ...')    
            loading_phe = dok_from_tsv(loading_phe_f).todense()
            logger_main.info('reading variant loading (correlation) ...')    
            loading_var = dok_from_tsv(loading_var_f).todense()
        else:
            loading_phe = np.ones((n_phes, n_PCs))
            loading_var = np.ones((n_vars, n_PCs))             
    else:
        # results from sparse SVD
        is_ssvd = True
        ssvd_res = read_ssvd_rds(ssvd_f)      
        
        # needs clean up...
        loading_phe = np.ones(ssvd_res['v'].shape)
        loading_var = np.ones(ssvd_res['u'].shape)             
        
        
    # convert to factor scores
    logger_main.info('computing scores ...') 
    if(is_ssvd):
        factor_phe = np.dot(ssvd_res['v'], np.transpose(ssvd_res['d']))
        factor_var = np.dot(ssvd_res['u'], ssvd_res['d'])
    else:
        factor_phe = compute_factor(eigen_phe, eigen_v)
        factor_var = compute_factor(eigen_var, eigen_v)

    # compute cosine scores & contribution scores
    contribution_phe = compute_contribution(factor_phe)
    contribution_var = compute_contribution(factor_var)
    cos_phe = compute_cos(factor_phe)
    cos_var = compute_cos(factor_var)
    
    contribution_gene, label_gene = compute_contribution_gene(
        var2gene_dict, label_var, contribution_var
    )    
    
    # compute data for stacked bar plots
    stackedbar_phe, label_phe_stackedbar = generate_data_mat_for_stacked_bar(
        contribution_phe, label_phe, stackedbar_threshold
    )
    stackedbar_gene, label_gene_stackedbar = generate_data_mat_for_stacked_bar(
        contribution_gene, label_gene, stackedbar_threshold
    )
    
    if(python_out):
        out_file = os.path.join(out_dir, '{}.npz'.format(dataset))
        # write to a python npz file
        logger_main.info('writing to npz file: {} ...'.format(out_file))
        if(is_ssvd):
            np.savez_compressed(
                out_file, 

		ssvd_u = ssvd_res['u'],
		ssvd_v = ssvd_res['v'],
		ssvd_d = ssvd_res['d'],
		ssvd_center = ssvd_res['center'],
		ssvd_scale = ssvd_res['scale'],
		ssvd_lambda = ssvd_res['lambda'],
		ssvd_iter = ssvd_res['iter'],
		ssvd_n = ssvd_res['n'],
		ssvd_alpha = ssvd_res['alpha'],

		total_inertia     = np.array([total_inertia]),
            	label_phe         = np.array(label_phe), 
            	label_var         = np.array(label_var),
            	label_phe_code    = np.array(label_phe_code),
            	label_gene        = np.array(label_gene),        
            	label_phe_stackedbar  = np.array(label_phe_stackedbar),
            	label_gene_stackedbar = np.array(label_gene_stackedbar),        
            	factor_phe        = np.array(factor_phe), 
            	factor_var        = np.array(factor_var),
            	contribution_phe  = np.array(contribution_phe),
            	contribution_var  = np.array(contribution_var),
            	contribution_gene = np.array(contribution_gene),
            	cos_phe           = np.array(cos_phe),
            	cos_var           = np.array(cos_var),
            	stackedbar_phe    = np.array(stackedbar_phe),
            	stackedbar_gene   = np.array(stackedbar_gene),
            	metadata          = np.array(metadata)
            )  
        else:
            np.savez_compressed(
                out_file, 
                eigen_v           = np.array(eigen_v),             
                eigen_phe         = np.array(eigen_phe),
                eigen_var         = np.array(eigen_var),
                loading_phe       = np.array(loading_phe),
                loading_var       = np.array(loading_var),                

		total_inertia     = np.array([total_inertia]),
            	label_phe         = np.array(label_phe), 
            	label_var         = np.array(label_var),
            	label_phe_code    = np.array(label_phe_code),
            	label_gene        = np.array(label_gene),        
            	label_phe_stackedbar  = np.array(label_phe_stackedbar),
            	label_gene_stackedbar = np.array(label_gene_stackedbar),        
            	factor_phe        = np.array(factor_phe), 
            	factor_var        = np.array(factor_var),
            	contribution_phe  = np.array(contribution_phe),
            	contribution_var  = np.array(contribution_var),
            	contribution_gene = np.array(contribution_gene),
            	cos_phe           = np.array(cos_phe),
            	cos_var           = np.array(cos_var),
            	stackedbar_phe    = np.array(stackedbar_phe),
            	stackedbar_gene   = np.array(stackedbar_gene),
            	metadata          = np.array(metadata)
            )
        
    if(json_out and (not is_ssvd)):
        # write to a JSON file
        logger_main.info('writing to JSON files ...')            
        write_json_misc(
            out_dir, dataset, metadata, n_PCs, total_inertia, eigen_v, 
            label_phe.tolist(), 
            label_var.tolist(), 
            label_phe_code.tolist(), 
            label_gene.tolist(),
            label_phe_stackedbar.tolist(),
            label_gene_stackedbar.tolist(),
            stackedbar_phe, 
            stackedbar_gene
        )

        # write to small json files
        write_json_data(
            out_dir, dataset, n_PCs, n_phes, n_vars, 
            label_phe, label_var, label_phe_code, gene2Ensembl_dict, 
            factor_phe, factor_var, contribution_phe, contribution_var, 
            cos_phe, cos_var, 
            label_phe_stackedbar, label_gene_stackedbar,
            stackedbar_phe, stackedbar_gene,
            loading_phe, loading_var
        )   
    

def main():
    logger_main = logging.getLogger('main')    
    
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=_README_
    )
    
    parser.add_argument('-i', metavar='i', required=True,
                        help='input dir')
    parser.add_argument('-o', metavar='o', required=True,
                        help='output dir')
    parser.add_argument('-n', metavar='n', required=True,
                        help='name of the dataset')
    parser.add_argument('-a', metavar='a', required=True,
                        help='variant annotation file (example: /home/ytanigaw/repos/rivas-lab/decomposition/private_data/variant_and_gene_labels.tsv.gz)')
    parser.add_argument('-t', metavar='t', default=0.1, type=float,
                        help='stacked bar plot cutoff threshold')        
    parser.add_argument('-p', action='store_true',
                        help='prepare npz file for python')
    parser.add_argument('-j', action='store_true',
                        help='prepare JSON files')
    parser.add_argument('-l', action='store_true',
                        help='loading scores')


    args = parser.parse_args()

    logger_main.info('  in       : {}'.format(args.i))
    logger_main.info('  out      : {}'.format(args.o))
    logger_main.info('  data     : {}'.format(args.n))
    logger_main.info('  variant  : {}'.format(args.a))
    logger_main.info('  bar plot : {}'.format(args.t))
    logger_main.info('  JSON     : {}'.format(args.j))
    logger_main.info('  python   : {}'.format(args.p))
    logger_main.info('  loading  : {}'.format(args.l))
    
    if(args.j or args.p): 
        data_prep_for_gbe_main(
            in_data_dir = args.i,
            out_dir     = args.o,
            dataset     = args.n,
            variant_tsv = args.a,
            stackedbar_threshold = args.t,
            json_out   = args.j,
            python_out = args.p,
            loading     = args.l 
        )
    
        
if __name__ == "__main__":
    main()            
            

# -*- coding: utf-8 -*-

from functools import reduce
import numpy as np
import pandas as pd
import itertools as it
import os
import re
import textwrap
import logging
import matplotlib

class decomposition:
    """a class to store & query npz file"""    
    # https://matplotlib.org/api/collections_api.html#matplotlib.collections.Collection.set_hatch
    def __init__(
        self, 
        npz_file,
        colors = matplotlib.pyplot.cm.tab20c((np.arange(13)).astype(int)),
        hatches = ['', '/', '\\', '|', '-', '+', 'x', 'o', 'O', '.', '*', '++', 'xx', 'oo', 'OO', '..', '**'],
    ):
        """read npz file and store it into a python dict and store"""
        self.logger = logging.getLogger('data_load_from_npz') 
        data = dict([])
        self.logger.info('reading data from {}'.format(npz_file))
        with np.load(npz_file) as f:
            for k in f.keys():
                data[k] = f[k]
        if('eigen_v' in data):
            data['n_PCs'] = len(data['eigen_v'])
        else: # for ssvd results 
            data['n_PCs'] = data['factor_phe'].shape[1]
        data['n_phes']  = data['contribution_phe'].shape[0]
        data['n_vars']  = data['contribution_var'].shape[0]  
        data['n_genes']  = data['contribution_gene'].shape[0]    
        if('eigen_v' in data):        
            data['variance_explained'] = (data['eigen_v'] ** 2) / data['total_inertia']
        self.d = data

        # dict to look up index given a label
        self.dict_inv_label = { 
            target: dict(zip(self.d['label_{}'.format(target)], 
                             np.arange(len(self.d['label_{}'.format(target)])))) 
            for target in ['phe', 'var', 'gene']} 
        
        self.colors = colors
        self.hatches = hatches
        self.color_len = len(self.colors)
        self.hatch_len = len(self.hatches)                
        
    def get_index_generic(self, target, label):
        return self.dict_inv_label[target][label]
    def get_phe_index(self, label):
        return self.get_index_generic('phe', label)
    def get_var_index(self, label):
        return self.get_index_generic('var', label)
    def get_gene_index(self, label):
        return self.get_index_generic('gene', label)
    def get_phe_str(self, index):
        return self.d['label_phe'][index]
    def get_cos_phe_by_index(self, index):
        return self.d['cos_phe'][index]
    def get_cos_phe_by_label(self, label):
        return self.d['cos_phe'][self.get_phe_index(label)]   
    def get_topk_pcs_for_phe_by_index(self, index, topk=None):
        if topk is None:
            topk = self.d['n_PCs']
        return np.argsort(-self.d['cos_phe'][index, :])[:topk]
    def get_topk_pcs_for_phe_by_label(self, label, topk=None):
        return self.get_topk_pcs_for_phe_by_index(self.get_phe_index(label), topk=topk)
    def get_topk_pcs_for_phe_with_scores_by_index(self, index, topk=None):
        sort_order = self.get_topk_pcs_for_phe_by_index(index, topk)
        labels = self.d['cos_phe'][index, sort_order]
        return sort_order, labels
    def get_topk_pcs_for_phe_with_scores_by_label(self, label, topk=None):
        return self.get_topk_pcs_for_phe_with_scores_by_index(
            self.get_phe_index(label), topk=topk)
    def get_kth_pc_for_phe_by_index(self, index, k):
        return self.get_topk_pcs_for_phe_by_index(index, topk=k)[-1]
    def get_kth_pc_for_phe_by_label(self, label, k):
        return self.get_kth_pc_for_phe_by_index(self.get_phe_index(label), k)
    def get_contribution_score_by_label_generic(self, target, pc_0index, label):
        return self.d['contribution_{}'.format(target)][self.get_index_generic(target, label), pc_0index]
    
    def plot_data_cos_phe_by_index(self, index):
        return {
            'x': np.arange(self.d['n_PCs']) + 1,
            'y': self.get_cos_phe_by_index(index),
            'title': 'Squared cosine score for {}'.format(self.get_phe_str(index)),
            'xlabel': 'component',
            'ylabel': 'Squared cosine score'
        }
    def plot_data_cos_phe_by_label(self, label):
        return self.plot_data_cos_phe_by_index(self.get_phe_index(label))
    
    def get_stacked_bar_color(self, target, labels):
        return np.array([
            self.colors[self.get_index_generic(target, label) % self.color_len]
            if label != 'others' 
            else '0.8'
            for label in labels
        ])
    
    def get_stacked_bar_hatch(self, target, labels):
        return np.array([
            self.hatches[self.get_index_generic(target, label) % self.hatch_len]
            if label != 'others' 
            else ''
            for label in labels
        ])   
    
    def plot_data_contribution_generic(
        self, target, pc_index, topk=None, sort=True, append_others = False,
        contribution_thr = None, contribution_cumsum_thr = None,
    ):
        assert(topk is None or sort)
        if target is 'phe':
            topk_default = self.d['n_phes']
            contribution_key = 'contribution_phe'
            xlabel = 'Phenotype'
            label_key = 'label_phe'
            label_key_sub = 'label_phe_code'
        elif target is 'gene':
            topk_default = self.d['n_genes']
            contribution_key = 'contribution_gene'            
            xlabel = 'Gene'
            label_key = 'label_gene'            
            label_key_sub = 'label_gene'
        elif target is 'var':
            topk_default = self.d['n_vars']
            contribution_key = 'contribution_var'            
            xlabel = 'Variant'
            label_key = 'label_var'            
            label_key_sub = 'label_var'
        else:
            raise(ValueError("target ({}) is not in ['phe', 'var', 'gene']".format(target)))
        if topk is None:
            topk = topk_default
        if (contribution_thr is not None) or (contribution_cumsum_thr is not None):
            sort = True
        if sort:
            sort_order = np.argsort(-self.d[contribution_key][:, pc_index])
        else:
            sort_order = np.arange(topk)
        if contribution_thr is not None:
            topk = np.sum(
                self.d[contribution_key][:, pc_index] >= contribution_thr)            
        if contribution_cumsum_thr is not None:
            topk = 1 + np.sum(
                np.cumsum(self.d[contribution_key][sort_order, pc_index]) < contribution_cumsum_thr)            
        pd_x = np.arange(topk) + 1
        pd_y = self.d[contribution_key][sort_order[:topk], pc_index]
        pd_title = '{} contribution score for PC{}'.format(xlabel, pc_index + 1)
        pd_xlabel = xlabel
        pd_ylabel = '{} contribution score'.format(xlabel)
        pd_xticklabels = self.d[label_key][sort_order[:topk]]
        pd_xticklabels_sub = self.d[label_key_sub][sort_order[:topk]]
        
        if(append_others):
            pd_y = np.append(pd_y, 1 - np.sum(pd_y))
            pd_xticklabels = np.append(pd_xticklabels, 'others')
            pd_xticklabels_sub = np.append(pd_xticklabels_sub, 'others')            
            pd_x = np.arange(topk + 1) + 1            
            
        stacked_bar_color = self.get_stacked_bar_color(target, pd_xticklabels)
        stacked_bar_hatch = self.get_stacked_bar_hatch(target, pd_xticklabels)

        return {
            'x': pd_x,
            'y': pd_y,
            'title': pd_title,
            'xlabel': pd_xlabel,
            'ylabel': pd_ylabel,
            'xticklabels': pd_xticklabels,
            'xticklabels_sub': pd_xticklabels_sub,
            'stacked_bar_color': stacked_bar_color,
            'stacked_bar_hatch': stacked_bar_hatch,
        }            
    def plot_data_contribution_phe(
        self, pc_index, topk=None, sort=True, append_others = False,
        contribution_thr = None, contribution_cumsum_thr = None,):
        return self.plot_data_contribution_generic(
            'phe', pc_index, topk, sort, append_others, 
            contribution_thr, contribution_cumsum_thr)
    def plot_data_contribution_var(
        self, pc_index, topk=None, sort=True, append_others = False,
        contribution_thr = None, contribution_cumsum_thr = None,):
        return self.plot_data_contribution_generic(
            'var', pc_index, topk, sort, append_others, 
            contribution_thr, contribution_cumsum_thr)
    def plot_data_contribution_gene(
        self, pc_index, topk=None, sort=True, append_others = False,
        contribution_thr = None, contribution_cumsum_thr = None,):
        return self.plot_data_contribution_generic(
            'gene', pc_index, topk, sort, append_others, 
            contribution_thr, contribution_cumsum_thr)
    def plot_data_contribution_phe_list(self, phe_list, pc_index=0):
        pd_title = 'Phenotype contribution score'    
        phe_set = set(phe_list)
        phe_slice = np.array([x in phe_set for x in self.d['label_phe']])
        labels = self.d['label_phe'][phe_slice]
        return {
            'x': np.arange(len(phe_set)),
            'y': self.d['contribution_phe'][phe_slice, pc_index],
            'title': pd_title,
            'xlabel': '',
            'ylabel': pd_title,
            'xticklabels': labels,
            'xticklabels_sub': self.d['label_phe_code'][phe_slice],
            'stacked_bar_color': self.get_stacked_bar_color('phe', labels),
            'stacked_bar_hatch': self.get_stacked_bar_hatch('phe', labels),
        }            
    def plot_data_pca_generic(self, target, pc_index1, pc_index2):
        assert ((target == 'phe') or (target == 'var'))
        if target == 'phe':            
            fig_title = 'Phenotype PCA'
        elif target == 'var':
            fig_title = 'Variant PCA'
        else:
            raise(ValueError("target ({}) is not in ['phe', 'var']".format(target)))
        return {
            'x': self.d['factor_{}'.format(target)][:, pc_index1],
            'y': self.d['factor_{}'.format(target)][:, pc_index2],            
            'title': fig_title,
            'xlabel': 'Component {}'.format(pc_index1 + 1),
            'ylabel': 'Component {}'.format(pc_index2 + 1)            
        }
    def plot_data_pca_phe(self, pc_index1, pc_index2):
        return self.plot_data_pca_generic('phe', pc_index1, pc_index2)
    def plot_data_pca_var(self, pc_index1, pc_index2):
        return self.plot_data_pca_generic('var', pc_index1, pc_index2)
    def plot_data_scree(self):
        return {
            'x': np.arange(self.d['n_PCs']) + 1,
            'y': self.d['variance_explained'],
            'title': 'Scree plot',
            'xlabel': 'Components',
            'ylabel': 'Variance explained score',
        }
    def bed_data_contribution_var(
        self, pc_index, topk=None, sort=True, append_others = False,
        contribution_thr = None, contribution_cumsum_thr = None,
    ):
        plotd = self.plot_data_contribution_generic(
            'var', pc_index, topk, sort, append_others, 
            contribution_thr, contribution_cumsum_thr)
        topk_var_chrs, topk_var_poss = zip(*[x.split('-') for x in plotd['xticklabels']])
        return pd.DataFrame({
            'chrom': ['chr{}'.format(x) for x in topk_var_chrs],
            'chromStart': topk_var_poss,
            'chromEnd': [str(int(x) + 1) for x in topk_var_poss],
            'name': ['{:.4e}'.format(x) for x in plotd['y']]
        })[['chrom', 'chromStart', 'chromEnd', 'name']]
    def text_data_contribution_gene(
        self, pc_index, topk=None, sort=True, append_others = False,
        contribution_thr = None, contribution_cumsum_thr = None,
    ):
        plotd = self.plot_data_contribution_generic(
            'gene', pc_index, topk, sort, append_others, 
            contribution_thr, contribution_cumsum_thr)
                
        mask_by_gene_name = np.array([len(x.split(',')) == 1 for x in plotd['xticklabels']])
        labels = plotd['xticklabels'][mask_by_gene_name]
        scores = plotd['y'][mask_by_gene_name]
        query_str = '\n'.join(['{},{:.4f}'.format(x[0], x[1]) 
                               for x in zip(labels, scores)])            
        return query_str, len(labels), np.sum(scores)
    def text_contribution_generic(        
        self, target, pc_index, topk=None, sort=True, append_others = False,
        contribution_thr = None, contribution_cumsum_thr = None, 
        labels='xticklabels', values='y', value_multiplier = 100, 
        string_format='{} ({:.2f}%)', delim=', '
    ):
        plot_d = self.plot_data_contribution_generic(
            target, pc_index, topk, sort, append_others, 
            contribution_thr, contribution_cumsum_thr
        )
        return delim.join([string_format.format(l, y*value_multiplier) for l, y in zip(plot_d[labels], plot_d[values])])
    def text_description_phe_by_index(
        self, index, num_of_top_pcs, target, 
        pcs_delim = '; ', pcs_string_format='[PC{}, {:.2f}]: {}',
        topk=None, sort=True, append_others = False,
        contribution_thr = None, contribution_cumsum_thr = None, 
        labels='xticklabels', values='y', value_multiplier = 100, 
        string_format='{} ({:.2f}%)', delim=', ', 
    ):
        pcs, squared_cos_scores = self.get_topk_pcs_for_phe_with_scores_by_index(index, topk=num_of_top_pcs)
        return 'Top {} components for {} is: '.format(
            num_of_top_pcs, self.get_phe_str(index)
        ) + pcs_delim.join([
            pcs_string_format.format(
                pc_index + 1, squared_cos_score, 
                self.text_contribution_generic(
                    target = target, pc_index = pc_index, 
                    topk = topk, sort = sort, append_others = append_others,  
                    contribution_thr = contribution_thr, 
                    contribution_cumsum_thr = contribution_cumsum_thr, 
                    labels = labels, values = values, 
                    value_multiplier = value_multiplier, 
                    string_format=string_format, delim=delim)) 
            for pc_index, squared_cos_score 
            in zip(pcs, squared_cos_scores)
        ])
    def text_description_phe_by_label(
        self, label, num_of_top_pcs, target, 
        pcs_delim = '; ', pcs_string_format='[PC{}, {:.2f}]: {}',
        topk=None, sort=True, append_others = False,
        contribution_thr = None, contribution_cumsum_thr = None, 
        labels='xticklabels', values='y', value_multiplier = 100, 
        string_format='{} ({:.2f}%)', delim=', ', 
    ):
        return self.text_description_phe_by_index(
            self.get_phe_index(label), num_of_top_pcs, target, 
            pcs_delim, pcs_string_format,
            topk, sort, append_others,
            contribution_thr, contribution_cumsum_thr, 
            labels, values, value_multiplier, 
            string_format, delim) 
    def get_biplot_arrow_by_phenotypes(self, pc_list, phenotypes):
        matrix_key = 'eigen_phe'
        return np.array([
            [self.d['eigen_phe'][phe_idx][pc_idx] for pc_idx in pc_list]
            for phe_idx in [self.get_phe_index(x) for x in phenotypes]
        ])
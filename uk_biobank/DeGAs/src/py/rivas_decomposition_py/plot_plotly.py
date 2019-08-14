import os, logging, collections
import numpy as np
import pandas as pd

from logging.config import dictConfig
from logging import getLogger

import plotly
import plotly.plotly as py
import plotly.graph_objs as go


from .label_txt_formatter import label_txt_formatter_gene


dictConfig(dict(
    version = 1,
    formatters = {'f': {'format': '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'}},
    handlers = {
        'h': {'class': 'logging.StreamHandler','formatter': 'f',
              'level': logging.DEBUG}},
    root = {'handlers': ['h'], 'level': logging.DEBUG,},
))

logger = getLogger('plot_plotly')


def plotly_generic(d, pc_index1, pc_index2, score_type, target, var2gene_dict=None, var_percentile_filter=99):    
    data_mat_key = '{}_{}'.format(score_type, target)
    xs = d.d[data_mat_key][:, pc_index1]
    ys = d.d[data_mat_key][:, pc_index2]
    if(target == 'phe'):        
        texts = ['{} ({})'.format(x[0], x[1]) for x in 
                 zip(d.d['label_phe'], d.d['label_phe_code'])]
    if(target == 'var'):
        dists = xs ** 2 + ys ** 2
        dists_filter = (dists >= np.percentile(dists, var_percentile_filter))
        xs = xs[dists_filter]
        ys = ys[dists_filter] 
        texts_raw = d.d['label_var'][dists_filter]
        if(var2gene_dict is not None):
            texts = [
                '{} ({})'.format(x[0], x[1]) for x in
                zip(
                    texts_raw, [
                        label_txt_formatter_gene(var2gene_dict[y]) if y in var2gene_dict else y 
                        for y in texts_raw
                    ]
                )
            ]
        else:
            texts = texts_raw
        
    return go.Figure(
        data=[
            go.Scatter(
                x = xs, 
                y = ys,
                text = texts,
                mode = 'markers',
                marker=dict(size=5, line=dict( width=0.5 ), opacity=1)
            )
        ], 
        layout=go.Layout(
            hovermode = 'closest',
            margin=dict(l=0, r=0, b=0, t=0),            
            yaxis = dict(
                  scaleanchor = "x",
            ),            
        )
    )

def plotly_factor_phe(d, pc_index1, pc_index2):
    return plotly_generic(d, pc_index1, pc_index2, 'factor', 'phe')
def plotly_factor_var(d, pc_index1, pc_index2, var2gene_dict=None, var_percentile_filter=99):
    return plotly_generic(d, pc_index1, pc_index2, 'factor', 'var', var2gene_dict, var_percentile_filter)
def plotly_eigen_phe(d, pc_index1, pc_index2):
    return plotly_generic(d, pc_index1, pc_index2, 'eigen', 'phe')
def plotly_eigen_var(d, pc_index1, pc_index2, var2gene_dict=None, var_percentile_filter=99):
    return plotly_generic(d, pc_index1, pc_index2, 'eigen', 'var', var2gene_dict, var_percentile_filter)
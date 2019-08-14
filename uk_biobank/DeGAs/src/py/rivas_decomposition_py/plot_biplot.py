import os, logging, collections
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

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

logger = getLogger('plot_biplot')


def plot_biplot(
    d, pc_index1, pc_index2, 
    biplot_phes, 
    variants_of_interest=None,
    arrow_max_scale=1.0,
    figsize=(12,12), 
    flip_xaxis=False, flip_yaxis=False,    
    save=None
):
    """scatter plot of variants with phenotype annotation (arrows)"""
        
    if (variants_of_interest is None):
        variants_set = set([])
    else:
        variants_set = set(variants_of_interest)

    # prepare data
    plot_d = d.plot_data_pca_var(pc_index1, pc_index2)
    plot_d['color'] = np.where([x in variants_set for x in d.d['label_var']], 'red', 'blue')  
    if(biplot_phes is not None and len(biplot_phes) > 0):
        biplot_arrow_2d = d.get_biplot_arrow_by_phenotypes([pc_index1, pc_index2], biplot_phes)    
    else:
        biplot_arrow_2d = np.array([])
        
    # prepare fig grid
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 1)
    fig_axs = [fig.add_subplot(sp) for sp in gs]
    

    # configure main and sub plots
    ax_main = fig_axs[0]
    ax_main.set_aspect('equal')   
    ax_main.set_adjustable('box')    
    ax_sub = ax_main.twinx().twiny()
    ax_sub.set_aspect('equal')  
    ax_sub.set_adjustable('datalim')
     
    # axis range
    scatter_max = 1.1 * np.max([plot_d['x'], -plot_d['x'], plot_d['y'], -plot_d['y']])    
    arrow_max = arrow_max_scale * np.max([biplot_arrow_2d, -biplot_arrow_2d])
        
    if(flip_xaxis):
        ax_main.set_xlim([scatter_max, -scatter_max])
        ax_sub.set_xlim([arrow_max, -arrow_max])        
    else:
        ax_main.set_xlim([-scatter_max, scatter_max])
        ax_sub.set_xlim([-arrow_max, arrow_max])        
    if(flip_yaxis):
        ax_main.set_ylim([scatter_max, -scatter_max])
        ax_sub.set_ylim([arrow_max, -arrow_max])        
    else:
        ax_main.set_ylim([-scatter_max, scatter_max])
        ax_sub.set_ylim([-arrow_max, arrow_max])        
    
        
    # plot arrows
    for ar in biplot_arrow_2d:
        if((ar[1]) ** 2 < (ar[0]) ** 2 ):            
            ax_sub.plot(
                np.linspace(-scatter_max, scatter_max, 1000),
                np.linspace(-scatter_max * ar[1] / ar[0], scatter_max * ar[1] / ar[0], 1000),
                linestyle='dashed',
                color='0.8'                
            )        
        else:
            ax_sub.plot(
                np.linspace(-scatter_max * ar[0] / ar[1], scatter_max * ar[0] / ar[1], 1000),
                np.linspace(-scatter_max, scatter_max, 1000),                
                linestyle='dashed',
                color='0.8'
            )            
        ax_sub.annotate(
            '', 
            xy=(ar[0], ar[1]), 
            xytext=(0, 0),
            arrowprops=dict(facecolor='red', shrinkA=0,shrinkB=0),
        )        

    # scatter plot    
    ax_main.scatter(
        plot_d['x'], plot_d['y'], 
        color=(plot_d['color'] if 'color' in plot_d else None),
        marker='x', s=(15**2)
    )        
    
    gs.tight_layout(fig, rect=[0, 0, 1, 1]) 

    if (biplot_phes is not None and len(biplot_phes) > 0):
        # construct a data frame of arrow coordinate for manual annotation
        df = pd.DataFrame(collections.OrderedDict([
                ('phe', biplot_phes),
                ('x', biplot_arrow_2d[:, 0]),
                ('y', biplot_arrow_2d[:, 1]),       
        ]))
        df['r'] = (df['y'] ** 2 + df['x'] ** 2) ** 0.5
        df['slope'] = df['y'] / df['x']
        df = df.sort_values(by='slope')
    else:
        df = None
    
    # save to file
    if save is not None:
        for ext in ['pdf', 'png']:
            tmp_save_name='{}.{}'.format(save, ext)            
            logger.info('saving the image to {}'.format(tmp_save_name))
            if(not os.path.exists(os.path.dirname(tmp_save_name))):
                os.makedirs(os.path.dirname(tmp_save_name))                    
            fig.savefig(tmp_save_name, bbox_inches="tight", pad_inches=0.0)
        if(df is not None):
            tmp_save_name='{}.tsv'.format(save)
            logger.info('saving the table to {}'.format(tmp_save_name))
            df.to_csv(tmp_save_name, sep='\t', index=False)
            
    return df
import os, logging, collections
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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

logger = getLogger('plot_pca')

def plot_pca_phe(
    d, pc_index1, pc_index2, 
    figsize=(12,12), 
    flip_xaxis=False, flip_yaxis=False,
    save=None, save_exts=['pdf', 'png'],
):
    """scatter plot of variants with phenotype annotation (arrows)"""
    # prepare data
    plot_d_phe = d.plot_data_pca_phe(pc_index1, pc_index2)        
    plot_ds=[plot_d_phe]
    
    # prepare fig grid
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 1)
    fig_axs = [fig.add_subplot(sp) for sp in gs]
    
    for ax, plot_d in zip(fig_axs, plot_ds):
        ax.set_aspect('equal')
        ax.scatter(
            plot_d['x'], plot_d['y'], 
            marker='x', s=(15**2),
            color='blue'
        )    
        ax_max = 1.1 * np.max([plot_d['x'], -plot_d['x'], plot_d['y'], -plot_d['y']])
        if(flip_xaxis):
            ax.set_xlim([ax_max, -ax_max])            
        else:
            ax.set_xlim([-ax_max, ax_max])
        if(flip_yaxis):
            ax.set_ylim([ax_max, -ax_max])
        else:
            ax.set_ylim([-ax_max, ax_max])    
        

    gs.tight_layout(fig, rect=[0, 0, 1, 1]) 

    
    # save to file
    if save is not None:
        for ext in save_exts:
            tmp_save_file='{}.{}'.format(save, ext)
            if(not os.path.exists(os.path.dirname(tmp_save_file))):
                os.makedirs(os.path.dirname(tmp_save_file))            
            logger.info('saving the image to {}'.format(tmp_save_file))
            fig.savefig(tmp_save_file, bbox_inches="tight", pad_inches=0.0)
            

def plot_pca(
    d, pc_index1, pc_index2, 
    figsize=(12,6), 
    flip_xaxis=False, flip_yaxis=False,    
    save=None, save_exts=['pdf', 'png'],
):
    """scatter plot of variants with phenotype annotation (arrows)"""
    # prepare data
    plot_d_phe = d.plot_data_pca_phe(pc_index1, pc_index2)        
    plot_d_var = d.plot_data_pca_var(pc_index1, pc_index2)
    plot_ds = [plot_d_phe, plot_d_var]
    
        
    # prepare fig grid
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 2)
    fig_axs = [fig.add_subplot(sp) for sp in gs]
    
    for ax, plot_d in zip(fig_axs, plot_ds):
        ax.set_aspect('equal')           
        ax.scatter(
            plot_d['x'], plot_d['y'], 
            marker='x', s=(15**2),
            color='blue',
        )    
        ax.set_title(plot_d['title'])
        ax.set_xlabel(plot_d['xlabel'])
        ax.set_ylabel(plot_d['ylabel'])   
        
        ax_max = 1.1 * np.max([plot_d['x'], -plot_d['x'], plot_d['y'], -plot_d['y']])
        if(flip_xaxis):
            ax.set_xlim([ax_max, -ax_max])            
        else:
            ax.set_xlim([-ax_max, ax_max])
        if(flip_yaxis):
            ax.set_ylim([ax_max, -ax_max])
        else:
            ax.set_ylim([-ax_max, ax_max])    
        

    gs.tight_layout(fig, rect=[0, 0, 1, 1]) 

    
    # save to file
    if save is not None:
        for ext in save_exts:
            tmp_save_file='{}.{}'.format(save, ext)
            if(not os.path.exists(os.path.dirname(tmp_save_file))):
                os.makedirs(os.path.dirname(tmp_save_file))            
            logger.info('saving the image to {}'.format(tmp_save_file))
            fig.savefig(tmp_save_file, bbox_inches="tight", pad_inches=0.0)
            
import os, logging
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

from logging.config import dictConfig
from logging import getLogger

from functools import reduce

dictConfig(dict(
    version = 1,
    formatters = {'f': {'format': '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'}},
    handlers = {
        'h': {'class': 'logging.StreamHandler','formatter': 'f',
              'level': logging.DEBUG}},
    root = {'handlers': ['h'], 'level': logging.DEBUG,},
))

logger = getLogger('plot_contribution')

from .label_txt_formatter import label_txt_formatter, label_txt_formatter_gene


def plot_ax_stacked_bar(
    ax, y, label = None, title = None, xlabel = None, ylabel = None,
    color = None, hatch = None,
    bar_width = 1, show_legend=False, 
    show_xticklabels = False, show_yticklabels = False,
    legend_order_rev = False, edgecolor=None,
):
    '''
    Given an Axes class object, plot a stacked bar plot on that object
    https://matplotlib.org/api/axes_api.html
    - ax: Axes object
    - y: data
    - label: legend
    - title: sub plot title
    - xlabel: x-axis label
    - ylabel: y-axis label
    - bar_width: width of the stacked bar plot
    - show_legend: Boolean, whether we will show a legend or not
    - show_xticklabels: Boolean, whether we will show xticklabels  
    - show_yticklabels: Boolean, whether we will show yticklabels
    '''
    assert label is None or len(y) == len(label)
    data_len = len(y)
    bottom = np.append(np.zeros(1), np.cumsum(y)[:-1])
    p = [None] * data_len
    if color is None:
        colors = color_pallet_with_gray(data_len)
    else:
        colors = color
    if hatch is None:
        hatches = [''] * len(colors)
    else:
        hatches = hatch
    for i in range(data_len):        
        p[i] = ax.bar(
            0, y[i], 
            bottom=bottom[i], 
            width=(2 * bar_width), 
            color=colors[i],
            hatch=hatches[i],
            linewidth=1,
            edgecolor=edgecolor,
        )
    ax.set_xlim((0, 1))
    ax.set_ylim((0, 1))
    ax.set_xticks([])
    if(show_legend and label is not None):
        if(legend_order_rev):
            ax.legend(
                (p[data_len - 1 - i] for i in range(data_len)),
                tuple(reversed(label)),
                bbox_to_anchor=(1, 1)
            )
        else:
            ax.legend(
                (p[i] for i in range(data_len)),
                tuple(label),
                bbox_to_anchor=(1, 1)
            )
    if(not show_xticklabels):            
        ax.set_xticklabels([]) 
    if(not show_yticklabels):            
        ax.set_yticklabels([]) 
    if(title is not None):
        ax.set_title(title)
    if(xlabel is not None):
        ax.set_xlabel(xlabel)        
    if(ylabel is not None):
        ax.set_ylabel(ylabel)        
    return ax


def sort_plot_d_by_phe_gs_find_new_sort_order(plot_d, phe_gs, include_groups=True, additional_topk_filter=None):
    def flatten_list(l):
        return [item for sublist in l for item in sublist]    
    
    plot_d_gs_slice_vec = {
        k2:v2 for k2, v2 in {
            k1: np.array([x in v1 for x in plot_d['xticklabels']]) 
            for k1, v1 in phe_gs.dict.items()
        }.items() if np.any(v2)
    }

    plot_d_gs_keys = np.array(list(plot_d_gs_slice_vec.keys()))
    plot_d_gs_weights = np.array([
        np.sum(np.array(plot_d['y'])[plot_d_gs_slice_vec[k]]) for k in plot_d_gs_keys
    ])

    plot_d_gs_sort_order = np.argsort(-plot_d_gs_weights)

    plot_d_gs_weights_sorted = plot_d_gs_weights[plot_d_gs_sort_order]
    plot_d_gs_keys_sorted = plot_d_gs_keys[plot_d_gs_sort_order]

    logger.info('; '.join(['{}: {:.4f}'.format(x, y) for x, y in zip(plot_d_gs_keys_sorted, plot_d_gs_weights_sorted)]))

    if(len(plot_d_gs_slice_vec.values()) > 0):
        slice_vec_non_group = list(np.where(np.logical_not(
            reduce(lambda x, y: x | y, plot_d_gs_slice_vec.values())
        ))[0])
    else:
        slice_vec_non_group = np.arange(len(plot_d['xticklabels']))
    
    if(include_groups):
        plot_d_new_sort_order = np.array(flatten_list(
            [list(np.where(plot_d_gs_slice_vec[x])[0]) for x in plot_d_gs_keys_sorted] + 
            [slice_vec_non_group]
        ))
    elif(additional_topk_filter is None):
        plot_d_new_sort_order = np.array(slice_vec_non_group)
    else:
        plot_d_new_sort_order = np.array(slice_vec_non_group[:additional_topk_filter] + [slice_vec_non_group[-1]])        

    return plot_d_new_sort_order


def sort_plot_d_by_phe_gs(plot_d, phe_gs, include_groups=True, additional_topk_filter=None):
    sort_order = sort_plot_d_by_phe_gs_find_new_sort_order(plot_d, phe_gs, include_groups, additional_topk_filter)
    new_plot_d = dict([])    
    copy_key = set(['title', 'x', 'xlabel', 'ylabel'])
    for k, v in plot_d.items():
        if(k in copy_key):
            new_plot_d[k] = v
        else:
            new_plot_d[k] = [v[x] for x in sort_order]
    return new_plot_d


def contribution_score_plot(
    d, target, figsize=(32, 16), 
    pc_list = None, bar_width = 1,
    num_pcs_per_row=50, save=None,
    topk=None, sort=True, append_others=True,
    contribution_thr=None, contribution_cumsum_thr = None,
    yaxis_label_size = None,
    phe_gs=None, 
    include_groups=True,
    additional_topk_filter=None,
):
    if ((contribution_thr is None) and
        (contribution_cumsum_thr is None) and
        (topk is None)):
        contribution_thr = 0.01    
    if pc_list is None:
        pc_list = range(0, d.d['n_PCs'])
    if(num_pcs_per_row > len(pc_list)):
        num_pcs_per_row = len(pc_list)        
        
    num_pcs = len(pc_list)
    num_rows = int((num_pcs - 1)/num_pcs_per_row) + 1    
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(num_rows, num_pcs_per_row)
    fig_axs = [fig.add_subplot(sp) for sp in gs]
        
    for i in pc_list:
        plotd_raw = d.plot_data_contribution_generic(
            target, i,
            topk=topk, 
            contribution_thr=contribution_thr, 
            contribution_cumsum_thr=contribution_cumsum_thr,
            sort=sort, 
            append_others=append_others,
        )
        if(phe_gs is None):
            plotd = plotd_raw
        else:
            plotd = sort_plot_d_by_phe_gs(plotd_raw, phe_gs, include_groups, additional_topk_filter)            
        
        if(target == 'phe'):
            edgecolor='black'
        else:
            edgecolor=None
            
        fig_axs[i] = plot_ax_stacked_bar(
            fig_axs[i], plotd['y'], 
            plotd['xticklabels'], 
            ylabel = plotd['ylabel'] if i % num_pcs_per_row == 0 else None, 
            xlabel = 'PC{}'.format(i + 1) if i % num_pcs_per_row == 0 else '{}'.format(i + 1),
            show_xticklabels = False, 
            show_yticklabels = i % num_pcs_per_row == 0,
            color = plotd['stacked_bar_color'],
            hatch = plotd['stacked_bar_hatch'],
            bar_width = bar_width,
            edgecolor=edgecolor,
        )
        if(i % num_pcs_per_row == 0 and yaxis_label_size is not None):
            fig_axs[i].yaxis.label.set_size(yaxis_label_size)

    gs.tight_layout(fig, rect=[0, 0, 1, 1], w_pad=-.1)
    
    for ext in ['pdf', 'png']:
        if save is not None:
            fig.savefig(
                '{}.{}'.format(save, ext), 
                bbox_inches="tight", pad_inches=0.0
            )



def plot_contribution_and_save(
    d, phe_or_gene, labels, topk, 
    out_dir, fig_title, 
    fig_spacing, 
    contribution_thr=0.005, 
    pc_y_max = 1, 
    figsize=(28,21),
    phe_gs=None, 
    include_groups=True,
    additional_topk_filter=None,
):
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(len(labels), fig_spacing * topk)
    fig_axs = [fig.add_subplot(sp) for sp in gs]

    for i_label, label in enumerate(labels):
        pcs_0based, cos_scores = d.get_topk_pcs_for_phe_with_scores_by_label(label=label, topk=topk)
        for i_pc, pc in enumerate(pcs_0based):
            ax_idx = i_label * fig_spacing * topk + i_pc * fig_spacing
            ax = fig_axs[ax_idx]
            
            plotd_raw = d.plot_data_contribution_generic(
                phe_or_gene, 
                pc_index = pc, 
                contribution_thr=contribution_thr, 
                append_others=True,
            )
            if(phe_gs is None):
                plot_d = plotd_raw
            else:
                plot_d = sort_plot_d_by_phe_gs(plotd_raw, phe_gs, include_groups, additional_topk_filter)            
                        
            pc_y = plot_d['y']
            pc_y_offset = np.cumsum(plot_d['y']) - plot_d['y']
            pc_label = plot_d['xticklabels']
            pc_color = plot_d['stacked_bar_color']
            pc_hatch = plot_d['stacked_bar_hatch']
            pc_y_max = pc_y_max

            for i in range(len(pc_label)):
                if(phe_or_gene == 'phe'):
                    edgecolor='black'
                else:
                    edgecolor=None
                
                ax.bar(
                    0, pc_y[i], bottom=pc_y_offset[i], 
                    color=pc_color[i], 
                    hatch=pc_hatch[i],
                    width=1, 
                    label=pc_label[i],
                    linewidth=1,
                    edgecolor=edgecolor,                    
                )
            ax.set_ylim(0, pc_y_max)
            ax.get_xaxis().set_visible(False)        
            ax.set_title('PC{} ({:.2f})'.format(pc + 1, cos_scores[i_pc]), size=20)
            if(i_pc == 0):
                if(phe_or_gene == 'phe'):
                    ax.set_ylabel('Phenotype contribution score')
                else:
                    ax.set_ylabel('Gene contribution score')
            for axis_off_idx in range(fig_spacing - 1):
                fig_axs[ax_idx + axis_off_idx + 1].axis("off")

    for ext in ['pdf', 'png']:
        tmp_save_name = os.path.join(out_dir, '{}.{}'.format(fig_title, ext))
        logging.info(tmp_save_name)
        if(not os.path.exists(os.path.dirname(tmp_save_name))):
            os.makedirs(os.path.dirname(tmp_save_name))
        fig.savefig(
            tmp_save_name,
            bbox_inches="tight", pad_inches=0.0,
        )

        
def plot_contribution_legend_phe(
    d, pc_index, 
    bar_width=0.5, label_txt_formatter_len=24, 
    phe_gs=None, include_groups=True, topk=10,
    phe_list=None,
    save=None, 
    save_exts=['png', 'pdf'],
):
    fig = plt.figure(figsize=(6,6))
    gs = gridspec.GridSpec(1, 1)
    fig_axs = [fig.add_subplot(sp) for sp in gs]
    ## here we get the data    
    if(phe_list is not None):
        plotd = d.plot_data_contribution_phe_list(phe_list, pc_index=pc_index)    
    elif(phe_gs is None):
        plotd = d.plot_data_contribution_generic('phe', pc_index, append_others=True, topk=topk)
    elif(include_groups):
        plotd_raw = d.plot_data_contribution_generic('phe', pc_index, append_others=True, topk=topk)
        plotd = sort_plot_d_by_phe_gs(plotd_raw, phe_gs, include_groups)
    else:
        plotd_raw = d.plot_data_contribution_generic('phe', pc_index, append_others=True, topk=topk + phe_gs.get_size())        
        plotd = sort_plot_d_by_phe_gs(plotd_raw, phe_gs, include_groups, additional_topk_filter=topk)

    tmp_labels = [
        label_txt_formatter(x, label_txt_formatter_len)
        for x in plotd['xticklabels']
    ]
    fig_axs[0] = plot_ax_stacked_bar(
        fig_axs[0], 
        plotd['y'], 
        tmp_labels, 
        title = plotd['title'],
        ylabel = plotd['ylabel'], 
        show_yticklabels = True,
        color = plotd['stacked_bar_color'],
        hatch = plotd['stacked_bar_hatch'],        
        show_legend=True, bar_width=bar_width,
        legend_order_rev = True,
        edgecolor='black',
    )
    
    for t, o in zip(tmp_labels, plotd['xticklabels']):
        if(t != o):
            logger.debug('{}\t{}'.format(t, o))    
    if save is not None:
        for ext in save_exts:
            tmp_save_name = '{}.{}'.format(save, ext)
            logger.info('saving to {}'.format(tmp_save_name))            
            if(not os.path.exists(os.path.dirname(tmp_save_name))):
                os.makedirs(os.path.dirname(tmp_save_name))
            fig.savefig(tmp_save_name, bbox_inches="tight", pad_inches=0.0)

            
def plot_contribution_legend_phe_batch_fig(
    phe_labels_dict, out_dir, fig_title, d, topk, phe_gs=None, topk_in_pc=10,
):
    for phe, phe_short in phe_labels_dict.items():
        pcs_0based, cos_scores = d.get_topk_pcs_for_phe_with_scores_by_label(label=phe, topk=topk)
        for i_pc, pc in enumerate(pcs_0based):
            savefilename = os.path.join(out_dir, '{}_legend_{}_{}_PC{}'.format(
                fig_title, phe_short, i_pc + 1, pc + 1,
            ))
            logger.debug(savefilename)
            if(phe_gs is not None):
                plot_contribution_legend_phe(d, pc, topk=topk_in_pc, phe_gs=phe_gs, include_groups=False, save=savefilename)
            else:
                plot_contribution_legend_phe(d, pc, topk=topk_in_pc, save=savefilename)
            

def plot_contribution_legend_gene(
    d, pc_index, gene_dict=None, bar_width=0.5, topk=10, 
    save=None, 
    save_exts=['png', 'pdf'],    
):
    fig = plt.figure(figsize=(6,6))
    gs = gridspec.GridSpec(1, 1)
    fig_axs = [fig.add_subplot(sp) for sp in gs]
    ## here we get the data
    plotd = d.plot_data_contribution_generic('gene', pc_index, append_others=True, topk=topk)
    tmp_labels = [
        label_txt_formatter_gene(x, gene_dict)
        for x in plotd['xticklabels']
    ]
    
    fig_axs[0] = plot_ax_stacked_bar(
        fig_axs[0], 
        plotd['y'], 
        tmp_labels, 
        title = plotd['title'],
        ylabel = plotd['ylabel'], 
        show_yticklabels = True,
        color = plotd['stacked_bar_color'],
        hatch = plotd['stacked_bar_hatch'],        
        show_legend=True, bar_width=bar_width,
        legend_order_rev = True,
    )
    for t, o in zip(tmp_labels, plotd['xticklabels']):
        if(t != o):
            logger.debug('{}\t{}'.format(t, o))    
    if save is not None:
        for ext in save_exts:
            tmp_save_name = '{}.{}'.format(save, ext)
            logger.info('saving to {}'.format(tmp_save_name))
            if(not os.path.exists(os.path.dirname(tmp_save_name))):
                os.makedirs(os.path.dirname(tmp_save_name))            
            fig.savefig(tmp_save_name, bbox_inches="tight", pad_inches=0.0)

            
def plot_contribution_legend_gene_batch_fig(
    phe_labels_dict, gene_dict, out_dir, fig_title, d, topk, 
):
    for phe, phe_short in phe_labels_dict.items():
        pcs_0based, cos_scores = d.get_topk_pcs_for_phe_with_scores_by_label(label=phe, topk=topk)
        for i_pc, pc in enumerate(pcs_0based):
            savefilename = os.path.join(out_dir, '{}_legend_{}_{}_PC{}'.format(
                fig_title, phe_short, i_pc + 1, pc + 1,
            ))
            logger.debug(savefilename)
            plot_contribution_legend_gene(d, pc, gene_dict, save=savefilename)

def plot_contribution_legend_batch_topk_fig(
    gene_dict, out_dir, fig_title, d, topk, phe_gs=None, 
):
    for pc in range(topk):
        for phe_var_gene in ['phe', 'gene']:
            savefilename = os.path.join(
                out_dir, '{}_legend_{}_PC{}'.format(
                fig_title, phe_var_gene, pc+1
            ))
            logger.debug(savefilename)
            if(phe_var_gene == 'phe'):
                plot_contribution_legend_phe(d, pc, phe_gs=phe_gs, include_groups=False, save=savefilename)
            elif(phe_var_gene == 'gene'):
                plot_contribution_legend_gene(d, pc, gene_dict, save = savefilename)
      
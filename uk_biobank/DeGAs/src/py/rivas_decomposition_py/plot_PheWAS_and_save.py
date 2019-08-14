import os, logging
import numpy as np
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

logger = getLogger('plot_PheWAS_and_save')


from .label_txt_formatter import label_txt_formatter

def plot_PheWAS_and_save(df, out_dir, save_name, color, effect_size_label, figsize=(18,12)):
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 2)
    fig_axs = [fig.add_subplot(sp) for sp in gs]

    vertical_pos = np.arange(len(df)) + 0.5
    pheno_labels = df['Name'].map(lambda x: label_txt_formatter(x, 30))

    fig_axs[0].barh(
        vertical_pos, df['l10pval'], 
        color=color
    )
    fig_axs[1].errorbar(
        df['LOR'], vertical_pos,
        xerr=df['196SE'], fmt='o', 
        color=color
    )

    for ax in fig_axs:
        ax.invert_yaxis()
    fig_axs[0].set_yticklabels(pheno_labels, fontsize=16)
    fig_axs[0].set_yticks(vertical_pos)
    fig_axs[0].set_xlabel('-log10(p-value)')
    fig_axs[1].set_xlabel(effect_size_label)
    fig_axs[1].get_yaxis().set_visible(False)
    fig_axs[1].plot(
        np.zeros(len(df)+1), 
        np.arange(len(df)+1), 
        color='black')

    effect_size_plot_max = max(
        np.max((df['LOR'] + df['196SE']).map(lambda x: abs(x))), 
        np.max((df['LOR'] - df['196SE']).map(lambda x: abs(x)))
    )
    fig_axs[1].set_xlim((-effect_size_plot_max*1.1, effect_size_plot_max*1.1))

    for ext in ['pdf', 'png']:
        tmp_save_name = os.path.join(out_dir, '{}.{}'.format(save_name, ext))        
        logger.info('saving to {}'.format(tmp_save_name))
        if(not os.path.exists(os.path.dirname(tmp_save_name))):
            os.makedirs(os.path.dirname(tmp_save_name))        
        fig.savefig(
            tmp_save_name,
            bbox_inches="tight", pad_inches=0.0,
        )

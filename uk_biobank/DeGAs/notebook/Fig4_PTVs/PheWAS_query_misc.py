import os, sys, collections
import numpy as np
import pandas
gbe_repo = os.path.join('/', 'opt', 'biobankengine', 'GlobalBioBankEngineRepo')
sys.path.append(os.path.join(gbe_repo, 'gbe_browser'))
import config, lookups


def variant_page_data_prep_sub(icdstats, sort_key='log10pvalue'):
    plot_d_raw = collections.defaultdict(list)
    keys = icdstats[0].keys()
    for key in keys:
        plot_d_raw[key] = np.array([x[key] for x in icdstats])
    plot_df = pandas.DataFrame(plot_d_raw)
    #.sort_values(
    #    by=['Group', sort_key], ascending=[True, False]
    #)
    plot_d_dict = collections.defaultdict(collections.defaultdict)
    
    groups = sorted(set(plot_df['Group']))
    for group in groups:
        for key in keys:
            plot_d_dict[group][key] = list(plot_df[plot_df['Group'] == group][key])
    for group in groups:
        for key in ['OR', 'LOR', 'L95OR', 'U95OR', 'pvalue', 'SE', 'log10pvalue']:
            plot_d_dict[group][key] = [float(x) for x in plot_d_dict[group][key]]
    for group in groups:
        #error_bar = {'L95OR': -1, 'U95OR': 1}
        #for key in error_bar.keys():
        #    diff = np.array(plot_d_dict[group][key]) - np.array(plot_d_dict[group]['LOR'])
        #    plot_d_dict[group]['d{}'.format(key)] = [0 if np.isnan(x) else np.abs(x) for x in diff]
        plot_d_dict[group]['196SE'] = list( 1.96 * np.array(plot_d_dict[group]['SE']) )

    for group in groups:
        if group in set(['INI', 'INI_FC', 'BROADQT']):
            beta_or_lor = 'BETA'
            beta_or_lor_val = plot_d_dict[group]['LOR']
            beta_or_lor_l95 = np.array(plot_d_dict[group]['LOR']) - 1.96 * np.array(plot_d_dict[group]['SE'])
            beta_or_lor_u95 = np.array(plot_d_dict[group]['LOR']) + 1.96 * np.array(plot_d_dict[group]['SE'])

        else:
            beta_or_lor = 'OR'
            beta_or_lor_val = np.exp(np.array(plot_d_dict[group]['LOR']))
            beta_or_lor_l95 = np.exp(np.array(plot_d_dict[group]['LOR']) - 1.96 * np.array(plot_d_dict[group]['SE']))
            beta_or_lor_u95 = np.exp(np.array(plot_d_dict[group]['LOR']) + 1.96 * np.array(plot_d_dict[group]['SE']))

        group_len = len(plot_d_dict[group]['icd'])
        plot_d_dict[group]['text'] = [
            '{}. Case: {}, P-value: {:.3e}, {} = {:.5f} (95% [{:.5f}, {:.5f}]), SE = {:.5f}'.format(
                ''.join([c if c != '_' else ' ' for c in x[0]]), x[1], x[2], x[3], x[4], x[5], x[6], x[7]
            ) for x in zip(
                plot_d_dict[group]['Name'],
                plot_d_dict[group]['Case'],
                plot_d_dict[group]['pvalue'],
                [beta_or_lor] * group_len,
                beta_or_lor_val,
                beta_or_lor_l95,
                beta_or_lor_u95,
                plot_d_dict[group]['SE'],
                #plot_d_dict[group]['L95OR'],
                #plot_d_dict[group]['U95OR'],
            )
        ]
    return plot_d_dict

def get_data(variant_str, db):

    chrom, pos = variant_str.split('-')
    icdstats = lookups.get_icd_by_chrom_pos(db, chrom, pos)

    indexes = []
    seend = {}

    for idx in range(0,len(icdstats)):
        # ICD10=T81/Name=Complications of procedures, not elsewhere classified/Chapter=T/L95OR=0.97/U95OR=2.04/OR=1.40/pvalue=0.0756463/l10pval=1.12/Case=1229
        item = icdstats[idx]
        icd10 = item['icd']
        item['Code'] = icd10
        # icd10info = lookups.get_icd_info(db, icd10)
        if 'Name' not in item:
            item['Name'] = 'NA'
            item['Group'] = 'NA'
            item['OR'] = 1
            item['LOR'] = 0
            item['L95OR'] = 1
            item['U95OR'] = 1
            item['pvalue'] = 1
            item['l10pval'] = 0
            item['Case'] = 'NA'
            item['SE'] = 0
            indexes.append(idx)
        else:
            # item['Name'] = icd10info[0]['Name']
            item['Group'] = icd10[0] # default value
            groups = ['RH', 'FH', 'HC', 'cancer', 'ADD', 'INI', 'MED', 'BIN', 'BRMRI', 'BROADBIN', 'BROADQT', 'INI_FC', 'BIN_FC']
            for group in groups:
                if icd10.startswith(group):
                    item['Group'] = group
                    break
            item['OR'] = format(float(item['or_val']), '.4g')
            item['LOR'] = format(float(item['lor']), '.4g')
            item['L95OR'] = format(float(item['l95or']), '.4g')
            item['U95OR'] = format(float(item['u95or']), '.4g')
            item['pvalue'] = format(float(item['pvalue']), '.4g')
            item['l10pval'] = format(float(item['log10pvalue']), '.4g')
            item['SE'] = format(float(item['se']), '.4g')
            if float(item['pvalue']) == 0:
                item['pvalue'] = numpy.finfo(float).eps
                item['pvalue'] = format(float(item['pvalue']),'.4g')
                item['l10pval'] = 250
            # item['Case'] = icd10info[0]['Case']
            se =  format(float(item['se']), '.4g')
            if float(item['l10pval']) < 1 or float(se) >= .5 or (float(se) >= .08 and item['OR'] == item['LOR']) or int(item['Case']) <= 100  or item['Code'] == "HC67" or icd10 in seend:
                indexes.append(idx)
            seend[icd10] = icd10


    d = variant_page_data_prep_sub(icdstats)
    return d
    

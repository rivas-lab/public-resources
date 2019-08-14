import numpy as np
import pandas as pd
import collections    

def read_PheWAS_npz(npz_file):    
    d = np.load(npz_file)
    return collections.OrderedDict([
        (k, d[k].item()) for k in d.keys()
    ])


def filter_PheWAS_by_p_value(data_dict, p_val_thr):     
    new_d = collections.OrderedDict()
    for category, cat_data in data_dict.items():
        data_filter = np.array(
            [x <= p_val_thr for x in cat_data['pvalue']]
        )
        new_d[category] = collections.OrderedDict([
            (k, np.array(v)[data_filter]) 
            for k, v in zip(cat_data.keys(), cat_data.values())
        ])
    return new_d        


def df_PheWAS(data_dict):    
    if(len(data_dict) == 0):
        return pd.DataFrame()
    else:
        float_cols = [
            'l10pval', 'log10pvalue',
            '196SE',
            'lor', 'LOR', 'or_val',
            'l95or', 'u95or',             
            'L95OR', 'U95OR',            
        ]
        df_concat = pd.DataFrame()
        for data_per_category in data_dict.values():
            df_concat = pd.concat(
                [
                    df_concat, 
                    pd.DataFrame(data_per_category)
                ], 
                ignore_index=True
            )
        for col in float_cols:
            df_concat[col] = np.array([
                float(x) for x in df_concat[col]
            ])
        df_concat['is_binary'] = df_concat['Group'].map(
            lambda x: x not in set(['INI', 'INI_FC', 'BROADQT'])
        )            
        return df_concat.sort_values(
            by=['is_binary', 'l10pval'], 
            ascending=[False, False],
        )

    
def df_PheWAS_filter_by_icd(df, icd_set):
    return df[df['icd'].map(lambda x: x in icd_set)]


def df_PheWAS_filter_by_case_count(df, min_case_count):
    return df[df['Case'].map(lambda x: x >= min_case_count)]


def PheWAS_data_loader(npz_file, icd_set, p_value_thr, min_case_count):
    
    npz_data = read_PheWAS_npz(npz_file)
    return collections.OrderedDict([
        (k, df_PheWAS_filter_by_case_count(df_PheWAS_filter_by_icd(df_PheWAS(filter_PheWAS_by_p_value(v, p_value_thr)), icd_set), min_case_count))
        for k, v in npz_data.items()
    ])


import re

import mygene

def label_txt_formatter(label, max_len = None):
    '''
    Given a label text, return an abbreviated text 
    (for figure)
    '''

    replace_strs = [
        ('_', ' '),
        (r'\(right\)$', '(R)'),
        (r'\(left\)$',  '(L)'),
        (r'percentage', '%'),
        (r'90th percentile', '90PCTL'),
        (r'number', '#'),
        (r'Number', '#'),
        (r'Volume of', 'Vol.'),
        (r'predicted', 'pred.'),
        (r'blood pressure', 'BP'),
        (r'Average', 'Ave.'),
        (r'average', 'ave.'),
        (r', automated reading', ' (AR)'),   
        (r'cholelithiasis/gall stones', 'Gallstones'),
        (r'Body mass index \(BMI\)', 'BMI'),
        (r'Creatinine \(enzymatic\) in urine', 'Creatinine in urine'),
        (r'Weighted-mean', 'WA'),
        (r'treatments/medications', 'medications'),
        (r'Peak expiratory flow \(PEF\)', 'PEF'),
        (r'Forced expiratory volume in 1-second \(FEV1\)', 'FEV1'),
        (r'Forced vital capacity \(FVC\)', 'FVC'),
        (r'statistic', 'stat.'),
        (r"Alzheimer's disease/dementia", "Alzheimer's/dementia"),
        (r'Time spend outdoors in', 'Outdoor time,'),
        (r'Time spent outdoors in', 'Outdoor time,'),
        (r'Nitrogen dioxide', 'NO2'),
        (r'Particulate matter ', ''),
        (r'sound level of noise pollution', 'noise lvl.'),
        (r'platelet \(thrombocyte\)', 'platelet'),
        (r'White blood cell \(leukocyte\)', 'White blood cell'),
        (r'Red blood cell \(erythrocyte\)', 'Red blood cell'),
        (r'Age at menopause \(last menstrual period\)', 'Age at menopause'),
        (r'difficulty/problems', 'problems'),
        (r'Nucleated red blood cell', 'Nuc. red blood cell'),
        (r'night-time', 'nighttime'),
        (r'air pollution', 'air poll.'),
        (r';', ''),          
        (r'heart attack/myocardial infarction', 'heart attack (MI)'),
        (r'Childhood sunburn occasions', 'Childhood sunburn'),
        (r'deep venous thrombosis \(dvt\)', 'DVT'),
        (r'dvt', 'DVT'),
        (r'pulmonary embolism', 'PE'),
        (r'hereditary/genetic', 'genetic'),
        
    ]
    formatted = label
    for replace in replace_strs:
        formatted = re.sub(replace[0], replace[1], formatted)
    if(max_len is not None and len(formatted) > max_len):
        formatted = '{}..'.format(formatted[:max_len])
        
    return formatted[0].capitalize() + formatted[1:]

def mygene_conv(mg, ens):
    ginfo = mg.querymany(
        ens, 
        scopes='ensembl.gene',
        fields='symbol', 
        species='human', 
        as_dataframe=True,
        df_index=False,
    ).dropna()
    if('symbol' in ginfo.keys()):
        conv_dict = dict(zip(ginfo['query'], ginfo['symbol']))    
    else:
        conv_dict = dict([])
    return [conv_dict[x] for x in ens if x in conv_dict]


def label_txt_formatter_gene(label, gene_dict=None, italic=True):
    def to_italic(str):
        return r"$\it{" + str + "}$"
    if(label.startswith('ENSG')):
        if(gene_dict is not None):
            label_remap = [
                gene_dict[x] for x in label.split(',') 
                if x in gene_dict
            ]
        else:
            label_remap = []
        if(len(label_remap) > 0):
            if(italic):
                return to_italic(','.join(label_remap))
            else:
                return ','.join(label_remap)
        else:
            mg = mygene.MyGeneInfo()
            label_mygene_map = mygene_conv(mg, label.split(','))
            if(len(label_mygene_map) > 0):
                if(italic):
                    return to_italic(','.join(label_mygene_map))
                else:
                    return ','.join(label_mygene_map)
            else:
                if(italic):
                    return to_italic(label)
                else:            
                    return label
    elif(label == 'others'):
        return 'Others'
    else:
        if(italic):
            return to_italic(label)
        else:            
            return label
        
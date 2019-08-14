import os, logging, collections
import subprocess
import numpy as np
import pandas as pd

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

logger = getLogger('generate_great_file')


def generate_great_file(d, out_dir, topk=5000, generate_targz=True):
    logger = logging.getLogger('generate_great_file')
    logger.info('writing bed files for GREAT')
    great_bed_dir = os.path.join(out_dir, 'great', 'bed')
    if(not os.path.exists(os.path.join(out_dir, 'great'))):
        os.makedirs(os.path.join(out_dir, 'great'))    
    if(not os.path.exists(great_bed_dir)):
        os.makedirs(great_bed_dir)
    
    if(not os.path.exists(os.path.join(out_dir, 'great', 'bed.tar.gz'))):    
        # write bed files
        for i in range(d.d['n_PCs']):
            d.bed_data_contribution_var(pc_index=i, topk=topk).to_csv(
                os.path.join(great_bed_dir, 'PC{}.bed'.format(i)),
                sep='\t', index=False, header=False)

        # generate tar.gz file
        if(generate_targz):
            subprocess.run(
                ('tar', '-czf', 'bed.tar.gz', 'bed'), 
                cwd=os.path.join(out_dir, 'great'), check=True)
            print(os.path.join(out_dir, 'great', 'bed.tar.gz'))
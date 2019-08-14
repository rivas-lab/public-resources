import collections
import re
from functools import reduce

class Phe_label_groups():    
    '''
    Class to group related phenotype labels using regEx
    '''    

    def __init__(self, labels):        
        '''
        labels: iterable
        '''
        self.labels = labels        
        self.dict = collections.defaultdict(set)
        
    def add(self, name, regex_key = None):                
        if (regex_key is None):
            regex_key = name
        self.dict[name] = set(sorted(filter(
            lambda x: re.search(regex_key, x.lower()) is not None, 
            self.labels
        )))
        
    def get_dict(self):
        return self.dict
    
    def get_size(self):
        return reduce(lambda x, y: x + y, [len(l) for l in self.dict.values()])
    
    def len(self):
        return len(self.dict)
    
    def __str__(self):
        return '\n\n'.join(['{}\n{}'.format(
            k, 
            '\n'.join(['  {}'.format(x) for x in sorted(v)])
        ) for k, v in self.dict.items()])
    
from scidbpy import connect
import numpy as np
from PheWAS_query_misc import get_data

def main():
    db = connect('http://localhost:8080')


    np.savez(
        'PheWAS.npz', 
        GPR151 = get_data('5-145895394',  db),
        PDE3B =  get_data('11-14865399', db)
    )

if __name__ == '__main__':
    main()


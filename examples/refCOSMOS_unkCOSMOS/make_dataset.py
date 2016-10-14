"""
Go from input catalog to hdf5 file because I didn't put in ability to handle csv files
"""

import pandas as pd

from wizard.config import check_make

###############################################################################
# run
###############################################################################

def main(out_directory, cosmos_path, cosmos_randoms_path, max_objects=50000000):

    # load up the files
    cosmos = pd.read_csv(cosmos_path, sep='\t', header=None, names=['RA', 'DEC', 'Z_MEAN_DES', 'Z_MC_DES', 'Z_COSMOS', 'PDFL68', 'PDFH68'], skiprows=[0], index_col=False)

    # save
    cosmos.to_hdf(out_directory + 'cosmos_galaxies.h5', 'catalog')

    # repeat with randoms
    cosmos_randoms = pd.read_csv(cosmos_randoms_path, sep='\t', header=None, names=['RA', 'DEC', 'Z_COSMOS'], skiprows=[0], index_col=False)
    cosmos_randoms['Z_MEAN_DES'] = cosmos_randoms['Z_COSMOS']
    cosmos_randoms['Z_MC_DES'] = cosmos_randoms['Z_COSMOS']
    # save
    cosmos_randoms.to_hdf(out_directory + 'cosmos_randoms.h5', 'catalog')

if __name__ == '__main__':
    out_directory = './dataset/'
    check_make(out_directory)
    cosmos_path = 'des_matched_05arcsec.cat'
    cosmos_randoms_path = 'randoms_matched_0_00002.cat'
    cosmos_randoms_path = 'final_cat_randoms20.txt'
    main(out_directory, cosmos_path, cosmos_randoms_path)

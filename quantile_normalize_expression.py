#####################################################################################################################
### Sara Geraghty ###
### July 2022 ###
#####################################################################################################################

# QUANTILE-NORMALIZATION
# Use scikit learn's implementation of the quantile normalization function (quantile_transform)
# https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.quantile_transform.html#sklearn.preprocessing.quantile_transform

import sys
import os
import csv
import re
import numpy as np
import pandas as pd
from sklearn.preprocessing import quantile_transform

CURRENT_PATH = os.path.dirname(os.path.realpath(__file__))

# Command Line Arguments:
# Argument 1 is local path to the input gene expression count matrix(es). This
# has already been filtered using edgeR to remove lowely expressed genes. We will
# also store the output files here.

def quantile_normalize(fn, outpath, name):

    exp_df = pd.read_csv(fn, header = 'infer', skiprows = None)
    col_names = exp_df.columns
    n_cols = len(col_names)
    gene_names = pd.read_csv(fn, index_col = False, usecols = [0], skiprows = 0)
    gene_names = list(gene_names.iloc[:,0])
    #print(gene_names[:5])
    #print(len(gene_names))

    # Import gene expression data table in read mode
    try:
        gene_expr_df = np.loadtxt(open(fn, "rb"), 
            delimiter = ",", skiprows = 1, usecols = np.arange(1, n_cols))
        #print(gene_expr_df[:5])

        # Call sklearn quantile_transform function with output distribution
        # set to normal
        quantile_transform(gene_expr_df, n_quantiles = (n_cols-1), 
            copy = False, output_distribution = "normal")
        #print(gene_expr_df.shape)

        gene_expr_df_final = pd.DataFrame(gene_expr_df, columns = col_names[1:], 
            index = gene_names)

        # Save to a file
        outfn = CURRENT_PATH + outpath + name + "expression_quantile_normalized_sklearn.csv"
        gene_expr_df_final.to_csv(outfn)

        #gene_expr_df.close()

    except Exception as e:
        print("an exception occurred:" + e)


def main():

    local_path = sys.argv[1]

    # Get the names of the files matching the given pattern in the given path
    filenames = next(os.walk(CURRENT_PATH + local_path), (None, None, []))[2]
    print(filenames[0:15])
    r = re.compile("_expression_counts_DF_edgeRfilt_TO.csv")
    filenames_sub = list(filter(r.match, filenames))

    # Call the helper function to get a quantile-normalized version of each of these files
    for f in filenames:
        if(re.search('_expression_counts_DF_edgeRfilt_TO.csv', f) != None):
            full_f = CURRENT_PATH + local_path + f
            name = f.split("_")[0] + "_"
            quantile_normalize(full_f, local_path, name)

    #f = "expression_counts_DF_edgeRfilt_TO.csv"
    #filename = CURRENT_PATH + local_path + f
    #quantile_normalize(filename, local_path, "")

if __name__ == "__main__":
    # execute only if run as a script
    main()
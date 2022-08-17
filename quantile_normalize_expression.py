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
import numpy as np
import pandas as pd
from sklearn.preprocessing import quantile_transform

CURRENT_PATH = os.path.dirname(os.path.realpath(__file__))

def main():
    fn = CURRENT_PATH + sys.argv[1]

    # Determine number of columns (samples) from the first line of text
    with open(fn) as f:
            n_cols = len(f.readline().split(","))
            col_names = f.readline().split(",")
            #row_names = f.

    # Import gene expression data table in read mode
    gene_expr_df = np.loadtxt(open(fn, "rb"), 
        delimiter = ",", skiprows = 1, usecols = np.arange(1, n_cols))
    print(gene_expr_df)

    # Get the number of samples (columns), as this will be the number of quantiles
    #num_samples = np.shape(gene_expr_df)[1]
    #print(num_samples)

    # Call sklearn quantile_transform function with output distribution
    # set to normal
    quantile_transform(gene_expr_df, n_quantiles = n_cols, 
        copy = False, output_distribution = "normal")

    print(gene_expr_df)

    # Save to a file
    outfn = CURRENT_PATH + sys.argv[2] + "expression_quantile_normalized_sklearn.csv"
    #outdf = pd.DataFrame(gene_expr_df, columns = col_names)
    np.savetxt(outfn, gene_expr_df, delimiter = ",")

    #gene_expr_df.close()

if __name__ == "__main__":
    # execute only if run as a script
    main()
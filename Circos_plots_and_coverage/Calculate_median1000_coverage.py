'''
Title: Merging bed files containing median1000 coverage
Author: Lewis Fisher
Date: 01/08/2022
'''

import pandas as pd
import os
import sys
import numpy as np

# Working directory
work_dir = sys.argv[1]
# bed file containing coverage
file_in = sys.argv[2]
os.chdir(work_dir)
os.listdir()

# import bed file
cov_df = pd.read_csv(file_in, header=None, index_col = None, sep="\t")

new_chr = []
new_start = []
new_end = []
new_median = []

# Using a for loop to take the median coverage at every 1000bp
# Creating lists to host the new information
# The remaining coverage that's not a multiple of 1000  is taken between 1 and <1000
for i in range(0,len(cov_df),1000):
    if i +999 > len(cov_df):
        median_1000 = np.median(cov_df.loc[i:len(cov_df)-1, 2])
        new_median.append(median_1000)
        new_chr.append(cov_df.loc[len(cov_df)-1, 0])
        new_start.append(cov_df.loc[i, 1])
        new_end.append(cov_df.loc[len(cov_df)-1, 1])
    else:
        median_1000 = np.median(cov_df.loc[i:i+999, 2])
        new_median.append(median_1000)
        new_chr.append(cov_df.loc[i+999, 0])
        new_start.append(cov_df.loc[i, 1])
        new_end.append(cov_df.loc[i+999, 1])

# a new dataframe is created with the median coverage
cov_bed = {
    'chr': new_chr,
    'start': new_start,
    'end': new_end,
    'median_cov':new_median
    }

# The same file name is used but with the addition of ""_median
file_out = file_in.split(".")
file_out = file_out[0] + "_median." + file_out[1]
new_df = pd.DataFrame(cov_bed)
new_df.to_csv(file_out, header = False, index=False, sep="\t")


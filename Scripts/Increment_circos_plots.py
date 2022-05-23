
import pandas as pd
import os
import sys
import numpy as np


work_dir = sys.argv[1]
file_in = sys.argv[2]
os.chdir(work_dir)
os.listdir()

cov_df = pd.read_csv(file_in, header=None, index_col = None, sep="\t")



new_chr = []
new_start = []
new_end = []
new_median = []


for i in range(0,len(cov_df),1000):
    # print(cov_df.loc[i:i+999,])
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




cov_bed = {
    'chr': new_chr,
    'start': new_start,
    'end': new_end,
    'median_cov':new_median
    }
file_out = file_in.split(".")
file_out = file_out[0] + "_median." + file_out[1]
new_df = pd.DataFrame(cov_bed)
new_df.to_csv(file_out, header = False, index=False, sep="\t")


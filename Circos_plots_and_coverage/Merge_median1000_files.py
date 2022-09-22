'''
Title: Merging bed files containing median1000 coverage
Author: Lewis Fisher
Date: 01/08/2022
'''


import pandas as pd
import sys
import os


# Working directory
work_dir = sys.argv[1]

# output file name
outfile = sys.argv[2]
# file extention/ what the files endwith e.g. "_median.bed"
file_ends = sys.argv[3]

os.chdir(work_dir)
dir_files = os.listdir(work_dir)

# create empty dataframe
new_df = pd.DataFrame()
count = 0
# Merging a series of bed files containing median coverage every 1000bp
for file in dir_files:
    if file.endswith(file_ends):
        file_split = file.split("_")
        sample = " ".join(file_split[0:2])
        coverage = pd.read_csv(file, header= None, index_col= None, sep="\t")
        sample_list = [sample] * len(coverage.index)
        coverage[0] = sample_list
        if count == 0:
            new_df = coverage
            count =1
        else:
            new_df = pd.concat([new_df, coverage])
new_df.to_csv(outfile, header=None, index=None, sep="\t")




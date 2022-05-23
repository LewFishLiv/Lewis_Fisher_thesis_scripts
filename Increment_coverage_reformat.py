import pandas as pd
import sys
import os
import re


# Working directory
work_dir = sys.argv[1]

# output file name
outfile = sys.argv[2]
# file extention/ what the files endwith
file_ends = sys.argv[3]

os.chdir(work_dir)
dir_files = os.listdir(work_dir)

# create empty dataframe
new_df = pd.DataFrame()
count = 0
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
print(new_df)



'''
Author: Lewis Fisher
Date: 19/05/2022
Apply functional annotations and data to breseq output dataframe
'''
import pandas as pd
import os
import sys
import numpy as np

# COG annotation file
enog_input = sys.argv[1]
# Breseq input file
query_input = sys.argv[2]
# COG table with corresponding descriptions
cog_trans_table = sys.argv[3]
# output directory
out_dir = sys.argv[4]

df_egg = pd.read_csv(enog_input, header = 0, skiprows= 4,skipfooter= 3,  sep= "\t", engine= 'python')
df_query = pd.read_csv(query_input, header = 0, sep= ',', index_col=0)
df_cog_trans = pd.read_csv(cog_trans_table, header = 0, sep = ',')

# print(df_egg.columns)

# Creating names of lists to make maps to substitute for cog annotations and functional annotations
gene_names = df_egg['Preferred_name'].tolist()
locus_tags = df_egg['#query'].tolist()
cog_names = df_egg['COG_category']
cog_descriptions = df_cog_trans['functional_annotation'].tolist()
# Translations of cog letters and respective function
cog_trans = dict(zip(df_cog_trans['cog'].tolist(),df_cog_trans['functional_annotation'].tolist()))

# gene_dict, creating a dictionary to map to the query database
gene_dict = dict(zip(locus_tags, gene_names))

# creating a dictionary of locus tags and cogs to match cogs to genes from the breseq query.
cog_dict = dict(zip(locus_tags, cog_names))
# print(df_query.columns)

new_list = []
egg_query = df_egg['#query'].tolist()


# Adding locus tags from intergenic loci to the dictionary with the appropriate data for mapping
# Two flanking genes are delimited by a "/", output is also delimited in the same way
# Genes not found in the eggnog dataset are called "unknown" for convenience, it is probably these missing values
# are non-protein coding genes.
for i in df_query['locus_tag'].tolist():
    intergenic_snp = i.split('/')
    if len(intergenic_snp) == 2:
        if not (intergenic_snp[0] in egg_query):
            intergenic_gene1 = "unknown"
            cog_1 = ""
        else:
            for e in range(0,len(egg_query)):
                if intergenic_snp[0] == egg_query[e]:
                    intergenic_gene1 = gene_names[e]
                    cog_1 = cog_names[e]
        if not (intergenic_snp[1] in egg_query):
            intergenic_gene2 = "unknown"
            cog_2 = ""
            gene_dict[i] = intergenic_gene1 + "/" + intergenic_gene2
        else:
            for f in range(0,len(egg_query)):
                if intergenic_snp[1] == egg_query[f]:
                    intergenic_gene2 = gene_names[f]
                    cog_2 = cog_names[f]
                    gene_dict[i] =intergenic_gene1+ "/"+intergenic_gene2
                    cog_dict[i] = cog_1 + "/" + cog_2
    else:
        if not (i in egg_query):
            gene_dict[i] = "unknown"

# Mapping the cogs and eggnog annotations by locus tag (the key value)
df_query['eggnog_gene'] = df_query['locus_tag'].tolist()
df_query['eggnog_gene'] = df_query['eggnog_gene'].map(gene_dict)
df_query['cog_id'] = df_query['locus_tag'].tolist()
df_query['cog_id'] = df_query['cog_id'].map(cog_dict)

# cog IDs were convoluted, some genes posess multiple functional cogs so these were added and delimited by ", "
# intergenic mutations were treated the same with "/" delimited outputs below
for i in df_query['cog_id'].tolist():

    if i is np.nan:
        continue
    elif len(i.split('/')) == 2:
        cog_split = i.split('/')
        if cog_split[0] == "":
            replaced_list1 = ["unknown"]
        elif cog_split[0] == "-":
            replaced_list1 = ["-"]
        elif len(cog_split[0]) == 1:
            replaced_list1 = [x if x not in cog_trans else x + ":" + cog_trans[x] for x in cog_split[0]]
        elif len(cog_split[0]) > 1:
            replaced_list1 = [x if x not in cog_trans else x + ":" + cog_trans[x] for x in cog_split[0]]
        if cog_split[1] == "":
            replaced_list2 = ["unknown"]
        elif cog_split[1] == "-":
            replaced_list2 = ["-"]
        elif len(cog_split[1]) == 1:
            replaced_list2 = [x if x not in cog_trans else x + ":" + cog_trans[x] for x in cog_split[1]]
        elif len(cog_split[1]) > 1:
            replaced_list2 = [x if x not in cog_trans else x + ":" + cog_trans[x] for x in cog_split[1]]
        cog_trans[i] = "/".join([", ".join(replaced_list1), ", ".join(replaced_list2)])

    elif len(i.split("/")) == 1:
        if len(i) == 1:
            continue
        elif len(i) > 1:
            relpaced_list =[x if x not in cog_trans else x + ":"+cog_trans[x] for x in [i[0], i[1]]]
            cog_trans[i] = ", ".join(relpaced_list)

print(cog_trans)
df_query['cog_description'] = df_query['cog_id'].tolist()
df_query['cog_description'] = df_query['cog_description'].map(cog_trans)
df_query.to_csv(os.path.join(out_dir, "Increment_breseq_sub_reannotation.csv"), sep= "," , index=None)
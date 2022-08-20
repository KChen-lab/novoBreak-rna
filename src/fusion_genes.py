import pandas as pd
import numpy as np
import sys

path = sys.argv[1]

count_dict = {}
df = pd.read_csv(path+"/novoBreak.rna.pass.vcf", sep='\t')
for index, row in df.iterrows():
	if (row.iloc[11], row.iloc[13]) in count_dict:
		count_dict[(row.iloc[11], row.iloc[13])] += 1
	else:
		count_dict[(row.iloc[11], row.iloc[13])] = 1
sorted_filtered = sorted(count_dict.items(), key=lambda item: item[1], reverse=True)
result = {k: v for k, v in sorted_filtered}

dicts = {}
for key in result:
	if (key in dicts):
		dicts[key] += int(result[key])
	elif (key[::-1] in dicts):
		dicts[key[::-1]] += int(result[key])
	else:
		dicts[key] = int(result[key])

df = pd.DataFrame(data=dicts, index=[0])
df = (df.T)
df.columns = ['The number of isoforms fusion']

df.to_excel('fusion_genes.xlsx')

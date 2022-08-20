#!/usr/bin/python

import numpy as np
import pandas as pd
import sys
import re

inputvcf = sys.argv[1]
position = sys.argv[2]

ribosomal_protein = pd.read_csv(position+"/source/ribosomal_proteins.txt", sep="\t")
rib_p = {}
for idx, r in ribosomal_protein.iterrows():
	if r.iloc[1] in rib_p:
		pass
	else:
		rib_p[r.iloc[1]] = 1

rRNAs = pd.read_csv(position+"/source/rRNAs.txt", sep="\t")
rrna = {}
for idx, r in rRNAs.iterrows():
	if r.iloc[1] in rrna:
		pass
	else:
		rrna[r.iloc[1]] = 1

gene_pair = pd.read_csv(position+"/source/hsapiens.Pairs.Relaxed.2R.txt", sep="\t")
pair_dic = {}
for idx, r in gene_pair.iterrows():
	if (r.iloc[2], r.iloc[3]) in pair_dic:
		pass
	elif (r.iloc[3], r.iloc[2]) in pair_dic:
		pass
	else:
		pair_dic[(r.iloc[2], r.iloc[3])] = 1

unchar = pd.read_csv(position+"/source/uncharacterized.txt", sep="\t")
unch = {}
for idx, r in unchar.iterrows():
        if r.iloc[5] in unch:
                pass
        else:
                unch[r.iloc[5]] = 1
        for i in str(r.iloc[6]).split(','):
                i = i.replace(" ","")
                if i in unch:
                        pass
                else:
                        unch[i] = 1

immun = pd.read_csv(position+"/source/immunoglobulin.txt", sep="\t")
imm = {}
for idx, r in immun.iterrows():
        if r.iloc[5] in imm:
                pass
        else:
                imm[r.iloc[5]] = 1
        for i in str(r.iloc[6]).split(','):
                i = i.replace(" ","")
                if i in imm:
                        pass
                else:
                        imm[i] = 1

paralog = pd.read_csv(position+"/source/paralog_clusters.dat", sep="\t")
para = {}
for idx, r in paralog.iterrows():
        newlist = [x for x in r if str(x) != 'nan']
        if len(newlist) == 1:
                continue
        for i in range(len(newlist)-1):
                for j in range(i+1, len(newlist)):
                        if (r[i], r[j]) in para:
                                pass
                        else:
                                para[(r[i], r[j])] = 1

homo_gene = pd.read_csv(position+"/source/homo_human.tsv", sep="\t")
homo = {}
length = len(homo_gene)
for idx, r in homo_gene.iterrows():
        if idx == length-1:
                break
        count = 0
        while r[0] == homo_gene.iloc[idx+1+count, 0]:
                count += 1
                if idx+1+count > length-1:
                        break
        if count == 0:
                continue
        else:
                for i in range(count):
                        if (r[3], homo_gene.iloc[idx+i+1, 3]) in homo:
                                pass
                        else:
                                homo[(r[3], homo_gene.iloc[idx+i+1, 3])] = 1


df = pd.read_csv(inputvcf, sep='\t', header=None)
print('aaa',len(df))

count_dict = {}
for index, row in df.iterrows():
	row[11] = re.split(r",|;|/", row.iloc[11])
	row[13] = re.split(r",|;|/", row.iloc[13])
	if set(row[11]) & set(row[13]):
		df.drop(index, inplace=True)
		continue
	else:
		row.iloc[11] = row[11][0]
		row.iloc[13] = row[13][0]

	if row.iloc[11] in rib_p or row.iloc[13] in rib_p:
		df.drop(index, inplace=True)
		continue
	if row.iloc[11] in rrna or row.iloc[13] in rrna:
		df.drop(index, inplace=True)
		continue
	if (row.iloc[11], row.iloc[13]) in pair_dic:
		df.drop(index, inplace=True)
		continue
	elif (row.iloc[13], row.iloc[11]) in pair_dic:
		df.drop(index, inplace=True)
		continue
	#if row.iloc[11] in unch or row.iloc[13] in unch:
	#	df.drop(index, inplace=True)
	#	continue
	if row.iloc[11] in imm and row.iloc[13] in imm:
		df.drop(index, inplace=True)
		continue
	if (row.iloc[11], row.iloc[13]) in homo:
		df.drop(index, inplace=True)
		continue
	elif (row.iloc[13], row.iloc[11]) in homo:
		df.drop(index, inplace=True)
		continue
	if (row.iloc[11], row.iloc[13]) in para:
		df.drop(index, inplace=True)
		continue
	elif (row.iloc[13], row.iloc[11]) in para:
		df.drop(index, inplace=True)
		continue
	if row.iloc[11].split("-AS")[0] == row.iloc[13].split("-AS")[0]:
		df.drop(index, inplace=True)
		continue
	if row.iloc[11].split("-IT")[0] == row.iloc[13].split("-IT")[0]:
		df.drop(index, inplace=True)
		continue
	if row.iloc[11].startswith("MIR") or row.iloc[13].startswith("MIR"):
		df.drop(index, inplace=True)
		continue
	if row.iloc[11].startswith("LOC") or row.iloc[13].startswith("LOC"):
		df.drop(index, inplace=True)
		continue
	if row.iloc[11].startswith("LIN") or row.iloc[13].startswith("LIN"):
		df.drop(index, inplace=True)
		continue
	if re.search(r'P\d$', row.iloc[11]):
		if row.iloc[11][:-2] == row.iloc[13] or row.iloc[11][:-2] == row.iloc[13][:-2]:
			df.drop(index, inplace=True)
			continue
	elif re.search(r'P\d$', row.iloc[13]):
		if row.iloc[13][:-2] == row.iloc[11] or row.iloc[11][:-2] == row.iloc[13][:-2]:
			df.drop(index, inplace=True)
			continue
	if (row.iloc[11], row.iloc[13]) in count_dict:
		df.loc[index, 'mark'] = 2
		df.loc[count_dict[(row.iloc[11], row.iloc[13])], 'mark'] = 2
	elif (row.iloc[13], row.iloc[11]) in count_dict:
		df.loc[index, 'mark'] = 2
		df.loc[count_dict[(row.iloc[13], row.iloc[11])], 'mark'] = 2
	else:
		count_dict[(row.iloc[11], row.iloc[13])] = index
		df.loc[index, 'mark'] = 1

	df.at[index, 11] = row.iloc[11]
	df.at[index, 13] = row.iloc[13]

df.to_csv('mark.vcf', sep='\t', index=False)

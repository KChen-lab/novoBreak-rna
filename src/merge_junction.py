#!/usr/bin/python

import pandas as pd
import glob
import re
import sys
import csv

position = sys.argv[1]

df = pd.read_csv(position+"/temp_splice_pass.tsv", sep="\t", header=None)
total_row = len(df.index) - 1

table = {}
drop_idx = []
for index, row in df.iterrows():
	if index in drop_idx:
		continue
	temp = []
	temp.append([index, row[6], row[17], row[18]])
	n = 1
	if index < total_row:
		while row[0] == df.iloc[index+n][0] and row[1] >= df.iloc[index+n][1] - 5:
			if row[2] >= df.iloc[index+n][2] - 5 and row[2] <= df.iloc[index+n][2] + 5:
				temp.append([index+n, df.iloc[index+n][6], df.iloc[index+n][17], df.iloc[index+n][18]])
			if index+n != total_row:
				n += 1
			else:
				break
	if len(temp) == 1:
		table[index] = row
	else:
		fst = []
		sec = []
		thd = []	
		rdm = []
		for i in range(len(temp)):
			drop_idx.append(temp[i][0])
			if re.match("GT-AG", temp[i][1], re.IGNORECASE):
				if not fst:
					fst.append(temp[i][0])
					fst.append(temp[i][2]) # quality
					fst.append(temp[i][3]) # number of reads
				else:
					if fst[1] < temp[i][2] and fst[2] < temp[i][2]:
						fst[0] = temp[i][0]
						fst[1] = temp[i][2]	
						fst[2] = temp[i][3]
			elif not fst and re.match("GC-AG", temp[i][1], re.IGNORECASE):
				if not sec:
					sec.append(temp[i][0])
					sec.append(temp[i][2]) # quality
					sec.append(temp[i][3]) # number of reads
				else:
					if sec[1] < temp[i][2] and sec[2] < temp[i][2]:
						sec[0] = temp[i][0]
						sec[1] = temp[i][2]	
						sec[2] = temp[i][3]
			elif not fst and not sec and re.match("AT-AC", temp[i][1], re.IGNORECASE):
				if not thd:
					thd.append(temp[i][0])
					thd.append(temp[i][2]) # quality
					thd.append(temp[i][3]) # number of reads
				else:
					if thd[1] < temp[i][2] and thd[2] < temp[i][2]:
						thd[0] = temp[i][0]
						thd[1] = temp[i][2]	
						thd[2] = temp[i][3]
			elif not fst and not sec and not thd:
				if not rdm:
					rdm.append(temp[i][0])
					rdm.append(temp[i][2]) # quality
					rdm.append(temp[i][3]) # number of reads
				else:
					if rdm[1] < temp[i][2] and rdm[2] < temp[i][2]:
						rdm[0] = temp[i][0]
						rdm[1] = temp[i][2]	
						rdm[2] = temp[i][3]
				
		if fst:
			table[fst[0]] = df.iloc[fst[0]]
		elif sec:
			table[sec[0]] = df.iloc[sec[0]]
		elif thd:
			table[thd[0]] = df.iloc[thd[0]]
		elif rdm:
			table[rdm[0]] = df.iloc[rdm[0]]
		
			
df1 = pd.DataFrame.from_dict(table, orient = 'index')
df1 = df1.drop(df1.columns[[3,4,7,8,9,11,12,13]], axis=1)
df1.to_csv(position+"/splice_novoRNABreak_pass.tsv", sep='\t', index = False, header = False)


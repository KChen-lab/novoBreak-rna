#!/usr/bin/python

import pandas as pd
import glob
import re
import sys
import csv

position = sys.argv[1]

table = {}
known = {}
for name in glob.glob(position+"/sj_new_*_annotated_novobreak.tsv"):
	with open(name) as file:
		tsv_file = csv.reader(file, delimiter="\t")
		for index, row in enumerate(tsv_file):
			#if index%10000 == 0:
			#	print(name, index)
			#if index not in table and (re.match("[GN][TN]-[AN][GN]", row[6], re.IGNORECASE) or re.match("[GN][CN]-[AN][GN]", row[6], re.IGNORECASE) or re.match("[AN][TN]-[AN][CN]", row[6], re.IGNORECASE)):
			if index not in table and (re.match("GT-AG", row[6], re.IGNORECASE) or re.match("GC-AG", row[6], re.IGNORECASE) or re.match("AT-AC", row[6], re.IGNORECASE)):
				table[index] = row	
			#elif index in table and (re.match("[GN][TN]-[AN][GN]", row[6], re.IGNORECASE) or re.match("[GN][CN]-[AN][GN]", row[6], re.IGNORECASE) or re.match("[AN][TN]-[AN][CN]", row[6], re.IGNORECASE)):
			elif index in table and (re.match("GT-AG", row[6], re.IGNORECASE) or re.match("GC-AG", row[6], re.IGNORECASE) or re.match("AT-AC", row[6], re.IGNORECASE)):
				if row[-1] < table[index][-1]:
					table[index] = row
				elif row[-1] == table[index][-1]:
					#if re.match("[GN][TN]-[AN][GN]", row[6], re.IGNORECASE):
					if re.match("GT-AG", row[6], re.IGNORECASE):
						table[index] = row
					#elif re.match("[GN][CN]-[AN][GN]", row[6], re.IGNORECASE) and re.match("[AN][TN]-[AN][CN]", table[index][6], re.IGNORECASE):
					elif re.match("GC-AG", row[6], re.IGNORECASE) and re.match("AT-AC", table[index][6], re.IGNORECASE):
						table[index] = row
			if re.match("DA", row[10], re.IGNORECASE):
				known[index] = 1

for key in known:
	table.pop(key, None)
		
df1 = pd.DataFrame.from_dict(table, orient = 'index')
df1.to_csv(position+"/merged_total.csv", sep='\t', index = False)

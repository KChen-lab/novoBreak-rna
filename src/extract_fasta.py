#!/usr/bin/python
import sys

position = sys.argv[1]

with open(position+'/ssake.ctg.fa', 'r') as f:
	lines = f.readlines()

with open(position+'/unmapped.txt', 'r') as id_file:
	id_list = id_file.readlines()
	

output = []
i = 0
while i < len(lines):
	if lines[i].startswith(">") and lines[i].split("\t")[0].split(">")[1] in id_list:
		output.append(lines[i].strip())
		i += 1
		output.append(lines[i].strip())
	else:
		i += 1

with open('ssake.ctg.unmap.fa', 'w') as f:
	f.write('\n'.join(output))

#!/usr/bin/python
import sys

position = sys.argv[1]

with open(position+'/ssake.ctg.fa', 'r') as f:
	lines = f.readlines()

output1 = []
output2 = []
i = 0
while i < len(lines):
	if lines[i].startswith(">") and int(lines[i].split("|")[2].split("size")[1]) >= 300:
		output1.append(lines[i].strip())
		i += 1
		output1.append(lines[i].strip())
	elif lines[i].startswith(">") and int(lines[i].split("|")[2].split("size")[1]) < 300:
		output2.append(lines[i].strip())
		i += 1
		output2.append(lines[i].strip())
	else:
		i += 1

with open('ssake.ctg.long.fa', 'w') as f:
	f.write('\n'.join(output1))
with open('ssake.ctg.short.fa', 'w') as f:
	f.write('\n'.join(output2))

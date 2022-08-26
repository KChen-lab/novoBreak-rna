#!/usr/bin/python

import sys
import subprocess
import os
import re

inputVCF = sys.argv[1]
position = sys.argv[2]
mode = sys.argv[3]

outputAnnoPath = sys.argv[1].replace("vcf","annovar")

outputAnno = open(outputAnnoPath, "w")
with open(inputVCF, 'r') as reader:
    for line in reader.readlines():
        if not line.startswith("#"):
            cols = line.split("\t")
            outputAnno.write(cols[0]+"\t"+str(cols[1])+"\t" +
                            str(cols[1])+"\t"+str(0)+"\t"+str(0)+"\t"+line)
reader.close()
outputAnno.close()
subprocess.run([position+"/app/annotate_variation.pl", "-geneanno", "-dbtype", "refGene",
            "-buildver", "hg38", outputAnnoPath,  position+"/humandb_hg38"])

inputAnnoPath = outputAnnoPath+".variant_function"
outputAnno_second_Path = sys.argv[1].replace("vcf", "second.annovar")
with open(inputAnnoPath, 'r') as reader:
    outputAnno_second = open(outputAnno_second_Path, "w")
    for line in reader.readlines():
        line = line.strip()
        cols = line.split("\t")
        annoA = cols[0]
        geneA = re.sub("[\(].*?[\)]", "", cols[1])
        scols = cols[11].split(";")
        chr = ""
        pos = 0
        for scol in scols:
            if scol.startswith("CHR2"):
                chr = scol.split("=")[1]
            if scol.startswith("END"):
                pos = scol.split("=")[1]
        originalString = "\t".join(str(x) for x in cols[7:])
        outputAnno_second.write(chr+"\t"+str(pos)+"\t"+str(pos)+
                                "\t"+str(0)+"\t"+str(0)+"\t"+originalString+"\t"+annoA+"\t"+geneA+"\n")
outputAnno_second.close()
reader.close()
subprocess.run([position+"/app/annotate_variation.pl", "-geneanno", "-dbtype", "refGene",
                "-buildver", "hg38",  outputAnno_second_Path,  position+"/humandb_hg38"])

outputAnno_all_Path = sys.argv[1].replace("vcf", "all.annovar")
inputAllAnnoPath = outputAnno_second_Path+".variant_function"
with open(inputAllAnnoPath, 'r') as reader:
    outputAnno_all = open(outputAnno_all_Path, "w")
    for line in reader.readlines():
        line = line.strip()
        cols = line.split("\t")
        cols = line.split("\t")
        annoB = cols[0]
        geneB = re.sub("[\(].*?[\)]", "", cols[1])
        annoA = cols[-2]
        geneA = cols[-1]
        if mode == 'fusion':
                if geneA != geneB and annoA!='intergenic' and annoB!='intergenic':
                    originalString = "\t".join(str(x) for x in cols[7:])
                    outputAnno_all.write(originalString+"\t"+annoB+"\t"+geneB+"\n")
        elif mode == 'splice':
                if geneA == geneB and annoA!='intergenic' and annoB!='intergenic':
                    originalString = "\t".join(str(x) for x in cols[7:])
                    outputAnno_all.write(originalString+"\t"+annoB+"\t"+geneB+"\n")


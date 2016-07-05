# -*- coding: utf-8 -*-
"""
Created on Tue Jul 05 11:22:58 2016

@author: wuji
"""

################################
## Generate kabat numbering 
################################
import csv
import re

kabat_v_input = open("C:\Software\igblast-1.5.0\data\internal_data\human\human.ndm.kabat", "r")
outfile = open("C:\Software\igblast-1.5.0\data\germline_kabat.csv",'wb')
writer = csv.writer(outfile)

DELIM=","
newLine ="\n"

gene_pattern = '(^IGHV\d{1,2}-\d{1,3}|^IGHJ\d{1,2})'
def extractGeneName(x):
    m = re.search(gene_pattern, x)
    found = 'NA'        
    if m:
        found=m.group()
    return found
        
for line in kabat_v_input.readlines():
    line=line.rstrip()
    tokens=line.split("\t")
    chain_type = tokens[-2]
    germline = tokens[0]
    gene = extractGeneName(germline)
    fr1_start = int(tokens[1])
    fr1_end = int(tokens[2])
    cdr1_start = int(tokens[3])
    cdr1_end = int(tokens[4])
    fr2_start = int(tokens[5])
    fr2_end = int(tokens[6])
    cdr2_start = int(tokens[7])
    cdr2_end = int(tokens[8])
    fr3_start = int(tokens[9])
    fr3_end = int(tokens[10])
    if chain_type=="VH":
        for p in range(fr1_start, fr3_end+1):
            #FR1
            if p in range(fr1_start, fr1_end+1):
                kabat = (p-fr1_start+3)/3
            #CDR1
            elif p in range(cdr1_start, cdr1_end+1):
                kabat = (p-cdr1_start+3)/3 + 30
                if kabat == 36:
                    kabat = "35A"
                elif kabat == 37:
                    kabat = "35B"
            #FR2
            elif p in range(fr2_start, fr2_end+1):
                kabat = (p-fr2_start+3)/3 + 35
            #CDR2
            elif p in range(cdr2_start, cdr2_end+1):
                if p - cdr2_start < 9:
                    kabat = (p-cdr2_start+3)/3 + 49
                elif cdr2_end - p < 39:     
                    kabat = 66 - (cdr2_end - p + 3)/3
                else:
                    kabat = (p-cdr2_start+3)/3 + 49
                    if kabat == 53:
                        kabat = "52A"
                    elif kabat ==54:
                        kabat = "52B"
                    elif kabat ==55:
                        kabat = "52C"
                    else:
                        sys.stderr.write("HCDR2 exceeds maximum length: " + gene)   
            #FR3
            elif p in range(fr3_start, fr3_end+1):
                kabat = (p-fr3_start+3)/3 + 65
                if kabat == 83:
                    kabat = "82A"
                elif kabat == 84:
                    kabat = "82B"
                elif kabat == 85:
                    kabat = "83C"
                elif kabat > 85:
                    kabat = kabat - 3
                    
        
            kabat="H" + str(kabat)
            record = [germline, gene, p, kabat]
            writer.writerow(record)

outfile.close()

############################################
## Generate kabat numbering for IGHJ
############################################
from Bio import SeqIO
import csv


outfile = open("C:\Software\igblast-1.5.0\data\germline_kabat.csv",'ab')
writer = csv.writer(outfile)
jh_fasta = "C:/Software/igblast-1.5.0/data/imgt/fasta/human_IGHJ"

if os.path.isfile(jh_fasta):
    fasta_dict=SeqIO.index(jh_fasta,"fasta")

kabat_j_input = open("C:\Software\igblast-1.5.0\data\optional_file\human_gl.aux", "r")
count = 0 
kabatj = {}
for line in kabat_j_input.readlines():
    count = count +  1
    if count > 3:
        line=line.rstrip()
        tokens=re.split(r"\t|\s+", line)
        germline = tokens[0]
        gene = extractGeneName(germline)
        frame = tokens[1]
        kabatj[germline]=frame

for germline in fasta_dict:
    size = len(fasta_dict[germline].seq)
    frame = int(kabatj[germline])
    bf = (size - frame) % 3
    start = frame
    gene = extractGeneName(germline)
    for i in range(size, 0, -1):
        if i > end:
            kabat = "NA"
        elif i <= start:
            kabat = "NA"
        else:
            kabat = "H" + str(113 - (end - i + 2) / 3)
        record = [germline, gene, str(i), kabat]    
        writer.writerow(record)
    
outfile.close()    

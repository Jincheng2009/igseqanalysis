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
import sys
import os
from Bio import SeqIO

kabat_v_input = open("C:\Software\igblast-1.5.0\data\internal_data\human\human.ndm.kabat", "r")
outfile = open("C:\Software\igblast-1.5.0\data\germline_kabatVL.csv",'wb')
writer = csv.writer(outfile)

igkv_fasta = "C:/Software/igblast-1.5.0/data/imgt/fasta/human_IGKV"
iglv_fasta = "C:/Software/igblast-1.5.0/data/imgt/fasta/human_IGLV"

if os.path.isfile(igkv_fasta):
    temp1=SeqIO.to_dict(SeqIO.parse(igkv_fasta,"fasta"))
if os.path.isfile(iglv_fasta):
    temp2=SeqIO.to_dict(SeqIO.parse(iglv_fasta,"fasta"))

fasta_dict = dict(temp1.items() + temp2.items())    
    
DELIM=","
newLine ="\n"

        
for line in kabat_v_input.readlines():
    line=line.rstrip()
    tokens=line.split("\t")
    chain_type = tokens[-2]
    germline = tokens[0]
    gene = germline.split("*")[0]
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
    if (chain_type=="VL" or chain_type=="VK") and re.match("^IG[KL]", germline):
        if not fasta_dict.has_key(germline):
            print(germline + " does not exist")
            continue
        nsize = len(fasta_dict[germline].seq)
        for p in range(fr1_start, nsize+1):
            #FR1
            if p in range(fr1_start, fr1_end+1):
                kabat = (p-fr1_start+3)/3
            #CDR1
            elif p in range(cdr1_start, cdr1_end+1):
                if p - cdr1_start < 12:
                    kabat = (p-cdr1_start+3)/3 + 23
                elif cdr1_end - p < 21:     
                    kabat = 35 - (cdr1_end - p + 3)/3
                else:
                    kabat = (p-cdr1_start+3)/3 + 23
                    if kabat == 28:
                        kabat = "27A"
                    elif kabat == 29:
                        kabat = "27B"
                    elif kabat == 30:
                        kabat = "27C"
                    elif kabat == 31:
                        kabat = "27D"
                    elif kabat == 32:
                        kabat = "27E"
                    elif kabat == 33:
                        kabat = "27F"
                    else:
                        sys.stderr.write("LCDR2 exceeds maximum length: " + gene)   
            #FR2
            elif p in range(fr2_start, fr2_end+1):
                kabat = (p-fr2_start+3)/3 + 34
            #CDR2
            elif p in range(cdr2_start, cdr2_end+1):
                kabat = (p-cdr2_start+3)/3 + 49
                
            #FR3
            elif p in range(fr3_start, fr3_end+1):
                kabat = (p-fr3_start+3)/3 + 56
            
            #CDR3
            elif p in range(fr3_end+1, nsize+1):
                kabat = (p-fr3_end+3)/3 + 88
                if kabat == 96:
                    kabat = "95A"
                elif kabat == 97:
                    kabat = "95B"
                elif kabat == 98:
                    kabat = "95C"
                elif kabat == 99:
                    kabat = "95D"
                elif kabat == 100:
                    kabat = "95E"
                elif kabat == 101:
                    kabat = "95F"
                       
            kabat="L" + str(kabat)
            record = [germline, gene, p, kabat]
            writer.writerow(record)
            
outfile.close()

############################################
## Generate kabat numbering for IGHJ
############################################

outfile = open("C:/Software/igblast-1.5.0/data/germline_kabatVL.csv",'ab')
writer = csv.writer(outfile)
igkj_fasta = "C:/Software/igblast-1.5.0/data/imgt/fasta/human_IGKJ"
iglj_fasta = "C:/Software/igblast-1.5.0/data/imgt/fasta/human_IGLJ"

if os.path.isfile(igkj_fasta):
    temp1=SeqIO.to_dict(SeqIO.parse(igkj_fasta,"fasta"))
if os.path.isfile(iglj_fasta):
    temp2=SeqIO.to_dict(SeqIO.parse(iglj_fasta,"fasta"))

fasta_dict = dict(temp1.items() + temp2.items())  

kabat_j_input = open("C:\Software\igblast-1.5.0\data\optional_file\human_gl.aux", "r")
count = 0 
kabatj = {}
for line in kabat_j_input.readlines():
    count = count +  1
    if count > 3:
        line=line.rstrip()
        tokens=re.split(r"\t|\s+", line)
        germline = tokens[0]
        gene = germline.split("*")[0]
        frame = tokens[1]
        kabatj[germline]=frame

for germline in fasta_dict:
    size = len(fasta_dict[germline].seq)
    frame = int(kabatj[germline])
    bf = (size - frame) % 3
    start = frame
    gene = germline.split("*")[0]
    for i in range(start, size + 1):
        kabat = (i - start + 3) / 3 + 95
        if (size - frame)/3 == 13 and kabat == 107:
            kabat = "106A"
        if (size - frame)/3 == 13 and kabat > 107:
            kabat = kabat - 1
        kabat = "L" + str(kabat)
        record = [germline, gene, str(i+1), kabat]    
        writer.writerow(record)
    
outfile.close()    

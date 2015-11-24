import os
import pandas as pd
from Bio import SeqIO
import re
from Bio.Seq import Seq

datapath="C:/Java/data/genewiz/"
file="sample7.tophit.tsv"
filepattern="cdr.h3.pattern.txt"

## Load reference data for extracting CDR.H3
blast_table = pd.read_table(datapath+file, sep='\t',header=None)
blast_table.head()

table = blast_table.drop_duplicates([0])

cdr3pattern= pd.read_table(datapath+filepattern, sep='\t',header=None)

cdr3pattern.columns = ['ref', 'pattern', 'hcdr3']
table.columns = ['id', 'ref', 'mismatch', 'align_length','evalue']

seqref = pd.merge(table, cdr3pattern, on='ref')

seqrefDict = {}
cdr3patternDict = {}

for index, row in seqref.iterrows():
    seqrefDict[row['id']]=row['ref']

for index, row in cdr3pattern.iterrows():
    cdr3patternDict[row['ref']]=row['pattern']
    
## Load sample Fasta
samplefile = "sample7.fasta"
count = 0 
printcount = 0
output = "sample7.CDR.H3.fasta"
fileout = open(datapath+output, "w")

with open(datapath+samplefile) as filein:
    for line in filein:
        line = line.rstrip()
        if line.startswith('>'):
            fastaId = line.strip()
            count +=1
        else:
            if not seqrefDict.has_key(fastaId[1:]):
                continue
            ref = seqrefDict[fastaId[1:]]
            pattern = cdr3patternDict[ref]
            seq = Seq(line.strip())
            # Forward search
            found = re.search(pattern, str(seq))
            # Reverse search
            if not found:
                found = re.search(pattern, str(seq.reverse_complement()))
            if found:
                fileout.write(fastaId + "|" + ref + "\n")
                ##trim off the regex framework region
                sequence = found.group(0)[12:len(found.group(0))-12]
                fileout.write(sequence + "\n")
            if count >= printcount:
                print "processed: " + str(count)
                printcount += 10000
            
fileout.close()        
        
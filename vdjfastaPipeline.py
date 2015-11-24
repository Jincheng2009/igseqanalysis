import subprocess
import os, sys
import re
import json
from seqanalysis import countUnique
from seqanalysis import extractCDR3
from seqanalysis import extractGermlineCount
from seqanalysis import pairByFastaID
from seqanalysis import clusterUnique
from seqanalysis import formatCluster

pypath="C:/java/uk.co.catplc.deepseq.analysis/pycode/"
deepseqpath="C:/java/uk.co.catplc.deepseq.analysis/build/classes"
datapath="C:/Java/data/sarav2/"
vdjfasta="vdjfasta/"

cdrfiles=[ f for f in os.listdir(datapath+vdjfasta) if os.path.splitext(f)[1]==".txt"]

readchain_map={}
readchain_map["R1"]="VL"
readchain_map["R2"]="VH"

# Connect to metadata file to update metadata information
metafname=datapath+'metadata.json'
if os.path.isfile(metafname) and os.path.getsize(metafname)>0:
    metafile = open(metafname, 'r+')
    samples = json.load(metafile)
    metafile.close()
else:
    samples=[]

# Parse the experiment conditions from file names and return a dict of metainformation
def parseSampleName(name, info):
    tokens = re.split("-|_", name)
#   For Sarav's dataset  (vdjfasta pipeline)  
    info["donor"]=tokens[0]
    info["time"]=re.findall("\d+", tokens[1])[0]
    info["protocol"]=re.findall("[a-zA-Z]+", tokens[1])[0]
    info["sampleID"]=tokens[3]
    info["region"]=readchain_map[tokens[5]]    
    return info

def vdjCDR2fasta(file_in, file_out):
    if os.path.isfile(file_out) and os.path.getsize(file_out)>0:
       print file_out + " exists, skip"
       return
    inputfile = open(file_in)
    output = open(file_out, "w")
    
    for line in inputfile:
        tokens=re.split('\t| ', line)
        seqid = tokens[0]
        seq = tokens[1]
        output.write('>' + seqid + '\n')
        output.write(seq)

    inputfile.close()
    output.close()

vdjsamples=[]
for f in cdrfiles:
    vdjsample = parseSampleName(f, {})
    for sample in samples:
        if sample["sampleID"] == vdjsample["sampleID"] and \
           sample["region"] == vdjsample["region"]:
                sample["vdjCDR3"] = vdjfasta + f
                vdjsamples.append(sample)
                break
             
  
for sample in samples:
    key = "vdjCDR3"
    if sample.has_key(key):
        fname = sample["basename"] + ".CDR3.fasta"
        filein = datapath + sample[key]
        fileout = datapath + vdjfasta + fname
        vdjCDR2fasta(filein, fileout)
        sample["vdjCDR3_fasta"] = vdjfasta + fname
    
    
# Pair the VH and VL CDR3 sequence
sampleVL={}
sampleVH={}
for sample in samples:
    if sample["region"]=="VH" and sample.has_key("vdjCDR3_fasta"):
        sampleVH[sample["sampleID"]]=sample
    elif sample["region"]=="VL" and sample.has_key("vdjCDR3_fasta"):
        sampleVL[sample["sampleID"]]=sample
        
for sid in sampleVH:
    if sampleVL.has_key(sid):
        key="vdjCDR3_fasta"
        if not sampleVL[sid].has_key(key) or not sampleVH[sid].has_key(key):
            print sample["annotation"] + " does not have associated file for: " + key 
            continue
        # 4.1 Pair VH and VL CDR3 into new fasta file
        fileVH=datapath + sampleVH[sid][key]
        fileVL=datapath + sampleVL[sid][key]
        fname=sid + ".paired.CDR3.vdjfasta.prot"
        file_out=datapath+vdjfasta+fname
        pairByFastaID(fileVH, fileVL, file_out)   
        sampleVH[sid]["vdj_paired_CDR3"] = vdjfasta + fname
        sampleVL[sid]["vdj_paired_CDR3"] = vdjfasta + fname
        # Count Unique
        file_in=file_out
        fname=sid + ".paired.CDR3.vdjfasta.count"
        file_out=datapath+vdjfasta+fname
        countUnique(file_in, file_out, sid)
        sampleVH[sid]["vdj_paired_CDR3_count"] = vdjfasta + fname
        sampleVL[sid]["vdj_paired_CDR3_count"] = vdjfasta + fname
    
# Overwrite the metadata file to update
# Close and reopen for writing
metafile=open(metafname, "w")
json.dump(samples, metafile, indent=2, separators=(',', ': '))
metafile.flush()
metafile.close()
    
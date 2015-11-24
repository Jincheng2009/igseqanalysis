import os
import sys
import subprocess
import json
import re
from seqanalysis import countUnique
from seqanalysis import extractCDR3
from seqanalysis import extractGermlineCount
from seqanalysis import pairByFastaID
from seqanalysis import clusterUnique
from seqanalysis import formatCluster

pypath="C:/java/uk.co.catplc.deepseq.analysis/pycode/"
deepseqpath="C:/java/uk.co.catplc.deepseq.analysis/build/classes"
datapath="C:/Java/data/sarav2/"

prot="prot/"
germline="germline/"
CDR="CDR/"
cluster="cluster/"

if not os.path.exists(datapath+germline):
    os.mkdir(datapath+germline)

if not os.path.exists(datapath+CDR):
    os.mkdir(datapath+CDR)
    
if not os.path.exists(datapath+cluster):
    os.mkdir(datapath+cluster)

# Get the prot fasta files in the folder
protfiles=[ f for f in os.listdir(datapath+prot) if os.path.splitext(f)[1]==".prot"]

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
#   For Sarav's dataset    
    info["donor"]=tokens[0]
    info["time"]=re.findall("\d+", tokens[1])[0]
    info["protocol"]=re.findall("[a-zA-Z]+", tokens[1])[0]
    info["sampleID"]=tokens[2]
    info["region"]=readchain_map[tokens[3]]    
#     For george's dataset from NCBI
#     info["donor"]=tokens[0]
#     info["time"]="NA"
#     info["protocol"]="em"
#     info["sampleID"]=tokens[0]
#     info["region"]=tokens[1]   
    return info

# 1. Parse the metadata information from the data file names 
for f in protfiles:
    sample={}
    sample["annotation"]=prot+f
    sample["basename"]=str.split(f,".")[0]
    sample = parseSampleName(sample["basename"], sample)
    # Replace the metadata for the file, if the file is not recorded in 
    # metadatafile, append a new record to the metadata
    exist=False
    for index, item in enumerate(samples):
        if item["annotation"]==sample["annotation"]:
            for k in sample:
                samples[index][k]=sample[k]
            exist=True
    if not exist:
        samples.append(sample)

# 2. Extract germline count
for sample in samples:
    file_in=datapath+sample["annotation"]  
    fname=sample["basename"]+".germline.count"
    file_out=datapath+germline+fname
    extractGermlineCount(file_in, file_out)
    sample["germline_count"]=germline+fname
    
# 3. Extract CDR3 sequence
for sample in samples:
    file_in=datapath+sample["annotation"]
    fname=sample["basename"] + "_" + sample["region"] + ".CDR3.fasta.prot"
    file_out=datapath+CDR+fname
    extractCDR3(file_in, file_out, sample["region"])
    sample["CDR3"]=CDR+fname

# 4. Pair the VH and VL CDR3 sequence
sampleVL={}
sampleVH={}
for sample in samples:
    if sample["region"]=="VH":
        sampleVH[sample["sampleID"]]=sample
    elif sample["region"]=="VL":
        sampleVL[sample["sampleID"]]=sample
        
for sid in sampleVH:
    if sampleVL.has_key(sid):
        key="CDR3"
        if not sampleVL[sid].has_key(key) or not sampleVH[sid].has_key(key):
            print sample["annotation"] + " does not have associated file for: " + key 
            continue
        # 4.1 Pair VH and VL CDR3 into new fasta file
        fileVH=datapath + sampleVH[sid][key]
        fileVL=datapath + sampleVL[sid][key]
        fname=sid + ".paired.CDR3.fasta.prot"
        file_out=datapath+CDR+fname
        pairByFastaID(fileVH, fileVL, file_out)   
        sampleVH[sid]["paired_CDR3"] = CDR + fname
        sampleVL[sid]["paired_CDR3"] = CDR + fname

# 5. Count unique paired CDR3 sequences        
for sample in samples:
    key="paired_CDR3"
    if not sample.has_key(key):
        print sample["annotation"] + " does not have associated file for: " + key 
        continue
    file_in = datapath + sample[key]
    fname=sample["sampleID"] + ".paired.CDR3.count"
    file_out = datapath+CDR+fname
    countUnique(file_in, file_out, sample["sampleID"])
    sample["paired_CDR3_count"] = CDR + fname

# 6. Cluster paired CDR3 sequences by CDR.H3
for sample in samples:
    key="paired_CDR3"
    if not sample.has_key(key):
        print sample["annotation"] + " does not have associated file for: " + key 
        continue
    file_in = datapath + sample[key]
    fname=sample["sampleID"] + ".paired.CDR.H3.uc"
    file_out = datapath+cluster+fname
    clusterUnique(file_in, file_out, sample["sampleID"], 0, 25)
    sample["CDR_H3_uc"] = cluster + fname

# 7. Cluster paired CDR3 sequences by CDR.H3
for sample in samples:
    key="CDR_H3_uc"
    if not sample.has_key(key):
        print sample["annotation"] + " does not have associated file for: " + key 
        continue
    file_in = datapath + sample[key]
    fname=sample["sampleID"] + ".paired.CDR3.H3.cluster"
    file_out = datapath+cluster+fname
    formatCluster(file_in, file_out)
    sample["CDR_H3_cluster"] = cluster + fname


# 8. Cluster paired CDR3 sequences by CDR.L3
for sample in samples:
    key="paired_CDR3"
    if not sample.has_key(key):
        print sample["annotation"] + " does not have associated file for: " + key 
        continue
    file_in = datapath + sample[key]
    fname=sample["sampleID"] + ".paired.CDR.L3.uc"
    file_out = datapath+cluster+fname
    clusterUnique(file_in, file_out, sample["sampleID"], 26, 43)
    sample["CDR_L3_uc"] = cluster + fname

# 9. Cluster paired CDR3 sequences by CDR.L3
for sample in samples:
    key="CDR_L3_uc"
    if not sample.has_key(key):
        print sample["annotation"] + " does not have associated file for: " + key 
        continue
    file_in = datapath + sample[key]
    fname=sample["sampleID"] + ".paired.CDR3.L3.cluster"
    file_out = datapath+cluster+fname
    formatCluster(file_in, file_out)
    sample["CDR_L3_cluster"] = cluster + fname

# Overwrite the metadata file to update
# Close and reopen for writing
metafile=open(metafname, "w")
json.dump(samples, metafile, indent=2, separators=(',', ': '))
metafile.flush()
metafile.close()





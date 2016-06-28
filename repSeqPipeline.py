import os
import json
import re
from seqanalysis import countUnique
from seqanalysis import extractCDR3
from seqanalysis import extractGermlineCount
from seqanalysis import pairByFastaID
from seqanalysis import clusterUnique
from seqanalysis import formatCluster
from seqanalysis import extractGermlineCDRCount

# Point datapath to the folder that contain the fasta file
datapath="C:/Java/data/sarav3/"
# Input data files containing fasta format files after antibody annotation 
# by uk.co.catplc.seqanalysis.blaze2.FastaAlignedProt or FastaAlignedProtVL
prot="prot/"
# Output folder for germline counting files
germline="germline/"
# Output folder for CDR sequence unique counting file
CDR="CDR/"
# Output folder for clustered CDR sequences and their counting files
cluster="cluster/"
# Output folder for unpaired CDR sequences
unpaired="unpaired/"

# Create output folders if they do not exist
if not os.path.exists(datapath+germline):
    os.mkdir(datapath+germline)
if not os.path.exists(datapath+CDR):
    os.mkdir(datapath+CDR)
if not os.path.exists(datapath+cluster):
    os.mkdir(datapath+cluster)
if not os.path.exists(datapath+unpaired):
    os.mkdir(datapath+unpaired)

# Get the prot fasta files in the folder
protfiles=[ f for f in os.listdir(datapath+prot) if os.path.splitext(f)[1]==".prot"]

readchain_map={}
# R1 is Light Chain
readchain_map["R1"]="VL"
# R2 is Heavy Chain
readchain_map["R2"]="VH"

# If metadata file exists, load it and keep updating metadata information into it along the pipeline analysis
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
    info["time"]="NA" #re.findall("\d+", tokens[1])[0]
    info["protocol"]=tokens[1]
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
    extractGermlineCount(file_in, file_out, sample["region"])
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

# 6. Cluster paired CDR3 sequences by CDR-H3:CDR-L3
for sample in samples:
    key="paired_CDR3"
    if not sample.has_key(key):
        print sample["annotation"] + " does not have associated file for: " + key 
        continue
    file_in = datapath + sample[key]
    fname=sample["sampleID"] + ".paired.CDR3.HL.uc"
    file_out = datapath+cluster+fname
    clusterUnique(file_in, file_out, sample["sampleID"], 0, 43)
    sample["CDR3_HL_uc"] = cluster + fname

# 7. Format cluster results
for sample in samples:
    key="CDR3_HL_uc"
    if not sample.has_key(key):
        print sample["annotation"] + " does not have associated file for: " + key 
        continue
    file_in = datapath + sample[key]
    fname=sample["sampleID"] + ".paired.CDR3.HL.cluster"
    file_out = datapath+cluster+fname
    formatCluster(file_in, file_out)
    sample["CDR3_HL_cluster"] = cluster + fname

# 8. Count unique CDR3 without pairing R1 and R2 by fasta ID
for sample in samples:
    key="CDR3"
    if not sample.has_key(key):
        print sample["annotation"] + " does not have associated file for: " + key 
        continue
    file_in = datapath + sample[key]
    fname=sample["sampleID"] + "." + sample["region"] + ".CDR3.count"
    file_out = datapath+unpaired+fname
    countUnique(file_in, file_out, sample["sampleID"])
    sample["CDR3_count"] = unpaired + fname

# 9. Cluster unpaired CDR3 sequences
for sample in samples:
    key="CDR3"
    if not sample.has_key(key):
        print sample["annotation"] + " does not have associated file for: " + key 
        continue
    if sample["region"]=="VH":
        start=0
        end=25
    elif sample["region"]=="VL":
        start=0
        end=18
    file_in = datapath + sample[key]
    fname=sample["sampleID"] + "." + sample["region"] + ".CDR3.uc"
    file_out = datapath+unpaired+fname
    clusterUnique(file_in, file_out, sample["sampleID"], start, end)
    sample["unpaired_CDR3_uc"] = unpaired + fname

# 10. Format cluster results
for sample in samples:
    key="unpaired_CDR3_uc"
    if not sample.has_key(key):
        print sample["annotation"] + " does not have associated file for: " + key 
        continue
    file_in = datapath + sample[key]
    fname=sample["sampleID"] + "." + sample["region"] + ".CDR3.cluster"
    file_out = datapath+unpaired+fname
    formatCluster(file_in, file_out)
    sample["unpaired_CDR3_cluster"] = unpaired + fname

# 11. Get CDR3 with germline information
for sample in samples:
    key="CDR3"
    if not sample.has_key(key):
        print sample["annotation"] + " does not have associated file for: " + key 
        continue
    file_in = datapath + sample[key]
    fname=sample["sampleID"] + "." + sample["region"] + ".germline.CDR3.count"
    file_out = datapath+germline+fname
    extractGermlineCDRCount(file_in, file_out)
    sample["germline_CDR3_count"] = germline + fname

# 12. Get paired CDR3 with heavy and light chain germline information
for sample in samples:
    key="paired_CDR3"
    if not sample.has_key(key):
        print sample["annotation"] + " does not have associated file for: " + key 
        continue
    file_in = datapath + sample[key]
    fname=sample["sampleID"] + ".paired.germline.CDR3.count"
    file_out = datapath+germline+fname
    extractGermlineCDRCount(file_in, file_out, True)
    sample["paired_germline_CDR3_count"] = germline + fname


# Overwrite the metadata file to update
# Close and reopen for writing
metafile=open(metafname, "w")
json.dump(samples, metafile, indent=2, separators=(',', ': '))
metafile.flush()
metafile.close()





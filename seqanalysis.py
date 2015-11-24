import subprocess
import os, sys

pypath="C:/java/uk.co.catplc.deepseq.analysis/pycode/"
deepseqpath="C:/java/uk.co.catplc.deepseq.analysis/build/classes"

def extractGermlineCount(file_in, file_out):
    if os.path.isfile(file_out) and os.path.getsize(file_out)>0:
        print file_out + " exists, skip"
        return
    print "Extracting germline counts for " + file_in 
    data_input = open(file_in)
    output = open(file_out, "w")
    print "Output to: " + file_out
    # proc = subprocess.Popen(["head"], stdin=input, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc = subprocess.Popen(["python", pypath+"countGermline.py", "--igh", "--igl", "--igk"], stdin=data_input, stdout=output, stderr=subprocess.PIPE)
    for line in proc.stderr:
        sys.stdout.write(line)
    proc.wait()
    output.flush()

def extractCDR3(file_in, file_out, atype):
    if os.path.isfile(file_out) and os.path.getsize(file_out)>0:
        print file_out + " exists, skip"
        return
    data_input = open(file_in)
    output = open(file_out, "w")
    print "Extracting CDR3 sequence for " + file_in 
    print "Output to: " + file_out
    start = 0
    end = 0
    if atype == "VH":
        start = 109
        end = 11
    elif atype == "VL":
        start = 103
        end = 13
    proc1 = subprocess.Popen(["java", "-cp", deepseqpath, "uk.co.catplc.deepseq.analysis.TrimFasta", str(start), str(end)], 
                             stdin=data_input, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    proc2 = subprocess.Popen(["java", "-cp", deepseqpath, "uk.co.catplc.deepseq.analysis.FilterByFeature", "."], 
                              stdin=proc1.stdout, stdout=output, stderr=subprocess.PIPE, shell=True)
    for line in proc1.stderr:
        print line
    for line in proc2.stderr:
        print line
    proc1.wait()
    proc2.wait()
    output.flush()
    output.close()

def pairByFastaID(fileVH, fileVL, file_out):
    if os.path.isfile(file_out) and os.path.getsize(file_out)>0:
        print file_out + " exists, skip"
        return
    output = open(file_out, "w")
    print "Pairing CDR3 sequence of " + fileVH + " and " + fileVL 
    print "Output to: " + file_out
    proc = subprocess.Popen(["java", "-cp", deepseqpath, "uk.co.catplc.deepseq.analysis.ConcatenateFasta", 
                             "-l", fileVH, "-r", fileVL, "-germline", "-delim", ":"], 
                             stdout=output, stderr=subprocess.PIPE, shell=True)
    for line in proc.stderr:
        print line

def countUnique(file_in, file_out, sampleID):
    if os.path.isfile(file_out) and os.path.getsize(file_out)>0:
        print file_out + " exists, skip"
        return    
    data_input = open(file_in)
    output = open(file_out, "w")
    print "Counting unique sequence " + file_in 
    print "Output to: " + file_out
    proc1 = subprocess.Popen(["java", "-cp", deepseqpath, "uk.co.catplc.deepseq.analysis.Unique", sampleID],
                            stdin=data_input, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    proc2 = subprocess.Popen(["java", "-cp", deepseqpath, "uk.co.catplc.deepseq.analysis.Cluster2Count"],
                            stdin=proc1.stdout, stdout=output, stderr=subprocess.PIPE, shell=True)
    for line in proc1.stderr:
        print line
    for line in proc2.stderr:
        print line
    proc1.wait()
    proc2.wait()
    output.flush()
    output.close()
    
def clusterUnique(file_in, file_out, sampleID, start, end, identity=0.92):
    if os.path.isfile(file_out) and os.path.getsize(file_out)>0:
        print file_out + " exists, skip"
        return    
    data_input = open(file_in)
    output = open(file_out, "w")
    print "Cluster unique sequence " + file_in + " by identity " + str(identity)
    print "Output to: " + file_out
    proc1 = subprocess.Popen(["java", "-cp", deepseqpath, "uk.co.catplc.deepseq.analysis.Unique", sampleID],
                            stdin=data_input, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    clusterCmds = ["java", "-cp", deepseqpath, "uk.co.catplc.deepseq.analysis.FindClusterByCount",
                   "-c", "0", "-id", str(identity), "-v", "-s", str(start), "-e", str(end)]
    proc2 = subprocess.Popen(clusterCmds, stdin=proc1.stdout, stdout=output, stderr=subprocess.PIPE, shell=True)
    for line in iter(proc1.stderr.readline, b''):
        print line.rstrip()
    for line in iter(proc2.stderr.readline, b''):
        print line.rstrip()
    proc1.wait()
    proc2.wait()
    output.flush()
    output.close()
   
def formatCluster(file_in, file_out):
    if os.path.isfile(file_out) and os.path.getsize(file_out)>0:
        print file_out + " exists, skip"
        return    
    data_input = open(file_in)
    output = open(file_out, "w")
    print "Format the uc file " + file_in 
    print "Output to: " + file_out
    proc = subprocess.Popen(["python", pypath+"ucFormat.py"], stdin=data_input, stdout=output, stderr=subprocess.PIPE)
    for line in iter(proc.stderr.readline, b''):
        print line.rstrip()
    proc.wait()
    output.flush()
    output.close()    
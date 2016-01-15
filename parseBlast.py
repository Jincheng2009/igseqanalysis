import sys
import re
import os
from Bio import SeqIO
import csv
import getopt


def analyzeMutation(qstart, qseq, refstart, refseq, fastaid, ref, strandforward):
    if len(qseq) != len(refseq):
        sys.stderr.write("Error :" + fastaid + " does not have same length as reference " + ref)
    mutations =[]
    querygap = 0
    refgap=0
    for i in range(0,len(qseq)-1):
        details = []
        if refseq[i] == '-':
            refgap += 1
        if qseq[i] == '-':
            querygap += 1
        if qseq[i]!=refseq[i] and refseq[i]!='-' and qseq[i]!='-':
            details.append(fastaid)
            details.append(ref)
            details.append(str(qstart+i-querygap))
            if strandforward:
                details.append(str(refstart+i-refgap))
            else:
                details.append(str(refstart-i+refgap))
            details.append(qseq[i])
            details.append(refseq[i])
            details.append(getNeighborBase(qseq, i, -3))
            details.append(getNeighborBase(qseq, i, -2))
            details.append(getNeighborBase(qseq, i, -1))
            details.append(getNeighborBase(qseq, i, 1))
            details.append(getNeighborBase(qseq, i, 2))
            details.append(getNeighborBase(qseq, i, 3))
            mutations.append(details)
    return mutations

def getNeighborBase(seq, idx, shift):
    seq = re.sub("-","",seq)
    if (idx + shift) < 0 or (idx+shift) >=len(seq):
        return "NA"
    
    return seq[idx+shift].upper()

def main(argv):
    extractFastq=False
    extractCoverage=False
    try:
        opts, args = getopt.getopt(argv,"h", ["blast=","fastq=","mutation=","coverage="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            usage()
            sys.exit()
        elif opt == "--fastq":
            extractFastq=True
            fastqfile = arg
        elif opt == "--blast":
            blastfile = arg
        elif opt == "--mutation":
            fileout = arg
        elif opt == "--coverage":
            extractCoverage=True
            coverage_out = arg
    
    fileout = open(fileout, "wb")
    writer = csv.writer(fileout)
    # writer.writerow('id, germline, query_pos, ref_pos, query_base, ref_base, b3, b2, b1, a1, a2, a3, mismatch, length, gap, phred')
    
    if extractCoverage:
        coveragefile = open(coverage_out, "w")
                            
    if extractFastq and os.path.isfile(fastqfile):
        fastq_dict=SeqIO.index(fastqfile,"fastq")
    
    count = 0
    printNext = 0
    coverage = {}
    strandPlus=True
    with open(blastfile) as filein:
        
        for line in filein:
            line=line.rstrip()
            line=line.strip()
            if line.startswith("Query="): 
                nextRefLength = False
                fastaid=line[7:].split(' ')[0]
                quality=[]
                count += 1
                if extractFastq and fastq_dict.has_key(fastaid):
                    quality = fastq_dict[fastaid].letter_annotations["phred_quality"]
                if count >= printNext:
                    sys.stderr.write("Processed: " + str(count) + "\n")
                    printNext = printNext + 10000
    
            if line.startswith(">"):
                ref = line[1:] 
                nextRefLength = True
            if line.startswith("Length=") and nextRefLength:
                nextRefLength = False
                try:
                    found = re.search('Length=(\d+)', line).group(1)
                except AttributeError:
                    # AAA, ZZZ not found in the original string
                    found = 'NA'
                refLength = int(found)
                if not coverage.has_key(ref):
                    refInfo = {}
                    for i in range(1,refLength):
                        refInfo[i]=0
                    coverage[ref] = refInfo
                
                
            if line.startswith("Identities ="):
                try:
                    found = re.search('Identities = (\d+/\d+).*', line).group(1)
                except AttributeError:
                    # AAA, ZZZ not found in the original string
                    found = 'NA'
                tokens = found.split('/')
                mismatch = str(int(tokens[1]) - int(tokens[0]))
                length = tokens[1]
                try:
                    found = re.search('Gaps = (\d+/\d+).*', line).group(1)
                except AttributeError:
                    # AAA, ZZZ not found in the original string
                    found = 'NA'
                gap = found.split('/')[0]
                
            if line.startswith("Strand="):
                if line=="Strand=Plus/Minus":
                    strandPlus=False
                else:
                    strandPlus=True
    
            if line.startswith("Query") and not line.startswith("Query="):
                qtokens = line.split(" ")
                qtokens = filter(lambda a: a!="", qtokens)
                qstart = int(qtokens[1])
                qseq = qtokens[2].upper()
    
            if line.startswith("Sbjct"):
                stokens = line.split(" ")
                stokens = filter(lambda a: a!="", stokens)
                sstart = int(stokens[1])
                sseq = stokens[2].upper()
                send = int(stokens[3])
                # Update coverage based on alignment
                if sstart > send:
                    a1 = send
                    a2 = sstart
                else:
                    a1 = sstart
                    a2 = send
                for i in range(a1,a2):
                    coverage.get(ref)[i] += 1
                records = analyzeMutation(qstart, qseq, sstart, sseq, fastaid, ref, strandPlus)
                if len(records) > 0:
                    for mutation in records:
                        mutation.append(mismatch)
                        mutation.append(length)
                        mutation.append(gap)
                        if extractFastq:
                            mutation.append(quality[int(mutation[2])])
                        if strandPlus:
                            mutation.append(str(qstart))
                            mutation.append(str(qstart + length-1))
                        else:
                            mutation.append(str(qstart-length+1))
                            mutation.append(str(qstart))
                        writer.writerow(mutation)
    
    # Output the coverage report
    if extractCoverage:
        for key in coverage:
            for idx in coverage[key]:
                coveragefile.write(key+","+str(idx)+","+str(coverage.get(key)[idx])+"\n");
    coveragefile.close()
    
    fileout.close()


def usage():
    print 'python parseBlast.py [--fastq fastq_file] [--blast blastn_input_file] [--mutation mutation_info_outputfile] [--coverage coverage_report_file]'

if __name__ == "__main__":
    main(sys.argv[1:])

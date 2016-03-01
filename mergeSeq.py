import os
import sys
import getopt

datapath="C:/Java/data/sarav4/fasta/"
protfile="testR2.fasta.prot"
cdrfile="testR2_CDRL3.fasta"

cdr3dict = {}

def checkMatch(template, query):
    count = 0
    for i in range(0,len(query)):
        if query[i] != "." and template[i] != ".":
            count = count + 1
            if query[i]!=template[i]:
                return False
    if count == 0:
        return False
    else:
        return True        

def main(argv):
    try:
        opts, args = getopt.getopt(argv,"hp:c:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            usage()
            sys.exit()
        elif opt == "-p":
            protfile = arg
        elif opt == "-c":
            cdrfile = arg
    with open(cdrfile) as filein:
        for line in filein:
            line = line.rstrip()
            if line.startswith(">"):
                seqid=line
            else:
                seq=line
                cdr3dict[seqid]=seq 
    with open(protfile) as filein:
        for line in filein:
            line = line.rstrip()
            if line.startswith(">"):
                fullid=line
                seqid = line.split("|")[0].strip();
            else:
                seq=line
                if cdr3dict.has_key(seqid):
                    cdr3 = cdr3dict[seqid]
                    seqList = list(seq)
                    if checkMatch(seq[102:], cdr3):                 
                        for i in range(0,len(cdr3)):
                            seqList[102+i]=cdr3[i]
                        seq = "".join(seqList)
                    print fullid
                    print seq

if __name__ == "__main__":
    main(sys.argv[1:])
    
def usage():
    print 'python mergeSeq.py -p full.fasta.prot -c cdr.fasta > output.txt'    
    
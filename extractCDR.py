import sys
import re
import csv
import getopt
import collections

complement={}
complement['A'] = 'T'
complement['T'] = 'A'
complement['G'] = 'C'
complement['C'] = 'G'

CDRH1start = 115 + 1
CDRH1end = 130 + 1
CDRH3start = 319 + 1
CDRH3end = 334 + 1
CDRL1start = 517 + 1
CDRL1end = 532 + 1

def extractBlock(qstart, qend, qseq, refstart, refend, refseq, cdr, strandforward, region="CDRH1"):
    if len(qseq) != len(refseq):
        sys.stderr.write("Error : does not have same length as reference ")
    region_start = CDRH1start
    region_end = CDRH1end
    if region=="CDRH3":
        region_start = CDRH3start
        region_end = CDRH3end
    elif region=="CDRL1":
        region_start = CDRL1start
        region_end = CDRL1end
    
    # if no overlap, simply return the original dict
    if strandforward:
        if refstart > region_end or refend < region_start:
            return cdr
    else:
        if refend > region_end or refstart < region_start:
            return cdr     
    # Has overlap
    ref_position = refstart    
    for i in range(0,len(qseq)):
        if refseq[i] != '-' and strandforward:
            ref_position = refstart + i
        elif refseq[i]!='-' and not strandforward:
            ref_position = refstart - i
        if ref_position>=region_start and ref_position<region_end:
            if strandforward:
                cdr[ref_position]=cdr[ref_position].join(qseq[i])
            else:
                cdr[ref_position]=cdr[ref_position].join(complement[qseq[i]])
    return cdr

def getSequence(cdr):
    cdr = collections.OrderedDict(sorted(cdr.items()))
    outseq=""
    for k, v in cdr.iteritems():
        outseq+=str(v)
    outseq = outseq.replace("-","")
    return outseq
    

def main(argv):
    readFromFile=False
    try:
        opts, args = getopt.getopt(argv,"h", ["blast=","fastq=","output=","coverage="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            usage()
            sys.exit()
        elif opt == "--blast":
            blastfile = arg
            readFromFile=True
        elif opt == "--output":
            fileout = arg
    
    fileout = open(fileout, "wb")
    writer = csv.writer(fileout)
    # writer.writerow('id, germline, query_pos, ref_pos, query_base, ref_base, b3, b2, b1, a1, a2, a3, mismatch, length, gap, phred')
    
    count = 0
    printNext = 0
    coverage = {}
    strandPlus=True
    refstart=1000
    refend=0
    block1={}
    for i in range(CDRH1start, CDRH1end):
        block1[i]=""
    block2={}
    for i in range(CDRH3start, CDRH3end):
        block2[i]=""
    block3={}
    for i in range(CDRL1start, CDRL1end):
        block3[i]=""
                
    if readFromFile:
        filein = open(blastfile)
    else:
        filein = sys.stdin

    previousid="" 
    block1seq=""
    block2seq=""
    block3seq=""
       
    for line in filein:
        line=line.rstrip()
        line=line.strip()
        if line.startswith("Query="):
            nextRefLength = False
            fastaid=line[7:].split(' ')[0]
            count += 1
            if count >= printNext:
                sys.stderr.write("Processed: " + str(count) + "\n")
                printNext = printNext + 10000
        if count==1:
            previousid=fastaid
            
        if count>1 and fastaid!=previousid:
            # Write out     
            record=[]
            record.append(previousid)
            block1seq = getSequence(block1)
            block2seq = getSequence(block2)
            block3seq = getSequence(block3)
            record.append(block1seq)
            record.append(block2seq)
            record.append(block3seq)
            writer.writerow(record)
            # Initialize again
            refstart=1000
            refend=0
            block1={}
            for i in range(CDRH1start, CDRH1end):
                block1[i]=""
            block2={}
            for i in range(CDRH3start, CDRH3end):
                block2[i]=""
            block3={}
            for i in range(CDRL1start, CDRL1end):
                block3[i]=""
            previousid=fastaid

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
            qend = int(qtokens[3])

        if line.startswith("Sbjct"):
            stokens = line.split(" ")
            stokens = filter(lambda a: a!="", stokens)
            sstart = int(stokens[1])
            sseq = stokens[2].upper()
            send = int(stokens[3])
            # Update coverage based on alignment
            a1=0
            a2=0
            if sstart > send:
                a1 = send
                a2 = sstart
            else:
                a1 = sstart
                a2 = send
            if a1 < refstart:
                refstart = a1
            if a2 > refend:
                refend = a2
            block1 = extractBlock(qstart, qend, qseq, sstart, send, sseq, block1, strandPlus, "CDRH1")
            block2 = extractBlock(qstart, qend, qseq, sstart, send, sseq, block2, strandPlus, "CDRH3")
            block3 = extractBlock(qstart, qend, qseq, sstart, send, sseq, block3, strandPlus, "CDRL1")
    
    
    # Write out the last  
    record=[]
    record.append(previousid)
    block1seq = getSequence(block1)
    block2seq = getSequence(block2)
    block3seq = getSequence(block3)
    record.append(block1seq)
    record.append(block2seq)
    record.append(block3seq)
    writer.writerow(record)
    
    fileout.close()


def usage():
    print 'python extractCDR.py [--blast blastn_input_file] [--output outputfile]'

if __name__ == "__main__":
    main(sys.argv[1:])
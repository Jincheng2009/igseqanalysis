import sys
import re
import os
from Bio import SeqIO
import getopt
from sequtility import Alignment
from sequtility import Sequence
import csv 

MUTATION="mutation"
CDR="CDR"

def usage():
    print 'Parse the output from IgBlastn to extract CDRs (CDR1, CDR2, CDR3), mutations and coverage'
    print 'Support standard input and file input. Optional coverage output is written in a file specified by --coverage option'
    print 'Annotation region or mutation is output into standard output'
    print 'Use --type option to choose output type'
    print 'cat igblastn_outputfile.txt | python parseIgBlast.py -t mutation [--fastq fastq_file] [--coverage coverage_report_file] > mutation.csv'
    print '-b --blast\t input file from output of ncbi-igblast (default: stdin)'
    print '-f --fastq\t fastq file (optional to extract Phred score of mutation)'
    print '-t --type \t Output type (mutation or CDR), output file is CSV format'
    print '-c --coverage\t Outpufile (optional for coverage for each alignment)'

    
def main(argv):
    extractFastq=False
    extractCoverage=False
    extractID=False
    readFromFile=False
    analysisType = ""
    try:
        opts, args = getopt.getopt(argv,"hb:f:t:c:", ["blast=","fastq=","type=","coverage="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            usage()
            sys.exit()
        elif opt in ("-f", "--fastq"):
            extractFastq=True
            fastqfile = arg
        elif opt in ("-b", "--blast"):
            blastfile = arg
            readFromFile=True
        elif opt in ("c", "--coverage"):
            extractCoverage=True
            coverage_out = arg
        elif opt in ("-t", "--type"):
            if arg==MUTATION:
                analysisType = MUTATION
            elif arg==CDR:
                analysisType = CDR
            else:
                sys.stdout.write("Only mutation and CDR analysis type are allowed \n")
                usage()
                sys.exit()
                    
            
    count = 0
    printNext = 0
    strandPlus=True
    
    if readFromFile:
        filein = open(blastfile)
    else:
        filein = sys.stdin

    fastq_dict = {}
    extractGermline=False
    extractStrand=False
    extractAlign=False
    nafter=0
    previous_line=""
    previous_line2 = ""
    
    if extractCoverage:
        coveragefile = open(coverage_out, 'wb')
        coverage_writer = csv.writer(coveragefile)
                            
    if extractFastq and os.path.isfile(fastqfile):
        fastq_dict=SeqIO.index(fastqfile,"fastq")

    for line in filein:
        line = line.strip('\n')
        ## Do not trim the end when parsing the alignment section
        if not extractAlign:
            line=line.rstrip()
        if line.startswith("Query="):
            align_record=None
            inframe=True
            extractGermline=False
            extractStrand=False
            extractAlign=False
            nafter=0
            fastaid=line[7:].split(' ')[0]
            quality=[]
            count += 1
            strandPlus = True
            query_length = 0
            if extractFastq and fastq_dict.has_key(fastaid):
                quality = fastq_dict[fastaid].letter_annotations["phred_quality"]
            if count >= printNext:
                sys.stderr.write("Processed: " + str(count) + "\n")
                printNext = printNext + 10000
        
        if extractID and not line.startswith("Length="):
            line = line.rstrip()
            fastaid=fastaid + line
        
        if line.startswith("Length="):
            extractID=False
            query_length = int(line[7:])

        if line.startswith("Query="):        
            extractID=True
        # End of one alignment result
        # Write out alignment results        
        if line.startswith("Effective search space used:"):
            if align_record is not None:
                if analysisType == MUTATION:
                    for record in align_record.getMutations():
                        outline = ""
                        if len(quality) > 0:
                            phred = quality[int(record[2])]
                        else:
                            phred = "NA"
                        for item in record:
                            outline = outline + "," + str(item)
                        outline = outline[1:] + "," + str(phred) + "\n"
                        sys.stdout.write(outline)
                #align_record.printOut()
                if analysisType == CDR:
                    cdr1 = align_record.getRegion("CDR1")
                    cdr2 = align_record.getRegion("CDR2")
                    cdr3 = align_record.getRegion("CDR3")
                    sys.stdout.write(align_record.getQuery().getName() + "," + vgene + "," + jgene+ "," + cdr1 + "," + cdr2 + "," + cdr3 + "\n")
                if extractCoverage:
                    coverage_record = [align_record.getQuery().getName(), "NA", -1, -1, "NA", -1, -1]
                    if align_record.get_v_sequence() is not None:
                        vseq = align_record.get_v_sequence()
                        coverage_record[1] = vseq.getName()
                        coverage_record[2] = vseq.getRange()[0]
                        coverage_record[3] = vseq.getRange()[1]   
                    if align_record.get_j_sequence() is not None:
                        jseq = align_record.get_j_sequence()
                        coverage_record[4] = jseq.getName()
                        coverage_record[5] = jseq.getRange()[0]
                        coverage_record[6] = jseq.getRange()[1]
                    coverage_writer.writerow(coverage_record)
        ## Extract germline genes
        if line.startswith("Sequences producing significant alignments"):
            extractGermline=True
        if line.startswith("Domain classification requested"):
            extractGermline=False
    
        if extractGermline and line.startswith("IG"):
            gene = line.split(' ')[0]
            if gene.startswith("IGHV") or gene.startswith("IGLV") or gene.startswith("IGKV"):
                vgene = gene
            elif gene.startswith("IGHD"):
                dgene = gene
            elif gene.startswith("IGHJ") or gene.startswith("IGLJ") or gene.startswith("IGKJ"):
                jgene = gene
                
        ## Extract strand information
        if extractStrand:
            nafter =+ 1
        if line.startswith("V-(D)-J rearrangement summary for query sequence"):
            extractStrand=True
            nafter=0
        if nafter == 1 and extractStrand:
            extractStrand=False
            tokens = line.split('\t')
            strand=tokens[-1]
            frame=tokens[-3]
            if strand=="+":
                strandPlus=True
            else:
                strandPlus=False
            if frame!="In-frame":
                inframe=False
    
        ## Extract alignment information
        if extractAlign:
            nafter=+1
        if line=="Alignments" and inframe:
            extractAlign=True
            nafter=0
            align_record = Alignment(fastaid, query_length)
            query_seq = None
            vseq = None
            jseq = None
            astart = 0   
        if extractAlign:
            tokens=line.split(' ')
            # remove empty strings
            tokens=[x for x in tokens if x]
            if len(tokens)>0 and re.match("(lcl\|)?Query_", tokens[0]):
                if query_seq is None:
                    query_seq = Sequence(tokens[2], tokens[1], tokens[3], fastaid, strandPlus)
                    align_record.setQuery(query_seq)
                    astart = 0
                else:
                    astart = len(query_seq.getSequence())
                    query_seq.addSequence(tokens[2], tokens[1], tokens[3])
                offset = line.index(tokens[2])
                # Add translation row
                query_seq.addTranslation(previous_line[offset:])
                # Add annotation row
                align_record.addAnnoation(previous_line2[offset:])
            elif re.match("^V", line) and len(tokens) > 6:
                if vseq is None:
                    vseq = Sequence(tokens[5], tokens[4], tokens[6], tokens[3])
                    align_record.add_v_alignment(vseq, astart)
                else:
                    vseq.addSequence(tokens[5], tokens[4], tokens[6])    
            elif re.match("^J", line) and len(tokens) > 6:
                if jseq is None:
                    jseq = Sequence(tokens[5], tokens[4], tokens[6], tokens[3])
                    align_record.add_j_alignment(jseq, astart)
                else:
                    jseq.addSequence(tokens[5], tokens[4], tokens[6])    
        if extractAlign and line.startswith("Lambda"): 
            extractAlign=False
            
        # Cache previous two lines
        previous_line2 = previous_line
        previous_line = line    

    # Output the coverage report
    if extractCoverage:
        coveragefile.close()

if __name__ == "__main__":
    main(sys.argv[1:])

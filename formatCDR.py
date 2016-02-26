import sys
import getopt

def main(argv):
    try:
        opts, args = getopt.getopt(argv,"h")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            usage()
            sys.exit()
            
    for line in sys.stdin:
        line = line.rstrip()
        if not line.strip():
            continue
        elif line.startswith(">"):
            print line
        else:
            seq = line
            length = len(seq)
            index = length - 2
            fill= "-" * (18 - length)
            newSeq=seq[:index] + fill + seq[index:]
            print newSeq

    
if __name__ == "__main__":
    main(sys.argv[1:])
    
def usage():
    print 'cat CDRL3.fasta | python formatCDR.py'

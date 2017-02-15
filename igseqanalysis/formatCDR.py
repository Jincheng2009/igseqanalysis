import sys
import getopt

def main():
    argv = sys.argv[1:]
    position = 0
    trailing = False
    try:
        opts, args = getopt.getopt(argv,"hp:t:e")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            usage()
            sys.exit()
        elif opt == "-p":
            position = int(arg)
        elif opt == "-t":
            total = int(arg)
        elif opt == "-e":
            trailing = True
                                    
    for line in sys.stdin:
        line = line.rstrip()
        if not line.strip():
            continue
        elif line.startswith(">"):
            print line
        else:
            seq = line
            length = len(seq)
            if trailing:
                index = length - position
            else:
                index = position
            fill= "-" * (total - length)
            newSeq=seq[:index] + fill + seq[index:]
            print newSeq

def usage():
    print 'cat CDRL3.fasta | python formatCDR.py -p 6 -t 23 -e'
    
if __name__ == "__main__":
    main()
    


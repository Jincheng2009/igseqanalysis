import sys
import getopt
import sequtility

def main(argv):
    readFromFile = False
    try:
        opts, args = getopt.getopt(argv,"hi:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            usage()
            sys.exit()
        elif opt == "-i":
            readFromFile = True
            fastafile = arg

    if readFromFile:
        filein = open(fastafile)
    else:
        filein = sys.stdin

    for line in filein:
        line = line.rstrip()
        # Skip empty lines
        if not line:
            continue
        if line.startswith(">"):
            fastaid = line
        # Skip sequences not in frame
        elif len(line) % 3 == 0:
            i = 0
            prot = ""
            seq = line.upper()
            while i + 3 <= len(line):
                if line[i: i + 3] in sequtility.codontable:
                    prot = prot + sequtility.codontable[line[i: i + 3]]
                else:
                    prot = prot + "X"
                i += 3
            if "X" not in prot:
                sys.stdout.write(fastaid + "\n")
                sys.stdout.write(prot + "\n")

if __name__ == "__main__":
    main(sys.argv[1:])
    
def usage():
    print 'cat dna.fasta | python translate.py > prot.fasta'    
    
import sys
import getopt
import sequtility

def main():
    argv = sys.argv[1:]
    readFromFile = False
    idx = []
    try:
        opts, args = getopt.getopt(argv,"hi:p:")
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
        elif opt == "-p":
            idx_str = arg

    idx = map(int, idx_str.split(","))

    if len(idx) == 0:
        usage()
        sys.exit(2)

    if readFromFile:
        filein = open(fastafile)
    else:
        filein = sys.stdin

    for line in filein:
        line = line.rstrip()
        # Skip empty lines
        if not line:
            continue
        line = line.replace("\t",",")
        tokens = line.split(",")
        for i in idx:
            seq = tokens[i].upper()
            prot=""
            if len(seq) % 3 == 0:
                j = 0
                while j + 3 <= len(seq):
                    if seq[j: j + 3] in sequtility.codontable:
                        prot = prot + sequtility.codontable[seq[j: j + 3]]
                    else:
                        prot = prot + "X"
                    j += 3
                if "X" in prot:
                    prot = "."
            else:
                prot = "."
                
            tokens[i] = prot
        record = tokens[0]
        for i in range(1, len(tokens)):
            record += "\t" + tokens[i]
        sys.stdout.write(record + "\n")

def usage():
    print 'cat dna.fasta | python translateTable.py -p 3,4,5,8,9,10 > prot.fasta'
    print '-p \t comma-separated list of column index for translation (0-based index)'    

if __name__ == "__main__":
    main()
    

    
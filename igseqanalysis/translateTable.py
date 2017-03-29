#==============================================================================
#     Copyright (C) 2017  MedImmune, LLC
#     
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#     
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#     
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#==============================================================================

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
            prot="."
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
    

    
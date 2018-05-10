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
    print 'cat CDRL3.fasta | python format_CDR.py -p 6 -t 23 -e'
    
if __name__ == "__main__":
    main()
    


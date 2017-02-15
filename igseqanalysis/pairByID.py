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
import pandas as pd

def main():
    argv = sys.argv[1:]
    file1 = None
    file2 = None
    try:
        opts, args = getopt.getopt(argv,"hl:r:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            usage()
            sys.exit()
        elif opt == "-l":
            file1 = arg
        elif opt == "-r":
            file2 = arg
        elif opt == "-o":
            outfile = arg

    if file1 is None or file2 is None:
        sys.stderr.write("missing input file \n")

    df1 = pd.read_table(file1, sep="\t|,", header=None, engine='python')
    df2 = pd.read_table(file2, sep="\t|,", header=None, engine='python')

    df1_vh = df1[map(lambda x : x.startswith('IGHV'),df1[1])]
    df1_vl = df1[map(lambda x : x.startswith('IGKV') or x.startswith('IGLV'),df1[1])]

    df2_vh = df2[map(lambda x : x.startswith('IGHV'),df2[1])]
    df2_vl = df2[map(lambda x : x.startswith('IGKV') or x.startswith('IGLV'),df2[1])]

    part1 = pd.merge(df1_vh, df2_vl, on = 0)
    part2 = pd.merge(df2_vh, df1_vl, on = 0)

    df = pd.concat([part1, part2])
    df.to_csv(outfile, sep="\t", index=False, header=False, na_rep='.')

def usage():
    print 'python pairByID.py -l vh.csv -r vl.csv> paired.tsv'    

if __name__ == "__main__":
    main()


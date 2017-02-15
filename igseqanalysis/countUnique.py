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
"""
Created on Thu Jan 05 14:58:59 2017

@author: wuji
"""
import sys
import getopt
import pandas as pd

def main():
    argv = sys.argv[1:]
    idx = []
    try:
        opts, args = getopt.getopt(argv,"hi:p:d:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            usage()
            sys.exit()
        elif opt == "-p":
            idx_str = arg

    idx = map(int, idx_str.split(","))

    if len(idx) == 0:
        usage()
        sys.exit(2)

    df = pd.read_table(sys.stdin, sep="\t", header=None)
    if max(idx) > df.shape[1]:
    	sys.stderr.write("column index out of range \n")

    colnames = []
    for i in idx:
    	colnames.append(df.columns[i])

    count_df = df.groupby(colnames).size().reset_index()
    count_df.to_csv(sys.stdout, sep="\t", index=False, header=False, na_rep='.')
    
def usage():
    print 'cat CDR.csv | python countUnique.py -p 5,10 > out.tsv'
    print '-p \t comma-separated list of column index for counting (0-based index)' 

if __name__ == "__main__":
    main()


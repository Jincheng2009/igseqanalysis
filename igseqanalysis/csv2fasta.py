# -*- coding: utf-8 -*-
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
    idx = 5
    sizeout = False
    csvfile=sys.stdin
    try:
        opts, args = getopt.getopt(argv,"hi:p:s")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            usage()
            sys.exit()
        elif opt == "-i":
            csvfile = arg
        elif opt == "-p":
            idx = int(arg)
        elif opt == "-s":
            sizeout = True

    df = pd.read_table(csvfile, sep="\t|,", header=None, engine='python', index_col=False)
    if idx > df.shape[1]:
        sys.stderr.write("column index out of range \n")

    count_df = df.groupby(idx).size().reset_index()
    count_df.columns=[idx, "count"]
    id_df = df.groupby(idx).head(1)
    count_df = pd.merge(id_df, count_df, left_on=idx, right_on=idx)
    count_df.sort_values("count", ascending=False, inplace=True)

    for i, row in count_df.iterrows():
        tokens = row.values
        if len(tokens) <= idx:
            continue
        seq = tokens[idx]
        if pd.isnull(seq):
            continue
        fastaid = tokens[0]
        germline = tokens[1] + "+" + tokens[2]
        seq = seq.replace("-","")
        if len(seq) > 1:
            header = ">" + fastaid + ";" + germline
            if sizeout:
                count = row['count']
                header = header + ";size=" + str(count)
            header = header + "\n"
            sys.stdout.write(header)
            sys.stdout.write(seq + "\n")

def usage():
    print 'cat cdr.csv | python csv2fasta.py -p 5 full.fasta'    

if __name__ == "__main__":
    main()

    
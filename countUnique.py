# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 14:58:59 2017

@author: wuji
"""
import sys
import getopt
import pandas as pd

def main(argv):
    idx = 5
    delim = ","
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
            idx = int(arg)

    df = pd.read_table(sys.stdin, sep="\t", header=None)
    if idx > df.shape[1]:
    	sys.stderr.write("column index out of range \n")

    colnames = df.columns
    df = df[df[idx] != "."]
    count_df = df.groupby(colnames[idx]).size()

    count_df.to_csv(sys.stdout, sep="\t", index=True, header=False, na_rep='.')

if __name__ == "__main__":
    main(sys.argv[1:])
    
def usage():
    print 'cat CDR.csv | python countUnique.py > out.tsv'    


# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 14:01:32 2016

@author: wuji
"""

import pandas as pd
import numpy as np
import sys
import getopt

def usage():
    
    print 'sampleCount.py -i <inputfile.count> -o <outputfile.count> -e <error rate> [-s <seed>] [-d]'
    print '-i\t input file, count format'
    print '-o\t output file, count format'
    print '-m\t sampling size (>10 million may affect memory usage)'
    print '-s\t seed for random number generator'


def main(argv):
    try:
        opts, args = getopt.getopt(argv,"hdxi:o:m:s:",["ifile=","ofile="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-m"):
            nsample=int(arg)
        elif opt in ("-s"):
            seed = int(arg)
            np.random.seed(seed)
    print 'Input file is ', inputfile
    print 'Output file is ', outputfile
    print 'Sample size is ', str(nsample)

    count_df = pd.read_table(inputfile, header=None)
    count_df.columns = ["seq", "count"]
    count_df = count_df.sort_values("count", ascending=False)
    
    count_df.head(10)
    total = sum(count_df["count"])
    count_df["frac"] = count_df["count"] / total
    
    if nsample > total:
        print("sample size is larger than total size")
    
    count_df["sample_count"] = 0
    
    for i in range(count_df.shape[0]):
        print(i)
        count_df.loc[i, "frac"] = float(count_df.loc[i, "count"]) / sum(count_df["count"][i:])
        rands = np.random.rand(int(nsample))
        nseq = sum(rands < count_df["frac"][i])
        count_df.loc[i, "sample_count"] = nseq
        nsample = nsample - nseq
    
    sampled_df = count_df[count_df["sample_count"]>0][["seq", "sample_count"]]
    sampled_df.to_csv(outputfile, sep="\t", header=False, index=False)    

if __name__ == "__main__":
    main(sys.argv[1:])
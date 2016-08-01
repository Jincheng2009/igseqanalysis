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
    inputfile=None
    outputfile=None
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
    sys.stderr.write('\nSample size is ' + str(nsample) + "\n")

    if inputfile is not None:
        count_df = pd.read_table(inputfile, header=None)
    else:
        count_df = pd.read_table(sys.stdin, header=None)
        
    count_df.columns = ["seq", "count"]
    count_df = count_df.sort_values("count", ascending=False)

    total = sum(count_df["count"])
    count_df["frac"] = count_df["count"] / total
    
    if nsample > total:
        sys.stderr.write("sample size is larger than total size\n")
    
    count_df["sample_count"] = 0
    
    complete = 0.0
    for i in range(count_df.shape[0]):
        if float(i)/count_df.shape[0] > complete + 0.1:
            complete = float(i)/count_df.shape[0]
            sys.stderr.write("{0:.0f}%".format(complete * 100) + "\n")
        count_df.loc[i, "frac"] = float(count_df.loc[i, "count"]) / sum(count_df["count"][i:])
        rands = np.random.rand(int(nsample))
        nseq = sum(rands < count_df["frac"][i])
        count_df.loc[i, "sample_count"] = nseq
        nsample = nsample - nseq
    
    sampled_df = count_df[count_df["sample_count"]>0][["seq", "sample_count"]]
    sampled_df.sort_values("sample_count", ascending=False)
    if outputfile is not None:
        sampled_df.to_csv(outputfile, sep="\t", header=False, index=False)
    else:
        sampled_df.to_csv(sys.stdout, sep="\t", header=False, index=False)   
    sys.stderr.write("{0:.0f}%".format(100) + "\n")
    
if __name__ == "__main__":
    main(sys.argv[1:])
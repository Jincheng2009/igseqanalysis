# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 14:58:59 2017

@author: wuji
"""
import sys
import getopt

def main(argv):
    idx = 5
    readFromFile = False
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
            csvfile = arg
        elif opt == "-p":
            idx = int(arg)

    if readFromFile:
        filein = open(csvfile)
    else:
        filein = sys.stdin

    for line in filein:
        line = line.rstrip()
        line = line.replace("\t",",")
        tokens = line.split(",")
        if len(tokens) < 5 or len(tokens[5]) == 0:
            continue
        fastaid = tokens[0]
        germline = tokens[1] + "+" + tokens[2]
        cdr3 = tokens[idx]
        cdr3 = cdr3.replace("-","")
        if len(cdr3) > 1:
            sys.stdout.write(">" + fastaid + "||" + germline + "\n")
            sys.stdout.write(cdr3 + "\n")

if __name__ == "__main__":
    main(sys.argv[1:])
    
def usage():
    print 'python csv2Fasta.py -i full.csv'    
    
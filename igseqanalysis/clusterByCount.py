# -*- coding: utf-8 -*-
"""
Created on Thu Jun 01 14:15:24 2017

@author: wuji
"""
import sys
import getopt

def diff_1bp_more(seq1, seq2):
    diff = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            diff += 1
        if diff > 1:
            return True
    return False

def count_diff(seq1, seq2):
    n1 = len(seq1)
    n2 = len(seq2)
    if n1 != n2:
        return min(n1, n2)
    diff = sum(1 for a,b in zip(seq1, seq2) if a!=b)
    return diff

def find_parent(seq, seq_list, count_dict):
    idx = -1
    j = 0 #len(seq_list) - 1
    while j < len(seq_list) and idx < 0:
        seq2 = seq_list[j]
        ## 2-folder difference is the cutoff
        if not diff_1bp_more(seq, seq2) and count_dict[seq] * 2 < count_dict[seq2]:
            idx = j
        j += 1
    return idx

def main():
    argv = sys.argv[1:]
    readFromFile = False
    cutoff = 2.
    try:
        opts, args = getopt.getopt(argv,"hi:c:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            usage()
            sys.exit()
        elif opt == "-i":
            readFromFile = True
            infile = arg
        elif opt == "-c":
            cutoff = int(arg)
            
    if readFromFile:
        filein = open(infile)
    else:
        filein = sys.stdin
            
    seq_list = []
    seq_count = []
    count_dict = {}
    for line in filein:
        line = line.rstrip()
        # Skip empty lines
        if not line:
            continue
        line = line.replace("\t",",")
        tokens = line.split(",")
        seq_list.append(tokens[0])
        seq_count.append(int(tokens[1]))
        count_dict[tokens[0]] = int(tokens[1])
    
    # Sort the sequences by descending count value  
    count_seq = zip(seq_count, seq_list)
    count_seq.sort(reverse=True)
    seq_list = [seq for count, seq in count_seq]   
    seq_count = [count for count, seq in count_seq]
    filein.close()  
    ##
    ## Start the clustering algorithm
    ##    
    # List of cluster centroids
    c_list = [seq_list[0]]
    # Total count of each cluster
    c_count = [seq_count[0]]
    # Memembership of each cluster
    # Key is the cluster centroids
    # Value is the list of children sequences
    c_member = {}
    c_member[seq_list[0]] = []
    # Max depth of cluster - max hamming distance of child sequences to centroid 
    c_depth ={}
    c_depth[seq_list[0]] = 0
    
    i = 1
    nextp = 0.
    total = len(seq_list)
    
    while i < total:
        if float(i) / total > nextp:
            print("Complete " + str(float(i) / total * 100) + "%")
            nextp += 0.1
        seq = seq_list[i]
        count = seq_count[i]
        ## Start search
        idx = -1
        # Search in the centroids
        j = 0
        while j < len(c_list) and idx < 0:
            c = c_list[j]
            diff = count_diff(seq, c)
            if len(c) != len(seq):
                j += 1
                continue
            if count_dict[c] <= cutoff * count_dict[seq]:
                j += 1
                continue
            if  diff > c_depth[c] + 1:
                j += 1
                continue
            if diff <= 1:
                idx = j
            elif find_parent(seq, c_member[c], count_dict) >= 0:
                idx = j
            j += 1
        if idx >= 0:
            c = c_list[idx]
            c_member[c].append(seq)
            c_count[idx] += count
            if count_diff(seq, c) > c_depth[c]:
                c_depth[c] = count_diff(seq, c)
        else:
            c_list.append(seq)
            c_count.append(count)
            c_member[seq] = []
            c_depth[seq] = 0
        i += 1
    
    for c in c_list:
        record = c + "\t" + c + "\t" + str(count_dict[c]) + "\n"
        sys.stdout.write(record)
        for hit in c_member[c]:
            record = hit + "\t" + c + "\t" + str(count_dict[hit]) + "\n"
            sys.stdout.write(record)
        
    

def usage():
    print 'cat cdr.count.csv | python clusterByCount.py -c 2 > cluster.count'

if __name__ == "__main__":
    main()

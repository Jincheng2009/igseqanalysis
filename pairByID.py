import sys
import getopt
import pandas as pd
pd.set_option('display.max_colwidth', -1)

def main(argv):
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

    df1 = pd.read_table(file1, sep=",", header=None)
    df2 = pd.read_table(file2, sep=",", header=None)

    df1_vh = df1[map(lambda x : x.startswith('IGHV'),df1[1])]
    df1_vl = df1[map(lambda x : x.startswith('IGKV') or x.startswith('IGLV'),df1[1])]

    df2_vh = df2[map(lambda x : x.startswith('IGHV'),df2[1])]
    df2_vl = df2[map(lambda x : x.startswith('IGKV') or x.startswith('IGLV'),df2[1])]

    part1 = pd.merge(df1_vh, df2_vl, on = 0)
    part2 = pd.merge(df2_vh, df1_vl, on = 0)

    df = pd.concat([part1, part2])
    df.to_csv(outfile, sep="\t", index=False, header=False, na_rep='.')

if __name__ == "__main__":
    main(sys.argv[1:])
    
def usage():
    print 'python pairByID.py -l vh.csv -r vl.csv> out.csv'    

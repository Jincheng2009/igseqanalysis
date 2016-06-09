import pandas as pd
from natsort import natsorted
import sys
import getopt

def main(argv):
    append=False
    try:
        opts, args = getopt.getopt(argv,"ha", ["output="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            usage()
            sys.exit()
        elif opt == "--output":
            fileout = arg
        elif opt == "-a":
            append=True
    germline = pd.read_csv("/home/wuji/germlinedb/germline.csv")
    
    ##########Coverage
    df = pd.read_csv(sys.stdin, header=None)
    df.columns = ['gene', 'ref_pos', 'depth']
    
    df['gene'] = df['gene'].map(str.strip)
    
    # 2. Merge germline dataframe with coverge dataframe
    df = pd.merge(df, germline)
    df = df.sort(['gene','position'])
    
    # 3. Get dataframe with number of coverage for each kabat number for each sequence
    df_count = df.groupby(['gene','isotype','label'])['depth'].sum()
    df_count = df_count.reset_index()
    
    df_count['label'] = df_count['label'].map(lambda g : "H" + g)
    df_count['gene'] = df_count['gene'].map(lambda g : g.split("*")[0])
    
    df_count = df_count[['gene','label','depth','isotype']]

    # 4. Output data into csv file
    if append:
        df_count.to_csv(fileout, mode='a', index=False, header=True)
    else:
        df_count.to_csv(fileout, index=False)
    
def usage():
    print 'cat 0-IgM_S13.merged.coverage.csv | python calculateDepth.py [--output outputfile]'

if __name__ == "__main__":
    main(sys.argv[1:])
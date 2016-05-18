import pandas as pd
from natsort import natsorted
import sys
import getopt

def main(argv):
    germline = pd.read_csv("/home/wuji/germlinedb/germline.csv")
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
    ######################################
    ## Custom import data section
    ######################################  
    df = pd.read_csv(sys.stdin, header=None)

    df.columns = ['id', 'germline', 'query_pos', 'ref_pos', 'query_base', 'ref_base', 'b3', 'b2', 'b1', 'a1', 'a2', 'a3', 'mismatch', 'length', 'gap', 'phred', 'refstart','refend']
    
    # 1. Remove mutation with "N"
    n1 = df.shape[0]
    df = df[df['query_base']!="N"]
    df = df[df['phred']>=10]
    n2 = df.shape[0]
    print "removing " + str(n1 - n2) + " rows"
    df['germline'] = df['germline'].map(str.strip)
    
    # 2. Merge germline dataframe with mutation dataframe
    df = pd.merge(df, germline, left_on=['germline', 'ref_pos'], right_on = ['gene','position'] )
    df = df.sort(['id','query_pos'])
    
    # 3. Get dataframe with number of mutation for each kabat number for each sequence
    df_count = df[['id','label']]
    df_count['count'] = 1
    df_count = df_count.groupby(['id','label'])['count'].sum()
    df_count = df_count.reset_index()
    
    # 3.temp Get the germline alignment length data
    #df_germline = df[['id','isotype','germline','length']].drop_duplicates()
    #df_germline = df_germline[df_germline['germline'].str.startswith("IGHV")]
    #df_germline[df_germline['isotype']=='IgG']['length'].describe()
    #df_germline[df_germline['isotype']=='IgM']['length'].describe()
    
    # 4. Reorder columns
    count_table = df_count.pivot(index="id", columns="label", values="count")
    cols = count_table.columns.tolist()
    cols = natsorted(cols)
    count_table = count_table[cols]
    cols_final = map(lambda g : "H" + g, cols)
    count_table.columns = cols_final
    count_table = count_table.fillna(0)
    count_table['id'] = count_table.index
    
    # 5. Join kabat mutation table with sequence information
    df1 = df[['id', 'germline']].drop_duplicates()
    df1 = df1[df1['germline'].str.startswith("IGHV")]
    count_table = pd.merge(count_table, df1, on=['id'])
    count_table['gene'] = count_table['germline'].map(lambda g : g.split("*")[0])
    
    # 6. Output data into csv file
    if append:
        count_table.to_csv(fileout, mode='a', index=False, header=True)
    else:
        count_table.to_csv(fileout, index=False)
    

def usage():
    print 'python extractIgFeatures.py [--output outputfile]'

if __name__ == "__main__":
    main(sys.argv[1:])
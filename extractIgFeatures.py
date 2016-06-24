import pandas as pd
from natsort import natsorted
import sys
import getopt
import re
import numpy as np

def main(argv):
    germline = pd.read_csv("/home/wuji/tools/imgt/germline_kabat.csv")
    germline.columns = ['gene', 'position', 'label']
    append=False
    try:
        opts, args = getopt.getopt(argv,"ha", ["output=", "coverage="])
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
    col_names = ['id', 'germline', 'query_pos', 'ref_pos', 'query_base', 'ref_base', 'query_aa', 'ref_aa']
    col_names.extend(['query_b2', 'query_b1', 'query_a1', 'query_a2', 'ref_b2', 'ref_b1', 'ref_a1', 'ref_a2', 'phred'])
    df.columns = col_names
    
    # 1. Remove mutation with "N"
    n1 = df.shape[0]
    df = df[df['query_base']!="N"]
    df = df[df['phred']>=10]
    n2 = df.shape[0]
    print "removing " + str(n1 - n2) + " rows"
    df['germline'] = df['germline'].map(str.strip)
    
    # 2. Merge germline dataframe with mutation dataframe
    df = pd.merge(df, germline, left_on=['germline', 'ref_pos'], right_on = ['gene','position'] )
    df = df.sort_values(by=['id','query_pos'])
    
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
    # Fill in missing columns
    label_set = set(germline['label'])
    label_set = filter(lambda v: v==v, label_set)
    full_cols = natsorted(label_set)
    
    miss_cols = set(full_cols) - set(cols)
    if len(miss_cols) > 0:
        index = count_table.index
        miss_table = pd.DataFrame(index=index, columns=miss_cols)
        count_table = pd.concat([count_table, miss_table], axis=1, join='inner')
    
    count_table = count_table[full_cols]
    count_table = count_table.fillna(0)
    count_table['id'] = count_table.index
    
    # 5. Join kabat mutation table with sequence information
    df1 = df[['id', 'germline']].drop_duplicates()
    df1 = df1[df1['germline'].str.startswith("IGHV")]
    df1.columns = ['id', 'vgene']
    count_table = pd.merge(count_table, df1, on=['id'])
    df1 = df[['id', 'germline']].drop_duplicates()
    df1 = df1[df1['germline'].str.startswith("IGHJ")]
    df1.columns = ['id', 'jgene']
    count_table = pd.merge(count_table, df1, on=['id'])
    gene_pattern = '(^IGHV\d{1,2}-\d{1,3}|^IGHJ\d{1,2})'
    
    def extractGeneName(x):
        m = re.search(gene_pattern, x)
        found = 'NA'        
        if m:
            found=m.group()
        return found
    
    count_table['vgene'] = count_table['vgene'].map(lambda g : extractGeneName(g))
    count_table['jgene'] = count_table['jgene'].map(lambda g : extractGeneName(g))
    
    # 6. Count number of sense and nonsense mutation
    #     0 - nonsense
    #     1 - sense
    df2 = df[['id','query_aa','ref_aa','label']]
    df2['sense_mut'] = np.where(df2['query_aa']==df2['ref_aa'], 0, 1)
    mut_count = df2.groupby(['id'])['sense_mut'].sum()
    mut_count = mut_count.reset_index()
    temp = df2['id'].value_counts()
    temp = temp.reset_index()    
    temp.columns = ['id','total']
    mut_count = pd.merge(mut_count, temp)
    mut_count['missense_mut'] = mut_count.total - mut_count.sense_mut
    mut_count = mut_count.drop('total', 1)
    count_table = pd.merge(count_table, mut_count)
    
    # 6. Output data into csv file
    if append:
        count_table.to_csv(fileout, mode='a', index=False, header=True)
    else:
        count_table.to_csv(fileout, index=False)
    

def usage():
    print 'python extractIgFeatures.py [-a] [--output outputfile]'
    print '-a: append to the output file'

if __name__ == "__main__":
    main(sys.argv[1:])
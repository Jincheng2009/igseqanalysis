import sys
import getopt
import pandas as pd
from Bio import SeqIO

def main(argv):
    ucfile = None
    fafile = None
    outfile = None
    try:
        opts, args = getopt.getopt(argv,"hc:f:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            usage()
            sys.exit()
        elif opt == "-c":
            ucfile = arg
        elif opt == "-f":
            fafile = arg
        elif opt == "-o":
            outfile = arg

    if ucfile is None or fafile is None:
        sys.stderr.write("missing input file \n")

    with open(fafile) as fasta_file:  # Will close handle cleanly
        identifiers = []
        lengths = []
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
            identifiers.append(seq_record.id)
            lengths.append(''.join(seq_record.seq))
            
    #converting lists to pandas Series    
    s1 = pd.Series(identifiers, name='fastaid')
    s2 = pd.Series(lengths, name='seq')
    #Gathering Series into a pandas DataFrame and rename index as ID column
    seq_df = pd.DataFrame(dict(fastaid=s1, seq=s2))

    # Read in clustering result
    cluster = pd.read_table(ucfile, sep="\t", header=None)
    cluster = cluster[cluster[0].isin(['S','H'])]

    df = pd.merge(cluster, seq_df, left_on = 8, right_on = 'fastaid')
    df = pd.merge(df, seq_df, left_on = 9, right_on = 'fastaid', how='left')

    df.loc[pd.isnull(df['seq_y']), 'seq_y'] = df['seq_x']

    df = df[['seq_x', 'seq_y']]
    df = df.groupby(['seq_x','seq_y']).size()
    df = df.reset_index()
    df.columns = ['CDR3','CDR3c','count']

    if outfile is None:
        output = sys.stdout
    else:
        output = outfile
    df.sort_values(['CDR3c', 'count'], ascending=[True, False], inplace=True)
    df.to_csv(output, sep="\t", index=False, header=False, na_rep='.')

def usage():
    print 'python formatCluster.py -c usearch_cluster.uc -f input.fasta > out.tsv'    

if __name__ == "__main__":
    main(sys.argv[1:])


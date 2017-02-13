# Immunoglobulin Sequence Analysis (igseqanalysis)

Igseqanalysis is a package developed in Python that parse and process the NCBI-IgBlast alignment results of antibody sequences (NGS) to extract the CDR regions or mutations. 

## Prerequisite
* Python 2.7 [Installation instruction](https://www.python.org/download/releases/2.7/)
* Pandas [Installation instruction](http://pandas.pydata.org/)
* Biopython [Download site](http://biopython.org/wiki/Download)
* Standalone NCBI-IgBlast 1.5.0 [Installation FTP site](ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/)
* usearch [Installation site](http://www.drive5.com/usearch/) (Optional: needed if you want to cluster CDR sequence by identity)

*Note: For convenience, add `igblastn` and `usearch` to environment PATH*

## Workflow

Change directory to the igseqanalysis folder and start the workflow.

### 1. Annotate sequence by IgBlast

IMGT folder contains the VDJ database built by blast with the fasta files. If you want to use database other than the one provided in the package, you could build a new one following the instructions in NCBI-Igblast and specify the location using `-germline_db_V`, `-germline_db_D` and `-germline_db_J` options.

#### Annotate Sequences in Read 1

```
cat sample/sample_R1.fasta | igblastn -germline_db_V imgt/human_V \
                    -germline_db_D imgt/human_D \
                    -germline_db_J imgt/human_J \
                    -domain_system kabat -organism human \
                    -auxiliary_data optional_file/human_gl.aux \
                    -show_translation \
                    -num_alignments_V 1 \
                    -num_alignments_J 1 \
                    -num_alignments_D 1 \
| python parseIgBlast.py -t CDR > sample/sample_R1.csv
```    

#### Annotate Sequences in Read 2

```  
cat sample/sample_R2.fasta | igblastn -germline_db_V imgt/human_V \
                    -germline_db_D imgt/human_D \
                    -germline_db_J imgt/human_J \
                    -domain_system kabat -organism human \
                    -auxiliary_data optional_file/human_gl.aux \
                    -show_translation \
                    -num_alignments_V 1 \
                    -num_alignments_J 1 \
                    -num_alignments_D 1 \
| python parseIgBlast.py -t CDR > sample/sample_R2.csv
```

### 2. Pair VH with VL

Pair the VH annotation with the VL annotation by fasta ID. 
Each row in the output file will contain the following columns:
  1.  fasta id
  2.  VH V gene
  3.  VH J gene
  4.  CDR-H1 sequence
  5.  CDR-H2 sequence
  6.  CDR-H3 sequence
  7.  VL V gene
  8.  VL J gene
  9.  CDR-L1 sequence
  10. CDR-L2 sequence
  11. CDR-L3 sequence
```
python pairByID.py -l sample/sample_R1.csv -r sample/sample_R2.csv -o sample/sample.paired.tsv
```

### 3. Translate DNA into protein for CDR sequences

Output file has the same format as the input format. Only the DNA sequences has been translated into protein sequences.

	cat sample/sample.paired.tsv | python translateTable.py -p 3,4,5,8,9,10 > sample/sample.paired.prot.tsv

### 4. Count the unique paired CDR3

If you want to count the unique CDR3 in DNA sequences, you could provide the DNA sequences as the input `sample.paired.tsv` 

	cat sample/sample.paired.prot.tsv | python countUnique.py -p 5,10 > sample/sample.paired.CDR3.count

### 5. Convert CDR3 in CSV format to fasta format for usearch clustering

If you want to cluster the unique CDR3 in DNA sequences, you could provide the DNA sequences as the input `sample.paired.tsv` 

	cat sample/sample.paired.prot.tsv | python csv2Fasta.py -p 5  > sample/sample.VH.fasta
	cat sample/sample.paired.prot.tsv | python csv2Fasta.py -p 10 > sample/sample.VL.fasta

### 6. Clustering CDR3 by usearch

Clustering could efficiently reduce the effect of PCR and sequencing errors, but at expense of cluster a real unique VH/VL into another VH/VL. `-id 0.88` for protein sequences allows 1 amino acid difference when CDR length is between 9 and 16, and 2 amino acid difference when CDR length is between 17 and 24. If clustering DNA sequences, `-id 0.96` is similar to `-id 0.88` for protein sequences. `-sort size` will enable the most abundant sequence is considered as the centroid sequence of the cluster. `-fulldp -maxgaps 0` disallow any gaps.

	usearch -cluster_fast sample/sample.VH.fasta -id 0.88 -sort size -uc sample/sample.VH.uc -fulldp -maxgaps 0
	usearch -cluster_fast sample/sample.VL.fasta -id 0.88 -sort size -uc sample/sample.VL.uc -fulldp -maxgaps 0

### 7. Format the usearch result into tabular format

	python formatCluster.py -c sample/sample.VH.uc -f sample/sample.VH.fasta > sample/sample.R1.cluster.count
	python formatCluster.py -c sample/sample.VL.uc -f sample/sample.VL.fasta > sample/sample.R2.cluster.count

The output file contains 3 columns:

1. Unique CDR3 sequence
2. Centroid CDR3 sequences that this sequence belongs to
3. Count of the unique CDR3 sequence

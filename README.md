# Immunoglobulin Sequence Analysis (igseqanalysis)

Igseqanalysis is a package developed in Python that parse and process the NCBI-IgBlast alignment results of antibody sequences (NGS) to extract the CDR regions or mutations. 

## Prerequisite
* Python 2.7 [Installation instruction](https://www.python.org/download/releases/2.7/)
* Pandas [Installation instruction](http://pandas.pydata.org/)
* Biopython [Download site](http://biopython.org/wiki/Download)
* Standalone NCBI-IgBlast 1.5 [Installation FTP site](ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/)
* usearch [Installation site](http://www.drive5.com/usearch/) (Optional: needed if you want to cluster CDR sequence by identity)

## Workflow
1. IgBlast
```
cat sample/sample_R1.fasta | igblastn -germline_db_V ~/tools/imgt/human_V \
                           			  -germline_db_D ~/tools/imgt/human_D \
               						  -germline_db_J ~/tools/imgt/human_J \
               						  -domain_system kabat \
               						  -organism human -auxiliary_data ~/tools/ncbi-igblast-1.5.0/optional_file/human_gl.aux \
               						  -show_translation -num_alignments_V 1 \
               						  -num_alignments_J 1 \
               						  -num_alignments_D 1 \
| python igseqanalysis/parseIgBlast.py -t CDR > sample/sample_R1.csv
```
```
cat sample/sample_R2.fasta | igblastn -germline_db_V ~/tools/imgt/human_V \
                           			  -germline_db_D ~/tools/imgt/human_D \
               						  -germline_db_J ~/tools/imgt/human_J \
               						  -domain_system kabat \
               						  -organism human -auxiliary_data ~/tools/ncbi-igblast-1.5.0/optional_file/human_gl.aux \
               						  -show_translation -num_alignments_V 1 \
               						  -num_alignments_J 1 \
               						  -num_alignments_D 1 \
| python igseqanalysis/parseIgBlast.py -t CDR > sample/sample_R2.csv
```
2. Pair VH with VL
```
python igseqanalysis/pairByID.py -l sample/sample_R1.csv \
                                   -r sample/sample_R2.csv \
                                   -o sample/sample.paired.tsv
```
3. Translate DNA into protein for CDR sequences
```
cat sample/sample.paired.tsv | python ../igseqanalysis/translateTable.py -p 3,4,5,8,9,10 > sample/sample.paired.prot.tsv
```
4. Count the unique CDR3
```
cat sample/sample.paired.prot.tsv | python ../igseqanalysis/countUnique.py -p 5,10 > sample/sample.paired.CDR3.count
```
5. Convert CDR3 in CSV format to fasta format for usearch clustering
```
cat sample/sample.paired.prot.tsv | python ../igseqanalysis/csv2Fasta.py -p 5  > sample/sample.VH.fasta
cat sample/sample.paired.prot.tsv | python ../igseqanalysis/csv2Fasta.py -p 10 > sample/sample.VL.fasta
```
6. Clustering CDR3 by usearch
```
usearch -cluster_fast sample/sample.VH.fasta -id 0.88 -sort size -uc sample/sample.VH.uc -fulldp -maxgaps 0
usearch -cluster_fast sample/sample.VL.fasta -id 0.88 -sort size -uc sample/sample.VL.uc -fulldp -maxgaps 0
```
7. Format the usearch result into tabular format
```
python igseqanalysis/formatCluster.py -c sample/sample.VH.uc -f sample/sample.R2.fasta > sample/sample.R1.cluster.count
python igseqanalysis/formatCluster.py -c sample/sample.VL.uc -f sample/sample.R1.fasta > sample/sample.R1.cluster.count
```
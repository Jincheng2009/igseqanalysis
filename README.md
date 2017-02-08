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
cat sample/sample_R1.fasta | igblastn -germline_db_V ~/tools/imgt/human_V \
                           			  -germline_db_D ~/tools/imgt/human_D \
               						  -germline_db_J ~/tools/imgt/human_J \
               						  -domain_system kabat \
               						  -organism human -auxiliary_data ~/tools/ncbi-igblast-1.5.0/optional_file/human_gl.aux \
               						  -show_translation -num_alignments_V 1 \
               						  -num_alignments_J 1 \
               						  -num_alignments_D 1 \
| python igseqanalysis/parseIgBlast.py -t CDR > annotation/sample_R1.csv

cat sample/sample_R2.fasta | igblastn -germline_db_V ~/tools/imgt/human_V \
                           			  -germline_db_D ~/tools/imgt/human_D \
               						  -germline_db_J ~/tools/imgt/human_J \
               						  -domain_system kabat \
               						  -organism human -auxiliary_data ~/tools/ncbi-igblast-1.5.0/optional_file/human_gl.aux \
               						  -show_translation -num_alignments_V 1 \
               						  -num_alignments_J 1 \
               						  -num_alignments_D 1 \
| python igseqanalysis/parseIgBlast.py -t CDR > annotation/sample_R2.csv

2. Pair VH with VL
python ../repseqanalysis/pairByID.py -l ../annotation/638_em_S3_L001_R1_001.filter_VH.tsv \
                                   -r ../annotation/638_em_S3_L001_R2_001.filter_VL.tsv \
                                   -o ../annotation/S3.paired.tsv

3. Translate DNA into protein for CDR sequences
cat $f | python ../repseqanalysis/translateTable.py -p 3,4,5,8,9,10 > "../annotation/"$name".tsv"

4. Count the unique CDR3
cat $f | python ../repseqanalysis/countUnique.py -p 5,10 > "../CDR/"$name".CDR3.count"

5. Convert CDR3 in CSV format to fasta format for usearch clustering
cat $f | python ../repseqanalysis/csv2Fasta.py -p 5 > "../CDR/"$name".VH.prot.fasta"
cat $f | python ../repseqanalysis/csv2Fasta.py -p 10 > "../CDR/"$name".VL.prot.fasta"

6. Clustering CDR3 by usearch
usearch -cluster_fast $f -id 0.88 -sort size -uc "../cluster/"$name".uc" -fulldp -maxgaps 0

7. Format the usearch result into tabular format
python ../repseqanalysis/formatCluster.py -c $f -f "../CDR/"$name".fasta" > "../unpaired/"$name".cluster.count"
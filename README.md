# Immunoglobulin Sequence Analysis (igseqanalysis)

Igseqanalysis is a package developed in Python that parse and process the NCBI-IgBlast alignment results of antibody sequences (NGS) to extract the CDR regions or mutations. 

## Prerequisite
* Python 2.7 [Installation instruction](https://www.python.org/download/releases/2.7/)
* Pandas [Installation instruction](http://pandas.pydata.org/)
* Biopython [Download site](http://biopython.org/wiki/Download)
* Standalone NCBI-IgBlast 1.5 [Installation FTP site](ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/)
* usearch [Installation site](http://www.drive5.com/usearch/) (Optional: needed if you want to cluster CDR sequence by identity)

## Workflow
**1. Annotate sequence by IgBlast**

IMGT folder contains the VDJ database built by blast with the fasta files. If you want to use database other than the one provided in the package, you could build a new one following the instructions in NCBI-Igblast and specify the location using -germline_db_V, -germline_db_D and -germline_db_J options.
	
	cat sample/sample_R1.fasta | igblastn -germline_db_V imgt/human_V \
	                           			  -germline_db_D imgt/human_D \
	               						  -germline_db_J imgt/human_J \
	               						  -domain_system kabat \
	               						  -organism human -auxiliary_data optional_file/human_gl.aux \
	               						  -show_translation -num_alignments_V 1 \
	               						  -num_alignments_J 1 \
	               						  -num_alignments_D 1 \
	| python igseqanalysis/parseIgBlast.py -t CDR > sample/sample_R1.csv

	cat sample/sample_R2.fasta | igblastn -germline_db_V imgt/human_V \
	                           			  -germline_db_D imgt/human_D \
	               						  -germline_db_J imgt/human_J \
	               						  -domain_system kabat \
	               						  -organism human -auxiliary_data optional_file/human_gl.aux \
	               						  -show_translation -num_alignments_V 1 \
	               						  -num_alignments_J 1 \
	               						  -num_alignments_D 1 \
	| python igseqanalysis/parseIgBlast.py -t CDR > sample/sample_R2.csv

**2. Pair VH with VL**

Pair the VH annotation with the VL annotation by fasta ID. Each row will contain the following columns:
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
python igseqanalysis/pairByID.py -l sample/sample_R1.csv \
	                             -r sample/sample_R2.csv \
	                             -o sample/sample.paired.tsv
```

**3. Translate DNA into protein for CDR sequences**

	cat sample/sample.paired.tsv | python igseqanalysis/translateTable.py -p 3,4,5,8,9,10 > sample/sample.paired.prot.tsv

**4. Count the unique paired CDR3**

	cat sample/sample.paired.prot.tsv | python igseqanalysis/countUnique.py -p 5,10 > sample/sample.paired.CDR3.count

**5. Convert CDR3 in CSV format to fasta format for usearch clustering**

	cat sample/sample.paired.prot.tsv | python igseqanalysis/csv2Fasta.py -p 5  > sample/sample.VH.fasta
	cat sample/sample.paired.prot.tsv | python igseqanalysis/csv2Fasta.py -p 10 > sample/sample.VL.fasta

**6. Clustering CDR3 by usearch**

	usearch -cluster_fast sample/sample.VH.fasta -id 0.88 -sort size -uc sample/sample.VH.uc -fulldp -maxgaps 0
	usearch -cluster_fast sample/sample.VL.fasta -id 0.88 -sort size -uc sample/sample.VL.uc -fulldp -maxgaps 0

**7. Format the usearch result into tabular format**

	python igseqanalysis/formatCluster.py -c sample/sample.VH.uc -f sample/sample.R2.fasta > sample/sample.R1.cluster.count
	python igseqanalysis/formatCluster.py -c sample/sample.VL.uc -f sample/sample.R1.fasta > sample/sample.R1.cluster.count
	
import sys
import re
import sys, getopt
import operator

centroid = ""

# example use this script: cat CDR/Hs638-48op-K-L_S6_L001_R1.CDR.L3.prot-0.9.uc | python ucFormat.py > CDR/Hs638-48op_S6_L001_R1.CDR.L3.0.9.cluster.out
# cat CDR/Hs638-48op-K-L_S6_L001_R1.CDR.L3.prot-0.9.uc | python ucFormat.py > CDR/Hs638-48op-K-L_S6_L001_R1.CDR.L3.prot-0.9.cluster.out
# cat CDR/Hs717-0op-K-L_S3_L001_R1.CDR.L3.prot-0.9.uc | python ucFormat.py > CDR/Hs717-0op-K-L_S3_L001_R1.CDR.L3.prot-0.9.cluster.out
# cat CDR/Hs717-24op-K-L_S4_L001_R1.CDR.L3.prot-0.9.uc | python ucFormat.py > CDR/Hs717-24op-K-L_S4_L001_R1.CDR.L3.prot-0.9.cluster.out
# cat CDR/Hs717-48op-K-L_S5_L001_R1.CDR.L3.prot-0.9.uc | python ucFormat.py > CDR/Hs717-48op-K-L_S5_L001_R1.CDR.L3.prot-0.9.cluster.out
# cat CDR/Hs717-24em-K-L_S1_L001_R1.CDR.L3.prot-0.9.uc | python ucFormat.py > CDR/Hs717-24em-K-L_S1_L001_R1.CDR.L3.prot-0.9.cluster.out
# cat CDR/Hs717-48em-K-L_S2_L001_R1.CDR.L3.prot-0.9.uc | python ucFormat.py > CDR/Hs717-48em-K-L_S2_L001_R1.CDR.L3.prot-0.9.cluster.out

# cat CDR/Hs638-48op-K-L_S6_L001_R2.CDR.H3.prot-0.85.uc | python ucFormat.py > CDR/Hs638-48op-K-L_S6_L001_R2.CDR.H3.prot-0.85.cluster.out
# cat CDR/Hs717-0op-K-L_S3_L001_R2.CDR.H3.prot-0.85.uc | python ucFormat.py > CDR/Hs717-0op-K-L_S3_L001_R2.CDR.H3.prot-0.85.cluster.out
# cat CDR/Hs717-24op-K-L_S4_L001_R2.CDR.H3.prot-0.85.uc | python ucFormat.py > CDR/Hs717-24op-K-L_S4_L001_R2.CDR.H3.prot-0.85.cluster.out
# cat CDR/Hs717-48op-K-L_S5_L001_R2.CDR.H3.prot-0.85.uc | python ucFormat.py > CDR/Hs717-48op-K-L_S5_L001_R2.CDR.H3.prot-0.85.cluster.out
# cat CDR/Hs717-24em-K-L_S1_L001_R2.CDR.H3.prot-0.85.uc | python ucFormat.py > CDR/Hs717-24em-K-L_S1_L001_R2.CDR.H3.prot-0.85.cluster.out
# cat CDR/Hs717-48em-K-L_S2_L001_R2.CDR.H3.prot-0.85.uc | python ucFormat.py > CDR/Hs717-48em-K-L_S2_L001_R2.CDR.H3.prot-0.85.cluster.out

# cat cluster/Hs638-48op-K-L_S6_L001.CDR3.pair.prot.uc | python ucFormat.py > CDR/Hs638-48op-K-L_S6_L001.CDR3.pair.prot-0.92.cluster.out
# cat cluster/Hs717-0op-K-L_S3_L001.CDR3.pair.prot.uc | python ucFormat.py > CDR/Hs717-0op-K-L_S3_L001.CDR3.pair.prot-0.92.cluster.out
# cat cluster/Hs717-24op-K-L_S4_L001.CDR3.pair.prot.uc | python ucFormat.py > CDR/Hs717-24op-K-L_S4_L001.CDR3.pair.prot-0.92.cluster.out
# cat cluster/Hs717-48op-K-L_S5_L001.CDR3.pair.prot.uc | python ucFormat.py > CDR/Hs717-48op-K-L_S5_L001.CDR3.pair.prot-0.92.cluster.out
# cat cluster/Hs717-24em-K-L_S1_L001.CDR3.pair.prot.uc | python ucFormat.py > CDR/Hs717-24em-K-L_S1_L001.CDR3.pair.prot-0.92.cluster.out
# cat cluster/Hs717-48em-K-L_S2_L001.CDR3.pair.prot.uc | python ucFormat.py > CDR/Hs717-48em-K-L_S2_L001.CDR3.pair.prot-0.92.cluster.out

def printCluster(cluster, centroid):
	if len(cluster)<1:
		return

	sorted_cluster = sorted(cluster.items(), key=operator.itemgetter(1), reverse=True)
	for seq,number in sorted_cluster:
		print seq + "\t" + centroid + "\t" + str(number)
	
def main(argv):
	count=False
	try:
		opts, args = getopt.getopt(argv,"hc", ["count"])
	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt == "-h":
			usage()
			sys.exit()
		elif opt in ("-c", "--count"):
			count=True

	member={}
	centroid=""

	for line in sys.stdin:
		seq = re.split(r'\t+', line.rstrip())
		
		if seq[0] == 'S':
			# Print the previous cluster
			printCluster(member, centroid)
			# Initialize the new cluster
			centroid = seq[3]
			member={}
			member[centroid] = 1

		elif seq[0] == 'H':
			if seq[3] in member:
				member[seq[3]] = member[seq[3]] + 1
			else:
				member[seq[3]] = 1

	# Print the last cluster
	printCluster(member, centroid)


def usage():
	print 'ucFormat.py [-count]'

if __name__ == "__main__":
	main(sys.argv[1:])
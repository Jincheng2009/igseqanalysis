import os
import re
import sys, getopt
import operator

centroid = ""
# Point this path to the folder which contain the germline gene fasta files
# Used to extract all names of germline genes such as IGHV18-01, IGHV18-03
dbpath = os.path.join(os.path.dirname(os.path.realpath(__file__)),"germlinedb")

def printGermlineCount(germline):
	if len(germline)<1:
		return

	sorted_germlines = sorted(germline.items(), key=operator.itemgetter(0))
	for germline, number in sorted_germlines:
		print germline + "\t" + str(number)

def readGermlineName(filein):
	germline=[]
	for line in filein:
		if not line.startswith(">"):
			continue
		line=line.rstrip()
		germline.append(line[1:])
	return germline
	
def main(argv):
	pair=False
	try:
		opts, args = getopt.getopt(argv,"hp", ["pair","igh","igk","igl"])
	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt == "-h":
			usage()
			sys.exit()
		elif opt in ("-p", "--pair"):
			pair=True
		elif opt == "--igh":
			printIgh=True
		elif opt == "--igk":
			printIgk=True
		elif opt == "--igl":
			printIgl=True

	ighCount={}
	iglCount={}
	igkCount={}

	# Read IGH germline names
	f=open(os.path.join(dbpath,'HOMO SAPIENS_HEAVY_V-REGION.ft'),'r')
	ighnames = readGermlineName(f)
	f=open(os.path.join(dbpath,'HOMO SAPIENS_HEAVY_J-REGION.ft'),'r')
	ighnames = ighnames + readGermlineName(f)

	# Read IGK germline names
	f=open(os.path.join(dbpath,'HOMO SAPIENS_KAPPA_V-REGION.ft'),'r')
	igknames = readGermlineName(f)
	f=open(os.path.join(dbpath,'HOMO SAPIENS_KAPPA_J-REGION.ft'),'r')
	igknames = igknames + readGermlineName(f)
	
	# Read IGL germline names
	f=open(os.path.join(dbpath,'HOMO SAPIENS_LAMBDA_V-REGION.ft'),'r')
	iglnames = readGermlineName(f)
	f=open(os.path.join(dbpath,'HOMO SAPIENS_LAMBDA_J-REGION.ft'),'r')
	iglnames = iglnames + readGermlineName(f)

	for name in ighnames:
		ighCount[name]=0
	for name in igknames:
		igkCount[name]=0
	for name in iglnames:
		iglCount[name]=0

	for line in sys.stdin:
		if not line.startswith(">"):
			continue

		line=line.rstrip()
		seq = line.split("|");
		seqid=seq[0]
		vh=""
		jh=""
		vl=""
		jl=""

		info=seq[len(seq)-1]
		germline = re.split(r"\++", info)
		germline.remove('')

		# If paired fasta sequence, have both VH and VL germline. Append them together.
		if pair:
			germline1=re.split(r"\++", seq[len(seq)-2])
			germline1.remove('')
			germline = germline1+germline

		for item in germline:
			if item.startswith("IGHV"):
				vh = item
			elif item.startswith("IGHJ"):
				jh = item
			elif item.startswith("IGLV") or item.startswith("IGKV"):
				vl = item
			elif item.startswith("IGLJ") or item.startswith("IGKJ"):
				jl = item
		
		print seqid + "\t" + vh + "\t" + jh + "\t" + vl + "\t" + jl
				



def usage():
	print 'Convert the annotated fasta file into tabular format with germline information in separate columns'
	print 'cat CDR/Hs638-48op-K-L_S6_L001.CDR3.pair.prot.fasta | python formatGermline.py > germline.tsv'

if __name__ == "__main__":
	main(sys.argv[1:])
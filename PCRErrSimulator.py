# Simulator for PCR process considering the error rate of polymerase
# PCR error happans randomly and error rate is measured
# on probability per nucleotide per cycle
# Author: Jincheng Wu

import sys, getopt
import random
import math
import numpy

nameList=[]
seqDict={}
pNumOfErr = []

inputfile = ''
outputfile = ''
erate = 0
cycle = 0

def usage():
	print 'PCRErrSimulator.py -i <inputfile.fasta> -o <outputfile.count> -e <error rate> -c <cycle number> [-s <seed>] [-d]'
	print '-i\t input file, fasta format, template must have same length'
	print '-o\t output file, count format'
	print '-e\t error rate of PCR, per nucleotide per cycle'
	print '-c\t number of PCR cycle'
	print '-s\t seed for random number generator'
	print '-d\t show number of nucleotide difference to the most count sequence'

def main(argv):
	showDiff=False
	inputfile=None
	outputfile=None
	try:
		opts, args = getopt.getopt(argv,"hdi:o:e:c:s:",["ifile=","ofile="])
	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			usage()
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-o", "--ofile"):
			outputfile = arg
		elif opt in ("-e"):
			erate = float(arg)
		elif opt in ("-c"):
			cycle = int(arg)
		elif opt in ("-s"):
			seed = int(arg)
			random.seed(seed)
		elif opt in ("-d"):
			showDiff=True	
	sys.stderr.write('Error rate is ' + str(erate) + '\n')
	sys.stderr.write('PCR cycle is '+ str(cycle) + '\n')

	if inputfile is not None:
		fp = open(inputfile)
	else:
		fp = sys.stdin
		
	for line in fp:
		line = line.rstrip()
		if line.startswith(">"):
			nameList.append(line)
		else:
			# If sequence exists in the dict, just add the count by 1
			if line in seqDict:
				seqDict[line] += 1
			# If sequence does not exist in the dict, create a new entry with value = 1
			else:
				seqDict[line] = 1
				n = len(line)

	############################################################
	## Calculate cumulative probability function of how many 
	## errors occur
	for i in range(n):
		f = math.factorial
		pNumOfErr.append( f(n)/f(i)/f(n-i) * erate**i * (1-erate)**(n-i) )

	pAllErr=sum(pNumOfErr)
	for i in range(n):
		if i ==0:
			pNumOfErr[i]=pNumOfErr[i]/pAllErr
		else:
			pNumOfErr[i]=pNumOfErr[i]/pAllErr + pNumOfErr[i-1]
	#############################################################

	############################################################
	## PCR starts	
	for ct in range(0,cycle):
		newDict = {}
		sys.stderr.write("PCR Cycle: " + str(ct) + '\n')
		#for seq in seqDict:
		newDict = replicateSeq(seqDict, erate)
		for newSeq in newDict:
			if newSeq in seqDict:
				seqDict[newSeq] += newDict[newSeq]
			# If sequence does not exist in the dict, create a new entry with value = 1
			else:
				seqDict[newSeq] = newDict[newSeq]
		if sum(seqDict.itervalues()) > 1e8:
			# If size is larger than 100 million, down-sample half
			for seq in seqDict:
				count = seqDict[seq]
				count = numpy.random.binomial(count, 0.5)
				if count == 0:
					seqDict.pop(seq, None)
				else:
					seqDict[seq]=count						

	############################################################
	## Print out the results in descending count order
	if outputfile is not None:
		out = open(outputfile,'w')
	else:
		out = sys.stdout
		
	first=True

	for seq in sorted(seqDict, key=seqDict.get, reverse=True):
		if showDiff:
			if first:
				mostSeq=seq
				first=False
			diff = countDiff(mostSeq,seq)
			out.write(seq+'\t'+str(seqDict.get(seq))+'\t'+str(diff)+'\n')
		else:
			out.write(seq+'\t'+str(seqDict.get(seq))+'\n')

	if outputfile is not None:
		out.close()
	sys.stderr.write('Simulation completes')

def replicateSeq(seqDict, erate):
	newSeqDict = {}
	base = ["A","T","C","G"]

	for seq in seqDict:
		count = seqDict.get(seq)
		good_rep = 0
		n = len(seq)
		# Probability that this replication will have error
		p_Err = 1 - (1 - erate) ** n
		for i in range(0,count):
			r1 = random.random()
			# print r
			# Error will happen when random number is less than p_Err
			if r1>p_Err:
				good_rep = good_rep + 1

		seqDict[seq] = seqDict[seq] + good_rep

		
		for i in range(0,count - good_rep):
			nErr = calcNumOfErr()
			for j in random.sample(range(n),nErr):
				errBase = seq[j]
				newBase = random.sample(base,1)[0]
				while(newBase==errBase):
					newBase = random.sample(base,1)[0]
				newSeq = seq[:j] + newBase + seq[j+1:]
				if newSeq in newSeqDict:
					newSeqDict[newSeq] += 1
				# If sequence does not exist in the dict, create a new entry with value = 1
				else:
					newSeqDict[newSeq] = 1

	return newSeqDict

def calcNumOfErr():
	r2=random.random()
	for i in range(len(pNumOfErr)): 
		if(r2<pNumOfErr[i]):
			return i+1

def countDiff(seq1,seq2):
	count = 0
	for i in range(len(seq1)): 
		if seq1[i]!=seq2[i]:
			count += 1
	return count


if __name__ == "__main__":
	main(sys.argv[1:])

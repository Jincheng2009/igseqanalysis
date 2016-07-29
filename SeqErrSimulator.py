import sys, getopt
import random
import math
import numpy as np
# Simulator for sequencing process considering the error rate of sequencing
# sequencing error happans randomly on any nucleotide being sequenced
# Author: Jincheng Wu

seqDict={}
pNumOfErr = []
inputfile = ''
outputfile = ''
erate = 0

def usage():
	
	print 'SeqErrSimulator.py -i <inputfile.count> -o <outputfile.count> -e <error rate> [-s <seed>] [-d]'
	print '-i\t input file, count format'
	print '-o\t output file, count format'
	print '-e\t error rate of sequencing, based on per nucleotide'
	print '-s\t seed for random number generator'
	print '-d\t show number of nucleotide difference to the highest-count sequence'

def main(argv):
	showDiff=False
	biasedSub=False
	n=0
	try:
		opts, args = getopt.getopt(argv,"hdxi:o:e:s:",["ifile=","ofile="])
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
		elif opt in ("-s"):
			seed = int(arg)
			random.seed(seed)
		elif opt in ("-d"):
			showDiff=True
		elif opt in ("-x"):
			biasedSub=True	
	print 'Input file is ', inputfile
	print 'Output file is ', outputfile
	print 'Error rate is ', erate

	with open(inputfile) as fp:
		for line in fp:
			line = line.rstrip()
			record = line.split("\t")
			# If sequence exists in the dict, just add the count by its count
			if record[0] in seqDict:
				seqDict[record[0]] += int(record[1])
			# If sequence does not exist in the dict, create a new entry with value = count
			else:
				seqDict[record[0]] = int(record[1])
				n=len(record[0])

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
	## Sequencing error simulation
	newDict = {}
	#for seq in seqDict:
	newDict = readSeq(seqDict, erate, biasedSub)
	for newSeq in newDict:
		if newSeq in seqDict:
			seqDict[newSeq] += newDict[newSeq]
		# If sequence does not exist in the dict, create a new entry with value = 1
		else:
			seqDict[newSeq] = newDict[newSeq]				

	############################################################
	## Print out the results in descending count order
	out = open(outputfile,'w')
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

	out.close()
	print 'Simulation completes'

def readSeq(seqDict, erate, biasedSub):
	newSeqDict = {}
	

	for seq in seqDict:
		count = seqDict.get(seq)
		good_read = 0
		n = len(seq)
		# Probability that this replication will not have error
		p_Err = 1 - (1 - erate) ** n
		for i in range(0,count):
			r1 = random.random()
			# print r
			# Error will happen when random number is less than p_Err
			if r1>p_Err:
				good_read = good_read + 1

		seqDict[seq] = good_read
	
		for i in range(0,count - good_read):
			nErr = calcNumOfErr()
			for j in random.sample(range(n),nErr):
				base = seq[j]
				if biasedSub:
					newBase = getErrBaseBiased(base)
				else:
					newBase = getErrBase(base)
				newSeq = seq[:j] + newBase + seq[j+1:]
				if newSeq in newSeqDict:
					newSeqDict[newSeq] += 1
				# If sequence does not exist in the dict, create a new entry with value = 1
				else:
					newSeqDict[newSeq] = 1

	return newSeqDict

def getErrBase(base):
	allBase = ["A","T","C","G"]
	allBase.remove(base)
	return random.sample(allBase,1)[0]

def getErrBaseBiased(base):
	r3 = random.random()
	allBase = ["A","T","C","G"]
	# matrix with ATCG*ATCG to describe substitution biases
	pMatrix = np.matrix('0 0.1 0.77 0.13; 0.1 0 0.4 0.5; 0.57 0.31 0 0.12; 0.31 0.54 0.15 0')
	if base=="A":
		pArray=np.asarray(pMatrix[0,])
	elif base=="T":
		pArray=np.asarray(pMatrix[1,])
	elif base=="C":
		pArray=np.asarray(pMatrix[2,])
	else:
		pArray=np.asarray(pMatrix[3,])

	cp = 0
	for i in range(len(pArray[0])):
		cp += pArray[0][i]
		if r3 < cp:
			return allBase[i]

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

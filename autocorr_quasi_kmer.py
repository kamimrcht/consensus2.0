import scipy
from scipy import signal
import matplotlib
from matplotlib import pyplot
from matplotlib.pyplot import plot, draw, show, xcorr
import sys
import os
import shlex, subprocess
from subprocess import Popen, PIPE, STDOUT
import numpy as np
import peakutils
from peakutils.plot import plot as pplot
import numpy.fft
from scipy import fftpack
import compute_consensus_from_period
#~ sequence = "ATCAGAAGCTTTGTCGGGAGAGGGTTCCCAGAGTCAGAGGGTTTCTGGGCCGGAGATGAAACCATCGGTTGGGAAACTTTGACATGCACTACCAGAGGTCCTGAGTTAGAATCCAGCAACCTTAACTAGGCTGGTGGCTAGCAACCATCAGCATGAGAATCTAATGCCCTCTTTTTGGTGGGTGGAAGGCAGCCACGGTGTACTTACACATATAATAAATAAATCTTTTTTTTAAAAAGAAATCGAAATATGATTTTTCCCTTCCGTTTTCTTCACTGGCCTCCTCCCATGTACTCTCTCACGTTGATACCGCAGCATCTCCTCGCCAGCAACAGCAGCAACAACAACGGAGGAGGAAAGAGAGAGATAGGCAGTGGTATCAGTATGAAGAGAAGGAGGATGGAAGAAAAAAAAATGGGAAAATCTCTGGGGCATTCTACTAAAAAAGGCCTAACTATATATATGGTGCTACTATTAGCACACCAGCATCTCATTACAACACTGCTTTCCTGAGACAGAACCTCTGAAGAAAAGAATCTGATCTCTCAACCAATTCAAAATAAAAAAATAGCGACTGGTCGTGGTGGCACAGCACAATCCCAGCACTCAGGGGCTGGGGGGCAGGCGGTTTTTTGGTTAAATTGGCCTGATGCACAGTGCAGGTTCCGGGACAGCCAGGGGCTGCAGAAAAACCTGTGCTGAAAAAAAAAAAAAAAAAAAAAAAGTCAGCATTCTCCTGCGACAGTCCAGGGGATCAGGAAGGCGTCGTGTAAGATCTCTCTCAGCAACAACAGCGGAGGGAGGAGGAAGAGAGAGATCTACACGTATAAATCTCGAACTTGTCATCAGGAAAGAATAACTAGACATTTTTTTTTTCAAGACACTGCTATATGTGTAACCCCTATCCTCCCCAAATCCTACTAAGAGTGTCCCCTGCCTGACTCCCACTGCTGAGAGAAATATCTAATATGAGCTGGAGATGGCTCATCGGTTAAGGAGCACTGACTGCTCTTCCAGGGTCCTGAGTTCAAAATCCCAGCATGCCAACCAAATGGTGACTAGCAACCATCTATGGCAGGATCTAATGCCCTCTTTTGGTGGATGCAGAGAGGCAAGCCTGGTAGCACATATAATAAATAAATCTTTTTTAAAATCTGAAATAGTTGCTTTTCCCACATACCGTTTTTCTTCCTCCAACTCCTCCCATGTACTCTGCATTGTCCCCTTGCCTTATCTCTCTCAACAACATACCAACTTGGAGGAGGAGGAAAAGAAGAAGATAAGCAGTGGTATCAGACACTTGAGTACATGGGAGGAGTTGGAGGGGGGAAGGAAGAAAACGGAAGGGGAAAAATCATATTTCAATTTTTTTTTAAAAGTTGCTTTAACCTGCCACCATTACATAACATATAGACCCATACAATTCAATTCATTTGGTGGTGCTTGCTTCACACCATTTTAGTTGCTGGGATTTGAGAACTCAGGGACTCTGGAAAACAGTCAGTGCTCTTAACCGATGAACCCTCTCTCAGCCCATATTTCAGTTTCAAATAAAAACTCTTGAGCTGGTCGTGGTGGCTGCACCTTTAATCGGCACTCAGGAGGCAGAGGCAGGCGGATTTCTGGAGTTCAGGCCAGCCTGATCTACAAGGGAGTGAGTTCCAGGACAGCCAGGGCTGCAGGGGAAAGAGCCTGTCAGAAAAAAAAAAAAAAAAAAAAAGTCTAGCATTCTCACCGACGGTCCGAGATCAGGCGTCGTGGCAGATCTCTTCTCAGCAGCAACAACGGAGGGGAGGGGAGGGAAAGGGGAGAGAGATCTACACGACCTCTGTCTCAGACTATCATAGAAAGAAAGGCTCTTGAACTCAGGCTCAGGAGTCCCCTATCCTGGGCTCACTGCACCTGAATCGAGCTGGCCTGGCTCCGCCTGCCCTCTGAGAATGCTGGATGTCTCCCTACCTGACTCAGAAAGACCCCCAACCTGAAATCTAGCTCTCGAGCTGGAGAGATGGCTCTCATCGGTTGAGAGCACTGACTGCTCTTACAGAGGTCCTGAGTTCAGAATGTAGCCAGCAAATGGTGGCCTGTAACCATCTGTAATGAGATCTAATGCCCTCTTGTAGATTGACAACCACAGTGTGCCCTACATATAATAATAAATCTTTTAAAAAAATCTGAAATATGATTTTTGCCTTCGTTTTCTTCTCCAACTCCTCCCATGTATAGGTGATACCACTGCTTGTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGGAAGAGAAATAAGCAGTGGTATCAACAGTACATGGGAAGGAGTTGGAGGAAGAGAAAACGGAAGGAAAAATCATTATTTTACAGATTTCTTTTTTTTAAAAAGATTTATTTATTATATATGTAAGTACACTGTAACATGAATTTTTCAAAAGAGAGCATTAGATCTCATTACAGACGTTTGCCATACCATTTGGTTTGCTGGGATTTGAACTCAGCCCTCTGGAAGAGAGCAGACAAATTAGCACTGCTAACCAAGTGAGCCATCTCTCAGCCCAGCCTCAGGAACCTCTCACAGCTTGGAACCTCAGTTTGGCTTTGGCAGGCCTGTATGCATGTAACTTGT"

ANCHOR_SIZE = 4
K_SIZE = 7
THRESH_GAP = 50
THRESH2 = THRESH_GAP

def align(target, query):
	ntQ = 0
	ntT = 0
	nbError = 0
	while ntQ < len(query) - 1:
		#~ print(ntQ, len(query), ntT, len(target))
		if nbError > 1:
			return False
		if ntT >= len(target) - 1:
			break
		if query[ntQ] == target[ntT]:
			ntQ += 1
			ntT += 1
		elif query[ntQ+1] == target[ntT]:
			ntQ += 1
			nbError += 1
		elif query[ntQ] == target[ntT+1]:
			ntT += 1
			nbError += 1
		else:
			nbError += 1
	return True


def merge_clusters(clusters):
	#~ cluster = clusters.copy()
	cluster = []
	for element in clusters:
		cluster.append([element])
	ind1 = len(cluster) - 1
	ind2 = len(cluster) - 1
	while ind1 >= 0:
		while ind2 >= 0:
			if ind1 != ind2 :
				if not (len(cluster[ind1]) == 0 or len(cluster[ind2]) == 0):
					if align(cluster[ind1], cluster[ind2]):
						cluster[ind1] += cluster[ind2]
						cluster[ind2] = []
			ind2 -= 1
		ind1 -= 1
	clustersToReturn = []
	for c in cluster:
		if len(c) > 0:
			clustersToReturn.append(c)
	return clustersToReturn

def get_quasi_kmers(sequence, kSize, anchorSize):
	headToTails = dict()
	positionInSequence = 0
	while positionInSequence < len(sequence) - kSize + 1:
		kmer = sequence[positionInSequence:positionInSequence+kSize]
		kmer = getCanonical(kmer)
		anchor = kmer[:anchorSize]
		tail = kmer[anchorSize:]
		if anchor in headToTails.keys():
			headToTails[anchor].append(tail)
		else:
			headToTails[anchor] = [tail]
		positionInSequence += 1
	return headToTails

def cluster_tails(headToTails):
	headToClusters = dict()
	for head in headToTails.keys():
		tails = headToTails[head]
		if len(tails) > 1:
			clusters = merge_clusters(tails)
		else:
			clusters = tails
		headToClusters[head] = clusters
	return headToClusters

def get_position_quasi_kmer(qKToPosition, kmer, headToClusters, position, anchorSize):
	kmer = getCanonical(kmer)
	anchor = kmer[:anchorSize]
	tail = kmer[anchorSize:]
	for cluster in headToClusters[anchor]:
		if tail in cluster:
			kmerForPosition = anchor + max(cluster)
			if kmerForPosition not in qKToPosition.keys():
				qKToPosition[kmerForPosition] = position
				return position
			else:
				return qKToPosition[kmerForPosition]


def autocorrel(x,N,i,M):
    C = np.zeros(N)
    for k in range(i,i+M):
        for n in range(N):
            C[n] += x[k]*x[k-n]
    return C/M        
               
def getReverseComplement(sequence):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a':'t', 'c':'g','g':'c','t':'a'}
	return "".join(complement.get(base, base) for base in reversed(sequence))
def getReverseComplement(sequence):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a':'t', 'c':'g','g':'c','t':'a'}
	return "".join(complement.get(base, base) for base in reversed(sequence))

def getCanonical(sequence):
	rc = getReverseComplement(sequence)
	if rc > sequence:
		return rc
	else:
		return sequence


def split(periods, readList):
	allSplits = []
	for nb,read in enumerate(readList):
		period = periods[nb]
		splitSeqList = []
		if period is not None:
			period = int(period)
			position = 0
			while position < len(read):
				splitSeqList.append(read[position:position+period])
				position += period
		allSplits.append(splitSeqList)
	return allSplits

def writeFileForMsa(allSplits):
	#~ #TODO add a track of files that were given to MSA/not given in order to put them together for the consensus
	fileNb = 0
	files = []
	splitToWritten = []
	for splits in allSplits:
		totalFileNb = 0
		splitToWritten.append(dict())
		fileNb += 1
		if len(splits) > 2:
			ind = 0
			toWrite = ''
			for seq in splits:
				if len(seq) > 50:
					ind += 1
					toWrite += ">" + str(ind) + "\n" + seq.lower() + "\n"
					splitToWritten[-1][totalFileNb] = (True, fileNb)
				else:
					splitToWritten[-1][totalFileNb] = (False, seq)
				totalFileNb += 1
			if len(toWrite) > 1:
				toWrite = toWrite[:-1]
				fMSA = open("file_"+ str(fileNb)+".fa", 'w')
				fMSA.write(toWrite)
				fMSA.close()
				getPOA("file_"+ str(fileNb)+".fa", "msa_" + str(fileNb) +".fa")
				files.append(fileNb)
				#~ splitToWritten[-1][totalFileNb] = (True, fileNb)
			#~ else:
				#~ splitToWritten[-1][totalFileNb] = (False, seq)
		
	return files, splitToWritten
	#### test
	# fileNb = 1
	# getPOA("test_offset_poa.fa", "msa_" + str(fileNb) +".fa")

	# return [1]

def subprocessLauncher(cmd, argstdout=None, argstderr=None,	 argstdin=None):
		args = shlex.split(cmd)
		p = subprocess.call(args, stdin = argstdin, stdout = argstdout, stderr = argstderr)
		return p


def getPOA(reads,  outFileName, matrixFile="/home/marchet/detection-consensus-isoform/poa-graph/blosum80.mat"):
	cmdPOA = "/home/marchet/detection-consensus-isoform/poa-graph/poa -reads_fasta " + reads + " -pathMatrix " + matrixFile
	subprocessLauncher(cmdPOA)
	cmdMv = "mv default_output_msa.fasta " + outFileName 
	subprocess.check_output(['bash','-c', cmdMv])
					

def findGapStretches(msaLine, threshGap, thresh2):
	prev = None
	countGap = 0
	positionsStretch = []
	pos = 0
	#~ for pos,ntResult in enumerate(correctedSequence):   # look for gaps in splitted/trimmed corrected read
	for nt in msaLine:   # look for gaps in splitted/trimmed corrected read
		if prev == ".":
			if nt == "." and countGap > 0:  # gaps are dots in msa file
				countGap += 1
			if nt == "." and countGap == 0:
				countGap = 2
		if prev == None:
			if nt == ".":
				countGap += 1
		if nt != ".":
			if countGap > 0:
				positionsStretch.append([])
			countGap = 0
		if countGap >= threshGap:
				if len(positionsStretch) == 0:
					positionsStretch.append([pos-threshGap + 1, pos]) # start new stretch of gap with leftmost position
				else:
					if len(positionsStretch[-1]) == 0:
						positionsStretch[-1].extend((pos-threshGap + 1, pos))
					if len(positionsStretch[-1]) == 2:
						positionsStretch[-1][1] = pos # update position
		prev = nt
		pos += 1
	stretchTmp = []
	stretch = dict()
	for s in positionsStretch:
		if len(s) > 0:
			#~ print(s[1], s[0])
			if s[1] - s[0] > thresh2:
				stretchTmp.append([s[0], s[1]])

	if len(stretchTmp) > 1:
		newStretchTmp = [stretchTmp[0]]
		ind1 = 0
		ind2 = 1
		while ind1 < len(newStretchTmp):
			while ind2 < len(stretchTmp):
				if ind1 < ind2:
					if stretchTmp[ind2][0] - newStretchTmp[ind1][1] < thresh2:
						newStretchTmp[ind1][-1] = stretchTmp[ind2][1]
						ind2 += 1
					else:
						newStretchTmp.append(stretchTmp[ind2])
						ind1 += 1
						ind2 += 1
			ind1 += 1
		stretchTmp = newStretchTmp

	#beginning and end of sequence:
	if len(stretchTmp) > 0:
		if stretchTmp[0][0] < thresh2:
			stretchTmp[0][0] = 0
		if len(msaLine) - stretchTmp[-1][1] < thresh2:
			stretchTmp[-1][1] = len(msaLine) - 1
	#~ for s in stretchTmp:
		#~ stretch[s[0]] = s[1]
	return stretchTmp


def getExistingPositionsInLine(msaLine, stretches):
	existingPositions = [True] * len(msaLine)
	for s in stretches:
		rangePos = [x for x in range(s[0], s[1]+1)]
		for i in rangePos:
			existingPositions[i] = False
	return existingPositions
	
def consensusNt(ntList):
	stringList = str(ntList)
	val = 0
	ntToKeep = '.'
	for nt in ['a', 'c', 'g', 't', '.']:
		tmpVal = stringList.count(nt)
		if tmpVal > val:
			val = tmpVal
			ntToKeep = nt
	if ntToKeep == '.':
		return None
	return ntToKeep.upper()


def get_consensus(existingPositions, msaAllLines, splitToWritten, nb):
	wholeConsensus = ''
	listConsensus = []
	consensusNone = set()
	for indexLine,line in enumerate(msaAllLines):
		consensus = ''
		for position in range(len(line)):
			listNt = []
			#~ print(existingPositions[indexLine][position])
			if existingPositions[indexLine][position]:
				listNt.append(line[position])
				for otherLines in range(len(msaAllLines)):
					if otherLines != indexLine:
						if existingPositions[otherLines][position]:
							listNt.append(msaAllLines[otherLines][position])
				if len(listNt) > 2:
					nt = consensusNt(listNt)
					if nt is not None:
						consensus += consensusNt(listNt)
					else:
						consensusNone.add(position)
				else:
					if line[position] != '.':
						consensus += line[position].upper().rstrip()
					else:
						consensusNone.add(position)
						#~ print("---------------------->", consensus)
		#~ print("---------------------->", consensus)
		#~ wholeConsensus += consensus
		listConsensus.append(consensus)
	for i in range(len(splitToWritten[nb-1])):
		if splitToWritten[nb-1][i][0] is True :
			wholeConsensus += listConsensus[splitToWritten[nb-1][i][1]]
		else:
			#~ print(splitToWritten[i][0])
			wholeConsensus += splitToWritten[nb-1][i][1]
	
	
	return wholeConsensus

def consensusFromMSA(msaFiles, splitToWritten, readsList, periods):
	stretchesList = []
	existList = []
	for nb in msaFiles:
		stretchesList.append([])
		existList.append([])
		msa = open("msa_"+ str(nb)+".fa", 'r')
		msaAllLines = []
		for line in msa.readlines():
			if not ">" in line:
				#~ stretchesList[-1].append([])
				#~ existList[-1].append([])
				#~ print(line)
				msaAllLines.append(line)
				stretches = findGapStretches(line, THRESH_GAP, THRESH2)
				#~ print(stretches)
				stretchesList[-1].append(stretches)
				existingPositions = getExistingPositionsInLine(line, stretches)
				existList[-1].append(existingPositions)
		#~ print("*", stretchesList[-1])
		msa.close()
		#~ input("appel")
		consensus = get_consensus(existList[-1], msaAllLines, splitToWritten, nb)
		#~ print('llllllllllllllllllllllllll',len(readsList), nb)
		originalAndConsensus = alignConsensusToOriginal(readsList[nb-1], consensus)
		consensus = patchConsensusWithOriginal(originalAndConsensus)
		#~ print("################################")
		#~ print(consensus)
		#~ print("################################")
		#~ print(consensus)
		#~ print(len(consensus))
		consensusFile = open("consensus_"+ str(nb)+".fa", 'w')
		consensusFile.write(">consensus_len" + str(len(consensus)) + "|period_" + str(int(periods[nb-1])) +"\n" + consensus +"\n")
		#~ print("written file " + "consensus_"+ str(nb)+ "|period_" + str(periods[nb-1]) + ".fa")
		consensusFile.close()
		if len(consensus) > 0 and consensus is not None:
			compute_consensus_from_period.getFinalConsensus(periods[nb-1], consensus)
	#~ return stretchesList, existList


def alignConsensusToOriginal(readSeq, consensus):
	c = open("consensus_original.fa",'w')
	c.write(">o\n" + readSeq + "\n>c\n" + consensus)
	c.close()
	getPOA("consensus_original.fa", "msa_consensus_original.fa")
	msa = open("msa_consensus_original.fa", 'r')
	originalAndConsensus = []
	for line in msa.readlines():
		if '>' not in line:
			originalAndConsensus.append(line.rstrip())
	return originalAndConsensus
		
def patchConsensusWithOriginal(originalAndConsensus, thresh=THRESH_GAP, thresh2=THRESH2):
	stretchO = findGapStretches(originalAndConsensus[0], thresh, thresh2)
	stretchC = findGapStretches(originalAndConsensus[1], thresh, thresh2)
	toCorrect = [False] * len(originalAndConsensus[0])
	toRemove = [False] * len(originalAndConsensus[0])
	for s in stretchC:
		for i in range(s[0], s[1]+1):
			toCorrect[i] = True
	for s in stretchO:
		for i in range(s[0], s[1]+1):
			toRemove[i] = True
	index = 0
	newConsensus = ''
	for ntO, ntC in zip(originalAndConsensus[0], originalAndConsensus[1]):
		if toCorrect[index]:
			newConsensus += ntO
		elif not toRemove[index]:
			if ntC != '.' :
				newConsensus += ntC
		index += 1
	return newConsensus

def getPeriod():
	periodList = []
	readList = []
	reads = open(sys.argv[1], 'r')
	for line in reads.readlines():
		if ">" not in line:
			readList.append(line.rstrip())

	for read in readList:
		sequence = read
		count = 0
		positionKmers = []
		positionInSequence = 0
		kSize = K_SIZE
		anchorSize = ANCHOR_SIZE
		qKToPosition = dict()
		headToTails = get_quasi_kmers(sequence, kSize, anchorSize)
		headToClusters = cluster_tails(headToTails)
		while positionInSequence < len(sequence) - kSize + 1:
			kmer = sequence[positionInSequence:positionInSequence+kSize]
			position = get_position_quasi_kmer(qKToPosition, kmer, headToClusters, positionInSequence, anchorSize)
			positionKmers.append(position)
			positionInSequence += 1

		tfd = numpy.fft.fft(positionKmers)
		tfd = fftpack.fft(positionKmers)
		S = numpy.square(numpy.absolute(tfd))/len(positionKmers)
		autocorrelation = numpy.fft.ifft(S)

	
		cb = np.array(list(autocorrelation))
		peakind = signal.find_peaks_cwt(cb,np.arange(1,1500))
		#~ peakind = signal.find_peaks_cwt(cb,np.arange(1,2000))
		intervals = list((peakind))
		print(intervals)
		# #print(intervals)
		if len(intervals)> 1:
			i = 0
			newIntervals = []
			diffSmall = False
			while i < len(intervals) - 1:
				diff = intervals[i+1] - intervals[i]
				if diff < 50:
					newIntervals.append(np.mean([intervals[i],intervals[i+1]]))
					diffSmall = True
				else:
					if not diffSmall:
						newIntervals.append(intervals[i])
					else:
						diffSmall = False
				i+=1
			if intervals[-1] - intervals[-2] > 50:
				newIntervals.append(intervals[-1])
			i = 0
			diffForMean = []
			# print(newIntervals)
			while i < len(newIntervals) - 1:
				diff = newIntervals[i+1] - newIntervals[i]
				if diff > 50:
					diffForMean.append(diff)
				i += 1
			if len(diffForMean) > 1:
				period = np.mean(diffForMean)
			# elif diffForMean == 1:
				periodList.append(period)
				print(period)
			else:
				periodList.append(None)
		else:
			periodList.append(None)
		pyplot.figure(figsize=(10,6))
		positionsSolid = [x for x in range(len(autocorrelation))]
		pplot(np.array(positionsSolid),np.array(autocorrelation), peakind)
		#~ # print(indexes)
		#~ #plot(positions,autocorrelation)
		#~ draw()
		#~ show()
	return periodList, readList

periods, reads = getPeriod()
print(periods)
allSplits = split(periods, reads)
msaFiles, splitToWritten = writeFileForMsa(allSplits)
consensusFromMSA(msaFiles, splitToWritten, reads, periods)





#### TODO split sequences according to coordinates given in periods

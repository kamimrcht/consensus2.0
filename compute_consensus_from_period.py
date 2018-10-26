import sys
import os
import shlex, subprocess
from subprocess import Popen, PIPE, STDOUT
import numpy as np

THRESH_GAP = 8
THRESH2 = THRESH_GAP

A = 'AAGCAGTGGTATCAACGCAGAGTACAT'
B = 'ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT'

def getReverseComplement(sequence):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a':'t', 'c':'g','g':'c','t':'a'}
	return "".join(complement.get(base, base) for base in reversed(sequence))
	
def splitConsensus(period, consensus):
	period = int(period)
	splitSeqList = []
	position = 0
	while position < len(consensus):
		splitSeqList.append(consensus[position:position+period])
		position += period
	#~ print(splitSeqList)
	#split in two
	readsToWrite = ''
	#~ reads = []
	if len(splitSeqList) > 1:
		halfT = int(period/2.0)
		for i,periodSplit in enumerate(splitSeqList):
			#~ if not (i==0 or i == len(splitSeqList) - 1): #first and last can be broken
				if len(periodSplit) == period:
						readsToWrite += ">" + str(i) + "_1\n" + periodSplit[:halfT].lower() +"\n"
						readsToWrite += ">" + str(i) + "_2\n" + getReverseComplement(periodSplit[halfT:]).lower() +"\n"
				else:
					if halfT >= len(periodSplit):
						readsToWrite += ">" + str(i) + "_1\n" + periodSplit[:halfT].lower() +"\n"
					else:
						readsToWrite += ">" + str(i) + "_1\n" + periodSplit[:halfT].lower() +"\n"
						readsToWrite += ">" + str(i) + "_2\n" + getReverseComplement(periodSplit[halfT:]).lower() +"\n"
						#~ reads.append(periodSplit[:halfT])
						#~ reads.append(getReverseComplement(periodSplit[halfT:]))
		#~ print(readsToWrite)
		readsToWrite += ">bab\n" + getReverseComplement(B).lower() + A.lower() + B.lower() + "\n"
		readsListFile = open("list_splits_reads_msa.fa",'w')
		readsListFile.write(readsToWrite)
		readsListFile.close()
		getPOA("list_splits_reads_msa.fa", "msa_for_consensus_final_read.fa")
		return True
	else:
		return False

def subprocessLauncher(cmd, argstdout=None, argstderr=None,	 argstdin=None):
		args = shlex.split(cmd)
		p = subprocess.call(args, stdin = argstdin, stdout = argstdout, stderr = argstderr)
		return p


def getPOA(reads,  outFileName, matrixFile="/home/marchet/detection-consensus-isoform/poa-graph/blosum80.mat"):
	cmdPOA = "/home/marchet/detection-consensus-isoform/poa-graph/poa -reads_fasta " + reads + " -pathMatrix " + matrixFile
	subprocessLauncher(cmdPOA)
	cmdMv = "mv default_output_msa.fasta " + outFileName 
	subprocess.check_output(['bash','-c', cmdMv])


def getMsaLines(msaFileName = "msa_for_consensus_final_read.fa"):
	num_lines = sum(1 for line in open(msaFileName))

	msaFile = open(msaFileName, 'r')
	forward = []
	reverse = []
	index = 0
	for ind,line in enumerate(msaFile.readlines()):
		if not ">" in line:
			if not ind == num_lines - 1:
				if index % 2:
					#~ reverse.append(line.rstrip())
					forward.append(line.rstrip())
				else:
					#~ forward.append(line.rstrip())
					reverse.append(line.rstrip())
				index += 1
			else:
				adaptAndBell = [line.rstrip()]
	return forward, reverse, adaptAndBell


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



def getStretchesFwd(forward, stretchesFwd, thresh=THRESH_GAP, thresh2=THRESH2):
	#~ stretchesFwd = []
	for fwd in forward:
		stretches = findGapStretches(fwd, thresh, thresh2)
		if len(stretches) > 0:
			if stretches[0][0] == 0:
				stretchesFwd.append(stretches[0][1])
	if len(stretchesFwd) > 1:
		return int(np.mean(stretchesFwd))
	elif len(stretchesFwd) == 1:
			return stretchesFwd[0]
	else:
		return 0
	

def getStretchesRev(reverse, stretchesRev, thresh=THRESH_GAP, thresh2=THRESH2):
	#~ stretchesRev = []
	for rev in reverse:
		stretches = findGapStretches(rev, thresh, thresh2)
		if len(stretches) > 0:
			if stretches[-1][1] == len(reverse[0]) - 1:
				stretchesRev.append(stretches[-1][0])
	if len(stretchesRev) > 1:
		return int(np.mean(stretchesRev))
	elif len(stretchesRev) == 1:
			return stretchesRev[0]
	else:
		return 0

def consensusNt(ntList):
	stringList = str(ntList)
	val = 0
	ntToKeep = '.'
	for nt in ['a', 'c', 'g', 't', '.']:
		tmpVal = stringList.count(nt)
		if tmpVal > val:
			val = tmpVal
			ntToKeep = nt
	return ntToKeep.upper()


def getStretchesAdaptAndBell(adaptAndBell,thresh=THRESH_GAP, thresh2=THRESH2):
	stretchesFwd = []
	stretchesRev = []
	for ab in adaptAndBell:
		stretches = findGapStretches(ab, thresh, thresh2)
		if len(stretches) > 0:
			maxS = []
			diff = 0
			for s in stretches:
				if s[1]-s[0] >= diff:
					maxS = [s[0],s[1]]
					diff = s[1]-s[0]
			stretchesFwd.append(maxS[0] - 1)
			stretchesRev.append(maxS[-1] + 1)
			#TODO find longest
			#~ if stretches[0][0] != 0:
				#~ stretchesFwd.append(stretches[0][0] - 1)
			#~ if stretches[-1][-1] != len(ab) - 1:
				#~ stretchesRev.append(stretches[-1][-1] - 1)
	return stretchesFwd,stretchesRev
	

def computeConsensus(forward, reverse, limL, limR):
	allLines = forward + reverse
	consensusTmp = ''
	for position in range(len(allLines[0])):
		ntList = []
		for seq in allLines:
			ntList.append(seq[position])
		consensusTmp += consensusNt(ntList)
	consensusACRegion = ''
	
	for i in range(limL, limR):
		nt =  consensusTmp[i]
		if nt!= '.':
			consensusACRegion += nt
	return consensusACRegion

#~ consensusFile = open("consensus_1.fa", 'r')
#~ period = int(consensusFile.readline().rstrip().split('_')[2])
#~ print("period:" ,period)
#~ consensus = consensusFile.readline().rstrip()
#~ period = 2143

def getFinalConsensus(period, consensus):
	continu = splitConsensus(period, consensus)
	if continu:
		forward, reverse, adaptAndBell = getMsaLines()
		stretchF, stretchR = getStretchesAdaptAndBell(adaptAndBell)
		if len(stretchF) > 0:
			limitLeft = stretchF[0]
		else:
			limitLeft = getStretchesFwd(forward, stretchF)
		if len(stretchR) > 0:
			limitRight = stretchR[0]
		else:
			limitRight = getStretchesRev(reverse, stretchR)
		consensus = computeConsensus(forward, reverse, limitLeft, limitRight)
		#~ print(consensus)
	else:
		consensus = '.'
	final = open("final_consensus.fa", 'w')
	final.write(">c\n" + consensus )

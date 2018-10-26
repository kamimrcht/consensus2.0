import sys
import os
import shlex, subprocess




def subprocessLauncher(cmd, argstdout=None, argstderr=None,	 argstdin=None):
		args = shlex.split(cmd)
		p = subprocess.call(args, stdin = argstdin, stdout = argstdout, stderr = argstderr)
		return p


def getPOA(reads,  outFileName, matrixFile="/home/marchet/detection-consensus-isoform/poa-graph/blosum80.mat"):
	cmdPOA = "/home/marchet/detection-consensus-isoform/poa-graph/poa -reads_fasta " + reads + " -pathMatrix " + matrixFile
	subprocessLauncher(cmdPOA)
	cmdMv = "mv default_output_msa.fasta " + outFileName 
	subprocess.check_output(['bash','-c', cmdMv])


consensusFile = open("final_consensus.fa",'r')
perfectFile = open("consensus_test.fa", 'r')

w = perfectFile.readline() + perfectFile.readline().lower() + "\n" +  consensusFile.readline() + consensusFile.readline().lower().rstrip()

verif = open("verif_msa_final.fa",'w')
verif.write(w)
verif.close()
getPOA("verif_msa_final.fa", "verif_result_final.fa")

msa = open("verif_result_final.fa",'r')
msaLines = []
for line in msa.readlines():
	if ">" not in line:
		msaLines.append(list(line.rstrip()))
		

score = len(msaLines[0]) - msaLines[0].count('.')
for nt1,nt2 in zip(msaLines[0], msaLines[1]):
	if nt1 != nt2:
		#~ if nt1 == '.' or nt2 == '.':
			#~ score -= 1
		score -= 1
print("score: " + str(round(score *1.0/(len(msaLines[0]) - msaLines[0].count('.')) * 100,3)) + "%")

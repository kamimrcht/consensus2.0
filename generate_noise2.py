import random
import os
import sys
import argparse

LENGTH = 1000
ADAPTORS = ['AAGCAGTGGTATCAACGCAGAGTACAT'.lower(), 'CTACACGACGCTCTTCCGATCT'.lower(), 'ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT'.lower() ] #[0] 5', [1] 3', [2] bell
#~ ADAPTORS = ['ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGATAAGCAGTGGTATCAACGCAGAGTACAT'.lower(), 'CTACACGACGCTCTTCCGATCT'.lower(), 'ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT'.lower() ] #[0] bell+5', [1] 3', [2] bell
NUMBER_SEQ = 1
ERROR_RATE = 15

NUCLEOTIDES = ['A', 'C', 'G', 'T']


def reverseComplement(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
	return "".join(complement.get(base, base) for base in reversed(seq))



def subs(nt):
	if nt.lower() ==  'C':
		return random.choice(['A','G','T'])
	elif nt.lower() == 'G':
		return random.choice(['C','A','T'])
	elif nt.lower() == 'T':
		return random.choice(['C','A','G'])
	else:
		return random.choice(['C','G','T'])


#~ def getErrors(seq):
	#~ errorTypes = ["subs", "ins", "del"]
	#~ newSeq = ''
	#~ for nt in seq:
		#~ pc = random.randint(1, 100)
		#~ if pc < ERROR_RATE:
			#~ err = random.choice(errorTypes)
			#~ if err == "subs":
				#~ newSeq += subs(nt)
			#~ elif err == "ins":
				#~ newSeq += random.choice(NUCLEOTIDES)
		#~ else:
			#~ newSeq += nt
	#~ return newSeq



def putErrors(nt):
	errorTypes = ["subs", "ins", "del"]
	newSeq = ''
	pc = random.randint(1, 100)
	if pc < ERROR_RATE:
		err = random.choice(errorTypes)
		if err == "subs":
			return subs(nt)
		elif err == "ins":
			return nt + random.choice(NUCLEOTIDES)
		else:
			return ''
	else:
		return nt

if __name__ == "__main__":
	######## parse arguments #######
	currentDirectory = os.path.dirname(os.path.abspath(sys.argv[0]))
	installDirectory = os.path.dirname(os.path.realpath(__file__))
	parser = argparse.ArgumentParser()
	parser.add_argument('-e', nargs='?', type=float, action="store", dest="errorRate", help="Error rate", default=10)
	parser.add_argument('-n', nargs='?', type=int, action="store", dest="nbSeq", help="Number of sequences", default=5)
	parser.add_argument('-c', nargs='?', type=int, action="store", dest="cSize", help="Read Size", default=1000)
	args = parser.parse_args()

	#~ if len(sys.argv) < 2:
		#~ parser.print_help()
		#~ sys.exit()
	#~ else:
	if args.errorRate is not None:
		ERROR_RATE = args.errorRate
	if args.nbSeq is not None:
		NUMBER_REPET = args.nbSeq
	if args.cSize is not None:
		LENGTH = args.cSize
	sequences = []
	errSequences = []
	consensus = []

	a = 'AAGCAGTGGTATCAACGCAGAGTACAT'
	b = 'ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT'
	c = ''
	for i in range(LENGTH):
		c+= random.choice(NUCLEOTIDES)
	period = a + c + b + reverseComplement(c) + reverseComplement(a) + b
	seq = ''
	for i in range(NUMBER_REPET):
		seq += period
	seqE = ''
	seqP = ''
	seqs=[seq]
	for n in seqs: # generate n sequences
		seqE += ">s\n"
		seqP += ">s\n"
		for nt in n:
			seqE += putErrors(nt)
		seqE += "\n"
		seqP += n + "\n"
	seqE = seqE[:-1]
	seqP = seqP[:-1]

	#~ print(seqE)
	fileP = open('perfect_seq_test.fa', 'w')
	fileE = open('error_seq_test.fa', 'w')
	fileC = open('consensus_test.fa', 'w')
	fileP.write(seqP)
	fileE.write(seqE)
	fileC.write(">consensus\n"+c+"\n")
	#~ fileC.write(">consensus\n"+a+c+"\n")
	fileP.close()
	fileE.close()
	fileC.close()
	#~ for i,s in enumerate(consensus):
		#~ fileC.write(">" + str(i)+"\n")
		#~ fileC.write(s+"\n")

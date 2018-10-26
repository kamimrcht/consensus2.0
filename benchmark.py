import sys
import os
import shlex, subprocess
from subprocess import Popen, PIPE, STDOUT
import numpy as np

def subprocessLauncher(cmd, argstdout=None, argstderr=None,	 argstdin=None):
		args = shlex.split(cmd)
		p = subprocess.call(args, stdin = argstdin, stdout = argstdout, stderr = argstderr)
		return p

#~ for error in [10]:
	#~ for length in [500]:
		#~ for cycles in [3]:
for i in range(10):
	for error in [10,15,20]:
		for length in [100,500,1000]:
			for cycles in [2,3,4,5]:
				print('#error ' +str(error) + " length " +str(length) + " cycles " + str(cycles))
				cmd = "rm final_consensus.fa "
				subprocessLauncher(cmd)
				cmd = "python generate_noise2.py -e " +str(error) + " -c " + str(length) + " -n " + str(cycles)
				subprocessLauncher(cmd)
				cmd = "python3.6 autocorr_quasi_kmer.py error_seq_test.fa"
				subprocessLauncher(cmd)
				#~ cmd = "python validate_consensus.py"
				#~ subprocessLauncher(cmd)
				#~ cmd = "python compute_consensus_from_period.py"
				#~ subprocessLauncher(cmd)
				cmd = "python validate_final_consensus.py"
				subprocessLauncher(cmd)


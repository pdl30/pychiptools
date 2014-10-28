#!/usr/bin/python

########################################################################
# 1 August 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import os, sys, re
import subprocess
import argparse
from collections import defaultdict
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import IntVector, FloatVector, StrVector

def diffReps(inv_conds, cond1, cond2, chromsizes, method="nb", controls=None, counts=None):
	#Can't deal with reps yet
	if counts:
		sample1 = inv_conds[cond1][0]
		sample2 = inv_conds[cond2][0]
		c1 = read_gapdh_counts_file(counts[sample1])
		c2 = read_gapdh_counts_file(counts[sample2])
		c1 = float(c1)/1000
		c2 = float(c2)/1000
		output = open("norm.txt", "w")
		output.write("treatment\t{}\n".format(c1)),
		output.write("control\t{}\n".format(c2)),
		output.close()
		command = "diffReps.pl --treatment {0} --control {1} --report {2}_vs_{3}_diffReps.txt --chrlen {4} -me {5} --norm norm.txt --nproc 8".format(inv_conds[cond1][0], inv_conds[cond2][0], 
			cond1, cond2, chromsizes, method )
	elif controls:
		backt1 = []
		backt2 = []
		for sample in t1:
			backt1.append(controls[sample])
		for sample in t2:
			backt2.append(controls[sample])
		command = "diffReps.pl --treatment {0} --control {1} --btr {2} --bco {3} --report {4}_vs_{5}_diffReps.txt --chrlen {6} -me {7}".format(inv_conds[cond1], inv_conds[cond2], 
			backt1, backt2, cond1, cond2, chromsizes, method)
	else:
		command = "diffReps.pl --treatment {0} --control {1} --report {2}_vs_{3}_diffReps.txt --chrlen {4} -me {5}".format(inv_conds[cond1][0], inv_conds[cond2][0], cond1, cond2, 
			chromsizes, method)
	print command
	subprocess.call(command.split())

def read_gapdh_counts_file(ifile):
	with open(ifile) as f:
		header= next(f)
		header= header.rstrip()
		word = header.split("\t")
	return word[1]



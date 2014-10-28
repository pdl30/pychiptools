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

def read_flagstat(ibam):
	ifile = re.sub("_sort.bam", "_stats.txt", ibam)
	with open(ifile) as f:
		header= next(f)
		word = header.split(" ")
	return word[0]

def read_macs_output(ifile):
	tvalue1 = ""
	cvalue2 = ""
	with open(ifile) as f:
		for line in f:
			line = line.rstrip()
			m1 = re.search("filtering in treatment", line)
			m2 = re.search("filtering in control", line)
			if m1:
				word1 = line.split("treatment: ")
				tvalue1 = int(word1[1])
			if m2:
				word2 = line.split("control: ")
				cvalue1 = int(word2[1])
			value1 = min(tvalue1, cvalue1)
	return value1

def bgdiff(inv_conds, cond1, cond2, controls, chiptype, bedgraphs=None):
	#Can't deal with reps. Also consider options without controls
	#Can handle with custom inputs, use samtools flagstat for sizes. Will assume these are in same dir as bam and named "_stats.txt"
	#sizes from: https://github.com/taoliu/MACS/wiki/Call-differential-binding-events
	if chiptype=="histone":
		extsize=147
		g = 120
	else:
		extsize=200
		g = 100

	if bedgraphs:
		sample1 = inv_conds[cond1][0]
		sample2 = inv_conds[cond2][0]
		t1flag = read_flagstat(inv_conds[cond1][0])
		t2flag = read_flagstat(inv_conds[cond2][0])
		c1flag = read_flagstat(controls[sample1])
		c2flag = read_flagstat(controls[sample2])
		value1 = min(t1flag, c1flag)
		value2 = min(t2flag, c2flag)
		command = "macs2 bdgdiff --t1 {0} --c1 {1} --t2 {2} --c2 {3} --d1 {4} --d2 {5} -g {6} -l {7} --o-prefix {8}_vs_{9}".format(bedgraphs[sample1], bedgraphs[controls[sample1]], 
			bedgraphs[sample2], bedgraphs[controls[sample2]], value1, value2, g, extsize, cond1, cond2)
		print command
		subprocess.call(command.split())
	else:
		sample1 = inv_conds[cond1][0]
		sample2 = inv_conds[cond2][0]
		control1 = controls[sample1]
		control2 = controls[sample2]
		#Builds bedGraphs from bam files
		command1 = "macs2 callpeak -B -t {0} -c {1} -n {2} --nomodel --extsize {3}".format(sample1, control1, cond1, extsize)
		command2 = "macs2 callpeak -B -t {0} -c {1} -n {2} --nomodel --extsize {3}".format(sample2, control2, cond2, extsize)
		macs2_output = open("macs_output1.txt", "w")
		subprocess.call(command1.split(), stderr=macs2_output)
		macs2_output.close()
		macs2_output = open("macs_output2.txt", "w")
		subprocess.call(command2.split(), stderr=macs2_output)
		macs2_output.close()
		#Get library size from macs2 output and because we are scaling down, get the smaller value
		value1 = read_macs_output("macs_output1.txt")	
		value2 = read_macs_output("macs_output2.txt")

		command = "macs2 bdgdiff --t1 {0}_treat_pileup.bdg --c1 {0}_control_lambda.bdg --t2 {1}_treat_pileup.bdg --c2 {1}_control_lambda.bdg --d1 {2} --d2 {3} -g {4} \
		 	-l {5} --o-prefix {0}_vs_{1}".format(cond1, cond2, value1, value2, g, extsize)
		subprocess.call(command.split())


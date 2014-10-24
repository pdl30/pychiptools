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
#Adjusts peaks files to 3 column bed files for MANorm
def peak_adjust(ibed, obed):
	output = open(obed, "w")
	with open(ibed) as f:
		header = next(f)
		for line in f:
			line = line.rstrip()
			word= line.split("\t")
			output.write("{}\t{}\t{}\n".format(word[0], word[1], word[2])),
	output.close()

def manorm(inv_conds, peaks, cond1, cond2, download=False):
	#  ./MAnorm.sh   sample1_peakfile[BED]     sample2_peakfile[BED]     sample1_readfile[BED]    sample2_readfile[BED]     sample1_readshift_lentgh[INT]       sample2_readshift_length[INT]    
	#Not finished and can't deal with reps. Also need to adjust peak files
	if download:
		path = os.getcwd()
		wget_command = "wget http://bcb.dfci.harvard.edu/~gcyuan/MAnorm/MAnorm_data/MAnorm_Linux_R_Package.zip"
		subprocess.call(wget_command.split())
		subprocess.call(["unzip", "MAnorm_Linux_R_Package.zip"])
		subprocess.call(["mv", "MAnorm_Linux_R_Package/MAnorm.r", path])
		subprocess.call(["mv", "MAnorm_Linux_R_Package/MAnorm.sh", path])
	sample1 = inv_conds[cond1][0] 
	sample2 = inv_conds[cond2][0] 
	peak_adjust(peaks[sample1], "tmp1.bed")
	peak_adjust(peaks[sample2], "tmp2.bed")
	command = "./MAnorm.sh tmp1.bed tmp2.bed {} {} 150 150".format(sample1, sample2)
	subprocess.call(command, shell=True) 
	subprocess.call(["mv", "MAnorm_result.xls", "{}_{}_MAnorm.xls".format(cond1, cond2)])

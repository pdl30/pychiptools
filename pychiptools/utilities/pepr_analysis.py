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

#Pepr analysis, tutorial found at https://code.google.com/p/pepr-chip-seq/wiki/PePrMainPage?tm=6
pepr = "/usr/local/lib/python2.7/dist-packages/PePr/PePr.py"

#For TF samples peak finding
#python path/PePr.py -i inputFiles -c chipFiles -n ExperimentName -f format --peaktype=sharp --remove_artefacts
#For histones:
#python path/PePr.py -i inputFiles -c chipFiles -n ExperimentName -f format --peaktype=broad  --remove_artefacts
#For differential binding analysis of transcription factors:
#python path/PePr.py -i inputFiles_group1 --input2 inputFiles_group2 -c chipFiles_group1 --chip2 chipFiles_group2 -n ExperrrimentName -f format --peaktype=sharp --diff --remove_artefacts
#Example of exact command for a transcription factor with two input files and two ChIP files:
#python path/PePr.py -i input-1.bed,input-2.bed -c chip-1.bed,chip-2.bed -n myTF -f bed --peaktype=sharp --remove_artefacts 

#if using bam files, they must be sorted and indexed
def pepr_diff_analysis(inv_conds, controls, cond1, cond2):
	sample1 = inv_conds[cond1][0] 
	sample2 = inv_conds[cond2][0] 
	name =  "{}_vs_{}".format(cond1, cond2)
	#Will have to adjust for reps
	command = "python {} -i {} --input2 {} -c {} --chip2 {} -n {} -f bed --peaktype=broad --remove_artefacts".format(pepr, controls[sample1], controls[sample2], sample1, sample2, name)
	print command
	subprocess.call(command, shell=True)
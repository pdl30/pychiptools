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

#Meryem warning about minOverlap, only create input for conditions being compared!
def create_diffbind_input(inv_conds, peaks, cond1, cond2, controls=None):
	output = open("samplesheet.csv", "w")
	#SampleID,Tissue,Factor,Replicate,Condition,bamReads,ControlID,bamControl,Peaks,PeakCaller
	c1 = 1
	c2 = 1
	if controls:
		output.write("SampleID,Tissue,Factor,Replicate,Condition,bamReads,ControlID,bamControl,Peaks,PeakCaller\n"),
		for bam in inv_conds[cond1]:
			name = re.sub(".bam", "", bam)
			cont_id = re.sub(".bam", "", controls[bam])
			output.write("{0},ESC,{1},{2},{3},{4},{5},{6},{7},bed\n".format(name, cond1, c1, cond1, bam, cont_id, controls[bam], peaks[bam])),
			c1 += 1
		for bam in inv_conds[cond2]:
			name = re.sub(".bam", "", bam)
			cont_id = re.sub(".bam", "", controls[bam])
			output.write("{0},ESC,{1},{2},{3},{4},{5},{6},{7},bed\n".format(name, cond2, c2, cond2, bam, cont_id, controls[bam], peaks[bam])),
			c2 += 1
	else:
		output.write("SampleID,Tissue,Factor,Replicate,Condition,bamReads,Peaks,PeakCaller\n"),
		for bam in inv_conds[cond1]:
			name = re.sub(".bam", "", bam)
			output.write("{0},ESC,{1},{2},{3},{4},{5},bed\n".format(name, cond1, c1, cond1, bam, peaks[bam])),
			c1 += 1
		for bam in inv_conds[cond2]:
			name = re.sub(".bam", "", bam)
			output.write("{0},ESC,{1},{2},{3},{4},{5},bed\n".format(name, cond2, c2, cond2, bam, peaks[bam])),
			c2 += 1
	output.close()
	if c1==1 or c2==1:
		reps=False
	else:
		reps=True
	return reps

#Differences to diffbind seem to be mmdiff can detect differences in peak shape. Needs Replicates!!!!
def write_mmdiff(cond1, cond2, reps, method="MMD"):
	rscript = "library('MMDiff')\n"
	rscript += "data <- dba(sampleSheet='samplesheet.csv', minOverlap=1)\n"
	rscript += "Peaks <- dba.peakset(data,bRetrieve=TRUE)\n"
	rscript += "dataProfiles <- getPeakProfiles(data, Peaks, bin.length=50, save.files=FALSE,run.parallel=FALSE)\n" #Takes a while!
	rscript += "dataNorm <- getNormFactors(dataProfiles, method = 'DESeq')\n"
	rscript += "dataDists <- compHistDists(dataNorm, method='MMD', overWrite=FALSE, NormMethod='DESeq',run.parallel=FALSE)\n"#V. Resource intensive and slow
	rscript += "dataPvals <- detPeakPvals(dataDists, group1=dataDists$masks$'{0}', group2=dataDists$masks$'{1}', name1='{0}', name2='{1}', method='MMD')\n".format(cond1, cond2)
	return rscript

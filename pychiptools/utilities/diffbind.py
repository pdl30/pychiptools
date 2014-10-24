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
import rpy2.robjects.numpy2ri as rpyn

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

#Downstream analysis could be find nearest feature using closestBed
#Results seem to be unreliable without reps!
def write_diffbind(cond1, cond2, reps):
	rscript = "library('DiffBind')\n"
	rscript += "data <- dba(sampleSheet='samplesheet.csv', minOverlap=1)\n"
	rscript += "dataCount=dba.count(data)\n"
	rscript += "dataContrast = dba.contrast(dataCount, group1=dataCount$masks$'{0}', group2 =dataCount$masks$'{1}', name1='{0}', name2='{1}')\n".format(cond1, cond2)
	if reps==True:
		rscript += "dataAnalyze = dba.analyze(dataContrast, method=DBA_DESEQ)\n"
	else:
		rscript += "dataAnalyze = dba.analyze(dataContrast, bTagwise=FALSE, method=DBA_DESEQ)\n"
	#MA plot
	rscript += "pdf('{}_vs_{}_MA.pdf')\n".format(cond1, cond2)
	rscript += "dba.plotMA(dataAnalyze,method=DBA_DESEQ)\n"
	rscript += "dev.off()\n"
	#PCA plot
	rscript += "pdf('{}_vs_{}_PCA.pdf')\n".format(cond1, cond2)
	rscript += "dba.plotPCA(dataAnalyze,method=DBA_DESEQ)\n"
	rscript += "dev.off()\n"
	rscript += "data.db = dba.report(dataAnalyze, DataType = DBA_DATA_FRAME, method=DBA_DESEQ)\n"
	rscript += "write.table(file='{}_vs_{}_diffbind.tsv',data.db,sep='\\t',row.names=F,quote=F)\n".format(cond1, cond2)
	return rscript

#!/usr/bin/python

########################################################################
# 20 Oct 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import argparse
import subprocess
import sys, re, os

def paired_bowtie(fastq1, fastq2, name, index, outdir):
	sam1 = outdir + "/" + "tmp.sam"
	sam1_o = open(sam1, "wb")
	report = outdir+'/'+name+'_report.txt'
	report1_o = open(report, "wb")
	uniq = "bowtie -m 2 -v 1 --best --strata --seed 0 --sam {0} -1 {1} -2 {2}".format(index, fastq1, fastq2)
	p = subprocess.Popen(uniq.split(), stdout = sam1_o, stderr=report1_o)
	p.communicate()
	sam2 = outdir+"/"+name+".sam"
	grep_paired_unique(sam1, sam2)
	os.remove(sam1)

def single_bowtie(fastq, name, index, outdir):
	sam1 = outdir + "/" + "tmp.sam"
	sam1_o = open(sam1, "wb")
	report = outdir+'/'+name+'_report.txt'
	report1_o = open(report, "wb")
	uniq = "bowtie -m 2 -v 1 --best --strata --seed 0 --sam {0} {1}".format(index, fastq)
	p = subprocess.Popen(uniq.split(), stdout = sam1_o, stderr=report1_o)
	p.communicate()
	sam2 = outdir+"/"+name+".sam"
	grep_single_unique(sam1, sam2)
	os.remove(sam1)

def grep_paired_unique(samfile, outfile):
	output=  open(outfile, "w")
	with open(samfile) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if line.startswith("@"):
				output.write("{}\n".format(line)),
				continue
			if len(word) > 12:
				m = re.match("XS:i:", word[12])
				if not m:
					if int(word[1]) == 147 or int(word[1]) == 83 or int(word[1]) == 99 or int(word[1]) == 163 or int(word[1]) == 81 or int(word[1]) == 97 or int(word[1]) == 145 or int(word[1]) == 161:
						output.write("{}\n".format(line)),

def grep_single_unique(samfile, outfile):
	output=  open(outfile, "w")
	with open(samfile) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if line.startswith("@"):
				output.write("{}\n".format(line)),
				continue
			if len(word) > 12:
				m = re.match("XS:i:", word[12])
				if not m:
					if int(word[1]) == 0 or int(word[1]) == 16:
						output.write("{}\n".format(line)),


def paired_bowtie2(fastq1, fastq2, name, index, outdir, threads):
	report = outdir+'/'+name+'_report.txt'
	report1_o = open(report, "wb")
	uniq = "bowtie2 -p {4} -k 2 -N 1 --mm --no-mixed --no-discordant -x {0} -1 {1} -2 {2} -S {3}/tmp.sam".format(index, fastq1, fastq2, outdir, threads)
	p = subprocess.Popen(uniq.split(), stderr=report1_o)
	p.communicate()
	grep_paired_unique(outdir+"/tmp.sam", outdir+'/'+name+'.sam')
	os.remove(outdir+"/tmp.sam")

def single_bowtie2(fastq, name, index, outdir, threads):
	report = outdir+'/'+name+'_report.txt'
	report1_o = open(report, "wb")
	uniq = "bowtie2 -p {3} -k 2 -N 1 --mm  -x {0} -U {1} -S {2}/tmp.sam".format(index, fastq, outdir, threads)
	p = subprocess.Popen(uniq.split(), stderr=report1_o)
	p.communicate()
	grep_single_unique(outdir+"/tmp.sam", outdir+'/'+name+'.sam')
	os.remove(outdir+"/tmp.sam")

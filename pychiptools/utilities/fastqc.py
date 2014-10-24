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

def run_fastqc(fq1):
	devnull = open('/dev/null', 'w')
	command = ["fastqc", "{}".format(fq1)] #outdir must exist!
	subprocess.call(command,  stdout=devnull)

def find_adapters(fq):
	adapters = []
	idir = re.sub(".fastq","", fq)
	report = idir+"_fastqc/fastqc_data.txt"
	flist = open(report).readlines()
	parsing = False
	for line in flist:
		if line.startswith(">>Overrepresented sequences\tfail"):
			parsing = True
		elif line.startswith(">>END_MODULE"):
			parsing = False
		if parsing:
			if line.startswith(">>"):
				continue
			if line.startswith("#"):
				continue
			else:
				line = line.rstrip()
				word = line.split("\t")
				if word[3] != "No Hit":
					adapters.append(word[0])
	return adapters

def cut_adapters(adapters, fq1, outdir, rev_adapters=None, fq2=None):
	devnull = open('/dev/null', 'w')
	adapt1 = ""
	for i in adapters:
		adapters = "-a {} ".format(i)
		adapt1 = adapters+adapt1
	if rev_adapters:
		adapt2 = ""
		for i in rev_adapters:
			adapters = "-a {} ".format(i)
			adapt2 = adapters+adapt2

	if rev_adapters:
		command1 = "cutadapt -q 20 {0} --minimum-length=10 --paired-output {1}/tmp.2.fastq -o {1}/tmp.1.fastq {2} {3}".format(adapt1, outdir, fq1, fq2)
		p = subprocess.Popen(command1.split())
		p.communicate()
		command2 = "cutadapt -q 20 {0} --minimum-length=10 --paired-output {1}/trimmed_1.fastq -o {1}/trimmed_2.fastq {1}/tmp.2.fastq {1}/tmp.1.fastq".format(adapt2, outdir)
		p = subprocess.Popen(command2.split())
		p.communicate()
		cleanup = ["rm", "{0}/tmp.2.fastq".format(outdir), "{0}/tmp.1.fastq".format(outdir)]
		subprocess.call(cleanup, stdout=devnull)
	else:
		command1 = "cutadapt -q 20 {0} --minimum-length=10 -o {1}/trimmed.fastq {2}".format(adapt1, outdir, fq1)
		p = subprocess.Popen(command1.split())
		p.communicate()
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

def run_fastqc(fq1, outdir):
	command = "fastqc -q -o {0} {1}".format(outdir, fq1) #outdir must exist!
	subprocess.call(command.split())
	devnull = open('/dev/null', 'w')
	command = "fastqc -q {}".format(fq1) #outdir must exist!
	subprocess.call(command.split(),  stdout=devnull)
	
def find_adapters(outdir, fq1):
	adapters = []
	name = re.sub(".fastq", "", fq1)
	name = os.path.basename(name)
	command = "unzip -o -q {}/{}_fastqc.zip -d {}".format(outdir, name, outdir)
	subprocess.call(command.split())
	report = "{}/{}_fastqc/fastqc_data.txt".format(outdir, name)
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

def single_cut_adapters(adapters, fq1, outdir):
	trim = open('{}/trim_report.txt'.format(outdir), 'w')
	adapt1 = ""
	for i in adapters:
		adapters = "-a {} ".format(i)
		adapt1 = adapters+adapt1
	command1 = "cutadapt -q 20 {0} --minimum-length=10 -o {1}/trimmed.fastq {2}".format(adapt1, outdir, fq1)
	p = subprocess.Popen(command1.split(), stdout=trim)
	p.communicate()

def paired_cut_adapters(adapters, fq1, outdir, rev_adapters, fq2):
	devnull = open('/dev/null', 'w')
	trim = open('{}/trim_report.txt'.format(outdir), 'w')
	adapt1 = ""
	for i in adapters:
		adapters = "-a {} ".format(i)
		adapt1 = adapters+adapt1

	adapt2 = ""
	for i in rev_adapters:
		adapters = "-a {} ".format(i)
		adapt2 = adapters+adapt2

	command1 = "cutadapt -q 20 {0} --minimum-length=10 --paired-output {1}/tmp.2.fastq -o {1}/tmp.1.fastq {2} {3}".format(adapt1, outdir, fq1, fq2)
	p = subprocess.Popen(command1.split(), stdout=trim)
	p.communicate()
	command2 = "cutadapt -q 20 {0} --minimum-length=10 --paired-output {1}/trimmed_1.fastq -o {1}/trimmed_2.fastq {1}/tmp.2.fastq {1}/tmp.1.fastq".format(adapt2, outdir)
	p = subprocess.Popen(command2.split(), stdout=trim)
	p.communicate()
	cleanup = ["rm", "{0}/tmp.2.fastq".format(outdir), "{0}/tmp.1.fastq".format(outdir)]
	subprocess.call(cleanup, stdout=devnull)

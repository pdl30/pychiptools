#!/usr/bin/python

########################################################################
# 1 August 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import sys, os, re
import argparse

def run_macs2(sample, name, genome, outdir, pvalue, control, histone, qvalue=None):
	#In new version, could have multiple samples/controls??
	#Dealing with names could be tough for these samples
	samples = " ".join(sample)
	if histone==False:
		if not control:
			if pvalue:
				name = os.path.basename(name)
				name = "{}_p1e-{}".format(name, pvalue)
				command = "macs2 callpeak -t {} -g {} -n {} -f BED -p 1e-{} --nomodel --outdir {}".format(samples, genome, name, pvalue, outdir)
			else:
				name = os.path.basename(name)
				name = "{}_q1e-{}".format(name, qvalue)
				command = "macs2 callpeak -t {} -g {} -n {} -f BED -q {} --nomodel --outdir {}".format(samples, genome, name, qvalue, outdir)
		else:
			controls = " ".join(control)
			# WITH CONTROL
			if pvalue:
				name = os.path.basename(name)
				name = "{}_p1e-{}".format(name, pvalue)
				command = "macs2 callpeak -t {} -c {} -g {} -n {} -f BED -p 1e-{} --nomodel --outdir {}".format(samples, controls, genome, name, pvalue, outdir)
			else:
				name = os.path.basename(name)
				name = "{}_q1e-{}".format(name, qvalue)
				command = "macs2 callpeak -t {} -c {} -g {} -n {} -f BED -q {} --nomodel --outdir {}".format(samples, controls, genome, name, qvalue, outdir)
	else:
		if not control:
			if pvalue:
				name = os.path.basename(name)
				name = "{}_p1e-{}".format(name, pvalue)
				command = "macs2 callpeak -t {} -g {} -n {} -f BED -p 1e-{} --broad --nomodel --outdir {}".format(samples, genome, name, pvalue, outdir)
			else:
				name = os.path.basename(name)
				name = "{}_q1e-{}".format(name, qvalue)
				command = "macs2 callpeak -t {} -g {} -n {} -f BED -q {} --broad --nomodel --outdir {}".format(samples, genome, name, qvalue, outdir)
		else:
			controls = " ".join(control)
			# WITH CONTROL
			if pvalue:
				name = os.path.basename(name)
				name = "{}_p1e-{}".format(name, pvalue)
				command = "macs2 callpeak -t {} -c {} -g {} -n {} -f BED -p 1e-{} --broad --nomodel --outdir {}".format(samples, controls, genome, name, pvalue, outdir)
			else:
				name = os.path.basename(name)
				name = "{}_q1e-{}".format(name, qvalue)
				command = "macs2 callpeak -t {} -c {} -g {} -n {} -f BED -q {} --broad --nomodel --outdir {}".format(samples, controls, genome, name, qvalue, outdir)
	p = subprocess.Popen(command.split())
	p.communicate()


def convert_peaks_to_bed(name, outdir):
	width = 400
	output = open(outdir + '/' + name+'_tmp.bed', "w")
	with open(outdir + '/' + name+'_peaks.xls') as f:
		new_end = 0
		new_start = 0
		for line in f:
			line = line.rstrip()
			if line.startswith('#'):
				pass
			elif re.match('^chr\sstart', line):
				pass
			elif re.match('^\s', line):
				pass
			else:
				word=  line.split("\t")
				if len(word) > 9:
					chr1 = word[0]
					start = word[1]
					end = word[2]
					summit = int(word[4])
					height = word[5] # pileup
					new_start = summit - width/2;
					new_end = summit + width/2;
					if new_start <= 0:
						new_start = 0
					output.write("{}\t{}\t{}\n".format(chr1, new_start, new_end)),
					

def post_process_peaks_for_ucsc(name, chrom_sizes, outdir):
	command = ["bedClip", outdir + '/' + name+'_tmp.bed', chrom_sizes, outdir + '/' + name+'_400bp.bed']
	subprocess.call(command)
	command2 = ["bedToBigBed", outdir + '/' + name+"_400bp.bed", chrom_sizes, outdir + '/' + name+".bb"]
	subprocess.call(command2)

def post_process_histone_peaks(name, chrom_sizes, outdir):
	output = open(outdir + '/' + name+'_tmp.bed', "w")
	with open(outdir + '/' + name+'_peaks.gappedPeak') as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if int(word[4]) > 1000:
				score = 1000
			else:
				score = int(word[4])
			output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(word[0],word[1],word[2],word[3],
				score, word[5],word[6],word[7],word[8],word[9],word[10],word[11],word[12],word[13], word[14])),
	output.close()
	command2 = ["bedToBigBed", "type=bed12+3", outdir + '/' + name+"_tmp.bed", chrom_sizes, outdir + '/' + name+".bb"]
	subprocess.call(command2)
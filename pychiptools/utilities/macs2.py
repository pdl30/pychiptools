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

def run_macs2(sample, genome, pvalue, control, histone, qvalue=None):
	name = re.sub(".BED", "", sample, re.IGNORECASE)
	if histone==False:
		if not control:
			if pvalue:
				name = "{}_p1e-{}".format(name, pvalue)
				command = "macs2 callpeak -t {} -g {} -n {} -f BED -p 1e-{} --nomodel".format(sample, genome, name, pvalue)
			else:
				name = "{}_q1e-{}".format(name, qvalue)
				command = "macs2 callpeak -t {} -g {} -n {} -f BED -q {} --nomodel".format(sample, genome, name, qvalue)
		else:
			# WITH CONTROL
			if pvalue:
				name = "{}_p1e-{}".format(name, pvalue)
				command = "macs2 callpeak -t {} -c {} -g {} -n {} -f BED -p 1e-{} --nomodel".format(sample, control, genome, name, pvalue)
			else:
				name = "{}_q1e-{}".format(name, qvalue)
				command = "macs2 callpeak -t {} -c {} -g {} -n {} -f BED -q {} --nomodel".format(sample, control, genome, name, qvalue)
	else:
		if not control:
			if pvalue:
				name = "{}_p1e-{}".format(name, pvalue)
				command = "macs2 callpeak -t {} -g {} -n {} -f BED -p 1e-{} --broad --nomodel".format(sample, genome, name, pvalue)
			else:
				name = "{}_q1e-{}".format(name, qvalue)
				command = "macs2 callpeak -t {} -g {} -n {} -f BED -q {} --broad --nomodel".format(sample, genome, name, qvalue)
		else:
			# WITH CONTROL
			if pvalue:
				name = "{}_p1e-{}".format(name, pvalue)
				command = "macs2 callpeak -t {} -c {} -g {} -n {} -f BED -p 1e-{} --broad --nomodel".format(sample, control, genome, name, pvalue)
			else:
				name = "{}_q1e-{}".format(name, qvalue)
				command = "macs2 callpeak -t {} -c {} -g {} -n {} -f BED -q {} --broad --nomodel".format(sample, control, genome, name, qvalue)

	print command
	p = subprocess.Popen(command.split())
	p.communicate()


def convert_peaks_to_bed(name):
	print name 
	width = 400
	output = open(name+'_tmp.bed', "w")
	with open(name+'_peaks.xls') as f:
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
					

def post_process_peaks_for_ucsc(name, chrom_sizes):
	command = ["bedClip", name+'_tmp.bed', chrom_sizes, name+'_400bp.bed']
	subprocess.call(command)
	command2 = ["bedToBigBed", name+"_400bp.bed", chrom_sizes, name+".bb"]
	subprocess.call(command2)

def post_process_histone_peaks(name, chrom_sizes):
	output = open(name+'_tmp.bed', "w")
	with open(name+'_peaks.gappedPeak') as f:
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
	command2 = ["bedToBigBed", "type=bed12+3", name+"_tmp.bed", chrom_sizes, name+".bb"]
	subprocess.call(command2)
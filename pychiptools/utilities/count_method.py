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
import math

def read_gapdh_counts_file(ifile):
	with open(ifile) as f:
		header= next(f)
		header= header.rstrip()
		word = header.split("\t")
	return word[1]

def create_consensus_gtf_from_bed(peak1, peak2, outprefix):
	peak_data = {}
	command1 = "cat {} {} > {}_tmp1.bed".format(peak1, peak2, outprefix)
	command2 = "sortBed -i {}_tmp1.bed > {}_tmp2.bed".format(outprefix, outprefix)
	command3 = "mergeBed -i {}_tmp2.bed > {}.bed".format(outprefix, outprefix)
	subprocess.call(command1, shell=True)
	subprocess.call(command2, shell=True)
	subprocess.call(command3, shell=True)
	subprocess.call("rm {0}_tmp1.bed {0}_tmp2.bed".format(outprefix), shell=True)
	c = 0
	output = open(outprefix+".gtf", "w")
	with open(outprefix+".bed") as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			chromo = re.sub("chr", "", word[0])
			output.write("{}\tgene\texon\t{}\t{}\t.\t.\t.\tgene_id \"peak_{}\"; transcript_id \"peak_{}\";\n".format(chromo, word[1], word[2], c, c)),
			peak_data[c] = (chromo, word[1], word[2])
			c += 1
	output.close()
	return peak_data

def read_count_files(count):
	data = {}
	with open(count) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			#Get rid of the name from the GTF added in the previous function
			name = re.sub("peak_", "", word[0])
			data[name] = word[1]
	return data

def count(inv_conds, peaks, cond1, cond2, gapdh_counts):
	sample1 = inv_conds[cond1][0] 
	sample2 = inv_conds[cond2][0] 
	name =  "{}_vs_{}".format(cond1, cond2)
	#Combine, sort and merge peaks
	peak_data = create_consensus_gtf_from_bed(peaks[sample1], peaks[sample2], name)
	name1 = os.path.basename(sample1)
	name2 = os.path.basename(sample2)
	count1 = re.sub(".bam", ".count", name1)
	count2 = re.sub(".bam", ".count", name2)
	#command1 = "htseq-count -f bam -s no {} {} > {}".format(sample1, name+".gtf", count1)
	#command2 = "htseq-count -f bam -s no {} {} > {}".format(sample2, name+".gtf", count2)
	#subprocess.call(command1, shell=True)
	#subprocess.call(command2, shell=True)
	peak_count1 = read_count_files(count1)
	peak_count2 = read_count_files(count2)
	#GAPDH count normalisation
	if gapdh_counts:
		c1 = read_gapdh_counts_file(gapdh_counts[sample1])
		c2 = read_gapdh_counts_file(gapdh_counts[sample2])
		c1 = float(c1)/1000
		c2 = float(c2)/1000
		#Add 1 to values to prevent errors from 0 division
		for key in peak_count1:
			new_value = int(peak_count1[key])*float(c1)
			peak_count1[key] = new_value +1
		for key in peak_count2:
			new_value = int(peak_count2[key])*float(c2)
			peak_count2[key] = new_value+1
	result = {}
	output = open(name+".txt", "w")
	output.write("Chromosome\tStart\tEnd\tPeak ID\t{}\t{}\tLFC\n".format(cond1, cond2)),
	#Log fold change calculation
	for key in peak_count1:
		if key.startswith("_"):
			pass
		elif float(peak_count1[key]) > 200 or float(peak_count2[key]) > 200:
			fc = float(peak_count1[key])/float(peak_count2[key])
			lfc = math.log(fc, 2)
			if lfc > 1 or lfc < -1:
				output.write("chr{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(peak_data[int(key)][0], peak_data[int(key)][1],peak_data[int(key)][2],key, peak_count1[key], peak_count2[key], lfc)),
	output.close()
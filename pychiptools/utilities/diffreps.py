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

def diffReps(inv_conds, cond1, cond2, chromsizes, method="nb", controls=None, counts=None):
	#Can't deal with reps yet
	if counts:
		sample1 = inv_conds[cond1][0]
		sample2 = inv_conds[cond2][0]
		c1 = read_gapdh_counts_file(counts[sample1])
		c2 = read_gapdh_counts_file(counts[sample2])
		c1 = float(c1)/1000
		c2 = float(c2)/1000
		output = open("norm.txt", "w")
		output.write("treatment\t{}\n".format(c1)),
		output.write("control\t{}\n".format(c2)),
		output.close()
		command = "diffReps.pl --treatment {0} --control {1} --report {2}_vs_{3}_diffReps.txt --chrlen {4} -me {5} --norm norm.txt --nproc 8".format(inv_conds[cond1][0], inv_conds[cond2][0], 
			cond1, cond2, chromsizes, method )
	elif controls:
		backt1 = []
		backt2 = []
		for sample in t1:
			backt1.append(controls[sample])
		for sample in t2:
			backt2.append(controls[sample])
		command = "diffReps.pl --treatment {0} --control {1} --btr {2} --bco {3} --report {4}_vs_{5}_diffReps.txt --chrlen {6} -me {7}".format(inv_conds[cond1], inv_conds[cond2], 
			backt1, backt2, cond1, cond2, chromsizes, method)
	else:
		command = "diffReps.pl --treatment {0} --control {1} --report {2}_vs_{3}_diffReps.txt --chrlen {4} -me {5}".format(inv_conds[cond1][0], inv_conds[cond2][0], cond1, cond2, 
			chromsizes, method)
	print command
	subprocess.call(command.split())

def read_gapdh_counts_file(ifile):
	with open(ifile) as f:
		header= next(f)
		header= header.rstrip()
		word = header.split("\t")
	return word[1]

def annotate_ensembl(dict_obj):
	ens = importr("biomaRt")
	ensembl = ro.r.useMart("ensembl")
	genome="mmusculus_gene_ensembl"
	ensembl = ro.r.useDataset(genome, mart=ensembl)
	values = []
	for key1 in dict_obj.keys():
		values.append(key1)
	C1BM = ro.r.getBM(attributes=StrVector(["ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_id", "description", "gene_biotype"]), 
		filters="ensembl_gene_id", values=values, mart=ensembl)
	gene = list(C1BM.rx(True,1))
	chr1 = list(C1BM.rx(True,2))
	tss = list(C1BM.rx(True,3))
	end = list(C1BM.rx(True,4))
	st = list(C1BM.rx(True,5))
	name = list(C1BM.rx(True,6))
	des = list(C1BM.rx(True,7))
	bio = list(C1BM.rx(True,8))
	data = {}
	for index, g in enumerate(gene):
		data[g] = (chr1[index], tss[index], end[index], st[index], name[index], des[index], bio[index])
	return data

def postprocess_diffreps(cond1, cond2, gtf):
	resbed = "{}_vs_{}_diffReps.txt".format(cond1, cond2)
	output = open("tmp1.bed", "w")
	with open(resbed) as f:
		for line in f:
			if line.startswith("#") or line.startswith("Chrom"):
				pass
			else:
				line = line.rstrip()
				word = line.split("\t")
				if float(word[12]) < 1e-7:
			#		name = re.sub("chr", "", word[0])
					output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(word[0], word[1], word[2],word[6], word[7], word[10], word[11], word[12], word[13])),
	output.close()
	command = "windowBed -w 50000 -a {} -b {} > {}".format("tmp1.bed", gtf, "tmp2.bed")
	subprocess.call(command, shell=True)
	data = defaultdict(list)
	ens = {}
	with open("tmp2.bed") as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			data[(word[0], word[1], word[2])].append(line)
			ens[word[11]] = 1
	print ens
	anno = annotate_ensembl(ens)
	output1 = open("{}_vs_{}_diffReps_result.txt".format(cond1, cond2), "w")
	output1.write("Chrom\tStart\tEnd\tPvalue\tLFC\tNearest Genes\n")
	for key in data:
		output1.write("{}\t{}\t{}\t{}\t{}\t".format(word[0], word[1], word[2], data[key][0][7], data[key][0][6])),
		for v in data[key]:
			output1.write("{},".format(anno[data[key][0][11]][4])),
		output1.write("\n"),
	output1.close()

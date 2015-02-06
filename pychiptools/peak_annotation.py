#!/usr/bin/python

########################################################################
# 10 July 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import sys, re, os
import csv
from collections import defaultdict
import subprocess
from operator import itemgetter
import HTSeq
import argparse
import tempfile
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import IntVector, FloatVector, StrVector

def read_annotation(gtf):
	gene_info = {}
	chrom = {}
	gtffile = HTSeq.GFF_Reader( gtf )
	for feature in gtffile:
		#feature.name is gene name, useful for delving into GFF files again!
		if feature.type == "exon" and feature.attr["exon_number"] == "1":
			if feature.name not in gene_info:
				chrom[feature.name] = feature.iv.chrom
				gene_info[feature.name] = []
				gene_info[feature.name].append(feature.iv.start)
			else:
				gene_info[feature.name].append(feature.iv.start)
	return gene_info, chrom

def find_key(input_dict, value):
    return {k for k, v in input_dict.items() if v == value}

def extract(d, keys):
	return dict((k, d[k]) for k in keys if k in d)

def f7(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if not (x in seen or seen_add(x))]

def closest_tss(peak_data, starts, chrom, outname):
	tss_anno = {}
	for peak in sorted(peak_data):
		mid = int(peak_data[peak][2]) - int(peak_data[peak][1])
		mid = mid/2
		middle = int(peak_data[peak][1])+int(mid)
		min_dis = 1000000000
		status = ""
		for gene in starts:
			#print gene, starts[gene]
			if chrom[gene] == peak_data[peak][0]:
				unique_starts = f7(starts[gene])
				for start in unique_starts:
					distance = abs(starts[gene] - middle)
					if distance < min_dis:
						if distance < 0:
							status= "upstream"
						else:
							status = "downstream"
						tss_anno[key] = gene, distance, status
	output = open(outname, "w")
	for key in tss_anno:
		output.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(peak_data[key][0], peak_data[key][1], peak_data[key][2], tss_anno[key][0],tss_anno[key][1],tss_anno[key][2])),

def homer_annotation(peak, genome, output, gtf=None):
	if gtf:
		command = "annotatePeaks.pl {} {} -gtf {} > {}".format(peak, genome, gtf, output)
	else:
		command = "annotatePeaks.pl {} {} > {}".format(peak, genome, output)
	subprocess.call(command, shell=True)

def closestbed(peak, gtf):
	f = tempfile.NamedTemporaryFile(delete=False)
	f.close()
	command = "closestBed -D ref -a {} -b {} > {}".format(peak, gtf, f.name)
	subprocess.call(command, shell=True)
	return f.name

def parse_closest(gtf, outname):
	output = open(outname, "w")
	transcripts = []
	with open(gtf) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			info = word[11].split(";")
			trans_id = info[0].strip("^gene_id ")
			trans_id = trans_id.strip("\"")
			if trans_id not in transcripts:
				transcripts.append(trans_id)
				output.write("{}\t{}\t{}\t{}\t{}\n".format(word[0], word[1], word[2], trans_id, word[12])),

def read_peak_info(peak_file, ens):
	peak_data = {}
	count = 0
	with open(peak_file) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if ens:
				chrom = word[0].strip("chr")
				if chrom == "M":
					chrom = "MT"
				peak_data[count]=(chrom, word[1], word[2]) #For example should use pandas dataframe instead of dict here!
			else:
				peak_data[count]=(word[0], word[1], word[2])
			count +=1
	f = tempfile.NamedTemporaryFile(delete=False)
	for key in peak_data:
		f.write("{}\t{}\t{}\n".format(peak_data[key][0], peak_data[key][1],peak_data[key][2])),
	f.close()
	return peak_data, f.name

def main():
	parser = argparse.ArgumentParser(description='Annotation of Peaks.\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")
	custom_parser = subparsers.add_parser('custom', help="Runs custom Annotation")
	homer_parser = subparsers.add_parser('homer', help="Runs annotatePeaks.pl")

	custom_parser.add_argument('-p', '--peak', help='Peak file in bed format', required=True)
	custom_parser.add_argument('-g', '--gtf', help='GTF file', required=True)
	custom_parser.add_argument('-c', action='store_true', help='Find closest gene instead of closest TSS', required=False)
	custom_parser.add_argument('-e', action='store_true', help='If GTF supplied and ensembl, convert peaks to correct format', required=False)
	custom_parser.add_argument('-o', '--output', help='Output name.', required=True)
	
	homer_parser.add_argument('-p', '--peak', help='Peak file', required=True)
	homer_parser.add_argument('-g', '--genome', help='Genome, options are mm10/hg19', required=True)
	homer_parser.add_argument('-a', '--gtf', help='Optional GTF file', required=False)
	homer_parser.add_argument('-e', action='store_true', help='If GTF supplied and ensembl, convert peaks to correct format', required=False)
	homer_parser.add_argument('-o', '--out', help='Output name.', required=True)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())

	path1 = "/home/patrick/Reference_Genomes/pyngspipe_references/"

	if args["subparser_name"] == "custom":
		peak_data, peak_file = read_peak_info(args["peak"], args["e"])
		if args["c"]:
			tmpgtf = closestbed(peak_file, args["gtf"])
			parse_closest(tmpgtf, args["output"])
		else:
			starts, chrom = read_annotation(args["gtf"])
			tss_anno = closest_tss(peak_data, starts, chrom, args["output"])
		os.remove(peak_file)
	elif args["subparser_name"] == "homer":
		peak_data, peak_file = read_peak_info(args["peak"], args["e"])
		homer_annotation(peak_file, args["genome"], args["out"], args["gtf"])

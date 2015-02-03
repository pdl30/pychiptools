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

def read_peak_info(peak_file):
	peak_data = {}
	count = 0
	with open(peak_file) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			peak_data[count]=(word[0], word[1], word[2])
			count +=1
	return peak_data

def closest_tss(peak_data, starts, chrom):
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
	return tss_anno

#Find closest gene:
def closest_gene(peak_data, chromo, start, end):
	gene_anno = {}
	for key in peak_data:
		min_dis2 = 1000000000
		min_gene2 = ""
		start_dis = {}
		end_dis = {}	
		mid = int(peak_data[key][2]) - int(peak_data[key][1])
		middle = int(peak_data[key][1])+mid
		res = find_key(chromo, peak_data[key][0])
		res_starts = extract(start, res)
		res_ends = extract(end, res)
		for key2 in res_starts:
			distances = (abs(v - middle) for v in res_starts[key2])
			start_dis[key2] = min(distances)
			distances2 = (abs(v - middle) for v in res_ends[key2])
			end_dis[key2] = min(distances2)
		for dis in start_dis:
			if start_dis[dis] < min_dis2:
				min_gene2 = dis
				min_dis2 = start_dis[dis]
		for dis in end_dis:
			if end_dis[dis] < min_dis2:
				min_gene2 = dis
				min_dis2 = end_dis[dis]
		if res_starts[min_gene2][0]- middle <= 0:
			status= "upstream"
		else:
			status = "downstream"
		gene_anno[key] = min_dis2, min_gene2, status
	return gene_anno

def homer_annotation(peak, genome, output, gtf=None):
	if gtf:
		command = "annotatePeaks.pl {} {} -gtf {} > {}".format(peak, genome, gtf, output)
	else:
		command = "annotatePeaks.pl {} {} > {}".format(peak, genome, output)
	subprocess.call(command, shell=True)

def convert_to_ens(peak_file):
	name = peak_file.strip(".bed")
	output = open(name + "_tmp.bed", "w")
	with open(peak_file) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			chrom = word[0].strip("chr")
			if chrom == "M":
				chrom = "MT"
			output.write("{}\t{}\t{}\n".format(chrom, word[1], word[2])),
	return name

def annotate_ensembl(dict_obj, genome, transcript=False):
	ens = importr("biomaRt")
	ensembl = ro.r.useMart("ensembl")
	if genome == "mm10":
		genome="mmusculus_gene_ensembl"
	elif genome == "hg19":
		genome="hsapiens_gene_ensembl"
	ensembl = ro.r.useDataset(genome, mart=ensembl)
	values = []
	for key1 in dict_obj:
		values.append(key1)
	if transcript:
		C1BM = ro.r.getBM(attributes=StrVector(["ensembl_transcript_id", "external_gene_name", "description", "gene_biotype"]), filters="ensembl_transcript_id", 
			values=values, mart=ensembl)
	else:
		C1BM = ro.r.getBM(attributes=StrVector(["ensembl_gene_id", "external_gene_name", "description", "gene_biotype"]), filters="ensembl_gene_id", values=values, mart=ensembl)
	gene = list(C1BM.rx(True,1))
	name = list(C1BM.rx(True,2))
	des = list(C1BM.rx(True,3))
	bio = list(C1BM.rx(True,4))
	data = {}
	for index, g in enumerate(gene):
		data[g] = (name[index], des[index], bio[index])
	return data

def get_paths(path1, genome):
	if genome == "hg19":
		gtf = path1 + "hg19/hg19.gtf"
	elif genome == "mm10":
		gtf = path1 + "mm10/mm10.gtf"
	return gtf

def closestbed(peak, gtf):
	f = tempfile.NamedTemporaryFile(delete=False)
	f.close()
	command = "closestBed -D ref -a {} -b {} > {}".format(peak, gtf, f.name)
	subprocess.call(command, shell=True)
	os.remove(peak)
	return f.name

def parse_closest(gtf, genome, outname):
	transcripts = []
	output = open(outname, "w")
	print gtf
	with open(gtf) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			info = word[11].split(";")
			trans_id = info[0].strip("^gene_id ")
			trans_id = trans_id.strip("\"")
			transcripts.append(trans_id)
	anno = annotate_ensembl(transcripts, genome, transcript=True)
	with open(gtf) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			info = word[11].split(";")
			trans_id = info[0].strip("^gene_id ")
			trans_id = trans_id.strip("\"")
			if trans_id in anno:
				output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(word[0], word[1], word[2], anno[trans_id][0], anno[trans_id][1], anno[trans_id][2], word[12])),
			else:
				output.write("{}\t{}\t{}\n".format(word[0], word[1], word[2])),

def make_tmpbed(peak_data):
	f = tempfile.NamedTemporaryFile(delete=False)
	for key in peak_data:
		f.write("{}\t{}\t{}\n".format(peak_data[key][0], peak_data[key][1],peak_data[key][2])),
	f.close()
	return f.name #Remember to delete this!

def main():
	parser = argparse.ArgumentParser(description='Annotation of Peaks.\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")
	custom_parser = subparsers.add_parser('custom', help="Runs custom Annotation")
	homer_parser = subparsers.add_parser('homer', help="Runs annotatePeaks.pl")
	#chip_parser = subparsers.add_parser('chipanno', help="Runs chippeakanno, no longer recommended")
	#OK new strategy, use closestBed and custom!
	custom_parser.add_argument('-p', '--peak', help='Peak file in bed format', required=True)
	custom_parser.add_argument('-g', '--genome', help='genome to use, options are hg19/mm10', required=True)
	custom_parser.add_argument('-c', action='store_true', help='Find closest gene instead of closest TSS', required=False)
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
	gtf = get_paths(path1, args["genome"])

	if args["subparser_name"] == "custom":
		if args["c"]:
			peak_data = read_peak_info(args["peak"])
			tmpbed = make_tmpbed(peak_data)
			tmpgtf = closestbed(tmpbed, gtf)
			parse_closest(tmpgtf, args["genome"], args["output"])
		else:
			peak_data = read_peak_info(args["peak"])
			starts, chrom = read_annotation(args["gtf"])
			tss_anno = closest_tss(peak_data, starts, chrom)
			for key in tss_anno:
				print "{}\t{}\t{}\t{}\t{}\t{}\n".format(peak_data[key][0], peak_data[key][1], peak_data[key][2], tss_anno[key][0],tss_anno[key][1],tss_anno[key][2] ),
	elif args["subparser_name"] == "homer":
		if args["ens"]:
			name = convert_to_ens(args["peak"])
			homer_annotation(name + "_tmp.bed", args["genome"], args["out"], args["gtf"])
		else:
			homer_annotation(args["peak"], args["genome"], args["out"], args["gtf"])

main()
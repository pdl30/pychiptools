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

#Use GTF for this
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

def read_peak_info(peak_file, ens):
	#Since 
	peak_data = {}
	count = 0
	with open(peak_file) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if ens:
				chrom = word[0].strip("chr")
				peak_data[count]=(chrom, word[1], word[2])
			else:
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
	print command
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

def run_rcode(rscript, name):
	rcode = open(name, "w")
	rcode.write(rscript)
	rcode.close()
	try:
		subprocess.call(['Rscript', name])
	except:
		error("Error in running {}\n".format(name))
		error("Error: %s\n" % str(sys.exc_info()[1]))
		error( "[Exception type: %s, raised in %s:%d]\n" % ( sys.exc_info()[1].__class__.__name__, 
		os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]), 
		traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
		sys.exit(1)

def chippeakanno(peaks, output_prefix):
	rscript = "library('ChIPpeakAnno')\n"
	rscript += "library(biomaRt)\n"
	rscript += "ensembl = useMart('ensembl')\n"
	rscript += "ensembl = useDataset('mmusculus_gene_ensembl', mart=ensembl)\n"
	rscript += "TSS.mouse.NCBIM38 = getAnnotation(mart=ensembl, featureType='TSS')\n"
	rscript += "peaks <- BED2RangedData('{}')\n".format(peaks)
	rscript += "annotatedPeaktss = annotatePeakInBatch(peaks, AnnotationData=TSS.mouse.NCBIM38, PeakLocForDistance='middle')\n"
	rscript += "annotatedPeaktrans = annotatePeakInBatch(peaks, AnnotationData=TSS.mouse.NCBIM38, FeatureLocForDistance='middle', PeakLocForDistance='middle', output='shortestDistance')\n"
	rscript += "annotatedPeaktss <- as.data.frame(annotatedPeaktss)\n"
	rscript += "annotatedPeaktrans <- as.data.frame(annotatedPeaktrans)\n"
	rscript += "C1BM <- getBM(c('ensembl_gene_id', 'external_gene_name','description', 'gene_biotype'),filters = 'ensembl_gene_id', values = annotatedPeaktrans$feature, mart = ensembl)\n"
	rscript += "C2BM <- getBM(c('ensembl_gene_id', 'external_gene_name','description', 'gene_biotype'),filters = 'ensembl_gene_id', values = annotatedPeaktss$feature, mart = ensembl)\n"
	rscript += "annotrans <- cbind(annotatedPeaktrans, C1BM[match(annotatedPeaktrans$feature, C1BM[,1]), 2:4])\n"
	rscript += "annotss <- cbind(annotatedPeaktss, C2BM[match(annotatedPeaktss$feature, C2BM[,1]), 2:4])\n"
	rscript += "write.table(annotss, file='{}_nearest_tss.tsv', sep='\\t', quote=F, row.names=F)\n".format(output_prefix)
	rscript += "write.table(annotrans, file='{}_nearest_gene.tsv', sep='\\t', quote=F, row.names=F)\n".format(output_prefix)
	return rscript

def annotate_ensembl(dict_obj):
	ens = importr("biomaRt")
	ensembl = ro.r.useMart("ensembl")
	genome="mmusculus_gene_ensembl"
	ensembl = ro.r.useDataset(genome, mart=ensembl)
	values = []
	for key1 in dict_obj.keys():
		values.append(key1)
	C1BM = ro.r.getBM(attributes=StrVector(["external_gene_id", "description", "gene_biotype"]), filters="ensembl_gene_id", values=values, mart=ensembl)
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

def main():
	parser = argparse.ArgumentParser(description='Annotation of Peaks.\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")
	custom_parser = subparsers.add_parser('custom', help="Runs custom Annotation")
	homer_parser = subparsers.add_parser('homer', help="Runs annotatePeaks.pl")
	chip_parser = subparsers.add_parser('chipanno', help="Runs chippeakanno")
	custom_parser.add_argument('-p', '--peak', help='Peak file in ensembl format', required=False)
	custom_parser.add_argument('-g', '--gtf', help='GTF file', required=False)
	custom_parser.add_argument('-e', help='If using ensembl GTF, convert peaks to Ensembl Format', action="store_true", required=False)
	custom_parser.add_argument('-o', '--output', help='Output name.', required=False)
	
	homer_parser.add_argument('-p', '--peak', help='Peak file', required=True)
	homer_parser.add_argument('-g', '--genome', help='Genome, options are mm10/hg19', required=True)
	homer_parser.add_argument('-a', '--gtf', help='Optional GTF file', required=False)
	homer_parser.add_argument('-e', '--ens', action='store_true', help='If GTF supplied and ensembl, convert peaks to correct format', required=False)
	homer_parser.add_argument('-o', '--out', help='Output name.', required=False)

	chip_parser.add_argument('-p', '--peak', help='Peak file', required=True)
	chip_parser.add_argument('-o', '--output', help='Output prefix', required=True)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	if args["subparser_name"] == "custom":
		starts, chrom = read_annotation(args["gtf"])
		peak_data = read_peak_info(args["peak"], args["e"])
		tss_anno = closest_tss(peak_data, starts, chrom)
		for key in tss_anno:
			print "{}\t{}\t{}\t{}\t{}\t{}\n".format(peak_data[key][0], peak_data[key][1], peak_data[key][2], tss_anno[key][0],tss_anno[key][1],tss_anno[key][2] ),
	elif args["subparser_name"] == "homer":
		if args["ens"]:
			name = convert_to_ens(args["peak"])
			homer_annotation(name + "_tmp.bed", args["genome"], args["out"], args["gtf"])
		else:
			homer_annotation(args["peak"], args["genome"], args["out"], args["gtf"])
	elif args["subparser_name"] == "chipanno":
		rscript = chippeakanno(args["peak"], args["output"])
		run_rcode(rscript, "chipseqanno_rcode.R")
		
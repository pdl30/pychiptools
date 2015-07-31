#!/usr/bin/python

########################################################################
# 19 May 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import os
import sys
import itertools
import subprocess
import re
import argparse
from multiprocessing import Pool
import pychiptools
import pkg_resources
from pychiptools.utilities import macs2

#Sicer requires UCSC format peaks in BED6 format!!! and sample and control BED files must be in one directory!
def run_sicer(sample, genome, qvalue, control, outdir, window_size=200, fragsize=150, gapsize=400):
	file_dir = os.path.dirname(os.path.realpath(sample))
	root,ext = sample.split(".")
	if ext == "BED":
		subprocess.call(["mv", sample, root+".bed"])
		sample = root+".bed"
	if control:
		root,ext = control.split(".")
		if ext == "BED":
			subprocess.call(["mv", control, root+".bed"])
			control = root+".bed"
	if control:
		command = "sh /home/patrick/Programs/SICER_V1.1/SICER/SICER.sh {0} {1} {2} {8} {3} 1 {4} {5} 0.7 {6} {7}".format(file_dir, sample, control, genome, window_size, fragsize, gapsize, qvalue, outdir)
	else:
		command = "sh /home/patrick/Programs/SICER_V1.1/SICER/SICER-rb.sh {0} {1} {7} {2} 1 {3} {4} 0.7 {5} {6}".format(file_dir, sample, genome, window_size, fragsize, gapsize, qvalue, outdir)
	#Output is name with -W200-G400-FDR0.01-island.bed in BED4 format. Convert this to BED3 and run to bigbed
	subprocess.call(command.split())

def peak_ranger_ccat(sample, control, outname, qvalue, threads):
	command = "peakranger ccat -d {0} -c {1} -o {2} --format bed -t {4} -q {3}".format(sample, control, outname, qvalue, threads)

def post_process_peakranger(name, chrom_sizes):
	outname = re.sub("_peakout", "", name)
	command = "bedClip {} {} {}".format(name+"_region.bed", chrom_sizes, outname+"_tmp.bed")
	subprocess.call(command.split())
	command = "bedSort {} {}".format(outname+"_tmp.bed", outname+"_tmp.bed")
	subprocess.call(command.split())
	output = open(outname+"_peakranger.bed", "w")
	with open(outname+"_tmp.bed") as f:
		for line in f:
			line=  line.rstrip()
			word = line.split("\t")
			output.write("{}\t{}\t{}\n".format(word[0], word[1], word[2])),
	output.close()
	command = "bedToBigBed {} {} {}".format(outname+"_peakranger.bed", chrom_sizes, outname+"_peakranger.bb")
	subprocess.call(command.split())

def function1(args):
	return macs2.run_macs2(*args)

def function2(args):
	return macs2.convert_peaks_to_bed(*args)

def function3(args):
	return macs2.post_process_peaks_for_ucsc(*args)

def main():
	parser = argparse.ArgumentParser(description='Programs included are MACS2 and SICER\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")
	
	macs2_parser = subparsers.add_parser('macs2', help="Runs MACS2 on ChIPseq samples. If no pvalue specified, will run 1e-3,4,5,6,7,9,12,15 pvalues concurrently!")
	macs2_parser.add_argument('-s','--sample', help='Sample BED file', required=True, nargs="+")
	macs2_parser.add_argument('-c','--control', help='Control BED file - Optional', required=False, nargs="+")
	macs2_parser.add_argument('-g','--genome', help='Genome aligned to, options include mm10/mm9/hg19', required=True)
	macs2_parser.add_argument('-p','--pvalue', help='Pvalue to use, optional. Must be integer! e.g. specifying 5 will use 1e-5', required=False, type=int)
	macs2_parser.add_argument('-q','--qvalue', help='Qvalue to use, optional. Must be float! e.g. 0.05', required=False, type=float)
	macs2_parser.add_argument('-n','--name', help='Output name, just the name and not the directory!', required=True)
	macs2_parser.add_argument('-b', action='store_true', help='Are samples histones', required=False)
	macs2_parser.add_argument('-o','--outdir', help='Output file directory, default is current directory', required=False)
	
	sicer_parser = subparsers.add_parser('sicer', help="Runs SICER")
	sicer_parser.add_argument('-s','--sample', help='Sample BED file', required=True)
	sicer_parser.add_argument('-c','--control', help='Control BED file - Optional', required=False)
	sicer_parser.add_argument('-g','--genome', help='Genome aligned to, options include mm10/mm9/hg19', required=True)
	sicer_parser.add_argument('-q','--qvalue', help='Qvalue to use, optional. Must be float! e.g. 0.05', required=False, type=float)
	sicer_parser.add_argument('-w', '--window_size', help='SICER window_size, default=200', default=200, required=False)
	sicer_parser.add_argument('-f', '--fragsize', help='SICER fragsize, default=150', default=150, required=False)
	sicer_parser.add_argument('-gS', '--gapsize', help='SICER gapsize, default=400', default=400, required=False)
	sicer_parser.add_argument('-o','--outdir', help='Output file directory, default is current directory', required=False)

	peakrang_parser = subparsers.add_parser('peakranger', help="Runs PeakRanger CCAT")
	peakrang_parser.add_argument('-s','--sample', help='Sample BED file', required=True, nargs="+")
	peakrang_parser.add_argument('-c','--control', help='Control BED file - Optional', required=False, nargs="+")
	peakrang_parser.add_argument('-g','--genome', help='Genome aligned to, options include mm10/mm9/hg19', required=True)
	peakrang_parser.add_argument('-q','--qvalue', help='Qvalue to use. Must be float! e.g. 0.05', required=False, type=float)
	peakrang_parser.add_argument('-t' '--threads', help='threads, default=1', default=1, required=False)
	peakrang_parser.add_argument('-o','--outdir', help='Output file directory, default is current directory', required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	#Needs revision!!
#	head, tail = os.path.split(args["sample"])
#	name = re.sub(".BED", "", args["sample"], re.IGNORECASE)

	if args["genome"] == "mm10":
		genome = "mm"
	elif args["genome"] == "mm9":
		genome = "mm"
	elif args["genome"] == "hg19":
		genome = "hs"
	else:
		raise Exception("Unsupported Genome!")

	if args["outdir"]:
		outdir = args["outdir"]
	else:
		outdir = os.getcwd()

	#Must add support for converting histone peaks to bigBed
	chrom = pkg_resources.resource_filename('pychiptools', 'data/{}.chrom.sizes'.format(args["genome"]))

	if args["subparser_name"]=="macs2":
		if args["pvalue"]:
			name = os.path.basename(args["name"])
			name = "{}_p1e-{}".format(name, args["pvalue"])
			macs2.run_macs2(args["sample"], args["name"], genome, outdir, args["pvalue"], args["control"], args["b"], chrom)

		elif args["qvalue"]:
			name = os.path.basename(args["name"])
			name = "{}_q-{}".format(name, args["qvalue"])
			macs2.run_macs2(args["sample"], args["name"], genome, outdir, None, args["control"], args["b"], chrom, qvalue=args["qvalue"])
		else:
			#Run all Pvalues
			pvalues = [3,4,5,6,7,9,12,15]
			name_list = []
			for p in pvalues:
				name = os.path.basename(args["name"])
				new_name = "{}_p1e-{}".format(name, p)
				name_list.append(new_name) #Need to sort this!

			pool = Pool(8)
			pool.map(function1, itertools.izip(itertools.repeat(args["sample"]), itertools.repeat(args["name"]), itertools.repeat(genome), 
				itertools.repeat(outdir), pvalues, itertools.repeat(args["control"]), itertools.repeat(args["b"], itertools.repeat(chrom))))

			pool2 = Pool(8)
			pool2.map(function2, itertools.izip(name_list, itertools.repeat(outdir)))
			pool3 = Pool(8)
			pool3.map(function3, itertools.izip(name_list, itertools.repeat(chrom), itertools.repeat(outdir)))

	elif args["subparser_name"]=="sicer":
		run_sicer(args["sample"], args["genome"], args["qvalue"], args["control"], outdir, args["window_size"], args["fragsize"], args["gapsize"])
	elif args["subparser_name"]=="peakranger":
		outname = re.sub(".BED", "_peakout", args["sample"])
		peak_ranger_ccat(args["sample"], args["control"], outname, args["qvalue"], args["threads"])
		post_process_peakranger(outname, chrom)
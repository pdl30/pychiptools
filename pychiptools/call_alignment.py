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
import gzip
from pychiptools.utilities import fastqc, alignment

def zip_trimmed(adapt, outdir, paired):
	#Gzip trimmed files and then remove the originals. Must test!
	if adapt == True:
		if paired:
			f_in = open("{}/trimmed_1.fastq".format(outdir), "rb")
			f_out = gzip.open('{}/trimmed_1.fq.gz'.format(outdir), 'wb')
			f_out.writelines(f_in)
			f_out.close()
			f_in.close()
			f_in = open("{}/trimmed_2.fastq".format(outdir), "rb")
			f_out = gzip.open('{}/trimmed_2.fq.gz'.format(outdir), 'wb')
			f_out.writelines(f_in)
			f_out.close()
			f_in.close()
			os.remove("{}/trimmed_1.fastq".format(outdir))
			os.remove("{}/trimmed_2.fastq".format(outdir))
		else:
			f_in = open("{}/trimmed.fastq".format(outdir), "rb")
			f_out = gzip.open('{}/trimmed.fq.gz'.format(outdir), 'wb')
			f_out.writelines(f_in)
			f_out.close()
			f_in.close()
			os.remove("{}/trimmed.fastq".format(outdir))
			
def main():
	parser = argparse.ArgumentParser(description='ChIP-seq fastqc, trimmer and bowtie wrapper\n ')
	parser.add_argument('-f','--fastq', help='Single end fastq', required=False)
	parser.add_argument('-p','--pair', help='Paired end fastqs. Please put them in order!', required=False, nargs='+')
	parser.add_argument('-i','--index', help='Path to bowtie index', required=True)
	parser.add_argument('-v','--version', help='Which bowtie to use, options are 1/2. default=2', default="2")
	parser.add_argument('-t','--threads', help='For bowtie2, how many threads to use (i.e. -p option on bowtie2', default=1, required=False)
	parser.add_argument('-n','--samname', help='Name of output sam file', required=True)
	parser.add_argument('-o','--outdir', help='Name of results directory', required=True)
	args = vars(parser.parse_args())
	path = os.getcwd()
	if os.path.isdir(args["outdir"]):
		print("Results directory already exists!")
	else:
		subprocess.call(["mkdir", args["outdir"]])
	if args["pair"]:
		adapt = False
		fq1 = args["pair"][0]
		fq2 = args["pair"][1]
		print "==> Running FastQC...\n"
		fastqc.run_fastqc(fq1)
		fastqc.run_fastqc(fq2)
		fwd_adapt = fastqc.find_adapters(fq1)
		rev_adapt = fastqc.find_adapters(fq2)
		if fwd_adapt or rev_adapt:
			adapt = True
			print "==> Removing adapters...\n"
			fastqc.paired_cut_adapters(fwd_adapt, fq1, args["outdir"], rev_adapt, fq2)
			fq1 = args["outdir"]+"/trimmed_1.fastq" 
			fq2 = args["outdir"]+"/trimmed_2.fastq"
		print "==> Aligning fastq's...\n"
		if args["version"] == "1":
			alignment.paired_bowtie(fq1, fq2, args["samname"], args["index"], args["outdir"])
		elif args["version"] == "2":
			alignment.paired_bowtie2(fq1, fq2, args["samname"], args["index"], args["outdir"], args["threads"])
		zip_trimmed(adapt, args["outdir"], True)
	elif args["fastq"]:
		adapt = False
		fq1 = args["fastq"]
		print "==> Running FastQC...\n"
		fastqc.run_fastqc(fq1)
		adapt = fastqc.find_adapters(fq1)
		if adapt:
			adapt = True
			print "==> Removing adapters...\n"
			fastqc.single_cut_adapters(adapt, fq1, args["outdir"])
			fq1 = args["outdir"]+"/trimmed.fastq" 
		print "==> Aligning fastq's...\n"
		if args["version"] == "1":
			alignment.single_bowtie(fq1, args["samname"], args["index"], args["outdir"])
		elif args["version"] == "2":
			alignment.single_bowtie2(fq1, args["samname"], args["index"], args["outdir"], args["threads"])
		zip_trimmed(adapt, args["outdir"], False)
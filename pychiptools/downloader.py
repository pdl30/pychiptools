#!/usr/bin/python

########################################################################
# 09 Oct 2014
# Patrick Lombard, Centre for Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import os, re, sys
import argparse
import subprocess
import ConfigParser
from multiprocessing import Pool, Manager
import itertools
import shutil
from collections import defaultdict

def ConfigSectionMap(section, Config):
	dict1 = {}
	options = Config.options(section)
	for option in options:
		try:
			dict1[option] = Config.get(section, option)
			if dict1[option] == -1:
				DebugPrint("skip: %s" % option)
		except:
			print("exception on %s!" % option)
			dict1[option] = None
	return dict1

def gse_download(gse):
	gsm_samples = defaultdict(list)
	download2 = "wget -c -nv -q 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='{0}'&targ=self&view=quick&form=text&token=_blank' -O tmp/{0}.soft".format(gsm)
	subprocess.call(download2, shell=True)
	with open("tmp/{}.soft".format(gse)) as f:
		for line in f:
			line  = line.rstrip()
			#GSM Samples
			if line.startswith("!Series_sample_id = "):
				gsm = re.sub("!Series_sample_id = ", "", line)
				gsm_samples.append(gsm)
	return gsm_samples

def downloader(gsm):
	old_path = os.getcwd()
	download2 = "wget -c -nv -q 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='{0}'&targ=self&view=quick&form=text&token=_blank' -O tmp/{0}.soft".format(gsm)
	subprocess.call(download2, shell=True)
	if not os.path.isdir(gsm):
		os.mkdir(gsm)
	f = open("tmp/{}.soft".format(gsm), "r")
	lines = f.readlines()
	sra = None
	gsm_sra = []
	os.chdir(gsm)
	new_path = os.getcwd()
	for line in lines:
		line = line.rstrip()
		if line.startswith("!Sample_supplementary_file"):
			sra_path = line.split("ByExp/sra/")
			if len(sra_path) > 1:
				sra = sra_path[1]
				download3 = "wget -r -c --no-verbose -N -nd ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/{}".format(sra)
				subprocess.call(download3, shell=True)
	sras = [f for f in os.listdir(new_path) if f.endswith(".sra")]
	for sra in sras:
		command0 = "fastq-dump --split-3 {}".format(sra)
		subprocess.call(command0, shell=True)
		os.remove(sra)
	#Catting everything together
	fq1s = [f for f in os.listdir(new_path) if f.startswith("SRR") and f.endswith("_1.fastq")]
	if fqs:
		command = "cat"
		for fq in fqs:
			command += " {}".format(fq)
		command += " > {}_1.fastq".format(gsm)
		subprocess.call(command, shell=True)
		fq1s = [f for f in os.listdir(new_path) if f.startswith("SRR") and f.endswith("_2.fastq")]
		command = "cat"
		for fq in fqs:
			command += " {}".format(fq)
		command += " > {}_2.fastq".format(gsm)
		subprocess.call(command, shell=True)
	else:
		fq1s = [f for f in os.listdir(new_path) if f.startswith("SRR") and f.endswith(".fastq")]
		command = "cat"
		for fq in fqs:
			command += " {}".format(fq)
		command += " > {}.fastq".format(gsm)
		subprocess.call(command, shell=True)
	#Remove old SRR fastqs
	for fq in fqs:
		os.remove(fq)
	os.chdir(old_path)

def download_function(args):
	return downloader(*args)

def main():
	parser = argparse.ArgumentParser(description='Downloads data from GEO and Arrayexpress\n')
	parser.add_argument('-g','--GSE', help='GSE accession', required=False)
	parser.add_argument('-c','--config', help='Input ConfigParser file containing [Conditions] with GSMs as keys', required=False)
	parser.add_argument('-t','--threads', help='Number of threads, default=20', default=20, required=False)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())

	if not os.path.isdir("tmp"):
		os.mkdir("tmp/")
	
	if args["GSE"]:
		conditions = gse_download(args["GSE"])
	elif args["config"]:
		Config = ConfigParser.ConfigParser()
		Config.optionxform = str
		Config.read(args["config"])
		conditions = ConfigSectionMap("Conditions", Config)
		
	pool = Pool(int(args["threads"]))
	pool.map(download_function, itertools.izip(list(conditions.keys())))
	pool.close()
	pool.join()	
	shutil.rmtree('tmp/')
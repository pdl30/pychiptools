#!/usr/bin/python

########################################################################
# 15 May 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import sys, os, re
import pybedtools
import pysam
import argparse
import operator
import pkg_resources
import pychiptools
import ConfigParser

def combine_sam_files(list_of_sams, outname):
	count = 0
	outsam = outname + "/" + outname + ".sam"
	print "==> Combining sam files...\n"
	for sam in list_of_sams:
		original_file = pysam.Samfile(sam)
		if count == 0:
			new_file = pysam.Samfile(outsam, mode='wh', template=original_file)
			for read in original_file: 
				new_file.write(read)
		else:
			for read in original_file: 
				new_file.write(read)
		count += 1

def convert_sam_bed(name, paired):
	obed = name+"_tmp.BED"
	outbed = open(obed, "wb")
	samfile = pysam.Samfile(name+".sam", "r")
	data = {}
	count = 0
	print "==> Converting sam to bed...\n"
	for read in samfile.fetch():
		count += 1
		strand = '+'
		if read.is_reverse :
			strand = '-'
		if strand == '+':
			new_start = read.pos
			new_end = int(read.pos) + 200
		elif strand == '-':
			new_end = read.aend
			new_start = int(read.aend) - 200

		if new_start <= 0 :
			new_start = 1
			new_end = 200
		outbed.write("{}\t{}\t{}\t{}\t0\t{}\n".format(samfile.getrname(read.tid), new_start, new_end, read.qname, strand)),

	outbed.close()
	command = "sort -k1,1 -k2,2g -o {} {}".format(name+".BED", name+"_tmp.BED")
	subprocess.call(command.split())
	subprocess.call(["rm", name+"_tmp.BED"])
	if paired:
		count /= 2
	return count
	
def change_for_ucsc(name, chromsizes, ens=False):
	if ens:
		outbed2 = open(name+'_tmp1.BED', "w")
		print "==> Converting Ensembl to UCSC chromosomes...\n"
		with open(name+".BED") as f:
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				if re.match(r"^\d", word[0]):
					new_chr = "chr" + word[0]
				elif re.match(r"^X", word[0]):
					new_chr = "chrX"
				elif re.match(r"^Y", word[0]):
					new_chr = "chrY"
				elif word[0] == "MT":
					new_chr = "chrM"
				else:
					pass
				outbed2.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(new_chr, word[1], word[2], word[3], word[4], word[5])),
		outbed2.close()
		subprocess.call(["bedClip", name+"_tmp1.BED", chromsizes, name+"_ucsc.BED"])	
		subprocess.call(["rm", name+"_tmp1.BED"])
	else:
		subprocess.call(["bedClip", name+".BED", chromsizes, name+"_ucsc.BED"])	
		os.remove(name+".BED")


def genomeCoverage(name, genome, scale=None):
	if scale:
		outg2 = name+"_rpm.bedGraph"
	else:
		outg2 = name+"_ucsc.bedGraph"
	inbed = pybedtools.BedTool(name+"_ucsc.BED")
	print "==> Creating bedGraph...\n"
	if scale:
		outcov = inbed.genome_coverage(bg=True, genome=genome, scale=scale)
	else:
		outcov = inbed.genome_coverage(bg=True, genome=genome)
	outcov.saveas(outg2)

def bedgraphtobigwig(name, chrom, house=False, rpm=False):
	print "==> Converting bedGraph to bigWig...\n"
	if rpm:
		command = ["bedGraphToBigWig", name+"_rpm.bedGraph", chrom, name+"_rpm.bw"]
	else:
		command = ["bedGraphToBigWig", name+"_ucsc.bedGraph", chrom, name+".bw"]
	subprocess.call(command)

def ConfigSectionMap(Config, section):
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
	
def main():
	parser = argparse.ArgumentParser(description='Processes ChIP-seq samples to bigWig tracks. Use either input or config file\n')
	parser.add_argument('-c', '--config', help='Contains [Conditions] with bam files as keys.', required=False)
	parser.add_argument('-i','--input', help='Input sam file', required=False)
	parser.add_argument('-p', action='store_true', help='Use if samples are paired end. Required if using RPM normalisation', required=False)
	parser.add_argument('-g','--genome', help='Genome the samples are aligned to, options include mm10/mm9/hg19', required=True)
	parser.add_argument('-e', action='store_true', help='Are samples aligned to Ensembl genome?', required=False)
	parser.add_argument('-rpm', action='store_true', help='Scale resulting bigwig to RPM', required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	chrom = pkg_resources.resource_filename('pychiptools', 'data/{}.chrom.sizes'.format(args["genome"]))
	
	path0 = os.getcwd()

	if args["input"]:
		name = re.sub(".sam$", "", args["input"])
		count = convert_sam_bed(name, args["p"])
		scale = float(1000000)/int(count)

		change_for_ucsc(name, chrom, args["e"])

		if args["rpm"]:
			genomeCoverage(name, args["genome"], scale)
			bedgraphtobigwig(name, chrom, rpm=True)	
		else:
			genomeCoverage(name, args["genome"])
			bedgraphtobigwig(name, chrom)

	elif args["config"]:
		Config = ConfigParser.ConfigParser()
		Config.optionxform = str
		Config.read(args["config"])

		conditions = ConfigSectionMap(Config, "Conditions")
		for key in conditions:
			name = re.sub(".sam$", "", key)
			count = convert_sam_bed(name, args["p"])
			scale = float(1000000)/int(count)

			change_for_ucsc(name, chrom, args["e"])

			if args["rpm"]:
				genomeCoverage(name, args["genome"], scale)
				bedgraphtobigwig(name, chrom, rpm=True)	
			else:
				genomeCoverage(name, args["genome"])
				bedgraphtobigwig(name, chrom)
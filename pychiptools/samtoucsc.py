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
	#inbed = pybedtools.BedTool(name+"_tmp.BED")
	#New sort:
	command = "sort -k1,1 -k2,2g -o {} {}".format(name+".BED", name+"_tmp.BED")
	#outbed2 = inbed.sort()
	#outbed2.saveas(name+".BED")
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


def genomeCoverage(name, scale=None):
	if scale:
		outg2 = name+"_rpm.bedGraph"
	else:
		outg2 = name+"_ucsc.bedGraph"
	inbed = pybedtools.BedTool(name+"_ucsc.BED")
	print "==> Creating bedGraph...\n"
	if scale:
		outcov = inbed.genome_coverage(bg=True, genome='mm10', scale=scale)
	else:
		outcov = inbed.genome_coverage(bg=True, genome='mm10')
	outcov.saveas(outg2)

def bedgraphtobigwig(name, chrom, house=False, rpm=False):
	print "==> Converting bedGraph to bigWig...\n"
	if rpm:
		command = ["bedGraphToBigWig", name+"_rpm.bedGraph", chrom, name+"_rpm.bw"]
	else:
		command = ["bedGraphToBigWig", name+"_ucsc.bedGraph", chrom, name+".bw"]
	subprocess.call(command)

def main():
	parser = argparse.ArgumentParser(description='Processes ChIP-seq samples to bigWig tracks.\n')
	parser.add_argument('-i','--input', help='Input sam file', required=False)
	parser.add_argument('-p', action='store_true', help='Use if samples are paired end. Required if using RPM normalisation', required=False)
	parser.add_argument('-g','--genome', help='Genome the samples are aligned to, options include mm10/mm9/hg19', required=True)
	parser.add_argument('-e', action='store_true', help='Are samples aligned to Ensembl genome?', required=False)
	parser.add_argument('-rpm', action='store_true', help='Scale resulting bigwig to RPM', required=False)
	args = vars(parser.parse_args())
	chrom = pkg_resources.resource_filename('pychiptools', 'data/{}.chrom.sizes'.format(args["genome"]))
	
	path0 = os.getcwd()
	name = re.sub(".sam$", "", args["input"])
	count = convert_sam_bed(name, args["p"])
	scale = float(1000000)/int(count)

	change_for_ucsc(name, chrom, args["e"])

	if args["rpm"]:
		genomeCoverage(name, scale)
		bedgraphtobigwig(name, chrom, rpm=True)	
	else:
		genomeCoverage(name)
		bedgraphtobigwig(name, chrom)
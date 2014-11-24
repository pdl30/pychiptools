#!/usr/bin/python

########################################################################
# 28 July 2014
# Patrick Lombard, Centre for Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

#Not tested and requires further tweaking!
import os, re, sys
import argparse
import subprocess

def run_meme_finder(fasta, outdir, nmotifs=10):
	cmd1 = "meme {} -nmotifs {} -oc {} -minw 8 -maxw 14 -maxsize 500000000 -mod oops".format(fasta, nmotifs, outdir)
	subprocess.call(cmd1.split())

def run_fimo_scanner(meme_motifs, fasta, outdir):
	cmd2 = "fimo --oc {} {} {}".format(outdir, meme_motifs, fasta)
	subprocess.call(cmd2.split())

#Useful for Homer which takes bed files of a certain width
def adjust_peak_size(bed_file, outbed_file, size=400):
	output = open(outbed_file, "w")
	each_side = size/2
	with open(bed_file) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			data = word[:2]
			rest = word[3:]
			start = data[1]
			end = data[2]
			middle = int(start) + ((int(end) - int(start))/2)
			new_start = middle - each_side
			new_end = middle + each_side
			output.write("{}\t{}\t{}".format(word[0], new_start, new_end)),
			for k in rest:
				output.write("\t{}".format(k)),
			output.write("\n"),
	output.close()

def run_homer(bed, genome, nmotifs=25, lmotif=12):
	name = re.sub(".bed", "_homer", bed)
	command = "findMotifsGenome.pl {} {} {} -S {} -len {}".format(bed, genome, name, nmotifs, lmotif)
	subprocess.call(command.split())
	return name

def fasta_from_bed(bed, ifasta, outfa):
	command2 = "fastaFromBed -fi {} -bed {} -fo {}".format(ifasta, bed, outfa)
	subprocess.call(command2.split())

def process_homer_motifs(homer_dir, outdir):
	command = ["cp", "{}/homerMotifs.all.motifs", "{}/"]


def main():
	parser = argparse.ArgumentParser(description='''Complete motif analyis script.\n
		Use motif_analysis.py meme -h for example to see subcommands\n''')

	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")

	meme_parser = subparsers.add_parser('meme', help="Runs MEME")
	meme_parser.add_argument('-f','--FASTA', help='Input fasta file', required=False)
	meme_parser.add_argument('-o', '--OUTDIR', help='Output directory', required=False)
	meme_parser.add_argument('-n', '--NMOTIFS', help='Number of Motifs. Default is 10', required=False, default=10)
	meme_parser.add_argument('-fimo', help='Runs FIMO on Homer Motifs', required=False, action='store_true')

	homer_parser = subparsers.add_parser('homer', help='Runs Homer')
	homer_parser.add_argument('-b','--BED', help='Input BED file', required=False)
	homer_parser.add_argument('-n', '--NMOTIFS', help='Number of Motifs. Default is 25', required=False, default=25)
	homer_parser.add_argument('-g', '--GENOME', help='Depends on what genomes homer has installed e.g. mm10/hg19. Default is mm10', required=True, default="mm10")
	homer_parser.add_argument('-l', '--LMOTIFS', help='Length of Motifs. Default is 12', required=False, default=12)
	homer_parser.add_argument('-fimo', help='Runs FIMO on Meme Motifs', required=False, action='store_true')

	fimo_parser = subparsers.add_parser('fimo', help='Runs FIMO')
	fimo_parser.add_argument('-i','--INPUT', help='Input motifs in meme format', required=False)
	fimo_parser.add_argument('-f','--FASTA', help='Reference fasta', required=False)
	fimo_parser.add_argument('-o', '--OUTDIR', help='Output directory', required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())

	if args["subparser_name"] == "meme":
		run_meme_finder(args["FASTA"], args["OUTDIR"], args["NMOTIFS"])

	if args["subparser_name"] == "homer":
		run_homer(args["BED"], args["GENOME"], args["NMOTIFS"], args["LMOTIFS"])

	if args["subparser_name"] == "fimo":
		run_fimo_scanner(args["INPUT"], args["FASTA"], args["OUTDIR"])
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
import ConfigParser
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import IntVector, FloatVector, StrVector
import rpy2.robjects.numpy2ri as rpyn
from collections import defaultdict
import math
from pychiptools.utilities import mmdiff, count_method, manorm, macs2_bgdiff, diffbind, diffreps, pepr_analysis
from multiprocessing import Pool
import itertools

#Programs to consider are diffReps, diffbind, MMDiff and MACS2
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

def reverse_dict(idict):
	inv_map = {}
	for k, v in idict.iteritems():
		inv_map[v] = inv_map.get(v, [])
		inv_map[v].append(k)
	return inv_map

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

#All possibilities and then choose which ones you want
def get_config_args(args, arglist):
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["CONFIG"])
	dict_list = []
	for arg in arglist:
		dict_list.append(ConfigSectionMap(Config, arg))
	return dict_list

def count_function(args):
	return count_method.count(*args)

def run_parallel(comparisons, *argv):
	#just an example. Maybe use kwargs for named arguments to denote which ones are repeats
	#Order is also a problem
	new_comp = []
	for comp in comparisons:
		c = comparisons[comp].split(",")
		comps = [x.strip(' ') for x in c]
		new_comp.append((comps[0], comps[1]))
	print inv_conds
	pool = Pool(8)
	pool.map(count_function, itertools.izip(itertools.repeat(inv_conds), itertools.repeat(peaks), new_comp, itertools.repeat(counts)))
	pool.close()
	pool.join()

def main():
	parser = argparse.ArgumentParser(description='Differential ChIPseq Analysis.\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")

	diffbind_parser = subparsers.add_parser('diffbind', help="Runs DiffBind for Differential Binding!")
	diffbind_parser.add_argument('-c','--CONFIG', help='Config file containing [Conditions], [Peaks], [Controls] and [Comparisons]', required=True)
	diffbind_parser.add_argument('-a', action='store_true', help='Ignore controls from config', required=False)
	
	mmdiff_parser = subparsers.add_parser('mmdiff', help="Runs MMDiff for Differential Binding!")
	mmdiff_parser.add_argument('-c','--CONFIG', help='Config file containing [Conditions], [Peaks], [Controls] and [Comparisons]', required=True)
	mmdiff_parser.add_argument('-method', help='MANorm method, options are MMD, GMD or Pearson. MMD -slow and GMD -v.fast', required=False)

	macs2_parser = subparsers.add_parser('macs2', help="Runs MAC2 for Differential Binding!")
	macs2_parser.add_argument('-c','--CONFIG', help='Config file containing [Conditions], [Controls] and [Comparisons] using BAM files', required=True)
	macs2_parser.add_argument('-a', action='store_true', help='Are samples histones?', required=False)
	macs2_parser.add_argument('-b', action='store_true', help='Use [bedGraph] from config. This will also query the bam for useable reads value', required=False)
	macs2_parser.add_argument('-p', action='store_true', help='Postprocess macs2 output', required=False)
	macs2_parser.add_argument('-gtf', help='GTF for annotation of results, this assumes it is ensembl formatted!', required=False)

	diffReps_parser = subparsers.add_parser('diffReps', help="Runs diffReps for Differential Binding!")
	diffReps_parser.add_argument('-c','--CONFIG', help='Config file containing [Conditions], [Controls] and [Comparisons] using BED files', required=True)
	diffReps_parser.add_argument('-a', action='store_true', help='Ignore controls from config', required=False)
	diffReps_parser.add_argument('-chrom', help='Chomosome sizes', required=False)
	diffReps_parser.add_argument('-method', help='Diffrep method, options are nb=Negative binomial; gt=G-test; tt=T-test; cs=Chi-square test', required=False)
	diffReps_parser.add_argument('-n', action='store_true', help='Uses [Counts] from config to get normalisation constants for program', required=False)
	diffReps_parser.add_argument('-p', action='store_true', help='Postprocess macs2 output', required=False)
	diffReps_parser.add_argument('-gtf', help='GTF for annotation of results, this assumes it is ensembl formatted!', required=False)

	manorm_parser = subparsers.add_parser('MAnorm', help="Runs MAnorm for Differential Binding!")
	manorm_parser.add_argument('-c','--CONFIG', help='Config file containing [Conditions], [Peaks], and [Comparisons] using BED files', required=True)
	manorm_parser.add_argument('-d', action='store_true', help='Download the package to current directory', default=False, required=False)

	count_parser = subparsers.add_parser('count', help="Runs HTseq-count and also ranks changes in peaks")
	count_parser.add_argument('-c','--CONFIG', help='Config file containing [Conditions], [Peaks] and [Comparisons] using BAM files', required=True)
	count_parser.add_argument('-a', action='store_true', help='Use [Counts] from config', required=False)

	pepr_parser = subparsers.add_parser('pepr', help="Runs Pepr Differential binding analysis. Needs Reps and review of script")
	pepr_parser.add_argument('-c','--CONFIG', help='Config file containing [Conditions] and [Comparisons] using BED files', required=True)
	args = vars(parser.parse_args())

	#Most of these are slow so must figure out how to run everything together!
	if args["subparser_name"] == "diffbind":
		conditions, controls, peaks, comparisons = get_config_args(args, ["Conditions", "Controls", "Peaks", "Comparisons"])
		inv_conds = reverse_dict(conditions)
		for comp in comparisons:
			c = comparisons[comp].split(",")
			comps = [x.strip(' ') for x in c]
			if args["a"]:
				reps = diffbind.create_diffbind_input(inv_conds, peaks, comps[0], comps[1]) 
			else:
				reps = diffbind.create_diffbind_input(inv_conds, peaks, comps[0], comps[1], controls) 	
				
			rscript = diffbind.write_diffbind(comps[0], comps[1], reps)
			run_rcode(rscript, "diffbind_rcode.R")

	elif args["subparser_name"] == "mmdiff": ##Not finished, needed reps
		conditions, controls, peaks, comparisons = get_config_args(args, ["Conditions", "Controls", "Peaks", "Comparisons"])
		inv_conds = reverse_dict(conditions)
		for comp in comparisons:
			c = comparisons[comp].split(",")
			comps = [x.strip(' ') for x in c]
			reps = mmdiff.create_diffbind_input(inv_conds, controls, peaks, comps[0], comps[1]) #Same as diffbind so doesn't matter
			rscript = mmdiff.write_mmdiff(comps[0], comps[1], reps)
			run_rcode(rscript, "mmdiff_rcode.R")

	elif args["subparser_name"]=="macs2":
		conditions, controls, comparisons = get_config_args(args, ["Conditions", "Controls", "Comparisons"])
		inv_conds = reverse_dict(conditions)
		for comp in comparisons:
			c = comparisons[comp].split(",")
			comps = [x.strip(' ') for x in c]
			if args["p"]:
				macs2_bgdiff.postprocess_bgdiff(comps[0], comps[1], args["gtf"])
			elif args["b"]:
				bedgraphs= get_config_args(args, ["BedGraph"])[0]
				if args["a"]:
					macs2_bgdiff.bgdiff(inv_conds, comps[0], comps[1], controls, "histone", bedgraphs=bedgraphs)
				else:
					macs2_bgdiff.bgdiff(inv_conds, comps[0], comps[1], controls, "tf", bedGraphs=bedgraphs)
			else:
				if args["a"]:
					macs2_bgdiff.bgdiff(inv_conds, comps[0], comps[1], controls, "histone")
				else:
					macs2_bgdiff.bgdiff(inv_conds, comps[0], comps[1], controls, "tf")

	elif args["subparser_name"] == "diffReps":
		conditions, controls, comparisons = get_config_args(args, ["Conditions", "Controls", "Comparisons"])
		inv_conds = reverse_dict(conditions)
		for comp in comparisons:
			c = comparisons[comp].split(",")
			comps = [x.strip(' ') for x in c]
			if args["n"]:
				counts= get_config_args(args, ["Counts"])[0]
				if args["a"]:
					diffreps.diffReps(inv_conds, comps[0], comps[1], args["chrom"], method=args["method"], controls=None, counts=counts)
				else:
					diffreps.diffReps(inv_conds, comps[0], comps[1], args["chrom"], method=args["method"], controls=controls, counts=counts)
			elif args["p"]:
				diffreps.postprocess_diffreps(comps[0], comps[1], args["gtf"])
			else:
				if args["a"]:
					diffreps.diffReps(inv_conds, comps[0], comps[1], args["chrom"], method=args["method"], controls=None)
				else:
					diffreps.diffReps(inv_conds, comps[0], comps[1], args["chrom"], method=args["method"], controls=controls)
	elif args["subparser_name"] == "MAnorm":
		conditions, peaks, comparisons = get_config_args(args, ["Conditions", "Peaks", "Comparisons"])
		inv_conds = reverse_dict(conditions)
		for comp in comparisons:
			c = comparisons[comp].split(",")
			comps = [x.strip(' ') for x in c]
			manorm.manorm(inv_conds, peaks, comps[0], comps[1], download=args["d"])

	elif args["subparser_name"] == "count":
		#Not the best example for parallisation!
		conditions, peaks, comparisons = get_config_args(args, ["Conditions", "Peaks", "Comparisons"])
		if args["a"]:
			counts = get_config_args(args, ["Counts"])[0]
		else:
			counts = None
		inv_conds = reverse_dict(conditions)
		for comp in comparisons:
			c = comparisons[comp].split(",")
			comps = [x.strip(' ') for x in c]
			count_method.count(inv_conds, peaks, comps[0], comps[1], counts)
	elif args["subparser_name"] == "pepr":
		conditions, controls, comparisons = get_config_args(args, ["Conditions", "Controls", "Comparisons"])
		inv_conds = reverse_dict(conditions)
		for comp in comparisons:
			c = comparisons[comp].split(",")
			comps = [x.strip(' ') for x in c]
			pepr_analysis.pepr_diff_analysis(inv_conds, controls, comps[0], comps[1])
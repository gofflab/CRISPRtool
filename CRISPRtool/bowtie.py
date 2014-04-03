#!/usr/bin/env python

#Imports
from subprocess import Popen, PIPE, STDOUT

#Set defaults
bowtieBin = "bowtie"
bowtieIdx = "~/indexes/mm9/bowtie/mm9"

SAMFields = [
		'readIdx',
		'bitflag',
		'chrom',
		'start',
		'mapQual',
		'CIGAR',
		'mate', #Should be '*'
		'matepos', #Should always be 0
		'insertSize', #Should always be 0
		'readSequence',
		'readQual', # Defaults to string of "I" corresponding to PHRED 40 when no quality scores are available
		'optionalFields'
	]

BowtieFields = [
		'readIdx',
		'strand',
		'chrom',
		'start',
		'readSequence',
		'readQual',
		'numAlignSameMismatch',
		'mismatches'
]

#Take a Fasta file and run bowtie2
class BowtieAlignment:
	def __init__(self,):
		pass

def runBowtie(fasta,bin=bowtieBin,idx=bowtieIdx):
	bowtieCmd = "%s -k 50 --best -l 23 -v 3 -f %s -c '%s'" % (bowtieBin,bowtieIdx,fasta)
	p = Popen(bowtieCmd, shell=True,stdout=PIPE)
	res = p.communicate()
	return res[0]

def parseBowtie(bowtieRes,guides):
	# Ignore header lines (^@) and comment lines (^#)
	lines = bowtieRes.split("\n")
	alignments = []
	
	for line in lines:
		if line.startswith("#") or line == '':
			continue
		vals = line.split("\t")
		align = dict(zip(BowtieFields,vals))
		#mmpos = [0] * 23
		mmpos = []
		try:
			mmlist = align['mismatches'].split(",")
			align['numMM'] = len(mmlist)
			for x in mmlist:
				#mmpos[int(x.split(":")[0])]+=1
				mmpos.append(int(x.split(":")[0]))
		except:
			align['numMM'] = 0
		align['mmpos'] = mmpos
		guides[int(align['readIdx'])].alignments.append(align)
	return guides
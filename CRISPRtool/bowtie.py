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
	bowtieCmd = "%s -k 20 -l 23 -v 3 -f %s -c '%s'" % (bowtieBin,bowtieIdx,fasta)
	p = Popen(bowtieCmd, shell=True,stdout=PIPE)
	res = p.communicate()
	return res[0]

def parseBowtie(bowtieRes,guides=''):
	# Ignore header lines (^@) and comment lines (^#)
	lines = bowtieRes.split("\n")
	alignments = []
	for line in lines:
		vals = line.split("\t")
		alignments.append(dict(zip(BowtieFields,vals)))
	return alignments
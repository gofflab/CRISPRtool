#!/usr/bin/env python

#Imports
from subprocess import Popen, PIPE, STDOUT

#Set defaults
bowtieBin = "~/bin/bowtie2"
bowtieIdx = "/n/rinn_data1/indexes/mouse/mm9/bowtie2/mm9"


#Take a Fasta file and run bowtie2
class Bowtie:
	def __init__(self):
		pass


def runBowtie(bin=bowtieBin,idx=bowtieIdx,argDict):
	bowtieCmd = "%s -x %s"
	pipe = subprocess.Popen(bowtieCmd, shell=True, bufsize=bufsize, stdin=subprocess.PIPE).stdin
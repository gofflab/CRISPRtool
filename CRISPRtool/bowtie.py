#!/usr/bin/env python

#Imports
from subprocess import Popen, PIPE, STDOUT

#Set defaults
bowtieBin = "bowtie"
bowtieIdx = "~/indexes/mm9/bowtie2/mm9"



#Take a Fasta file and run bowtie2
class Bowtie:
	def __init__(self):
		pass


def runBowtie(fasta,bin=bowtieBin,idx=bowtieIdx):
	bowtieCmd = "%s -f %s -c '%s'" % (bowtieBin,bowtieIdx,fasta)
	p = Popen(bowtieCmd, shell=True)
	res = p.communicate()
	return res[0]

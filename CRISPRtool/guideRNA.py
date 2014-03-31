#!/usr/bin/env python
import sequencelib,getopt,sys,re
import tempfile
import bowtie

# Base class for a guide RNA
# Takes a 23mer 20 Guide + 3 PAM and processes it accordingly
class GuideError(Exception):
	def __init__(self,value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class GuideRNA:
	def __init__(self,sequence,startPos,seqName,strand='+'):
		self.sequence = sequence.upper()[:20]
		self.pam = sequence.upper()[-3:]
		self.start = startPos
		self.strand = strand
		self.seqName = seqName
		self.name = "%s:%s:%d" % (self.seqName,self.strand,self.start)
		self.validate()
	def validate(self):
		if self.pam[-1] != 'G' and self.pam[-2] not in ['A','G']:
			raise GuideError("%s: Not valid PAM sequence" % self.name)
		if (len(self.sequence)!=20):
			raise GuideError("%s: Guide not appropriate length" % self.name)
	def __repr__(self):
		return "%s - %s" % (self.sequence,self.pam)
	def toFasta(self):
		return ">%s\n%s" % (self.name,self.sequence+self.pam)


def usage():
	sys.stderr.write("")


#######################
# Variables
#######################
defaultPams = [
		'NGG',
		'NAG'
	]
guideSize = 20

#######################
# Scan input sequence #
#######################
def makePamRe(pamList=defaultPams):
	pamStrings = []
	for pam in pamList:
		tmp = ''
		for n in pam:
			tmp += sequencelib.iupacdict.get(n)
		pamStrings.append(tmp)
	pamSearchString = "|".join(pamStrings)
	pamRe = re.compile(pamSearchString)
	return pamRe

def scanSequence(sequence,seqName,pamList=defaultPams):
	pamRe = makePamRe(pamList)
	#instatiate matches
	guides = []
	#Find PAM hits on both strands
	for strand in ["+","-"]:
		searchSeq = sequence.upper()
		if strand == "-":
			searchSeq = sequencelib.rcomp(searchSeq)
		pamIter = pamRe.finditer(searchSeq,guideSize)
		for match in pamIter:
			if strand == "+":
				start = match.start()-guideSize
			elif strand == "-":
				start = len(sequence)-match.start()
			guides.append(GuideRNA(searchSeq[match.start()-guideSize:match.end()],start,seqName,strand))
	return guides

def alignGuides(guides):
	fname = open('guides.fasta','w')
	for g in guides:
		print >>fname, g.toFasta()
	return

def main():
	try:
		opts,args = getopt.getopt(sys.argv[1:],"ho:v",["help","output="])
	except getopt.GetoptError as err:
		print(err)
		usage()
		sys.exit(2)
	output = None
	verbose = False
	for o,a in opts:
		if o == "-v":
			verbose = True
		elif o in ("-h","--help"):
			usage()
			sys.exit()
		elif o in ("-o","--output"):
			output = a
		else:
			assert False, "Unhandled option"

	#Main workflow

def test():
	fname = "test/test.fasta"
	handle = open(fname,'r')
	fastaIter = sequencelib.FastaIterator(handle)

	mySeq = fastaIter.next()
	guides = scanSequence(mySeq['sequence'],mySeq['name'])
	
	alignRes = bowtie.runBowtie(",".join([x.sequence+x.pam for x in guides]))

	print alignRes

if __name__ == "__main__":
	test()
#!/usr/bin/env python
import sequencelib,getopt,sys,re
import tempfile
import bowtie
from itertools import tee,izip

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

		#Setup empty fields for later use
		self.alignments = []
		self.score = 100.0
		self.overlap = False

		#Validate guide
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

#Mismatch Matrix (M)
M = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.538]
PamM = [0,1,1]

#######################
# Helper functions
#######################
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)


def meanPairwiseDist(mmPos):
	if len(mmPos)==1:
		return 0
	else:
		count = 0
		dists = []
		for pair in pairwise(mmPos):
			count += 1
			dists.append(max(pair)-min(pair))
		return sum(dists)/float(count)

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

# def alignGuides(guides):
# 	fname = open('guides.fasta','w')
# 	for g in guides:
# 		print >>fname, g.toFasta()
# 	return

def scoreAlignments(guides,M=M,PamM=PamM):
	mismatchM = M + PamM
	for guide in guides:
		scores = []
		for alignment in guide.alignments:
			if len(alignment['mmpos'])>0:
				#Calculate 
				term1 = []
				for mm in alignment['mmpos']:
					term1.append(1-mismatchM[mm])
				term1 = reduce(lambda x, y: x*y, term1)
				#print term1
				term2 = 1.0/(((19.0-meanPairwiseDist(alignment['mmpos'])/19.0)*4.0)+1.0)
				#print term2
				if (len(alignment['mmpos'])!=0):
					#print len(alignment['mmpos'])
					term3 = 1.0/(len(alignment['mmpos'])**2)
				else:
					term3 = 1
				#print term3
				alignment['score'] = term1*term2*term3
			else:
				alignment['score'] = 0
	return guides

def summarizeGuideScores(guides):
	for guide in guides:
		sumStat = []
		for hit in guide.alignments:
			sumStat.append(hit['score'])
		sumStat = sum(sumStat)
		print sumStat
		guide.score = 100.0/(100.0+sumStat)
		print guide.score
	return guides

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
	#print guides[0].alignments

	alignRes = bowtie.runBowtie(",".join([x.sequence+x.pam for x in guides]))

	guides = bowtie.parseBowtie(alignRes,guides)
	
	#print guides[0].alignments
	#print len(guides[0].alignments)

	guides = scoreAlignments(guides)

	guides = summarizeGuideScores(guides)

	for guide in guides:
		print "%s\t%f\t%d" % (guide,guide.score,len(guide.alignments))

if __name__ == "__main__":
	test()
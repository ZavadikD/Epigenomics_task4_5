#!/usr/bin/env python


#************
# LIBRARIES *
#************

import sys
from optparse import OptionParser


#*****************
# OPTION PARSING *
#*****************

parser = OptionParser()
parser.add_option("-i", "--input", dest="input")
parser.add_option("-s", "--start", dest="start")
options, args = parser.parse_args()

open_input = open(options.input)
enhancer_start = int(options.start)
#print(enhancer_start)

#********
# BEGIN *
#********

x=1000000 # set maximum distance to 1 Mb
selectedGene="" # initialize the gene as empty
selectedGeneStart=0 # initialize the start coordinate of the gene as empty

for line in open_input.readlines(): # for each line in the input file
	y = line.strip().split('\t') # split the line into two columns based on a tab 
	position = int(y[1])
	# define a variable called position that correspond to the integer of the start of the gene
	value = abs(position-enhancer_start)
	# compute the absolute value of the difference between position and enhancer_start
	if value < x :
		if value != 0:
			x=value
			selectedGene=y[0]
			selectedGeneStart=position
	# if this absolute value is lower than x
		# this value will now be your current x
		# save gene as selectedGene
		# save position as selectedGeneStart

print ("\t".join([selectedGene, str(selectedGeneStart), str(x)]))

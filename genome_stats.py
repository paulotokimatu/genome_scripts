#!usr/bin/python

# Author: Paulo Tokimatu
# Last update: 2016/05/20
# This script is used to vizualize several useful statistics for genome avaliation.
# Requirements: Biopython (http://biopython.org/wiki/Biopython) installed.
# Usage: genome_stats.py [GENOME_FILE]
# In case of bugs, suggestions or questions, please contact me: paulomtf(at)lge.ibi.unicamp.br

import sys
from sys import argv
from Bio import SeqIO

def CountingLetters(sequence):
	for i in number_letters:
		number_letters[i] += sequence.count(i)

def CalcN50_N90(seq_size, genome):
	temp = 0
	n = 0
	switch = 0
	for i in sorted(seq_size.values(), reverse=True):
		n += 1
		temp += i
		if temp >= (genome/2) and switch == 0:
			print "# N50 = %s (L50 = %s)" % (i, n)
			switch = 1
		if temp >= (genome*0.9):
			print "# N90 = %s (L90 = %s)" % (i, n)
			break

size = {}
genome_size = 0
number_letters = {"A":0, "T":0, "C":0, "G":0, "N":0}
error_msg = "Fasta file not found. Please, check your command line.\nUsage: python genome_stats.py [GENOME_FILE]"

if len(argv) <= 1:
	print error_msg
	sys.exit()

# Reading the fasta file
try:
	for seq_record in SeqIO.parse(argv[1], "fasta"):
		CountingLetters(seq_record.seq)
		size[seq_record.id] = len(seq_record)
		genome_size += len(seq_record)
except:
	print error_msg
	sys.exit()

# Printing statistics
print "# FILE INPUT = %s" % argv[1]
print "# GENOME LENGTH = %s (%.2fMb)" % (genome_size, (float(genome_size)/1000000))
print "# NUMBER OF \"N\" BASES = %s (%.1f%%)" % (number_letters["N"], float(number_letters["N"])/genome_size*100)
print "# NUMBER OF SEQUENCES = %s" % len(size) 
CalcN50_N90(size, genome_size)
print "# A: %.1f%% T: %.1f%% C: %.1f%% G: %.1f%%" % ((float(number_letters["A"])/genome_size)*100, (float(number_letters["T"])/genome_size)*100, (float(number_letters["C"])/genome_size)*100, (float(number_letters["G"])/genome_size)*100)
print "#----------------------------------"
print "# Sequence\tLength"

# Printing each sequence in a descending bp length order
for key, value in sorted(size.iteritems(), key=lambda (k,v): (v,k), reverse=True):
	print "%s \t %s" % (key, value)



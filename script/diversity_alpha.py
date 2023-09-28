#!/usr/bin/env python
import os, sys, argparse
import math

def shannons_alpha(p):
	# shannons = -sum(pi ln(pi))
	h = []
	for i in p:
		h.append(i * math.log(i))
	print("Shannon's diversity: %s" %(-1 *sum(h)))
	return (-1 *sum(h))

def eqn_output(a):
	return a * np.log(1+N_f/a) - S_f

# Main method
def main():
	# get arguments
	parser = argparse.ArgumentParser(description='Pick an alpha diversity.')
	parser.add_argument('-i','--infile',dest='infile',help='abundance file')
	parser.add_argument('-o','--outfile',dest='outfile',default='Sh',type=str,  help='Shannon index')
	args = parser.parse_args()

	f = open(args.infile)
	
	f.readline()
	n = []
	# read in the file
	for line in f: 
		ind_abund = line.split('\t')[1] # finds the abundance estimate
		n.append(float(ind_abund)) 
	
	f.close()
	
	# calculations
	N = sum(n) # total number of individuals
	S = len(n) # total number of species
	# calculate all the pi's
	p = [] # store pi's 
	D = 0
	for i in n: # go through each species
		if i != 0: # there should not be any zeros
			p.append(i/N) # pi is the ni/N
			D += i*(i-1)

	D = D/(N*(N-1))
	index=shannons_alpha(p)

	fw = open(args.outfile,"w")
	fw.write("Shannon's diversity: %s" % index)
	fw.close()

if __name__ == "__main__":
	main()


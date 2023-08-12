#!/usr/bin/env python
import os, sys, argparse
import operator
from time import gmtime
from time import strftime
import numpy as np
#Main method
def main():
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', required=True, dest='infile',
        help='Input file where the first column is the category and the rest are counts')
    parser.add_argument('-o','--outfile', dest='outfile', help='Output file containing the bray diatance matrix')
    args=parser.parse_args()

    i2totals = {}
    i2counts = {}
    i2names = {}
    num_samples = 0
    num_categories = 0

    header = True
    infile = open(args.infile,'r')
    for line in infile:
        l_vals = line.strip().split("\t")
        #Read header
        if header:
            s_count = 0
            for i in range(1, len(l_vals)):
                i2names[s_count] = l_vals[i]
                i2counts[s_count] = {}
                i2totals[s_count] = 0
                s_count += 1
            num_samples = s_count
            header = False
        else:
            s_count = 0
            curr_categ = l_vals[0]
            for i in range(1, len(l_vals)):
                if float(l_vals[i]) > 0:
                    i2totals[s_count] += float(l_vals[i])
                    i2counts[s_count][curr_categ] = float(l_vals[i])
                s_count += 1
            num_categories += 1
    infile.close()

    bc = np.zeros((num_samples,num_samples))
    for i in range(0,num_samples):
        i_tot = i2totals[i]
        for j in range(i+1, num_samples):
            j_tot = i2totals[j]
            C_ij = 0.0
            for cat in i2counts[i]:
                if cat in i2counts[j]:
                    C_ij += min(i2counts[i][cat], i2counts[j][cat])
            
            bc_ij = 1.0 - ((2.0*C_ij)/float(i_tot+j_tot))
            
            bc[i][j] = bc_ij
            bc[j][i] = bc_ij

    #################################################
    out = open(args.outfile,'w')
    out.write("ID")
    for i in range(num_samples):
        out.write("\t%s" % i2names[i])
    out.write("\n")
    
    for i in range(num_samples):
        out.write("%s" % i2names[i])
        for j in range(num_samples):
            if i <= j:
                out.write("\t%0.3f" % bc[i][j])
            else:
                out.write("\tx.xxx")
        out.write("\n")

####################################################################
if __name__ == "__main__":
    main()

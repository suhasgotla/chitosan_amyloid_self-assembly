#!/usr/bin/python2.7
'''
    Generate a random charge profile based on given protonation fraction
    Input:
        n_chain        = number of chains
        n_monomer       = number of monomers
        prot_percentage = percentage of monomers to be protonated
    Output:
        proto.txt       = n_monomer x n_chain matrix of randomly placed charges 0 and 1. 
    Note that for 100% or 0% protonation, proto.txt just needs to be a matrix of 1s or 0s respectively                          
'''
import numpy as np

outdat = './proto.txt'

# User inputs
n_chain         = 2 
n_monomer       = 30
prot_percentage = 50 

# Generate a n_chain x n_monomer matrix of random numbers between 0 and 1
rds = np.random.rand(n_monomer, n_chain)

# User defined fraction of monomers to be protonated
rdsPercentile = np.percentile(rds, prot_percentage)
# if element_ij > rdsPercentile => False; elif element_ij < rdsPercentile => True
proto = np.less(rds, rdsPercentile).astype(int)

np.savetxt(outdat, proto, "%-2d"*n_chain)

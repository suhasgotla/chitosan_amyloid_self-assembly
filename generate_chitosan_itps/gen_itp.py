#!/usr/bin/python2.7

'''
    Generate itp files for each chitosan chain based on the charge profile in proto.txt
    Input: proto.txt must exist in current directory. See README.
    Output: 
        ./nglu/nglu_00**.itp = directory of itp files for chitosan chains
'''
import os
import numpy as np

if not os.path.isdir('./nglu/'):
	os.mkdir('./nglu')

with open('./chi_template.itp', 'r') as f:
	data = f.readlines()

charge = (np.loadtxt('./proto.txt')).T

for j in range(0, len(charge)):
	row = data[4].split()
	row[0] = 'NGLU%04d' % (j)
	data[4] = '%s %s\n' % (row[0], row[1])
	for i, ln in enumerate(range(10, 100, 3)):
		row = data[ln]
		row = row.split()
		if charge[j, i]==1:
			row[1] = 'Qd'
			row[4] = 'B3' 
			row[6] = 1.0
		elif charge[j, i]==0:
			row[1] = 'P5'
			row[4] = 'B3'
			row[6] = 0.0
		row = ' %-s %s  %-s %s %s %s %5.2f\n' % (row[0], row[1], row[2], row[3], row[4], row[5], row[6])
		data[ln] = row

	with open('./nglu/nglu_%04d.itp' % (j), 'w') as f:
		f.write(''.join(data))

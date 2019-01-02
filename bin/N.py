#!/usr/bin/env python3
#version 1.0

import os.path
import argparse
import csv

parser = argparse.ArgumentParser(description='Create a BED file of N blocks from FASTA files.')
parser.add_argument('-i','--infile', help='Comma separated input FASTA files.', action="store", type=str)
parser.add_argument('-o','--outfile', help='Output BED file.', action="store", type=str)
args = vars(parser.parse_args())

if not args['infile']:
    parser.print_help()
    raise SystemExit

if not args['outfile']:
    args['outfile'] = args['infile'] + '.N'

def read_fasta(fasta):
    name, seq = None, []
    for line in fasta:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def N(ID,sequence,writer):
	pos = 0
	switch = False
	N_start = ""
	N_end = ""
	N_switch = True
	written = True
	for nt in sequence:
		if nt == "N":
			written = False
			N_switch = True
			if not switch:
				N_start = pos
				N_end = pos + 1
				switch = True
			else:
				N_end = pos + 1
		else:
			switch = False
			if N_switch:
				if N_start != "":
					row = [ID,N_start,N_end]
					writer.writerow(row)
					written = True
				N_switch = False
		pos += 1
	if not written:
		row = [ID,N_start,N_end]
		writer.writerow(row)

inputfiles = args['infile'].split(',')
outfile = open(args['outfile'], 'w')
writer = csv.writer(outfile, delimiter='\t')

for fasta in inputfiles:
	handle = open(fasta, 'rU')
	for name, seq in read_fasta(handle):
		ID = name.replace(">","")
		sequence = seq.upper()
		N(ID,sequence,writer)

	handle.close()

outfile.close()

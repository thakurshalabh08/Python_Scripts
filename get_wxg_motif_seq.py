#!/usr/bin/env python

"""
Author:
Shalabh Thakur
Sycuro Lab, University of Calgary
02/01/2018

Description:
This is a script to extract sequences from a multi-fasta file based on list of sequence ids given as a file input

"""
import sys
import os
from Bio import SeqIO
import argparse
import re
from argparse import RawTextHelpFormatter


input_sequence_file=sys.argv[1]
output_sequence_file=sys.argv[2]



###### Read Sequence File #####
sequences_dict={}

output_seq=open(output_sequence_file,'w')

with open(input_sequence_file,"r") as in_seq_file:
	seq_records=SeqIO.parse(in_seq_file, "fasta")
	
	for record in seq_records:

		sequence=str(record.seq)

                if len(sequence)>=85 and len(sequence)>=135:

                	if re.search("W[A-Z]{1}G",sequence):

                		output_seq.write(">{}\n{}\n".format(record.id,sequence))
output_seq.close()

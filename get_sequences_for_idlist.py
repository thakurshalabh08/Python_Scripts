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
from argparse import RawTextHelpFormatter


gene_id_list=sys.argv[1]
input_sequence_file=sys.argv[2]
output_sequence_file=sys.argv[3]
add_annotation=sys.arg[4]

add_annotation=add_annotation.rstrip("\n")

###### Read Sequence File #####
sequences_dict={}

with open(input_sequence_file,"r") as in_seq_file:
	seq_records=SeqIO.parse(in_seq_file, "fasta")
 
	

	for record in seq_records:
		seq_record={}
		seq_record["id"]=record.id
		seq_record["desc"]=record.description
		seq_record["seq"]=record.seq

		sequences_dict[record.id]=seq_record





##### Read file with the list of sequence ids and write corresponding sequences in an output file ######

output_seq=open(output_sequence_file,'w')
seq_id_file=open(gene_id_list,"r")

for seq_id in seq_id_file:
	seq_id_trim=seq_id.rstrip('\n')

        if seq_id_trim in sequences_dict:

	    output_seq.write(">{} {} {}\n{}\n".format(seq_id_trim,sequences_dict[seq_id_trim]["desc"],add_annotation,sequences_dict[seq_id_trim]["seq"]))

seq_id_file.close()
output_seq.close()

			

















	









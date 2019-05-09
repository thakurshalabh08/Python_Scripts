#!/usr/bin/env python

"""
Author:
Shalabh Thakur
Sycuro Lab, University of Calgary
02/01/2018

Description:
This is a script to extract annotation description from a multi-fasta file and prepare a table

"""
import sys
import os
from Bio import SeqIO
import re
import argparse
from argparse import RawTextHelpFormatter

input_sequence_file=sys.argv[1]
output_file=sys.argv[2]


#output_table=open(output_file,'w')

###### Read Sequence File #####
sequences_dict={}

outfile=open(output_file,"w")

with open(input_sequence_file,"r") as in_seq_file:
	seq_records=SeqIO.parse(in_seq_file, "fasta")
	
	for record in seq_records:
		desc_line=record.description

		desc_column=re.split(r"[\[\]]",desc_line)

		#seq_id=desc_column[0]
		#db_xref=desc_column[1]
		#gene_name=desc_column[2]

		print desc_column

		if ' ' in desc_column:
			desc_column=filter(lambda a: a !=' ',desc_column)
		if '' in desc_column:
			desc_column=filter(lambda a: a !='',desc_column)
		
		print desc_column

		organism_name=desc_column.pop()
		family_desc=desc_column.pop()
		seq_desc=" ".join(desc_column)

		
		#### parse sequence description #######
		seq_desc_list=re.split(r"[\(\)]",seq_desc)

		if ' ' in seq_desc_list:
			seq_desc_list=filter(lambda a: a !=' ',seq_desc_list)
		if '' in seq_desc_list:
			seq_desc_list=filter(lambda a: a !='',seq_desc_list)

		seq_id=seq_desc_list.pop(0)
		dbxref=seq_desc_list.pop(0)
		gene_name=seq_desc_list.pop(0)
		gene_description=" ".join(seq_desc_list)

		#### parse family description ######

		family_desc_list=re.split(r"[\(\)]",family_desc)

		print family_desc

		if ' ' in family_desc_list:
			family_desc_list=filter(lambda a: a !=' ',family_desc_list)
		if '' in family_desc_list:
			family_desc_list=filter(lambda a: a !='',family_desc_list)

		family_id=family_desc_list.pop()
		family_desc=" ".join(family_desc_list)


		print "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(seq_id,dbxref,gene_name,gene_description,family_id,family_desc,organism_name)

		outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(seq_id,dbxref,gene_name,gene_description,family_id,family_desc,organism_name))




			

















	









#!/usr/bin/env python

"""
Author:
Shalabh Thakur
Sycuro Lab, University of Calgary
09/09/2018

Description:
This is a script to construct multiple alignment for clustered sequences

"""
import sys
import os
import glob
import argparse
import re
import shutil
import logging
import datetime
from argparse import RawTextHelpFormatter
from Bio import SeqIO
import subprocess
from collections import defaultdict


def main():

	"""
	Main Function
	"""
	now = datetime.datetime.now()
	
	c_args = parse_args(__file__)

	cluster_file=os.path.abspath(c_args['cluster_file'])
	sequence_file=os.path.abspath(c_args['sequence_file'])
	fasta_group_dir=os.path.abspath(c_args['fasta_group_dir'])
	align_group_dir=os.path.abspath(c_args['align_group_dir'])
	hmm_model_dir=os.path.abspath(c_args['hmm_model_dir'])


	cluster_dict=read_cluster_file(cluster_file)

	sequence_dict=read_sequence_file(sequence_file)

	create_fasta_cluster(cluster_dict,sequence_dict,fasta_group_dir)

	#create_align_cluster(cluster_dict,fasta_group_dir,align_group_dir)

	create_hmm_model(cluster_dict,align_group_dir,hmm_model_dir)



def parse_args(desc):

	"""
	Command line option parser
	
	:param desc: Short program description.
	:type desc: str
	
	:return arg: dict of command line arguments with key=command line
        argument and val=argument value
        
	:return_type: dict
	"""

	parser = argparse.ArgumentParser(description=desc, formatter_class=RawTextHelpFormatter)
	parser._optionals.title="Arguments"
	parser.add_argument("--cluster_file", help="Name and path of the tab-delimited cluster file", required=True)
	parser.add_argument("--sequence_file", help="Name of the sequence file", required=True)
	parser.add_argument("--fasta_group_dir", help="Name and path of the directory to store sequence files for each cluster in fasta format", required=True)
	parser.add_argument("--align_group_dir", help="Name and path of the directory to store multiple aligned sequence files for each cluster in fasta format", required=True)
	parser.add_argument("--hmm_model_dir", help="Name and path of the directory to store HMM models for multiple aligned sequence files for each cluster", required=True)
	
	args=parser.parse_args()

	return vars(args)


def read_cluster_file(cluster_file):

	"""
	Function to read sequence id for each cluster
	"""

	mycluster=open(cluster_file,"r")

	cluster_dict=defaultdict(list)

	for cluster_line in mycluster:

		cluster_line=cluster_line.rstrip("\n")
		cluster_row=cluster_line.split("\t")

		cluster_group_id=cluster_row.pop(0)
		cluster_group_id=re.sub(r"\:","",cluster_group_id)

		cluster_dict[cluster_group_id].append(cluster_row)

		#print cluster_group_id
		#print cluster_dict[cluster_group_id]
		
	return(cluster_dict)



def read_sequence_file(sequence_file):


	sequence_dict={}

	with open(sequence_file,"r") as in_seq_file:
		seq_records=SeqIO.parse(in_seq_file, "fasta")
	
		for record in seq_records:

			seq_id=record.id
			sequence=record.seq

			seq_id=re.sub(r"\(.+","",seq_id)

			sequence_dict[seq_id]=sequence

	return(sequence_dict)




def create_fasta_cluster(cluster_dict,sequence_dict,fasta_group_dir):

	
	if not os.path.exists(fasta_group_dir):
		os.makedirs(fasta_group_dir)

	
	for cluster_group_id in cluster_dict:

		print cluster_group_id

		cluster_fasta=open(os.path.join(fasta_group_dir,cluster_group_id+".fasta"),"w")

		cluster_seq_id_list=cluster_dict[cluster_group_id]

		for seq_id_list in cluster_seq_id_list:

			print seq_id_list

			for seq_id in seq_id_list:

				if sequence_dict[seq_id]:

					cluster_fasta.write(">{}\n{}\n".format(seq_id,sequence_dict[seq_id]))


def create_align_cluster(cluster_dict,fasta_group_dir,align_group_dir):

	
	if not os.path.exists(align_group_dir):
		os.makedirs(align_group_dir)
		os.makedirs(os.path.join(align_group_dir,"MULT_SEQ_ALIGN"))
		os.makedirs(os.path.join(align_group_dir,"SINGLETON"))


	for cluster_group_id in cluster_dict:

		print cluster_group_id

		cluster_fasta=os.path.join(fasta_group_dir,cluster_group_id+".fasta")
		
		num_cluster_seq=len(cluster_dict[cluster_group_id][0])

		print "NUM_SEQ:",num_cluster_seq

		if num_cluster_seq>1:

			os.system("muscle -in {} -out {}".format(cluster_fasta,os.path.join(align_group_dir,"MULT_SEQ_ALIGN",cluster_group_id+".aln.fasta")))

		else:

			os.system("cp {} {}".format(cluster_fasta,os.path.join(align_group_dir,"SINGLETON",cluster_group_id+".aln.fasta")))
			

def create_hmm_model(cluster_dict,align_group_dir,hmm_model_dir):

	if not os.path.exists(hmm_model_dir):
		os.makedirs(hmm_model_dir)

	for cluster_group_id in cluster_dict:

		print cluster_group_id
		
		num_cluster_seq=len(cluster_dict[cluster_group_id][0])

		if num_cluster_seq>1:

			align_fasta_file=os.path.join(align_group_dir,"MULT_SEQ_ALIGN",cluster_group_id+".aln.fasta")

			os.system("hmmbuild {} {}".format(os.path.join(hmm_model_dir,cluster_group_id+".hmm"),align_fasta_file))
	

if __name__ == '__main__':
	main()

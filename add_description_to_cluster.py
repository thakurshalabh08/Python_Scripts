#!/usr/bin/env python

"""
Author:
Shalabh Thakur
Sycuro Lab, University of Calgary
14/08/2018

Description:
This is a script to map seqeunce description to sequence cluster

"""
import sys
import os
import re
import argparse
from argparse import RawTextHelpFormatter


cluster_file_group=sys.argv[1]
sequence_description_file=sys.argv[2]
cluster_table=sys.argv[3]


sequence_dict={}

with open(sequence_description_file,"r") as desc_handel:

	for line in desc_handel:
		
		line=line.rstrip("\n")
		column=line.split("\t")

		seq_id=column[0]

		seq_id=seq_id.strip()

		sequence_dict[seq_id]=line



tab_cluster=open(cluster_table,"w")

with open(cluster_file_group,"r") as in_cluster:

	for cluster_line in in_cluster:

		cluster_line=cluster_line.rstrip("\n")

		cluster_column=cluster_line.split("\t")

		group_id=cluster_column.pop(0)
		group_id=re.sub(":","",group_id)

		for seq_id in cluster_column:

			tab_cluster.write("{}\t{}\n".format(group_id,sequence_dict[seq_id]))


		

		



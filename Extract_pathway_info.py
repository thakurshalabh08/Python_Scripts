#!/usr/bin/env python

"""
Author:
Shalabh Thakur
Sycuro Lab, University of Calgary
26/03/2018

Description:
This is a script to extract information from kegg pathway details
http://rest.kegg.jp/get/clo00010

"""
import sys
import os
import argparse
import re
import logging
from argparse import RawTextHelpFormatter


input_dir=sys.argv[1]
output_dir=sys.argv[2]


### Read Input Directory ####

for file_name in os.listdir(input_dir):

	input_file_path=os.path.join(input_dir,file_name)

	output_file_path=os.path.join(output_dir,file_name)
	
	pathway_info_file=open(input_file_path,"r")

	pathway_gene_file=open(output_file_path,"w")

	for line in pathway_info_file:

		line=str(line)
		if re.search("HMPREF",line):

			pathway_gene_file.write("{}\t{}".format(file_name,line))
               		

		







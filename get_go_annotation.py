#!/usr/bin/env python

"""
Author:
Shalabh Thakur
Sycuro Lab, University of Calgary
18/01/2018

Description:
This is a script to extract go annotation

"""
import sys
import os
import argparse
from argparse import RawTextHelpFormatter

input_file=sys.argv[1]
output_file=sys.argv[2]


annotation_file=open(input_file,"r")

result_file=open(output_file,'w')

for line in annotation_file:

    line=line.rstrip('\n')
    columns=line.split("\t")

    gene_id=""

    for index,value in enumerate(columns):

	if index==0:

           gene_id=value    
    
        elif index>=1:

           print "{}\t{}\n".format(gene_id,value)
           result_file.write("{}\t{}\n".format(gene_id,value))
    

    


#!/usr/bin/env python

"""
Author:
Shalabh Thakur
Sycuro Lab, University of Calgary
02/01/2018

Description:
This is a script to extract information from KEGG BlastKOALA result for KEGG BRITE.

A<ID1> <LEVEL1>
B  <ID2>  <LEVEL2>
C    <ID3>  <LEVEL3>
D     	<ID4>  <LEVEL4>

Output file format:

Gene_ID	KO_ID	KO_Description	Pathway_Name	Pathway_ID	Pathway_Class Pathway_Class_ID	Functional_Class Functional_Class_ID

"""
import sys
import os
from Bio import SeqIO
import argparse
import re
from argparse import RawTextHelpFormatter

### Read Input From CommandLine####
kegg_infile=sys.argv[1]
mapping_outfile=sys.argv[2]

kegg_file=open(kegg_infile,"r")

outfile=open(mapping_outfile,'w')

level_code=""
a_level_desc=""
b_level_desc=""
c_level_desc=""
d_level_desc=""
a_level_id=""
b_level_id=""
c_level_id=""
d_level_id=""

for line in kegg_file:

	line=line.rstrip('\n')

	if re.search("^#",line):
		continue
	
	if re.search("^A",line):
		match_pattern=re.match(r"(\w)(\d+)(\s+)(.+)",line)
		
		level_code=match_pattern.group(1)
		a_level_id=match_pattern.group(2)
		a_level_desc=match_pattern.group(4)

		if not a_level_desc:
			a_level_desc="unknown"

	if re.search("^B",line):
		match_pattern=re.match(r"(\w)(\s+)(\d+)(\s+)(.+)",line)

		level_code=match_pattern.group(1)
		b_level_id=match_pattern.group(3)
		b_level_desc=match_pattern.group(5)

		if not b_level_desc:
			b_level_desc="unknown"

	if re.search("^C",line):
		match_pattern=re.match(r"(\w)(\s+)(\d+)(\s+)(.+)",line)

		level_code=match_pattern.group(1)
		c_level_id=match_pattern.group(3)
		c_level_desc=match_pattern.group(5)

		if not c_level_desc:
			c_level_desc="unknown"

	if re.search("^D",line):

		match_pattern=re.match(r'(\w)(\s+)(\w+\|\w+\.\w+\.\w+\.\w+)(\s+)(.+)',line)

		level_code=match_pattern.group(1)
		d_level_id=match_pattern.group(3)
		ko_line=match_pattern.group(5)

		ko_line=re.sub(r'<a.*?>',"",ko_line)
		ko_line=re.sub(r'<\/a>',"",ko_line)

		ko_pattern=re.match(r'(\w+)(\s+)(.+)',ko_line)
		ko_id=ko_pattern.group(1)
		ko_desc=ko_pattern.group(3)

		print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(a_level_id,a_level_desc,b_level_id,b_level_desc,c_level_id,c_level_desc,d_level_id,ko_id,ko_desc)

		outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(a_level_id,a_level_desc,b_level_id,b_level_desc,c_level_id,c_level_desc,d_level_id,ko_id,ko_desc))


		
		
			
		




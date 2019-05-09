#!/usr/bin/env python

"""
Author:
Shalabh Thakur
Sycuro Lab, University of Calgary
02/01/2018

Description:
This is a script to extract information from EC to KO reverse link information from KEGG.
Input file format obtained from KEGG at the location: http://www.genome.jp/kegg-bin/get_htext?ko01000.keg

A<b>1. Oxidoreductases</b>
B  1.1  Acting on the CH-OH group of donors
C    1.1.1  With NAD+ or NADP+ as acceptor
D      1.1.1.1  alcohol dehydrogenase
E        K00001  E1.1.1.1, adh; alcohol dehydrogenase [EC:1.1.1.1]
E        K00121  frmA, ADH5, adhC; S-(hydroxymethyl)glutathione dehydrogenase / alcohol dehydrogenase [EC:1.1.1.284 1.1.1.1]

Output file format:
EC	Level_Code	Top_Level	Level_Desc				KO_ID	KO_Desc
1.-.-.-	A		-.-.-.-		Oxidoreductases				K00001  E1.1.1.1, adh; alcohol dehydrogenase
1.1.-.- B		1.-.-.-		Acting on the CH-OH group of donors 	K00001  E1.1.1.1, adh; alcohol dehydrogenase

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

keggec2ko=open(kegg_infile,"r")

outfile=open(mapping_outfile,'w')

a_level_ec=""
b_level_ec=""
c_level_ec=""
d_level_ec=""
e_level_ec=""
a_top_level=""
b_top_level=""
c_top_level=""
d_top_level=""
e_top_level=""
level_code=""
a_level_desc=""
b_level_desc=""
c_level_desc=""
d_level_desc=""
ko_id=""
ko_desc=""

for line in keggec2ko:

	line=line.rstrip('\n')
	
	if re.search("^A",line):
		line=re.sub("<b>","\t",line)
		line=re.sub("</b>","",line)
		line=re.sub("\s{1}","\t",line)
		columns=line.split("\t")
		
		level_code=columns[0]
		a_level_ec=columns[1]+"-.-.-"
		a_level_desc=columns[2]
		a_top_level="-.-.-.-"

		if not a_level_desc:
			a_level_desc="unknown"
		

		print "level_code:{}\na_level_ec:{}\na_top_level:{}\nlevel_desc:{}\n".format(level_code,a_level_ec,a_top_level,a_level_desc)

	if re.search("^B",line):
		line=re.sub("\s{2}","\t",line)
		columns=line.split("\t")

		level_code=columns[0]
		b_level_ec=columns[1]+".-.-"
		b_level_desc=columns[2]
		b_top_level=a_level_ec

		if not b_level_desc:
			b_level_desc="unknown"

		print "level_code:{}\nb_level_ec:{}\nb_top_level:{}\nlevel_desc:{}\n".format(level_code,b_level_ec,b_top_level,b_level_desc)

	if re.search("^C",line):
		line=re.sub("\s{4}","\t",line)
		line=re.sub("\s{2}","\t",line)
		columns=line.split("\t")

		level_code=columns[0]
		c_level_ec=columns[1]+".-"
		c_level_desc=columns[2]
		c_top_level=b_level_ec

		if not c_level_desc:
			c_level_desc="unknown"

		print "level_code:{}\nc_level_ec:{}\nc_top_level:{}\nlevel_desc:{}\n".format(level_code,c_level_ec,c_top_level,c_level_desc)

	if re.search("^D",line):
		line=re.sub("\s{6}","\t",line)
		line=re.sub("\s{2}","\t",line)
		columns=line.split("\t")

		level_code=columns[0]
		d_level_ec=columns[1]
		d_level_desc=columns[2]
		d_top_level=c_level_ec

		if not d_level_desc:
			d_level_desc="unknown"

		print "level_code:{}\nd_level_ec:{}\nd_top_level:{}\nlevel_desc:{}\n".format(level_code,d_level_ec,d_top_level,d_level_desc)

	if re.search("^E",line):
		line=re.sub("\s{8}","\t",line)
		line=re.sub("\s{2}","\t",line)
		columns=line.split("\t")

		level_code=columns[0]
		ko_id=columns[1]
		ko_desc=columns[2]
		e_top_level=d_level_ec

		print "level_code:{}\nko_id:{}\ne_top_level:{}\nko_desc:{}\n".format(level_code,ko_id,e_top_level,ko_desc)

		outfile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format("A",a_level_ec,a_top_level,a_level_desc,ko_id,ko_desc))
		outfile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format("B",b_level_ec,b_top_level,b_level_desc,ko_id,ko_desc))
		outfile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format("C",c_level_ec,c_top_level,c_level_desc,ko_id,ko_desc))
		outfile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format("D",d_level_ec,d_top_level,d_level_desc,ko_id,ko_desc))

		
		
			
		




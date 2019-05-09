#!/usr/bin/env python

"""
Author:
Shalabh Thakur
Sycuro Lab, University of Calgary
08/11/2018

Description:
This is a script to combine multi-hsp hit from tabular BLAST result

"""
import sys
import os
import glob
import argparse
import re
import shutil
import copy
from argparse import RawTextHelpFormatter
from Bio import SeqIO
from collections import defaultdict
from collections import OrderedDict


def main():

	"""
	Main Function
	"""
	
	c_args = parse_args(__file__)

	infile=os.path.abspath(c_args['in'])
	outfile=os.path.abspath(c_args['out'])

	header=list()

	if not (c_args['column_names'] and c_args['header']):
		
		header=["QSEQID","SSEQID","PIDENT","LENGTH","MISMATCH","GAPOPEN","QSTART","QEND","SSTART","SEND","EVALUE","BITSCORE"]

	elif (c_args['column_names']):
		c_args['column_names']=c_args['column_names'].rstrip("\n")
		header=c_args['column_names'].split(" ")
		header=map(lambda x:x.upper(),header)

	#### Parse Blast result ###
	formatted_blast_table,header=parse_blast(infile,header,c_args)

	####Fix multiple HSPs #####
	fixed_hsp_blast_table=fix_multi_hsp(formatted_blast_table)

	####Convert HSPs to HIT ######
	formatted_hit_blast_table,header=convert_hsp2hit(fixed_hsp_blast_table,header)

	#### Write Formatted Blast Table to the output file ####
	write_blast_table(formatted_hit_blast_table,header,outfile)


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
	parser.add_argument("--in", help="Name and path of the tab-delimited blast result file", required=True)
	parser.add_argument("--out", help="Name and path of the tab-delimited parsed blast result file", required=True)
	columns=parser.add_mutually_exclusive_group()
	columns.add_argument("--column_names", help="column order and names as provided by -outfmt parameter in blast. If not given, use default column order and names for tabular output", default=False)
	columns.add_argument("--header", help="Read header column names from top row of the blast output file", action="store_true", default=True)
	
	args=parser.parse_args()

	return vars(args)


def parse_blast(infile,header,c_args):

	"""
	Reads in tabular BLAST results, calculate percent identity for full length query sequence and create a dictionary for each row and column

	param blast_file: full path and name of blast output file

	return blast_table: dict of blast columns sorted by query_id/subject_id pair
		
	The method/formula of calculating full length pident was taken from https://github.com/DenefLab/EMIRGE/blob/master/calc_full_length_pident.R
	"""

	#### Read Blast Raw output file #####
	blast_result=open(infile,"r")

	formatted_blast_table=defaultdict(list)

	header2=list()
	add_nident=False
	add_qcovs=False

	for index, blast_row in enumerate(blast_result):

		blast_row=blast_row.rstrip("\n")
		blast_column=blast_row.split("\t")

		if(c_args['header'] and index==0):
			header=blast_column
			header=map(lambda x:x.upper(),header)
			

			if(len(header)!=len(blast_column)):
				print "Number columns in header and blast file does not match\n"
				print "{} columns found in blast result file but {} given in header\n".format(len(blast_column),len(header))
				sys.exit()
			continue

		hsp_hit={}
		header2=header
		
		for field_index, field_name in enumerate(header):
			
			hsp_hit[field_name]=blast_column[field_index]

		query_subject_pair="{}:{}".format(hsp_hit['QSEQID'],hsp_hit['SSEQID'])

		### check if NIDENT columns and QCOVS column present in raw blast output, if not add to the parsed result file ###
		if not 'NIDENT' in hsp_hit:
			hsp_hit['NIDENT']=((float(hsp_hit['PIDENT']) * int(hsp_hit['LENGTH']))/100)+1
			
			if add_nident==False:
				add_nident=True
			

		if 'QLEN' in hsp_hit:
			if not 'QCOVS' in hsp_hit:
				hsp_hit['QCOVS']=((((int(hsp_hit['QEND']) - int(hsp_hit['QSTART']))+1)/int(hsp_hit['QLEN']))*100)

				if add_qcovs==False:
					add_qcovs=True

		
		#print query_subject_pair,":",hsp_hit['PIDENT'],":",hsp_hit['LENGTH'],":",int(hsp_hit['NIDENT']),":",float(hsp_hit['QCOVS']),"\n"

		formatted_blast_table[query_subject_pair].append(hsp_hit)

	if add_nident==True:
		header2.append("NIDENT")
	if add_qcovs==True:
		header2.append("QCOVS")
	
	print header2	
			
	return(formatted_blast_table,header2)


def fix_multi_hsp(formatted_blast_table):

	"""
	Reads formatted blast result, find compatible hsps and re-calculate identities and bitscore for multiple hsp hits
	"""
	fixed_hsp_blast_table={}

	for query_subject_pair in formatted_blast_table:

		#sorted_hsp_list=sorted(formatted_blast_table[query_subject_pair], key=lambda k: k['qstart'], reverse=False)
	
		##### find compatible hsps with highest bitscore that do not cross and overlap each other ##
		compatible_hsp=list()

		compatible_hsp = compute_best_hsp(formatted_blast_table[query_subject_pair],hsp_list=list())	

		fixed_hsp_blast_table[query_subject_pair]=compatible_hsp

	return(fixed_hsp_blast_table)


# Return summed bitscores of best combo. To figure out which HSPs
# correspond to this optimal score, we must also store whether the optimal
# solution includes the last HSP in `hsps` at each step, then recursively
# reconstruct the optimal solution after the algorithm finishes. See the
# linked code for details.

def compute_best_hsp(hsps,hsp_list=list()):
  	"""
  	Memoize the function so that, if we have already computed the solution
  	for `hsps`, just fetch and return that value immediately without
  	performing below computations.

	Based on code found here:
    	http://jeff.wintersinger.org/posts/2014/07/designing-an-algorithm-to-
    	compute-the-optimal-set-of-blast-hits/
  	"""
	if len(hsps) == 0:
		return hsp_list
	if len(hsps) == 1:
    		# Trivial solution: one HSP, so optimal solution is just itself.
		hsp_list.append(hsps[0])
		return hsp_list

  	# Last HSP
  	last_hsp = hsps[-1]
  	# All HSPs except last
  	previous = hsps[:-1]

  	# Find subset of HSPs in `previous` that don't overlap `last_hsp`.
  	compatible = find_compatible(last_hsp, previous)

	without_list = copy.deepcopy(hsp_list)
    	with_list = copy.deepcopy(hsp_list)
    	with_list.append(last_hsp)

  	best_without_last = compute_best_hsp(previous,without_list)
	best_with_last    = compute_best_hsp(compatible, with_list)

	compatible_hsp=max([best_without_last, best_with_last],key=lambda (x): sum([float(y['BITSCORE']) for y in x]))
	
	return compatible_hsp

##### Find Compatible HSP #####
def find_compatible(target, hsps):
  	'''Find hsps amongst `hsps` that don't overlap or cross `target`. They
  	may not be mutually compatible, as they are only guaranteed to be
  	compatible with `target`.'''
	compatible = []

  	for hsp in hsps:
    	# Don't define target as being compatible with itself.
    		if hsp == target:
      			continue

   		if target['QSTART'] <= hsp['QSTART']:
      			first, second = target, hsp
    		else:
      			first, second = hsp, target

    		overlap = (second['QSTART'] <= first['QEND'] or
              		second['SSTART'] <= first['SEND'])
    		if not overlap:
      			compatible.append(hsp)

	return compatible

##### convert HSP to HIT format ####

def convert_hsp2hit(fixed_hsp_blast_table,header):

	formatted_hit_blast_table={}

	print "Converting HSPs into Hit format\n"

	for query_subject_pair in fixed_hsp_blast_table:

		multi_hsp_hit=list()
		hsp_hit_list=fixed_hsp_blast_table[query_subject_pair]

		#print hsp_hit_list

		total_pident=hsp_hit_list[0]['PIDENT']
		total_bitscore=hsp_hit_list[0]['BITSCORE']
		total_qcov=None
		total_query_len=None
		total_subject_len=None
		num_hsp=len(hsp_hit_list)

		if 'QCOVS' in hsp_hit_list[0]:
			total_qcov=hsp_hit_list[0]['QCOVS']

		if 'QLEN' in hsp_hit_list[0]:
			total_query_len=hsp_hit_list[0]['QLEN']

		if 'SLEN' in hsp_hit_list[0]:
			total_subject_len=hsp_hit_list[0]['SLEN']

		if num_hsp>1:

			sum_num_ident=sum(int(k['NIDENT']) for k in hsp_hit_list)
			sum_aln_length=sum(int(k['LENGTH']) for k in hsp_hit_list)
			total_bitscore=sum(float(k['BITSCORE']) for k in hsp_hit_list)
			total_pident=round(float(int(sum_num_ident)/int(sum_aln_length))*100,2)
			
			if not total_query_len==None:
				total_qcov=round(float(int(sum_aln_length) / int(total_query_len))*100,2)
		
 		for fieldname in header:

			if re.search(r"^(\d+)(\.)(\d+)$",str(hsp_hit_list[0][fieldname])):
				multi_hsp_hit.append('_'.join(map(str,(hsp_hit[fieldname] for hsp_hit in  hsp_hit_list))))
			elif re.search(r"^(\d+)$",str(hsp_hit_list[0][fieldname])):
				multi_hsp_hit.append('_'.join(map(str,(hsp_hit[fieldname] for hsp_hit in  hsp_hit_list))))
			else:
				multi_hsp_hit.append(hsp_hit_list[0][fieldname])
		
		multi_hsp_hit.append(total_bitscore)
		multi_hsp_hit.append(total_pident)
		multi_hsp_hit.append(total_qcov)		

		formatted_hit_blast_table[query_subject_pair]=multi_hsp_hit

	header.append("TOTAL_BITSCORE")
	header.append("TOTAL_PIDENT")
	header.append("TOTAL_QCOVS")

	return(formatted_hit_blast_table,header)


def write_blast_table(formatted_hit_blast_table,header,output_file):

	"""
	Reads formatted blast result dictonary and write it to output file
	"""	
	blast_report=open(output_file,"w")

	##### Write Column header ######

	#print "[PARSER_STEP 3]: Writing Formatted Blast Results to the output file {}\n".format(output_file)

	blast_report.write("\t".join(header))
	blast_report.write("\n")

	for query_subject_pair in formatted_hit_blast_table:

			hit_column=formatted_hit_blast_table[query_subject_pair]

			column_value="\t".join(map(str,(hit_column)))

			#print column_value

			blast_report.write("{}\n".format(column_value))



if __name__ == '__main__':
	main()

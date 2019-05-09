#!/usr/bin/env python

"""
Author:
Shalabh Thakur
Sycuro Lab, University of Calgary
14/08/2018

Description:
This is a script to add group id to sequence cluster

"""
import sys
import os
import re
import argparse
from argparse import RawTextHelpFormatter


cluster_file_no_group=sys.argv[1]
cluster_file_group=sys.argv[2]
group_prefix=sys.argv[3]

cluster_num=10000

out_cluster=open(cluster_file_group,"w")

with open(cluster_file_no_group,"r") as in_cluster:

	for cluster_line in in_cluster:

		group_id=group_prefix+"_"+str(cluster_num)

		out_cluster.write("{}:\t{}".format(group_id,cluster_line))

		cluster_num=cluster_num + 1

		

		



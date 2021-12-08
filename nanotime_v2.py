#!/usr/bin/env python
#coding=utf-8

from argparse import ArgumentParser
import os, gzip
# import numpy as np
# import pandas as pd
# from dateutil.parser import parse as dparse
# from dateutil.relativedelta import relativedelta

"""

In order to obtain the pool start time, the nanotime.py firstly scans all fastq/signal files, 
finds the earliest start_time among files.


"""

def get_args():
    parser = ArgumentParser(description="Nanotime v2. Nanopore reads extraction based on the start time of sequencing; The process on fast5(signal) file is still under development")
    parser.add_argument('-q', "--fqz", help="the path of compressed fastq file (e.g. *.fastq.gz)", required=True, type = str)
    parser.add_argument('-s', "--ses", help="the path to sequencing_summary.txt", required=True, type = str)
    parser.add_argument('-l', "--len", help="the length of run of sequencing (e.g. 1)", required=True, default=1, type=str)
    parser.add_argument('-o', "--out", help="the dir of output", required=True, type = 'str')
    return parser.parse_args()


def main():
    input_args = get_args()
    downsampling_fastq_creator(input_args)
    return 


def downsampling_fastq_creator(input_args):
    readids_dict = readids_dict_collector(input_args.ses) # {readids, }

    fq_dict = {}
    with gzip.open(input_args.fqz, 'rb') as file:
        for line in file:
            if line.startswith('@'):
                readid = line[1:].split(" ")[0]
                if readids_dict[readid] <= float(input_args.len) * 3600:  # convert hours to seconds
                    qname  = line
                    fq_dict[qname] = qname
                else:
                    qname  = None
            else:
                if qname is None:
                    continue
                else:
                    fq_dict[qname] += line.lstrip()
    fastq_saver(fq_dict, input_args)
    return 


def readids_dict_collector(summary_path):
    readids_time_dict = {}
    with open(summary_path, 'r') as file:
        for lcn, line in enumerate(file):
            line_lst = line.strip().split()

            if lcn == 0:
                for ind, i in enumerate(line_lst):
                    if i == 'start_time':
                        start_time_col = ind
                    elif i == 'read_id':
                        read_id_col    = ind
                    else:
                        next
            else:
                start_time = line_lst[start_time_col]
                read_id    = line_lst[read_id_col]
                readids_time_dict[read_id] = float(start_time)

    return readids_time_dict

def fastq_saver(fq_dict, input_args):
    output_dir = os.path.join(input_args.out, input_args.len)
    if not os.path.exists( output_dir ):
        os.makedirs(output_dir)
    output_name = os.path.basename(input_args.fqz)
    output_path = os.path.join( output_dir,  output_name)
    output = gzip.open( output_path, 'wb' )
    for readid in fq_dict:
        output.write(fq_dict[readid])
    output.close()
    return 


if __name__ == "__main__":
    main()
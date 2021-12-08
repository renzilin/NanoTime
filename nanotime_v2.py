#!/usr/bin/env python
#coding=utf-8

from argparse import ArgumentParser
import os, gzip
import numpy as np
import pandas as pd
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
    pool_start_time, readids_dict = readids_dict_collector(input_args.ses) # {readids, }

    fq_dict = {}
    with gzip.open(input_args.fqz, 'rb') as file:
        for line in file:
            if line.startswith('@'):
                readid = line[1:].split(" ")[0]
                if readids_dict[readid] - pool_start_time <= float(input_args.len) * 3600:  # convert hours to seconds
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
    dat = pd.read_csv(summary_path, sep="\t")
    pool_start_time = np.min( dat['start_time'].values )
    readids_time_dict = dat.loc[:, ["read_id", "start_time"]].set_index("read_id", )
    readids_time_dict = readids_time_dict.to_dict()['start_time']
    return pool_start_time, readids_time_dict

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











def fq_time_reader(fq_path):
    pool_start_time = None

    for record in SeqIO.parse( gzopen(fq_path, 'rt'), "fastq" ):
        time = dparse([i for i in record.description.split() if i.startswith('start_time')][0].split('=')[1])
        if pool_start_time == None:
            pool_start_time = time
            continue

        if time < pool_start_time:
            pool_start_time = time     
    return pool_start_time



def find_pool_start_time(pool_dir, ncpus=1):
    
    file_paths = pathlib.Path(pool_dir).glob("*.gz")
    
    
    if ncpus == 1:
        pbar = tqdm( list( file_paths ) )
        pool_start_time = None
        for char in pbar:
            pbar.set_description( "Processing %s" % char.name )

            time = fq_time_reader(char)        
            if pool_start_time == None:
                pool_start_time = time
                continue

            if time < pool_start_time:
                pool_start_time = time
                
    
    else:
        
        ## create multiprocessing
        pool     = mp.Pool(processes = ncpus)
        pool_lst = [ pool.apply_async( func=fq_time_reader, \
                                      args=(str(fq_path), ) ) for fq_path in file_paths]

        pool.close()
        pool.join()
        
        ## result fetch
        pool_start_time = None
        for time in pool_lst:
            time_get = time.get()
            
            if time_get != None:
                pool_start_time = time_get
                continue

            if time_get < pool_start_time:
                pool_start_time = time_get

    print( 'the start time of the pool is %s' % (pool_start_time) )            
    return pool_start_time

    
"""
outpath = 'time_intervval_fastqs/time_interval_0-XXhrs/barcodeXX.fq.gz'

"""


def generate_time_interval_fastq(fq_path, pool_start_time, time_interval = 1):
    
    def filter_time(descr, tfrom, tto):
        time = dparse([i for i in descr.split() if i.startswith('start_time')][0].split('=')[1])
        if tfrom <= time <= tto:
            return True
        else:
            return False
    
    
    outpath         = 'time_interval_fastqs/time_interval_0-%shrs' % time_interval
    pathlib.Path(outpath).mkdir(parents=True, exist_ok=True)
    
    
    pool_end_time   = pool_start_time - relativedelta( hours = -1 * time_interval )
    
    outfile = '%s/%s' % ( outpath, pathlib.Path(fq_path).name )
    fastq_out = gzopen(outfile, "wt")
    for record in SeqIO.parse( gzopen(fq_path, 'rt'), format="fastq" ):
        if filter_time(record.description, pool_start_time, pool_end_time ):
            # print(record.format("fastq"), end="\n")
            # break
            fastq_out.write(record.format("fastq"))
    
    fastq_out.close()
    
    return 1


def get_args():
    parser = ArgumentParser(description="Nanotime v0.1. Nanopore reads extraction based on the start time of sequencing; The process on fast5(signal) file is still under development")
    parser.add_argument("--fq",       help="the path of compressed fastq file (e.g. *.fq.gz or *.fastq.gz)", required=True, metavar = 'string')
    parser.add_argument("--pool_dir", help="The directory of pool where fastq files were stored (e.g. ./test/20210601-pool1)", required=True, metavar = 'string')
    parser.add_argument("--duration", help="The time duration (hour) of sequencing (e.g. 1)", required=False, default=1, metavar = 'INT')
    parser.add_argument("--ncpus",     help="The number of cores used in scan process for parallel acceleration (e.g. 3)", required=False, default=1, metavar = 'INT')

    return parser.parse_args()



if __name__ == '__main__':
    main()
    
#!/usr/bin/env python
import os
import sys
import argparse
import time
import gzip
import random as r
import itertools as it

des = """Paired or Single Read FASTQ Downsampler v0.1.0, 11/11/2017-11/23/2018 Timothy James Becker"""
parser = argparse.ArgumentParser(description=des)
parser.add_argument('-i', '--in_path',type=str,help='fastq input file(s)\t[None]')
parser.add_argument('-o', '--out_path',type=str,help='fastq output file(s)\t[None]')
parser.add_argument('-f', '--fraction',type=float,help='fraction [0.0,1.0] of the reads to sample\t[None]')
des = """comma seperated pattern to search for if paired reads are to be sampled\nthis pattern is concatenated to the in_path value\t[None]"""
parser.add_argument('-p','--paired_pattern',type=str,help=des)
args = parser.parse_args()

if args.in_path is not None:
    in_path = args.in_path
else:
    print('no in_path argument provided')
    raise IOError
if args.out_path is not None:
    out_path = args.out_path
    if not os.path.exists(out_path):
        os.mkdir(out_path)
else:
    raise IOError
if args.fraction is not None:
    fraction = args.fraction
    if fraction >= 1.0:
        print('fraction can not be set greater than 1.0')
        p = 1.0
    if fraction <= 0.0:
        print('fraction can not be set less than 0.0')
        p = 0.0
    if fraction > 0.0 and fraction < 1.0:
        p = fraction
if args.paired_pattern is not None:
    paired = True
    r1,r2 = args.paired_pattern.split(',')
    if not os.path.exists(in_path+r1) or not os.path.exists(in_path+r2):
        print('in_path and paired pattern not pointing to a valid file')
        raise IOError
else:
    paired = False
    r1,r2 = in_path,None

if __name__ == '__main__':
    x,y = 0,0
    start = time.time()
    if r2 is None:
        if r1.endswith('.gz'):
            with gzip.open(r1,'r') as q_in_1:
                with gzip.open(out_path,'w') as q_out_1:
                    for r1l1,r1l2,r1l3,r1l4 in it.izip_longest(*[q_in_1]*4):
                        y+=1
                        if r.random.random()<=p:
                            x+=1
                            q_out_1.write(''.join([r1l1,r1l2,r1l3,r1l4]))
        else:
            with open(in_path,'r') as q_in_1:
                with open(out_path,'w') as q_out_1:
                    for r1l1,r1l2,r1l3,r1l4 in it.izip_longest(*[q_in_1]*4):
                        y+=1
                        if r.random.random()<=p:
                            x+=1
                            q_out_1.write(''.join([r1l1,r1l2,r1l3,r1l4]))
    else:
        if r1.endswith('.gz') and r2.endswith('.gz'):
            with gzip.open(in_path+r1,'r') as q_in_1:
                with gzip.open(in_path+r2,'r') as q_in_2:
                    with gzip.open(out_path+r1,'w') as q_out_1:
                        with gzip.open(out_path+r2,'w') as q_out_2:
                            for r1l1,r2l1,r1l2,r2l2,r1l3,r2l3,r1l4,r2l4 in it.izip(*[q_in_1,q_in_2]*4):
                                y += 1
                                if r.random()<=p:
                                    x += 1
                                    q_out_1.write(''.join([r1l1,r1l2,r1l3,r1l4]))
                                    q_out_2.write(''.join([r2l1,r2l2,r2l3,r2l4]))
        else:
            with open(in_path+r1,'r') as q_in_1:
                with open(in_path+r2,'r') as q_in_2:
                    with open(out_path+r1,'w') as q_out_1:
                        with open(out_path+r2,'w') as q_out_2:
                            for r1l1,r2l1,r1l2,r2l2,r1l3,r2l3,r1l4,r2l4 in it.izip(*[q_in_1,q_in_2]*4):
                                y += 1
                                if r.random()<=p:
                                    x += 1
                                    q_out_1.write(''.join([r1l1,r1l2,r1l3,r1l4]))
                                    q_out_2.write(''.join([r2l1,r2l2,r2l3,r2l4]))
    stop = time.time()
    print('sampling results= %s reads sampled from %s total reads in %s sec'%(x,y,round(stop-start,2)))

#Timothy James Becker, PhD candidate, UCONN 01/10/2017-01/06/2018
#modifying this library to have the ability to skip the NNNNN regions of the regnome
#could also get into the alignibilty of a given read length onto a seq like genome strip...

import os
import json
try:
    import cPickle as pickle
except Exception as E:
    import pickle
    pass
import numpy as np
import pysam

#py_obj is the obj with a single name key
def dump_state(py_obj,name,path):
    if not os.path.exists(path): os.mkdir(path)
    dump_file = '%s.pickle'%name
    print('##DUMPING## %s to %s'%(name,path+'/'+dump_file))
    with open(path+'/'+dump_file, 'wb') as f:
        pickle.dump(py_obj,f,protocol=pickle.HIGHEST_PROTOCOL)
        return True
    return False

#takes in the multi-chrom fasta file and rads it by seq/chrom
#reads each sequence and then writes that portion to the chrom_fasta_dir
#using the chrom_base='' seq_name.fa
#write_fasta_by_chrom(path,path[0:path.rfind('/')],'')

def read_fasta_substring(fasta_path,chrom,pos,end):
    ss = ''
    with pysam.FastaFile(fasta_path) as f:
        ss = f.fetch(chrom)
    return ss[pos:end]
    
def read_fasta_chrom(fasta_path,chrom):
    ss = ''
    with pysam.FastaFile(fasta_path) as f:
        ss = f.fetch(chrom)
    return ss

#try this as alternate to reading all seqs at once for TRA
def read_fasta_chrom_pos(fasta_path,chrom,start,stop):
    ss =''
    with pysam.FastaFile(fasta_path) as f:
        ss = f.fetch(chrom,start,stop)
    return ss

def read_fasta(fasta_path):
    ss = {}
    with pysam.FastaFile(fasta_path) as f:
        names = f.references
        for chrom in names:
            ss[chrom] = f.fetch(chrom)
    return ss

def get_fasta_seq_names(fasta_path):
    ss = []
    with pysam.FastaFile(fasta_path) as f:
        ss = list(f.references)
    return ss

def get_fasta_seq_lens(fasta_path):
    ss = []
    with pysam.FastaFile(fasta_path) as f:
        ss = list(f.lengths)
    return ss

def get_fasta_seq_names_lens(fasta_path):
    ss = {}
    with pysam.FastaFile(fasta_path) as f:
        names = f.references
        lens  = f.lengths
    for i in range(len(names)):
        ss[names[i]]=lens[i]
    return ss   
    
def write_fasta(seqs,fasta_path,index=True):
    with open(fasta_path,'w') as fasta:
        for k in seqs:
            fasta.write('\n'.join(['>%s'%k]+[seqs[k][i:(i+80)] for i in range(0,len(seqs[k]),80)]+['\n']))
    if index: pysam.faidx(fasta_path) #reindex
    return True    
    
#ss is a HTSeq Sequence list?   
def write_fasta_by_chrom(seqs, chrom_fasta_dir, chrom_base=''):
    names = []
    for k in seqs:
        name = chrom_fasta_dir+'/'+chrom_base+k+'.fa'
        names += [name]
        with open(name, 'w') as fasta:
            fasta.write('\n'.join(['>%s'%k]+[seqs[k][i:(i+80)] for i in range(0,len(seqs[k]),80)]+['\n']))
    return names

def write_fasta_mask(M,json_path):
    with open(json_path,'w') as f:
        json.dump(M,f)
    return True

#compute an expectation given randomly distributed short reads for the RD windows (hist bins)
def expected_window(depth=20,length=100,target=100):
    return int(round((1.0*target/depth)*2.0*length,0))           

def get_coordinate_offsets(json_name):
    path = os.path.dirname(os.path.abspath(__file__))+'/data/'
    info,O = {},{}
    with open(path+json_name,'r') as f:
        info = json.load(f)
    for i in info:
        O[str(i)] = info[i]
    return O

def write_coordinate_offsets(fasta_path,json_path):
    #read in the fasta lengths and then sort them by length into a json offset map
    L = get_fasta_seq_names_lens(fasta_path)
    l_i = list(np.argsort([L[k] for k in L]))[::-1] #sort by max length
    O,offset = {},0 #starts at zero
    for i in l_i:
        O[L.keys()[i]] = offset
        offset = offset+L[L.keys()[i]] #old + new + 1
    with open(json_path,'w') as f:
        json.dump(O,f)
    return True

def get_chrom_dict(json_name):
    path = os.path.dirname(os.path.abspath(__file__))+'/data/'
    info,I = {},{}
    with open(path+json_name,'r') as f:
        info = json.load(f)
    for i in info:
        I[str(i)] = int(info[i])
    return I
    
#input is json data store for svmask regions and  offset map O
#output is a sorted by ref start pos list of list to filter on
def get_mask_regions(json_name,O):
    path = os.path.dirname(os.path.abspath(__file__))+'/data/'
    M = {}
    with open(path+json_name,'r') as f:
        M = json.load(f) #load the svmask per sequence
    N = []    
    for k in M:
        for i in range(len(M[k])):
            N += [[np.uint32(O[k])+np.uint32(M[k][i][0]),np.uint32(O[k])+np.uint32(M[k][i][1])]] #apply offsets
    N = sorted(N,key=lambda x: x[0])
    return N

def write_mask_regions(json_name):
    return True

def bed_mask_to_json_mask(bed_path,json_path):
    bed_data = []
    with open(bed_path, 'r') as f:
        bed_data = [i.split('\t') for i in f.read().split('\n')]
    data = {}
    for row in bed_data: #chrom,start,stop
        if len(row)==3:
            if data.has_key(row[0]): data[row[0]] += [[int(row[1]),int(row[2])]]
            else:                    data[row[0]]  = [[int(row[1]),int(row[2])]]
    for k in data:
        data[k] = sorted(data[k], key = lambda x: x[0])
    #no coordinate sort per contig
    with open(json_path,'w') as f:
        json.dump(data,f)
    return True
    
def get_offsets(chroms,order):
    x,m = 0,{k:0 for k in chroms}
    for k in order:
        m[k] = x
        x += chroms[k]
    return m

def get_reverse_complement(seq):
    m = {'A':'T','a':'t','T':'A','t':'a','C':'G','c':'g','G':'C','g':'c','N':'N'}
    return ''.join([m[k] for k in seq][::-1])

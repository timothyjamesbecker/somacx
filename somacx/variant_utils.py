#Timothy James Becker, PhD candidate, UCONN 01/10/2017-03/20/2020
import os
import sys
import re
import copy
import datetime
import itertools as it
import utils
import time
import gzip
import json
import hashlib
try:
    import cPickle as pickle
except Exception as E:
    import pickle
    pass
import numpy as np

import utils
import somacx.read_utils as ru

class VariantCall:
    def __init__(self,chrom=None,pos=None,identifier=None,
                 ref=None,alt=None,qual=None,flter=None,
                 info=None,frmat=None,end=None):
        self.chrom      = chrom
        self.pos        = pos
        self.identifier = identifier
        self.ref        = ref
        self.alt        = alt
        self.qual       = qual
        self.flter      = flter
        self.info       = info
        self.frmat      = frmat
        if self.info is not None:
            self.end    = get_info_end(self.info)
        else:
            self.end = self.pos+abs(len(self.ref)-len(self.alt))

    def __enter__(self):
        return self

    def __del__(self):
        return 0

    def __exit__(self,type,value,traceback):
        return 0

    def to_dict(self):
        return {"chrom":self.chrom,"pos":self.pos,"identifier":self.identifier,
                "ref":self.ref,"alt":self.alt,"qual":self.qual,"flter":self.flter,
                "info":self.info,"frmat":self.frmat,"end":self.end}

    def from_dict(self,D):
        self.chrom      = D['chrom']
        self.pos        = D['pos']
        self.identifier = D['identifier']
        self.ref        = D['ref']
        self.alt        = D['alt']
        self.qual       = D['qual']
        self.flter      = D['flter']
        self.info       = D['info']
        self.frmat      = D['frmat']
        self.end        = D['end']

def get_local_path(local_path=''):
    path = os.path.dirname(os.path.abspath(__file__)) + '/data/'
    return path + local_path

#json hook for conversion to str types as needed by fasta reader libs
if sys.version_info.major<3:
    def str_hook(obj):
        return {k.encode('utf-8') if isinstance(k,unicode) else k:
                v.encode('utf-8') if isinstance(v, unicode) else v for k,v in obj}
else:
    def str_hook(obj):
        return {k:v for k,v in obj}

#load in a JSON mutation map such as g1kp3_mut_map.json, somatic_mut_map.json,ect...
# mut_p = {layer(sorted integer 1,2,3,...):
#          type('SUB','INS','DEL',DUP',...:
#          'l:n':{size_center:rate, },'h':0.5
#          'TYPE':['PERFECT','COMPLEX'],TP:[0.75,0.25]
# defaults to g1kp3_mut_map for germline and somatic_mut_map for somatic
#converts all the json types to the native python used for distribution calculatoins
def read_json_mut_map(json_path):
    mut_map = {}
    if json_path.endswith('.gz'):
        with gzip.GzipFile(json_path,'rb') as f:
            raw_map = json.load(f)
    else:
        with open(json_path,'r') as f:
            raw_map = json.load(f)
    for l in raw_map:
        mut_map[int(l)] = {}
        for t in raw_map[l]:
            mut_map[int(l)][str(t)] = {}#'l:n' and 'h' are required
            mut_map[int(l)][str(t)]['l:n'] = {}
            for b in raw_map[l][t]['l:n']:
                mut_map[int(l)][t]['l:n'][float(b)] = raw_map[l][t]['l:n'][b]
            if 's:p' in raw_map[l][t]:
                mut_map[int(l)][str(t)]['s:p'] = {}
                for b in raw_map[l][t]['s:p']:
                    mut_map[int(l)][str(t)]['s:p'][float(b)] = raw_map[l][t]['s:p'][b]
            mut_map[int(l)][str(t)]['h'] = raw_map[l][t]['h']
            if sys.version_info.major<3:
                for e in set(raw_map[l][t]).difference(set(['l:n','h','s:p'])):
                    mut_map[int(l)][str(t)][str(e)] = [str(i) if type(i) is unicode else i for i in raw_map[l][t][e]]
            else:
                for e in set(raw_map[l][t]).difference(set(['l:n','h','s:p'])):
                    mut_map[int(l)][str(t)][str(e)] = [str(i) for i in raw_map[l][t][e]]
    return mut_map

def read_json_anueuploidy(json_path):
    raw,anueuploidy = {},{}
    if json_path.endswith('.gz'):
        with gzip.GzipFile(json_path,'rb') as f:
            raw = json.load(f)
    else:
        with open(json_path,'r') as f:
            raw = json.load(f)
    for k in raw:
        anueuploidy[str(k)] = {}
        for c in raw[k]: anueuploidy[str(k)][int(c)] = raw[k][c]
    return anueuploidy

#UCSC refGene table file, merge delta will merge overlaping cds regions by union
#obtained by direct mysql query:
#mysql --user=genome -N --host=genome-mysql.cse.ucsc.edu -A -D hg19 -e "select * from refGene" | gzip > refGene.txt.gz
def build_ucsc_gene_exon_map(refGene,swap=True,tx=False):
    G,E,C,W,raw = {},{},{},{},[]
    fields = {'bin':0,'name':1,'chrom':2,'strand':3,'txStart':4,'txEnd':5,'cdsStart':6,
              'cdsEnd':7,'exonCount':8,'exonStarts':9,'exonEnds':10,'name2':12}
    if refGene.endswith('.gz'):
        with gzip.GzipFile(refGene,'rb') as f:
            raw = f.readlines()
    else:
        with open(refGene,'r') as f:
            raw = f.readlines()
    for i in range(len(raw)):
        row = raw[i].split('\t')
        g,a1,x1,x2   = row[fields['name2']],row[fields['chrom']],int(row[fields['txStart']]),int(row[fields['txEnd']])
        if x1>x2: x1,x2 = x2,x1 #swap ordering if needed
        starts,ends = row[fields['exonStarts']].split(','),row[fields['exonEnds']].split(',')
        while len(starts)>=0 and starts[-1]=='': starts = starts[:-1]
        exons = [[min(starts[j],ends[j]),max(starts[j],ends[j])] for j in range(len(starts))]
        exons = sorted(exons,key=lambda x: x[0])
        if len(row)>=12 and (x2-x1)>0:
            if g in G:
                similiar = False
                for l in range(len(G[g])):
                    b1,y1,y2, = G[g][l]
                    if a1==b1 and intersect([x1,x2],[y1,y2]):
                        G[g][l]  = [a1,min(x1,y1),max(x2,y2)]
                        similiar = True
                        break
                if not similiar: G[g] += [[a1,x1,x2]]
                E[g] += [[a1,e[0],e[1]] for e in exons]
                E[g] = sorted(E[g],key=lambda x: (x[0],x[1]))
            else:
                G[g]  = [[a1,x1,x2]]
                E[g]  = [[a1,e[0],e[1]] for e in exons]
    for g in E:
        E[g] = sorted([list(j) for j in list(set([tuple(e) for e in E[g]]))],key=lambda x: (x[0],x[1]))
    if swap:  # swap from hg19->GCh37 or GCh37->hg19
        chroms = set([])
        for g in G:
            for l in range(len(G[g])):
                chroms.add(G[g][l][0])
        M = seq_name_map(chroms)
        for g in G:
            R,V = [],[]
            for l in range(len(G[g])):
                if G[g][l][0] in M: G[g][l][0] = M[G[g][l][0]]
                else:               R += [G[g][l]]
            for e in E[g]:
                if e[0] in M: e[0] = M[e[0]]
                else:               V += [e]
            for r in R: G[g].remove(r)
            for v in V: E[g].remove(v)
    for g in G:
        for l in range(len(G[g])):
            a1,x1,x2 = G[g][l]
            if a1 in C: C[a1] += [[g,x1,x2]]
            else:       C[a1]  = [[g,x1,x2]]
    for c in C:
        C[c] = sorted(C[c],key=lambda x: x[1])
    for c in C:
        W[c] = []
        for i in range(len(C[c])):
            W[c] += [[C[c][i][1],C[c][i][2],0.0,{'gene':set([C[c][i][0]])}]]
        W[c] = {'gene':sorted(W[c],key=lambda x: x[0])}
    return {'exon':E,'gene':G,'seq':C,'wcu':W}

#write a gzip compressed gene_map pickle
def write_gene_map(path,gene_map):
    with gzip.GzipFile(path,'wb') as f:
        pickle.dump(gene_map,f)
        return True
    return False

#read a gzip copressed gene map pickle
def read_gene_map(path):
    gene_map = {}
    with gzip.GzipFile(path,'rb') as f:
        gene_map = pickle.load(f)
    return gene_map

#given hg19 make GCh37 style names or
#given GCh37 style names make hg19 style
def seq_name_map(chroms):
    gc_set,gc  = set([str(s) for s in range(23)]+['X','Y']),0
    hg_set,hg  = set(['chr'+str(s) for s in list(gc_set)]),0
    for chrom in chroms:
        if chrom.upper().startswith('CHR'): hg += 1
        else:                               gc += 1
    C,D = {},{}
    if hg>=gc: #have hg take off the chr
        for chrom in chroms:
            k = chrom.upper().split('CHR')[1].split('_')[0]
            if chrom in hg_set:
                if k in C: C[k] += [chrom]
                else:      C[k]  = [chrom]
    else: #have 1 add the chr
        for chrom in chroms:
            k = 'chr'+chrom
            if chrom in gc_set:
                if k in C: C[k] += [chrom]
                else:      C[k]  = [chrom]
    for k in C:
        for i in range(len(C[k])):
            D[C[k][i]] = k
    return D

#given a gene list text file that is seperated by delim
def read_gene_list(path,delim='\n'):
    s,g,l = '',[],[]
    with open(path,'r') as f:
        s = ''.join(f.readlines())
    g = s.split(delim)
    for i in range(len(g)): #filter the empty genes
        if g[i]!='': l += [g[i]]
    return l

#implied gain or loss of function
#gene_symbol\tweight(+ or - fold change or something here...)
def read_de_gene_list(path,delim='\t'):
    l = []
    return l

#for each gene in the gene list lookup the position
#for every chrom, build a pos list and weight list w
#L = gene_list = [GS1,GS2,GS3,...]
#W = weight list or risk multiplier for each gene in the gene list
#G = gene_map = ['gene':{GS1:[c1,c2,...],GS2:[c1,c2,c3,...]}]
#P = wcu: {k:[[start,end,weight,{idx}],...]} -> idx = {'onco':set(['MYC2','BRAC1',...])}
#:::TO DO::: will need to write a js version of this conversion tool :::TO DO:::
def gene_list_to_wcu(L,W,G,label=''):
    P = {}
    for i in range(len(L)):
        if L[i] in G['gene']:
            for c in G['gene'][L[i]]:
                if c[0] in P: P[c[0]] += [[c[1],c[2],W[i],{label:set([L[i]])}]]
                else:         P[c[0]]  = [[c[1],c[2],W[i],{label:set([L[i]])}]]
    for k in P:
        P[k] = sorted(P[k],key=lambda x: x[0])
    return P

#pop off contigs that are not needed from the wcu
def remove_wcu_keys(rs,xs,wcu):
    R = {}
    for x in xs:
        if x in rs:
            R[x] = copy.deepcopy(wcu[x])
    return R

#class ids are not set type but are made into array for json
def wcu_to_json_wcu(wcu):
    J = {}
    for k in wcu:
        J[k] = []
        for i in range(len(wcu[k])):
            J[k] += [[wcu[k][i][0],wcu[k][i][1],wcu[k][i][2],
                      {l:sorted(list(wcu[k][i][3][l])) for l in wcu[k][i][3]}]]
    return J

#array type for class ids are array which become lists and then are converted to set
def json_wcu_to_wcu(raw,W=None):
    wcu = {}
    if equal_dim(raw,W):
        for k in raw:
            wcu[str(k)] = []
            for i in range(len(raw[k])):
                wcu[str(k)] += [[raw[k][i][0],raw[k][i][1],W[k][i],
                                 {str(l):set([str(x) for x in raw[k][i][3][l]]) for l in raw[k][i][3]}]]
    else:
        for k in raw:
            wcu[str(k)] = []
            for i in range(len(raw[k])):
                wcu[str(k)] += [[raw[k][i][0],raw[k][i][1],raw[k][i][2],
                                 {str(l):set([str(x) for x in raw[k][i][3][l]]) for l in raw[k][i][3]}]]
    return wcu

def equal_dim(X,Y):
    eq = True
    if X is dict and Y is dict:
        x,y = set(X.keys()),set(Y.keys)
        if len(x.difference(y))!=0:
            eq = False
        else:
            for k in X:
                if(len(X[k])!=len(Y[k])):
                    ed = False
                    break
    else:
        eq = False
    return eq

#reads json and then converts json wcu to wcu
def read_json_wcu(json_path,W=None):
    raw = {}
    with open(json_path,'r') as f:
        raw = json.load(f,object_pairs_hook=str_hook)
    wcu = json_wcu_to_wcu(raw,W)
    return wcu

#converts wcu to json wcu and then to disk
def write_json_wcu(json_path,wcu):
    with open(json_path,'w') as f:
        f.write(json.dumps(wcu_to_json_wcu(wcu)))
        return True
    return False

#will convert from the json for the 'wcu' section
def read_json_gene_map(json_path):
    gene_map = {}
    if json_path.endswith('.gz'):
        with gzip.open(json_path,'rb') as f:
            raw = json.load(f,object_pairs_hook=str_hook)
            gene_map['seq']  = raw['seq']
            gene_map['gene'] = raw['gene']
            gene_map['wcu']  = {k:json_wcu_to_wcu(raw['wcu'][k]) for k in raw['wcu']}
    else:
        with open(json_path, 'r') as f:
            raw = json.load(f,object_pairs_hook=str_hook)
            gene_map['seq']  = raw['seq']
            gene_map['gene'] = raw['gene']
            gene_map['wcu']  = {k:json_wcu_to_wcu(raw['wcu'][k]) for k in raw['wcu']}
    return gene_map

#will read a json based gene_map and convert to gene map
def write_json_gene_map(json_path,gene_map,gz=False):
    json_gene_map = {}
    json_gene_map['seq']  = gene_map['seq']
    json_gene_map['gene'] = gene_map['gene']
    json_gene_map['wcu']  = {k:wcu_to_json_wcu(gene_map['wcu'][k]) for k in gene_map['wcu']}
    if gz:
        if not json_path.endswith('.gz'): json_path += '.gz'
        with gzip.open(json_path,'wb') as f:
            f.write(json.dumps(json_gene_map))
            return True
        return False
    else:
        with open(json_path,'w') as f:
            f.write(json.dumps(json_gene_map))
            return True
        return False

def vca_to_json_vca(vca):
    json_vca = []
    for i in range(len(vca)):
        json_vca += [vca[i].to_dict()]
    return json_vca

def json_vca_to_vca(json_vca):
    vca = []
    for i in range(len(json_vca)):
        vca += [VariantCall()]
        vca[-1].from_dict(json_vca[i])
    return vca

def write_json_vca(json_path,vca,gz=True):
    if gz:
        if not json_path.endswith('.gz'): json_path += '.gz'
        with gzip.open(json_path,'wb') as f:
            f.write(json.dumps(vca_to_json_vca(vca)))
            return True
        return False
    else:
        with open(json_path,'w') as f:
            f.write(json.dumps(vca_to_json_vca(vca)))
            return True
        return False

def read_json_vca(json_path):
    vca = []
    if json_path.endswith('.gz'):
        with gzip.open(json_path,'rb') as f:
            vca = json_vca_to_vca(json.loads(f.read()))
    else:
        with open(json_path,'r') as f:
            vca = json_vca_to_vca(json.loads(f))
    return vca

#naive vcf reader that will pull the chrom and pos
#mainly intended for reading MELT VCF prior formated files
#to assist with MEI target locations in soMaCX workflows or in defining
#germline tollerant gene regions such as those published by 1000 genomes consortium
def vcf_to_wcu(path,delim='\t',skip='#',remove_dup=False,clean_coord=True,w=0.0,label='vcf',filter_size=None):
    s,g,l = '',[],{}
    if path.endswith('.gz'):
        with gzip.GzipFile(path,'rb') as f:
            s = ''.join(f.readlines())
    elif path.endswith('.vcf'):
        with open(path,'r') as f:
            s = ''.join(f.readlines())
    else: print('incorrect vcf file suffix...')
    has_end = s.find('END=')>0
    has_svlen = s.find('SVLEN=')>0
    has_svtype = s.find('SVTYPE=')>0
    has_sample = s.find('SAMPLE=')>0
    g = [t.split(delim) if not t.startswith(skip) else None for t in s.split('\n')]
    for i in g:
        if i is not None and len(i)>0 and i[0]!= '':
            if i[0] in l: l[i[0]] += [i]
            else:         l[i[0]]  = [i]
    if label is not None:
        s,g = '',{}
        if has_end and has_svtype and has_sample:
            S = set([])
            for k in l:
                g[k] = {label:[]}
                for i in range(len(l[k])):
                    if l[k][i][2] not in S:
                        S.add(l[k][i][2])
                        g[k][label] += [[int(l[k][i][1]),int(l[k][i][7].split('END=')[-1].split(';')[0]),
                                         w,{l[k][i][7].split('SVTYPE=')[-1].split(';')[0]:\
                                         set([l[k][i][7].split('SAMPLE=')[-1].split(';')[0]])}]]
        elif has_svlen and has_svtype:
            for k in l:
                g[k] = {label:[]}
                for i in range(len(l[k])):
                    g[k][label] += [[int(l[k][i][1]),int(l[k][i][1])+abs(int(l[k][i][7].split('SVLEN=')[-1].split(';')[0])),
                                     w,{l[k][i][7].split('SVTYPE=')[-1].split(';')[0]:set([l[k][i][2]])}]]
    else:
        s, g = '', {}
        if has_end and has_svtype and has_sample:
            S = set([])
            for k in l:
                g[k] = []
                for i in range(len(l[k])):
                    if l[k][i][2] not in S:
                        S.add(l[k][i][2])
                        g[k] += [[int(l[k][i][1]), int(l[k][i][7].split('END=')[-1].split(';')[0]),
                                  w, {l[k][i][7].split('SVTYPE=')[-1].split(';')[0]: \
                                      set([l[k][i][7].split('SAMPLE=')[-1].split(';')[0]])}]]
        elif has_svlen and has_svtype:
            for k in l:
                g[k] = []
                for i in range(len(l[k])):
                    g[k] += [[int(l[k][i][1]), int(l[k][i][1]) + abs(int(l[k][i][7].split('SVLEN=')[-1].split(';')[0])),
                              w, {l[k][i][7].split('SVTYPE=')[-1].split(';')[0]: set([l[k][i][2]])}]]
    if clean_coord:
        for k in g:
            if type(g[k]) is list:
                for i in range(len(g[k])):
                    r_min,r_max = min(g[k][i][0],g[k][i][1]),max(g[k][i][0],g[k][i][1])
                    g[k][i][0],g[k][i][1] = r_min,r_max
            elif type(g[k]) is dict:
                for l in g[k]:
                    for i in range(len(g[k][l])):
                        r_min, r_max = min(g[k][l][i][0],g[k][l][i][1]),max(g[k][l][i][0],g[k][l][i][1])
                        g[k][l][i][0],g[k][l][i][1] = r_min,r_max
    if remove_dup:
        u,h = {},{}
        for k in g:
            if type(g[k]) is list:
                u[k],h[k] = set([]),[]
                for i in range(len(g[k])):
                    v = (g[k][i][0],g[k][i][1],g[k][i][2],list(g[k][i][3].keys())[0])
                    if v not in u[k]:
                        u[k].add(v)
                        h[k] += [g[k][i]]
                h[k] = sorted(h[k])
                g[k] = h[k]
            elif type(g[k]) is dict and label in g[k]:
                u[k],h[k] = set([]),{label:[]}
                for i in range(len(g[k][label])):
                    v = (g[k][label][i][0],g[k][label][i][1],g[k][label][i][2],list(g[k][label][i][3].keys())[0])
                    if v not in u[k]:
                        u[k].add(v)
                        h[k][label] += [g[k][label][i]]
                h[k][label] = sorted(h[k][label])
                g[k][label] = h[k][label]
    if filter_size is not None:
        for k in g:
            if type(g[k]) is list:
                h = []
                for i in range(len(g[k])):
                    a = abs(g[k][i][1]-g[k][i][0]+1)
                    if a >= filter_size[0] and a <= filter_size[1]:
                        h += [g[k][i]]
                g[k] = h
            elif type(g[k]) is dict:
                h = {}
                for l in g[k]:
                    h[l] = []
                    for i in range(len(g[k][l])):
                        a = abs(g[k][l][i][1] - g[k][l][i][0] + 1)
                        if a >= filter_size[0] and a <= filter_size[1]:
                            h[l] += [g[k][l][i]]
                g[k] = h
    return g

#simple bed-3 region to wcu
def bed_to_wcu(path,delim='\t',w=0.0,label='bed'):
    s,g,x = '',{},-1
    with gzip.GzipFile(path,'rb') as f:
        s = ''.join(f.readlines())
    while s[x]=='\n': x -= 1
    s,y = s[:x+1],1
    for row in s.split('\n'):
        raw   =  row.split('\t')
        k,i,j = [raw[0],int(raw[1]),int(raw[2])]
        if k in g: g[k] += [[i,j,w,{label:set(['%s_%s'%(label,y)])}]]
        else:      g[k]  = [[i,j,w,{label:set(['%s_%s'%(label,y)])}]]
        y += 1
    return g

#chrom\tstart\tend\tp-value
#dump the wcu used for germline and somatic
#and check against the VCF files generated
#to see if the values are correctly distributed
def wcu_to_bedgraph(wcu,path,sign=1.0):
    s = ''
    for k in wcu:
        for c in wcu[k]:
            s += '\n'.join(['\t'.join([k,str(row[0]),str(row[1]),str(row[2]*sign)]) for row in wcu[k][c]])+'\n'
    with open(path, 'w') as f:
        f.write(s)
        return True
    return False

#place holder for uniprot tumor suppressor genes list
def uniprot_to_wcu(path):
    return True

def pos_to_wcu(pos,w=0.0,label=None):
    return [[pos[i][0],pos[i][1],w,{label:set([label+'_'+str(i+1)])}] for i in range(len(pos))]

#offers two way to reweight the wcu, using the counts per class or all counts
#mainly used to summarize distrubutions on regions reported by g1kp3, ect...
#min_w will filter entries that are less and max_w will trim those above,
#filter lesser will check the weights of each piled class and keep the greatest only
#perfect for building distrubtions that both capture some of the nature distributions
#while maintianing a infinite alleles style model where each location is unique
def weight_wcu_by_idx_count(wcu,count_all_classes=True,filter=[5,500,True]):
    W = {}
    for k in wcu:
        W[k] = {}
        G = weight_graph(wcu[k])
        P = scan_graph(G)
        if count_all_classes and filter is not None:
            W[k]['all'] = []
            for i in range(len(P)):
                if filter[2]: #implies majority rules
                    b = list(P[i][3].keys())[np.argmax([len(P[i][3][c]) for c in P[i][3]])]
                    P[i][3] = {b:P[i][3][b]}
                sums = {c:set([len(P[i][3][c])]) for c in P[i][3]}
                cs = list(sums.keys())
                for c in cs:
                    for s in sums[c]:
                        if s<filter[0] or s>filter[1]: sums.pop(c)
                if len(sums)>0:
                    W[k]['all'] += [[P[i][0],P[i][1],float(sum([len(P[i][3][c]) for c in P[i][3]])),sums]]
        elif count_all_classes:
            W[k]['all'] = []
            for i in range(len(P)):
                sums = {c: set(len(P[i][3][c])) for c in P[i][3]}
                W[k]['all'] += [[P[i][0],P[i][1],float(sum([len(P[i][3][c]) for c in P[i][3]])),sums]]
        else:
            for i in range(len(P)):
                for c in P[i][3]:
                    if t in W[k]:
                        W[k][c] += [[P[i][0],P[i][1],float(len(P[i][3][c]),P[i][3])]]
                    else:
                        W[k][c]  = [[P[i][0],P[i][1],float(len(P[i][3][c]),P[i][3])]]
    return W

#use the pregenerated wcu key in the gene_map
#to generate a list of genes, aka affected regions
def wcu_to_gene_list(wcu,gene_map):
    gain_genes,loss_genes,genes = set([]),set([]),set([])
    for k in wcu:
        for c in wcu[k]:
            if k in gene_map['wcu']:
                G = weight_graph({'gene':gene_map['wcu'][k]['gene'],c:wcu[k][c]})
                P = scan_graph(G)
                for i in range(len(P)):
                    if 'gene' in P[i][3] and len(P[i][3])>1:
                        for g in P[i][3]['gene']:
                            if 'DUP' in P[i][3]:   gain_genes.add(g)
                            elif 'DEL' in P[i][3]: loss_genes.add(g)
                            genes.add(g)
    return sorted(list(gain_genes)),sorted(list(loss_genes)),sorted(list(genes))

#special g1kp3 VCF format file to gene list tool
def g1kp3_to_gene_list(path,gene_map,filter=[int((1/50.0)*2500),int((1/2.0)*2500),True],write=None):
    vcf_wcu = vcf_to_wcu(path)
    idx_wcu = weight_wcu_by_idx_count(vcf_wcu,filter=filter)
    gain,loss,total = wcu_to_gene_list(idx_wcu,gene_map)
    if write is not None:
        with open(get_local_path(write),'w') as f:
            f.write('\n'.join(total))
    return gain,loss,total

# returns a dict of chroms and number of Variants
def vcf_chroms(path):
    C = {}
    with open(path,'r') as f:
        raw = [l.replace('\n','') for l in f.readlines()]
    for row in raw:
        if not row.startswith('#'):
            c = row.split('\t')[0]
            if c not in C: C[c]  = 1
            else:          C[c] += 1
    return C

def alter_lvcam_genotypes(lvcam,h=0.5):
    gs = [[1,0],[0,1],[1,1]]
    for l in lvcam:
        for k in lvcam[l]:
            i = 0
            while i < len(lvcam[l][k])-1:
                if lvcam[l][k][i].end>=lvcam[l][k][i+1].pos:
                    g = np.random.choice([0,1],p=[0.5,0.5])
                    lvcam[l][k][i].frmat   = '%s|%s'%(gs[g][0],gs[g][1])
                    lvcam[l][k][i+1].frmat = '%s|%s'%(gs[g][1],gs[g][0])
                    i += 1
                else:
                    g = np.random.choice([0,1,2],p=[0.5*h,0.5*h,1-h])
                    lvcam[l][k][i].frmat = '%s|%s'%(gs[g][0],gs[g][1])
                i += 1
    return lvcam

#given a single sample FusorSV VCF file, create a VCAM
#[1] read the FusorSV VCF file het? 0/1 => 1|0, 0|1, 1|1 randomly?
#[2] convert to VariantCall and store into vcam
#[3] map the types appropriately INS,DEL,DUP=> only tandem, INV=>only simple =>
#[4] this involves mining the INFO field correctly and reforming it again...
#[5] sort the vca for each chrom and save back again
def vcf_to_vcam(vcf_path,ref_path,chroms,skip='#',delim='\t',small=50):
    s,g,l,M = '',[],{},{}
    if vcf_path.endswith('.gz'):
        with gzip.GzipFile(vcf_path,'rb') as f:
            s = ''.join(f.readlines())
    elif vcf_path.endswith('.vcf'):
        with open(vcf_path,'r') as f:
            s = ''.join(f.readlines())
    else:
        print('incorrect vcf file suffix...')
    while s[-1]=='\n': s = s[:-1] #clean up any loose '\n'
    for t in s.split('\n'):
        if not t.startswith(skip): g += [t.split(delim)]

    #test for FusorSV G1k VCF formats
    if len(g)>0:
        if len(g[0])>9: #FusorSV/soMaCX vcf file
            for i in g:
                if i is not None and len(i) > 0 and i[0] in chroms:
                    #vcam example------------------------------------------------------------------------------------------
                    #{'alt': 'CCCAAC','chrom': '1','end': 23793569,'flter': None,'frmat': '1|1',
                    # 'identifier': '.','info': 'SVTYPE=INS;END=23793569;SVLEN=6;','pos': 23793569,'qual': None,'ref': 'N'}

                    #fusorsv vcf example------------------------------------------------------------------------------------------------
                    #[chrom,pos,id,ref,alt,qual,filter,info,format,GT,geno1]
                    #['1','2633901','fusorSV_1','C','<DEL>','.','PASS',
                    # 'SVTYPE=DEL;SVLEN=50300;END=2684201;CHR2=1;IMPRECISE;SVEX=0.30340265922;SVMETHOD=9:120,18:32,11:45,10:42;TARGET=0',
                    #'GT','0/1']
                    cdx = [i[0],int(i[1]),get_info_end(i[7]),get_info_type(i[7])]
                    if tuple(cdx) not in M:   #no duplicate entries
                        vc = VariantCall(chrom=i[0],pos=int(i[1]),identifier=i[2],ref=i[3],
                                         alt=i[4],qual=i[5],flter=i[6],info=i[7],frmat=i[9])
                        if vc.chrom in l: l[vc.chrom] += [vc]
                        else:             l[vc.chrom]  = [vc]
                        M[tuple(cdx)]  = 1
                    else:
                        M[tuple(cdx)] += 1
            for k in l: l[k] = sorted(l[k],key=lambda x: x.pos) #pos sorted by chrom
            for k in l:
                print('converting SVs to soMaCX form on chrom %s'%k)
                for i in range(len(l[k])):
                    svtype = get_info_type(l[k][i].info)
                    # DEL are fine.... nope, ref need to be the ref[vc.pos:vc.end]
                    if svtype == 'DEL':
                        l[k][i].ref = ru.read_fasta_substring(ref_path,k,l[k][i].pos,l[k][i].end)
                        l[k][i].alt = '<DEL>'
                        l[k][i].info = set_info_len(l[k][i].info,get_info_len(l[k][i].info)+1)
                    # 'SVTYPE=DUP;END=247916018;SVLEN=56793;DUP_POS=247859225;DUP_TYPE=TANDEM;ICN=4'
                    if svtype == 'DUP':  # CN, D_POS, ...
                        l[k][i].ref = ru.read_fasta_substring(ref_path,k,l[k][i].pos,l[k][i].end)
                        l[k][i].info += ';DUP_POS=%s;DUP_TYPE=TANDEM;ICN=4'%l[k][i].pos
                    if svtype == 'INV':  # TYPE == SIMPLE
                        l[k][i].info += ';INV_TYPE=PERFECT'
                        l[k][i].ref = ru.read_fasta_substring(ref_path,k,l[k][i].pos,l[k][i].end)
                        l[k][i].alt = utils.get_reverse_complement(l[k][i].ref)
                    if svtype == 'INS':
                        l[k][i].info = set_info_len(l[k][i].info,get_info_len(l[k][i].info)+1)
                    #if svtype == 'TRA' ...
                    #if svtype == 'INS' ? ...
        else: #g1kp3 vcf file
            for i in g:
                if i is not None and len(i) > 0 and i[0] in chroms:
                    # vcam example------------------------------------------------------------------------------------------
                    # {'alt': 'CCCAAC','chrom': '1','end': 23793569,'flter': None,'frmat': '1|1',
                    # 'identifier': '.','info': 'SVTYPE=INS;END=23793569;SVLEN=6;','pos': 23793569,'qual': None,'ref': 'N'}

                    # fusorsv vcf example------------------------------------------------------------------------------------------------
                    # [chrom,pos,id,ref,alt,qual,filter,info,format,GT,geno1]
                    # ['1','2633901','fusorSV_1','C','<DEL>','.','PASS',
                    # 'SVTYPE=DEL;SVLEN=50300;END=2684201;CHR2=1;IMPRECISE;']
                    cdx = [i[0],int(i[1]),get_info_end(i[7]),get_info_type(i[7])]
                    if tuple(cdx) not in M:  # no duplicate entries
                        if cdx[3]=='INS':
                            if abs(cdx[2]-cdx[1])<small:
                                cdx[2] = cdx[1]+int(np.random.uniform(0.5*small,1.5*small,1))
                                alt = gen_alt_seq(ru.read_fasta_substring(ref_path,cdx[0],cdx[1],cdx[2]))
                                ref = 'N'
                                i[7] = i[7].rsplit('END=')[0]+'END=%s;'%cdx[2]+';'.join(i[7].rsplit('END=')[-1].rsplit(';')[1:])
                        else:
                            alt = i[4]
                            ref = i[3]
                        i[7] += ';SVLEN=%s'%(cdx[2]-cdx[1]+1)
                        vc = VariantCall(chrom=cdx[0],pos=cdx[1],identifier=i[2],ref=ref,
                                         alt=alt,qual=i[5],flter=i[6],info=i[7],frmat='0|1')
                        if vc.chrom in l: l[vc.chrom] += [vc]
                        else:             l[vc.chrom]  = [vc]
                        M[tuple(cdx)]  = 1
                    else:
                        M[tuple(cdx)] += 1
            for k in l: l[k] = sorted(l[k], key=lambda x: x.pos)
            for k in l:
                print('converting SVs to soMaCX form on chrom %s'%k)
                for i in range(len(l[k])):
                    svtype = get_info_type(l[k][i].info)
                    # ref need to be the ref[vc.pos:vc.end]
                    if svtype == 'DEL':
                        l[k][i].ref = ru.read_fasta_substring(ref_path,k,l[k][i].pos,l[k][i].end)
                        l[k][i].alt = '<DEL>'
                    # 'SVTYPE=DUP;END=247916018;SVLEN=56793;DUP_POS=247859225;DUP_TYPE=TANDEM;ICN=4'
                    if svtype == 'DUP':  # CN, D_POS, ...
                        l[k][i].ref = ru.read_fasta_substring(ref_path,k,l[k][i].pos,l[k][i].end)
                        l[k][i].info += ';DUP_POS=%s;DUP_TYPE=TANDEM;ICN=4' % l[k][i].pos
                    if svtype == 'INV':  # TYPE == SIMPLE
                        l[k][i].ref = ru.read_fasta_substring(ref_path,k,l[k][i].pos,l[k][i].end)
                        l[k][i].alt = utils.get_reverse_complement(l[k][i].ref)
                    # if svtype == 'TRA' ...
                    # if svtype == 'INS' ? ...
    if 'Y' in l and 'X' in l:
        for i in range(len(l['Y'])):
            l['Y'][i].frmat = '1'
        for i in range(len(l['X'])):
            l['X'][i].frmat = '1'
    return l

#remove conflicting region vc by coordinate, sv-type frequency priority
def vcam_remove_conflicts(vcam):
    for k in vcam:
        F = {}
        for vc in vcam[k]: #freq analysis on that k
            svtype = get_info_type(vc.info)
            if svtype in F: F[svtype] += 1
            else:           F[svtype]  = 1
        C = {} #conflict graph
        for i,j in it.combinations(range(len(vcam[k])),2):
            a,b = [vcam[k][i].pos,vcam[k][i].end],[vcam[k][j].pos,vcam[k][j].end]
            if intersect(a,b):
                if i in C: C[i].add(j)
                else:      C[i]  = set([j])
                if j in C: C[j].add(i)
                else:      C[j]  = set([i])
        R = set([])
        cs = list(C.keys())
        for i in cs:
            if i in C:
                fwd,rev = 0,0
                for j in C[i]:
                    if j in C: #can be removed
                        if F[get_info_type(vcam[k][i].info)] >= F[get_info_type(vcam[k][j].info)]:
                            fwd += len(C[i].difference(C[j]))
                        if F[get_info_type(vcam[k][i].info)] <= F[get_info_type(vcam[k][j].info)]:
                            rev += len(C[j].difference(C[i]))
            if fwd >= rev and len(C[i])>0:
                R.add(i)
                for j in C[i]: C[j] = C[j].difference(set([i]))
                C.pop(i)
        print('removing %s VCF call conflicts for seq %s'%(len(R),k))
        vca = []
        for i in range(len(vcam[k])):
            if i not in R: vca += [vcam[k][i]]
        vcam[k] = vca
    return vcam

#returns the elements in vcam_b that are not part of vcam_a
#using the tuple (chrom, start, end, type)
def vcam_diff(vcam_a,vcam_b):
    D = {}
    for k in set(vcam_a.keys()).difference(set(vcam_b.keys())):
        D[k] = copy.deepcopy(vcam_a[k])
    for k in set(vcam_a.keys()).intersection(set(vcam_b.keys())):
        D[k] = []
        ks_a = {(vcam_a[k][i].chrom,vcam_a[k][i].pos,vcam_a[k][i].end):i for i in range(len(vcam_a[k]))}
        ks_b = {(vcam_b[k][i].chrom,vcam_b[k][i].pos,vcam_b[k][i].end):i for i in range(len(vcam_b[k]))}
        for v in set(ks_a.keys()).difference(set(ks_b.keys())):
            D[k] += [vcam_a[k][ks_a[v]]]
        D[k] = sorted(D[k],key=lambda x: x.pos)
    return D

#returns the unique elements in vcam_a and vcam_b
#using the tuple (chrom, start, end, type)
def vcam_union(vcam_a,vcam_b):
    D,A,B = {},vcam_a,vcam_b
    for k in set(A.keys()).union(set(B.keys())):
        D[k],ks_a,ks_b = [],{},{}
        if k in A:
            ks_a = {(A[k][i].chrom,A[k][i].pos,A[k][i].end):i for i in range(len(A[k]))}
            for v in set(ks_a.keys()):
                D[k] += [A[k][ks_a[v]]]
        if k in B:
            ks_b = {(B[k][i].chrom,B[k][i].pos,B[k][i].end):i for i in range(len(B[k]))}
            for v in set(ks_b.keys()).difference(set(ks_a.keys())):
                D[k] += [B[k][ks_b[v]]]
        D[k] = sorted(D[k],key=lambda x:x.pos)
    return D

#if aneuploidy causes a full chrom DEL, then we take the DEL in vcam_b and remove from vcam_a
def del_aneuploidy_vcam(vcam_a,vcam_b):
    D = {}
    for k in set(vcam_a.keys()).difference(set(vcam_b.keys())):
        D[k] = vcam_a[k] #copy all the vca from non affected k
    for k in set(vcam_a.keys()).intersection(set(vcam_b.keys())):
        D[k] = []
        #get the genotype
        #get_genotype_vca(vcam_a[k],g=,index=)
    return D

def merge_wcu(wcu):
    M = {}
    for c in wcu:
        for k in wcu[c]:
            if k in M:
                if c in M[k]: M[k][c] += wcu[c][k]
                else:         M[k][c]  = wcu[c][k]
            else:             M[k] = {c: wcu[c][k]}
    return M

def pretty_size(b,units='bp'):
    s,size_map = '',{0:'',1:'K',2:'M',3:'G',4:'T',5:'P',6:'E',7:'Z',8:'Y',9:'Z'}
    x,m = [b,0],0
    if x[0]/1000>0:
        while x[0]/1000>0:
            x[1] = x[0]%1000
            x[0],m = x[0]/1000,m+1
        d = ''
        if x[1]>0: d = '.'+str(x[1])[:1]
        s = str(x[0])+d+size_map[m]+units
    else:
        s = str(b)+units
    return s

#grab the exclude regions bed file and parse it into a seq dict
#with a exclude pos list inside
def read_exclude(exclude_path):
    exclude = {}
    with open(exclude_path, 'r') as f:
        for line in f:
            row = f.split(',')
            if len(row) > 0:
                if row[0] in exclude:
                    exclude[row[0]] += [[int(i) for i in row[1:]]]
                else:
                    exclude[row[0]]  = [[int(i) for i in row[1:]]]
    return exclude

#use this ratio to set the number of SNVs or SV according to the sorted length
def get_ratios(S):
    largest = max([S[k] for k in S])
    ratios = {k:1.0*S[k]/largest for k in S}
    return ratios

def get_mean_size(pos):
    if len(pos)>0:
        return np.mean([x[1]-x[0] for x in pos])
    else:
        return 0.0

def gen_random_seqs(seqs={'chr1':1E5},alphabet=['A','C','G','T'],
                   prob=[0.25,0.25,0.25,0.25],method='slow'):
    ss = []
    sigma = np.asarray(alphabet,dtype='S1') #do set difference for SNPs
    for k in seqs:
        if method=='slow': #slow method
            s = ''
            for i in range(0,seqs[k]):
                s += np.random.choice(sigma,size=1,p=prob).tostring()
        elif method=='fast':                                 #~100x faster
            s = np.random.choice(sigma,size=seqs[k],p=prob)  #size param
            s = s.tostring()
        else: raise KeyError
        ss += [s]
    return ss[::-1]

#[0] given a type and variant size interval VS --> {INS:1E0}
def gen_chrom_lens(L,chroms):
    c = len(chroms)
    ls = np.random.geometric(0.4,L) #partitioned lens of each chrom that will sum to L
    LS = np.histogram(ls,bins=list(set(ls)))  #bin the geometric dist
    chr_lens = list(LS[0][0:c])               #chop off largest c
    r = L-sum(chr_lens)                       #calculate leftovers
    chr_lens = [i+int(r/c) for i in chr_lens] #distribute leftovers
    if r%c!=0: chr_lens[0]+=L-sum(chr_lens)   #correct to sum to L
    return {chroms[i]:chr_lens[i] for i in range(0,c)}

def intersect(a,b):
    if (a[0]<=b[1] and a[1]>=b[0] or \
        b[0]<=a[0] and b[1]>=a[0] or \
        a[0]>=b[0] and a[1]<=b[1] or \
        b[0]>=a[1] and b[1]<=a[1]):
        return True
    else:
        return False

#given a mut_p and generated vca
#do a histogram and count events
#M is one sequence, m is its length
def mut_p_rate(mut_p,M,m,l=2):
    H = {t:{s:0.0 for s in mut_p[l][t]['l:n']} for t in mut_p[l]}
    for t in M:
        bins=sorted(mut_p[l][t]['l:n'])
        for i in range(len(M[t])):
            x = (M[t][i][1]-M[t][i][0])
            for s in bins:
                if x>=s*0.5001 and x<=s*1.5:
                    H[t][s] += 1.0
                    break
    for t in H:
        for s in H[t]:
            H[t][s] /= mut_p[l][t]['l:n'][s]*m
    return H

#given the type:pos that you make using the normal gen_map
#returns all the wcu that overlap each other therin violating IS model
def mut_overlap(M):
    O = {t:pos_to_wcu(M[t],1.0,t) for t in M}
    G = weight_graph(O)
    try:
        P = scan_graph(G,scale=False)
    except Exception as E:
        pad = ':'.join(['-' for i in range(20)])
        print(pad+'@@scan graph error@@'.upper()+pad)
        ru.dump_state({'O':O,'G':G},'mut_overlap-scan_graph',os.path.expanduser('~/'))
    Q = []
    for p in P:
        if sum([len(p[3][i]) for i in p[3]])>1:
            Q += [p]
    return Q

def mut_over_detected(over):
    over_det = False
    for i in over:
        idxs = set(i[3].keys())
        if len(idxs)>1: over_det = True
    return over_det

#update the probabilty to sample given the loss of key k
#transfering the k key equally back to the others and poping k
def update_sampling_prob(P,k):
    if k in P:
        p = P.pop(k)
        n = 1.0*len(P)
        P = {r:P[r]+p/n for r in P}
    return P

def class_hp_to_pos(class_p,cutoff=0.0):
    pos = []
    for t in class_p:
        for i in range(len(class_p[t][1])):
            if class_p[t][1][i]<=cutoff:
                pos += [[class_p[t][0][i][0],class_p[t][0][i][1]]]
    pos = sorted(pos,key=lambda x: x[0])
    return utils.merge_1D(pos)

    #start with just one level

def gen_class_mut_pos_map_slow(ref_seq,class_p,mut_p,l,y=10,
                               size_dist='uniform',size_prop=3.0,center=False,cutoff=0.0):
    #[1] generate candidates and probabilties distributions to sample from
    A,F,P,R,M,G,m,n = [],{'requested':0.0,'remaining':0.0},{},{},{},[],[],len(ref_seq) #x= total requested
    for t in mut_p[l]: #masked or areas that should not be utilized
        G = class_hp_to_pos(class_p,cutoff)
        for s in mut_p[l][t]['l:n']: m += [s]
    a,b,c,d = 1.0*np.min(m),1.0*np.max(m),1.0,1.0*size_prop
    for t in mut_p[l]: #requested is the desired prob of a type and size
        F[t]={}      #gene_p is a type and prob mapping that can raise or lower gene probabilities for types

        for s in mut_p[l][t]['l:n']:
            r = int(round(mut_p[l][t]['l:n'][int(round(s,0))]*n,0))
            x = ((d-c)/(b-a))*(s-a)+c
            r = int(round(r*x,0))

            if t in class_p: #mulitmodal distribution over the class ranges and weights
                s_class_p = copy.deepcopy(class_p[t])
                start_pos = utils.weighted_random(s_class_p[0],s_class_p[1],r*y)
                if center:
                    if size_dist == 'triangular':
                        start_pos = [[max(0,i-int(round(np.random.triangular(0.5001*s,s,1.5*s), 0))/2),
                                      min(n,i+int(round(np.random.triangular(0.5001*s,s,1.5*s), 0))/2)] for i in start_pos]
                    elif size_dist == 'uniform':
                        start_pos = [[max(0,i-int(round(np.random.uniform(0.5001*s,1.5*s),0))/2),
                                      min(n,i+int(round(np.random.uniform(0.5001*s,1.5*s),0))/2)] for i in start_pos]
                else:
                    if size_dist=='triangular':
                        start_pos=[[max(0,i),min(n,i+int(round(np.random.triangular(0.5001*s,s,1.5*s),0)))] for i in start_pos]
                    elif size_dist=='uniform':
                        start_pos=[[max(0,i),min(n,i+int(round(np.random.uniform(0.5001*s,1.5*s),0)))] for i in start_pos]
            else: #default will be uniform over the full range
                start_pos=np.array(np.random.uniform(0,n-s,r*y),dtype=int)
                if size_dist=='triangular':
                    start_pos=[[i,min(n,i+int(round(np.random.triangular(0.5001*s,s,1.5*s),0)))] for i in start_pos]
                elif size_dist=='uniform':
                    start_pos=[[i,min(n,i+int(round(np.random.uniform(0.5001*s,s,1.5*s),0)))] for i in start_pos]

            F[t][s],R[(t,s)] = {'r':r,'pos':start_pos},[]      #y*r requested are generated
            F['requested'],F['remaining'] = F['requested']+r,F['remaining']+len(start_pos)
    for t in mut_p[l]: #probabilties for smaple from each distribution are maintained
        for s in mut_p[l][t]['l:n']:
            if F['requested']>0: P[(t,s)] = F[t][s]['r']/F['requested']

    #[2] randomly sample from F using P into R and report results in F/R while P is maintained-------
    while F['requested']>0 and F['remaining']>0 and len(P)>0:        #select which pos list to add to
        t,s = list(P.keys())[np.random.choice(range(len(P)),p=[P[k] for k in P])]   #sample from (t,s) keys
        if len(R[(t,s)])<1 and len(F[t][s]['pos'])>0 and F[t][s]['r']>0:      #first sample only
            if len(utils.LRF_1D(G,[F[t][s]['pos'][0]])[0])<1:
                R[(t,s)] += [F[t][s]['pos'][0]] #get the first
                F[t][s]['r']   -= 1         #update the requested number
                F['requested'] -= 1         #update the total requested
                F['remaining'] -= 1         #update remaining
                F[t][s]['pos'].pop(0)       #update pos pool
            else:
                F[t][s]['pos'].pop(0)
                F['remaining'] -= 1
        else:
            while len(F[t][s]['pos'])>0 and F[t][s]['r']>0:
                if len(utils.LRF_1D(G,[F[t][s]['pos'][0]])[0])<1 and \
                   len(utils.LRF_1D(A,[F[t][s]['pos'][0]])[0])<1:
                    A += [F[t][s]['pos'][0]]
                    R[(t,s)] += [F[t][s]['pos'][0]]
                    F[t][s]['r']   -= 1
                    F['requested'] -= 1
                    F['remaining'] -= 1
                    F[t][s]['pos'].pop(0)
                    break
                else:
                    F[t][s]['pos'].pop(0)
                    F['remaining'] -= 1
            A = sorted(A)
        R[(t,s)]=sorted(R[(t,s)])
        if len(F[t][s]['pos'])<1 or F[t][s]['r']<1:
            F[t].pop(s)
            P = update_sampling_prob(P,(t,s))
    for k in R: #clean up and delever each type seperately as keys in M
        if k[0] in M:
            M[k[0]] += R[k]
        else:
            M[k[0]] = R[k]
    for t in M: M[t] = sorted(M[t])

    Q,H = mut_overlap(M),{}
    for q in Q:
        for t in q[3]:
            if t not in H: H[t] = {}
            for n in q[3][t]:
                i = int(n.split('_')[-1])-1
                if i in H[t]: H[t][i] += 1
                else:         H[t][i]  = 1
    for t in H:
        V = []
        for i in H[t]:
            if H[t][i]>1 or len(M[t][i])>1:
                V += [M[t][i]]
        for v in V:
            M[t].remove(v)
    #-------------------------------------------------------------------------------------------------
    return M

#alternate generation algorithm, that uses a coordinate_map distribution called class_p
#to set the joint distribution with mut_p. level is l and y is the multiple
#of the number of mutations that are requested that are provided as a fail-safe
#this alternate algorithm makes the SV pos list first and then does a merge and clean
#step instead of looking for collisions which is generally much faster...
def gen_class_mut_pos_map(ref_seq,class_p,mut_p,l,
                          size_dist='uniform',size_prop=2.0,
                          center=False,germ=True):
    A,R,Y,M,m,n = {},{},{},{},[],len(ref_seq)  # x= total requested
    for t in mut_p[l]:
        for s in mut_p[l][t]['l:n']: m += [s]
    a,b,c,d = 1.0*np.min(m),1.0*np.max(m),1.0,1.0*size_prop
    # [1] generation step--------------------------------------------------------------------------------------
    for t in mut_p[l]:
        for s in mut_p[l][t]['l:n']:
            R[(t,s)] = int(round(mut_p[l][t]['l:n'][s]*n,0)) #requested
            x = ((d-c)/(b-a+1))*(s-a)+c
            r = int(round(R[(t,s)]*x,0))
            if t in class_p: #mulitmodal distribution over the class ranges and weights
                s_class_p = copy.deepcopy(class_p[t])
                start_pos = utils.weighted_random(s_class_p[0],s_class_p[1],r)
                if center:
                    if size_dist == 'triangular':
                        local_size = int(round(np.random.triangular(0.5001*s,s,1.5*s), 0))/2
                        start_pos = [[max(0,i-local_size),min(n,i+local_size)] for i in start_pos]
                    elif size_dist == 'uniform':
                        local_size = int(round(np.random.uniform(0.5001*s,1.5*s),0))/2
                        start_pos = [[max(0,i-local_size),min(n,i+local_size)] for i in start_pos]
                else:
                    if size_dist=='triangular':
                        start_pos=[[max(0,i),min(n,i+int(round(np.random.triangular(0.5001*s,s,1.5*s),0)))] for i in start_pos]
                    elif size_dist=='uniform':
                        start_pos=[[max(0,i),min(n,i+int(round(np.random.uniform(0.5001*s,1.5*s),0)))] for i in start_pos]
            else: #default will be uniform over the full range
                start_pos=np.array(np.random.uniform(0,n-s,r),dtype=int)
                if size_dist=='triangular':
                    start_pos=[[max(0,i),#-int(round(np.random.triangular(0.5001*s,s,1.5*s),0))/2),
                                min(n,i+int(round(np.random.triangular(0.5001*s,s,1.5*s),0))/1)] for i in start_pos]
                elif size_dist=='uniform':
                    start_pos=[[max(0,i),#-int(round(np.random.uniform(0.5001*s,1.5*s),0))/2),
                                min(n,i+int(round(np.random.uniform(0.5001*s,1.5*s),0))/1)] for i in start_pos]
            Y[(t,s)] = [[start_pos[i][0],start_pos[i][1],1.0,{t:{s}}] \
                        for i in range(len(start_pos))] #try using a wcu here
    if not germ: #in somatic the loss_wcu and gain_wcu had the germ class put inside, so class_p.keys() might have germ
        C,t = {},list(class_p.keys())[0]
        for i in range(len(class_p[t][0])):
            if class_p[t][1][i]<=0.0:
                if ('germ') in C:
                    C[('germ')] += [[class_p[t][0][i][0],class_p[t][0][i][1],0.0,{'germ':set([0.0])}]]
                else:
                    C[('germ')]  = [[class_p[t][0][i][0],class_p[t][0][i][1],0.0,{'germ':set([0.0])}]]
        for k in C: Y[k] = sorted(C[k],key=lambda x: x[0]) #append into Y and scan
    # [2] merge and cleaning step algo--------------------------------------------------------------------------
    G = weight_graph(Y)
    P = scan_graph(G,scale=False)
    W = get_weights(G[0])
    P = scale_weights(P,W,scale=[0.0,10.0],method='add-mult') #special scaling
    for i in range(len(P)): #load conflict free, randomly sample
        idx = P[i][3]
        if sum([1 for s in [idx[k] for k in idx] for v in s])<=1:
            for t in idx:
                for s in idx[t]: #needs to be at least half as big as its bin
                    if 1.0*(P[i][1]-P[i][0])>=0.5*s:
                        if (t,s) in A: A[(t,s)] += [[P[i][0],P[i][1]]]
                        else:          A[(t,s)]  = [[P[i][0],P[i][1]]]
    for k in R:
        if k in A:
            I = np.random.choice(range(len(A[k])),min(R[k],len(A[k])),replace=False)
            if k[0] in M: M[k[0]] += [A[k][i] for i in I]
            else:         M[k[0]]  = [A[k][i] for i in I]
    for t in M: M[t] = sorted(M[t])
    return M

#generate inner SV regions
def gen_mut_pos_inside(mut_pos,mut_l,mut_n,break_ends=0.25,low_bp=25):
    if break_ends < 0.0 or break_ends > 1.0:
        print('break_ends not set to a proper probabilty [0.0,1.0]')
        raise AttributeError
    inside_pos = []
    if break_ends > 0.0: mut_n -= int(round(2.0*break_ends,0))
    if len(mut_pos)>0 and len(mut_pos[0])>1:
        for pos in mut_pos: #start with about 2x pos and filter collisions
            if break_ends<=0.0: #no breakend mut prob
                start_pos = np.array(np.random.uniform(pos[0],max(pos[0]+mut_l,pos[1]-mut_l),2*mut_n),dtype=int)
            else:               #leave room for breakend prob
                start_pos = np.array(np.random.uniform(pos[0]+mut_l,max(pos[0]+mut_l,pos[1]-2*mut_l),2*mut_n),dtype=int)
            start_pos = [[i,i+max(low_bp,int(round(np.random.uniform(0.25001*mut_l,1.75*mut_l),0)))] for i in start_pos]
            start_pos = [[j[0],j[1]] for j in list(set([(i[0],i[1]) for i in start_pos]))] #keep only unique
            inner_pos = []
            for i in start_pos:
                if len(inner_pos)<mut_n:
                    collision = False
                    for j in inner_pos:
                        if intersect(i,j): collision = True
                    if not collision: inner_pos += [i]
                else: break
            if break_ends > 0.0: #
                size = max(low_bp,int(round(np.random.uniform(0.25001*mut_l,1.75*mut_l),0)))
                if np.random.choice([True,False],p=[break_ends,1.0-break_ends]):
                    inner_pos += [[pos[0],pos[0]+size]]
                size = max(low_bp,int(round(np.random.uniform(0.25001*mut_l,1.75*mut_l),0)))
                if np.random.choice([True,False],p=[break_ends,1.0-break_ends]):
                    inner_pos += [[pos[1]-size,pos[1]]]
            inside_pos += inner_pos
        inside_pos = sorted(inside_pos)
    return utils.LRF_1D(mut_pos,inside_pos)[0] #intersection so they are clipped to boundries

#given one range outer_pos = [start,stop] and a list of inner_pos pos that fall inside
#that range, flip the coordinates as used for correcting INV coordinates
#if the corrdinates are not in range, print a warning and return empty
#always returns a sorted complementary pos list for downstream usage
def reverse_complement_pos(outer_pos,inner_pos,verbose=False):
    comp_pos = []
    for outer in outer_pos:
        pass_pos = []
        for pos in inner_pos:
            if len(utils.LRF_1D([outer],[pos])[0])>0:
                if verbose: print('inner_pos %s intersected outer_pos %s'%(pos, outer))
                pass_pos += [pos]
        for i in range(len(pass_pos))[::-1]:
            comp_pos += [sorted([outer[1]-(pass_pos[i][1]-outer[0]),
                                 outer[1]-(pass_pos[i][0]-outer[0])])]
    return sorted(comp_pos)

#pos is the disjoint pos list, limits is upper and lower bound
#usually will be set to = [0,len(seq)]
def complement_pos(pos,limits):
    return utils.LRF_1D(limits,pos)[2] #d1 is the sections of limits sliced away by the pos list

#given a mut_pos with generate a new pos with flanking probabilty mut_prob
def gen_mut_pos_flanking(mut_pos,mut_l,l_prob=0.25,r_prob=0.25,size_dist='uniform',lower_limit=10):
    if len(mut_pos)>0:
        flank_pos,minima,maxima = [],min([pos[0] for pos in mut_pos]),max([pos[1] for pos in mut_pos])
        for pos in mut_pos:
            mut_size = max(lower_limit,int(round(np.random.uniform(0.1*mut_l,10.1*mut_l),0)))
            if np.random.choice([True,False],p=[l_prob,1.0-l_prob]) and pos[0]>minima+mut_size:
                flank_pos += [[pos[0]-mut_size-1,pos[0]-1]]
            mut_size = max(lower_limit, int(round(np.random.uniform(0.1 * mut_l, 10.1 * mut_l), 0)))
            if np.random.choice([True,False],p=[r_prob,1.0-r_prob]) and pos[1]<maxima-mut_size:
                flank_pos += [[pos[1]+1,pos[1]+mut_size+1]]
        flank_pos = sorted(flank_pos)
        return utils.LRF_1D(mut_pos,flank_pos)[-1] #return only new area from mut_pos
    else:
        return []

def adjust_anueuploidy_effect(aneuploidy,ploidy,loss,loss_types=['DEL','INS','INV','TRA'],
                             driver='mitcp',dist='uniform'):
    if dist=='uniform':
        hits,total = {},{'effected':0.0,driver:0.0}
        for k in loss:
            if driver in loss[k]:
                hits[k] = {'effected':0.0,driver:0.0}
                for t in loss[k][driver][k]:
                    if t in loss_types:
                        hits[k]['effected'] = max(hits[k]['effected'],loss[k][driver][k][t][0])*1.0
                        hits[k][driver]     = max(hits[k][driver],loss[k][driver][k][t][1])*1.0
        for k in hits:
            total['effected'] += hits[k]['effected']
            total[driver]     += hits[k][driver]
        for k in aneuploidy:
            if k in ploidy:
                q   = aneuploidy[k][ploidy[k]]
                a_s = sorted(set(aneuploidy[k].keys()).difference(set([ploidy[k]])))
                if len(a_s)>0 and total[driver]>0:
                    p   = sum([aneuploidy[k][a] for a in a_s])
                    d   = min(1.0,pow(total['effected']/total[driver],0.5)*(1.0/(1.0*len(a_s))))
                    for a in a_s: aneuploidy[k][a] += d
                    aneuploidy[k][ploidy[k]]  = max(0.0,aneuploidy[k][ploidy[k]]-d*len(a_s))
                    x = 1.0*sum([aneuploidy[k][i] for i in aneuploidy[k]])
                    if x <= 0.0: x = 1.0
                    aneuploidy[k] = {i:aneuploidy[k][i]/x for i in aneuploidy[k]}
    return aneuploidy

#scane for
def adjust_del_effect(vca,driver='mmej'):
    return True

def adjust_ins_effect(vca,driver='nhej'):
    return True

#given a set of sequences and a sequency ploidy prob distribution
def gen_anueuploidy(seqs,ploidy,CT,anueuploidy): #aneuploidy: {k:{0:0.05,1:0.05,2:0.8,3:0.05,4:0.025,5:0.025}}
    A,clones = {},sorted(list(CT.freq.keys()),key=lambda x: int(x.rsplit('_')[-1]))
    for k in seqs:
        if k in anueuploidy:
            A[k] = []
            a = sorted(list(anueuploidy[k].keys()))
            p = [anueuploidy[k][s] for s in a]
            x = np.random.choice(a=a,p=p,size=1,replace=False)[0]
            start_clone = np.random.choice(a=clones)  #pick at random
            print('start clone %s'%start_clone)
            if x<ploidy[k]: #< 2 for 1-22 or < 1 for X,Y
                pos = [0,len(seqs[k])]
                ref = seqs[k][pos[0]]
                alt='<DEL>'
                if ploidy[k]<2: geno,base = '1','0'
                elif x==0:      geno,base = '1|1','0|0'
                else:           geno,base = np.random.choice(['1|0','0|1']),'0|0'
                genos = [base for i in range(len(clones))]
                #now update geno to include the CT tree structure given the starting clone id
                affected = CT.decendants(start_clone)+[start_clone]
                print('affected clones %s'%(affected,))
                for i in range(len(clones)):
                    if clones[i] in affected: genos[i] = geno
                print(genos)
                info='SVTYPE=DEL;END=%s;SVLEN=%s;'%(pos[1],pos[1]-pos[0])
                A[k] +=[gen_var_call(k,pos[0],ref,alt,info,'\t'.join(genos))]
            elif x>ploidy[k]: #> 2 for 1-22 or > 1 for X,Y
                pos = [0,len(seqs[k])]
                ref = seqs[k][pos[0]]
                alt = '<DUP>'
                if ploidy[k]<2: geno,base = '1','0'
                elif x>2:       geno,base = '1|1','0|0'
                else:           geno,base = np.random.choice(['1|0','0|1']),'0|0'
                genos = [base for i in range(len(clones))]
                #now update geno to include the CT tree structure given the starting clone id
                affected = CT.decendants(start_clone)+[start_clone]
                print('affected clones %s'%(affected,))
                for i in range(len(clones)):
                    if clones[i] in affected: genos[i] = geno
                print(genos)
                info='SVTYPE=DUP;END=%s;SVLEN=%s;DUP_POS=%s;DUP_TYPE=%s;ICN=%s'%\
                     (pos[1],pos[1]-pos[0],pos[0],'ANEUPLOIDY',x)
                A[k] += [gen_var_call(k,pos[0],ref,alt,info,'\t'.join(genos))]
    return A

#equally divide into p random partitions proportional to prop from a disjoint position list pos
def partition_pos(pos,p,prop=None,index_only=False):
    if prop is None:
        prop = [1.0/p for i in range(p)]
    elif sum(prop)>1.0:
        x = sum(prop)*1.0
        prop = [i/x for i in prop]
    S,P = {i:pos[i] for i in range(len(pos))},{}
    if p>0:
        if not index_only:
            if len(pos)>1:
                for x in range(p): #when this is finished pos will be the last partition
                    s = [i for i in S]
                    c = max(1,min(len(s),int(round(prop[x]*len(pos),0))))
                    if len(s)>0:
                        K = np.random.choice(s,c,replace=False)
                        P[x] = sorted([pos[k] for k in K])
                        for k in K: S.pop(k)
                if len(S)>0 and p>0:
                    ks,x = list(S.keys()),np.random.choice(range(p),1,replace=False)[0]
                    for k in ks:
                        P[x] += [S[k]]
                        S.pop(k)
            else:
                P[0] = [pos[i] for i in S]
        else:
            if len(pos)>1:
                for x in range(p): #when this is finished pos will be the last partition
                    s = [i for i in S]
                    c = max(1,min(len(s),int(round(prop[x]*len(pos),0))))
                    if len(s)>0:
                        K = np.random.choice(s,c,replace=False)
                        P[x] = sorted([k for k in K])
                        for k in K: S.pop(k)
                if len(S)>0 and p>0:
                    ks,x = list(S.keys()),np.random.choice(range(p),1,replace=False)[0]
                    for k in ks:
                        P[x] += [k]
                        S.pop(k)
            else:
                P[0] = [i for i in S]
    return P

#pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
def join_idx(C,j,p):
    I = {}
    for idx in [C[i][p] for i in j]:
        for k in idx:
            for i in idx[k]:
                if k in I: I[k].add(i)
                else:      I[k] =  {i}
    return I

#given two idx merge them into one
def merge_idx(idx1,idx2):
    I = {}
    for k in idx1:
        for i in idx1[k]:
            if k in I: I[k].add(i)
            else:      I[k] =  {i}
    for k in idx2:
        for i in idx2[k]:
            if k in I: I[k].add(i)
            else:      I[k] =  {i}
    return I

#returns G[g][d][C[g]=>i, C[g][wx]]
#as well as a flattened combination list for use
#key f G are the coordinates x1,x2
def weight_graph(C):
    G,X = {},[]
    c_i,i_c = group_indecies(C)
    if len(i_c)>0: #load all unique verticies
        for c in sorted(list(C.keys())): #for each classifier
            for i in range(len(C[c])): #for every entry in every classifier
                for j in range(2):     #for the start and stop point in the pos (range)
                    if C[c][i][j] in G:
                        if c in G[C[c][i][j]]:
                            if j in G[C[c][i][j]][c]: G[C[c][i][j]][c][j] += [[i+c_i[c],C[c][i][2]]]
                            else:                     G[C[c][i][j]][c][j]  = [[i+c_i[c],C[c][i][2]]]
                        else:                         G[C[c][i][j]][c]  = {j:[[i+c_i[c],C[c][i][2]]]}
                    else:                             G[C[c][i][j]]  = {c:{j:[[i+c_i[c],C[c][i][2]]]}}
            X += C[c] #load the entry into X which will then get sorted by position
        last = sorted(list(G.keys()))[-1]+2          #get last x2 value in G
        G[last] = {(None,):{None:[None,None]}} #to pad out a terminal vertex
        #check all indecies are paired
        V = {}
        for i in sorted(G)[:-1]:
            for c in G[i]:
                for e in G[i][c]:
                    for l in G[i][c][e]:
                        if l[0] in V: V[l[0]] += 1
                        else:         V[l[0]]  = 1
        #print('all vertices and edges added : %s'%all([V[i]==2 for i in V]))
    return X,G,c_i,i_c

#given the sorted keys of C
#calculate the offset you need to add to each groups
#index to resolve the correct row entry
def group_indecies(C):
    c_i,i_c,i = {},{},0
    for g in sorted(list(C.keys())):
        c_i[g]=i
        for j in range(i,i+len(C[g])): i_c[j]=g
        i+=len(C[g])
    return c_i,i_c

#return a list of indecies into X, given an edge in G
def get_x_i(A):
    x_i = set([])
    for c in A:
        for e in A[c]:
            for l in A[c][e]: x_i.add(l[0])
    return sorted(list(x_i))

#clear off terminal edges of A and clear keys as needed
def del_edges(A):
    k = list(A.keys())
    for i in range(len(k)): #take value 1 keys for each class
        if 1 in A[k[i]]:
            while len(A[k[i]][1]) > 0:
                A[k[i]][0].remove(A[k[i]][1][0])
                A[k[i]][1].remove(A[k[i]][1][0])
            if len(A[k[i]][0])<1: A[k[i]].pop(0)
            if len(A[k[i]][1])<1: A[k[i]].pop(1)
        if len(A[k[i]])<1: A.pop(k[i])

#active edges A, next vertext edges B
#add new edges into A, allowing a single class to overlap itself
#add class memebers to overlap other class members as well
def add_edges(A,B):
    for k in B:
        if k in A:
            for e in B[k]:
                if e in A[k]: A[k][e] += [l for l in B[k][e]]
                else:         A[k][e]  = [l for l in B[k][e]]
        else:                 A[k] = {e:[l for l in B[k][e]] for e in B[k]}

#return the set of weights from X still in coordinate sorted order
def get_weights(X):
    W = {}
    if len(X)>0:
        for i in range(len(X)):   #generate disjoint weight and length (area)
            c = list(X[i][3].keys())[0] #assume single classes to start
            if c in W:
                for g in X[i][3][c]:
                    W[c][g] = [X[i][2],1.0*(X[i][1]-X[i][0]+1.0)] #static positions here
            else:
                for g in X[i][3][c]:
                    W[c] = {g:[X[i][2],1.0*(X[i][1]-X[i][0]+1.0)]}
    return W

#scaling routines to apply tranformed values into the projected segments in P
def scale_weights(P,W,scale=[0.0,10.0],method='add-mult'):
    a,b,x,U,scale = 1.0E100,0.0,0.0,[],[min(scale),max(scale)]
    if method=='add':
        for i in range(len(P)):
            U += [sum([sum([W[c][g][0] for g in P[i][3][c]]) for c in P[i][3]])]
            if U[i]<a: a = U[i] #min
            if U[i]>b: b = U[i] #max
        x = float(scale[1]-scale[0])/float(b-a+1.0)
        for i in range(len(P)):
            P[i][2] = x*(U[i]-a)+scale[0]
    elif method=='len':
        for i in range(len(P)):
            U += [sum([sum([W[c][g][0]*W[c][g][1]/(P[i][1]-P[i][0]+1.0) for g in P[i][3][c]]) for c in P[i][3]])]
            if U[i]<a: a = U[i] #min
            if U[i]>b: b = U[i] #max
        x = float(scale[1]-scale[0])/float(b-a+1.0)
        for i in range(len(P)):
            P[i][2] = x*(U[i]-a)+scale[0]
    elif method=='log':
        for i in range(len(P)):
            U += [sum([sum([W[c][g][0] for g in P[i][3][c]]) for c in P[i][3]])]
            if U[i]<a: a = U[i] #min
            if U[i]>b: b = U[i] #max
        x = float(scale[1]-scale[0])
        for i in range(len(P)):
            P[i][2] = x*(np.log10(1.0+(U[i]/a))/np.log10(1.0+(b-a)))+float(scale[0])
    elif method=='mult':
        for i in range(len(P)):
            U += [np.product([np.product([W[c][g][0] for g in P[i][3][c]]) for c in P[i][3]])]
            if U[i]<a: a = U[i] #min
            if U[i]>b: b = U[i] #max
        if b-a!=0.0: x = float(scale[1]-scale[0])/float(b-a)
        else:        x = float(scale[1]-scale[0])
        for i in range(len(P)):
            P[i][2] = x*(U[i]-a)+scale[0]
    elif method=='add-mult':
        for i in range(len(P)):
            U += [sum([sum([W[c][g][0] for g in P[i][3][c]]) for c in P[i][3]])]
            if U[i]<=a: a = U[i] #min
            if U[i]>=b: b = U[i] #max
        x = float(scale[1]-scale[0])/float(b-a+1.0)
        for i in range(len(P)):
            P[i][2] = x*(U[i]-a)+scale[0]
        A = []
        for i in range(len(P)):
            A += [np.product([np.product([W[c][g][0] for g in P[i][3][c]]) for c in P[i][3]])]
            if A[i]<=a: a = U[i] #min
            if A[i]>=b: b = U[i] #max
        if b-a!=0.0: x = float(scale[1]-scale[0])/float(b-a)
        else:        x = float(scale[1]-scale[0])
        for i in range(len(P)):
            P[i][2] *= x*(A[i]-a)+scale[0]
    return P

#given a graph G = X,C,c_i,i_c where X is a coordinate sorted concat of all original rows
#and C is the graph edges that correspond to each unique vertex in X with c_i and i_c
#mappings from each set of indecies into the other, scan and project into 1D the unique
#segments or areas that overlap, while maintaining the idx = class:{name1,name2,..}
#on each segment the weights here are cleared out in the resulting solution P
#so that a sperate scaled weighting can be completed on the projected segments
def scan_graph(G,scale=True,filter=None): #while offer additional weighting corrections
    X,C,c_i,i_c = G
    P,Q,V,A,B,D = [],[],[],{},{},set([]) #A is alive edge stack, B is the current edge set
    if len(X)>0:
        V = sorted(list(C.keys()))  #V is sorted vertices, where disjoint => |V|-1 = |X|*2
        for i in range(0,len(V)-1):   #scan i-1,i,i+1 up to padding
            B = C[V[i]]               #get edges for v[i+1]
            D = {m for l in [C[V[i]][k] for k in C[V[i]]] for m in l}  #check edge direction for i
            if len(A) <= 0:         #[1] len(a) <= 0 (f has 0 edge)   #section starting, start new  p+=[]
                add_edges(A,B)
                x_i =  get_x_i(A)
                P += [[V[i],V[i],0.0,join_idx(X,x_i,3)]]
                if D == set([0,1]):   #check for singles
                    del_edges(A)      #clean the singles
                    if len(A)>0:      #with new sections that are not singles
                        x_i =  get_x_i(A)
                        P += [[V[i]+1,V[i]+1,0.0,join_idx(X,x_i,3)]]
            else:
                if D == set([0]):   #[2] len(a) > 0 and f has 0 edge  #close subsection, start new  p[-1],p+=[]
                    x_i =  get_x_i(A)
                    P[-1][1] = V[i]-1
                    P[-1][3] = merge_idx(P[-1][3],join_idx(X,x_i,3))
                    del_edges(A)      #clean the closed edges
                    add_edges(A,B)    #start the new section
                    x_i =  get_x_i(A)
                    P += [[V[i],V[i],0.0,join_idx(X,x_i,3)]]
                if D == set([0,1]): #[3] len(a) > 0 and f has 0 and 1 #close subsection, set single p[-1],p+=[]
                    x_i =  get_x_i(A)
                    P[-1][1] = V[i]-1 #close the last open section
                    P[-1][3] = merge_idx(P[-1][3],join_idx(X,x_i,3))
                    add_edges(A,B)
                    x_i =  get_x_i(A)
                    P += [[V[i],V[i],0.0,join_idx(X,x_i,3)]]
                    del_edges(A)
                    if len(A)>0:
                        x_i =  get_x_i(A)
                        P += [[V[i]+1,V[i]+1,0.0,join_idx(X,x_i,3)]]
                if D == set([1]):   #[4] len(a) > 0 and f has 1       #section closing,  fix last   p[-1]
                    x_i =  get_x_i(A)
                    P[-1][1] = V[i]   #close and clean the last section
                    P[-1][3] = merge_idx(P[-1][3],join_idx(X,x_i,3))
                    add_edges(A,B)
                    del_edges(A)
                    if len(A)>0:      #find any remaining open sections
                        x_i =  get_x_i(A)
                        P += [[V[i]+1,V[i]+1,0.0,join_idx(X,x_i,3)]]
    for p in P: #clean up dangling sections
        if p[1]<p[0]: P.remove(p)
    if scale:
        P = scale_weights(P,get_weights(X))
    if filter is not None:
        for p in P:
            if p[2]>filter: Q += [p]
    return P
#pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp

#given a position list weight w and label
def vcam_to_wcu(vcam,w=0.0,label=None):
    wcu = {}
    for k in vcam:
        wcu[k] = []
        for vc in vcam[k]:
            wcu[k] += [[vc.pos,get_info_end(vc.info),w,
                        {label:set([get_info_type(vc.info)])}]]
    return wcu

def wcu_coord_clean(wcu):
    for k in wcu:
        for i in range(len(wcu[k])):
            wcu[k][i][0],wcu[k][i][1] = sorted([wcu[k][i][0],wcu[k][i][1]])
        wcu[k] = sorted(wcu[k],key=lambda x: x[0])
    return wcu

#convert a weight graph above to a pos list and weight list
def wcu_to_pos_w(wcu,l,scale=[0.0,10.0],method='mult',background=1.0):
    wcu = wcu_coord_clean(wcu)
    G = weight_graph(wcu)
    try:
        P = scan_graph(G,scale=False)
    except Exception as E:
        pad = ':'.join(['-' for i in range(20)])
        print(pad+'@@scan graph error@@'.upper()+pad)
        ru.dump_state({'wcu':wcu,'l':l,'G':G},'wcu_to_pos_w-scan_graph',os.path.expanduser('~/'))

    if background>0.0:
        wcu['background'] = [[0,l,background,{'background':set(['background'])}]]
        G = weight_graph(wcu)
        P = scan_graph(G,scale=False)
    W = get_weights(G[0])
    P = scale_weights(P,W,scale=scale,method=method)
    pos,w = [],[]
    for i in range(len(P)):
        pos += [[P[i][0],P[i][1]]]
        w   += [P[i][2]]
    return np.array(pos,dtype=int),np.array(w,dtype=float)

#input is a vcam of positions
#output is a pos_mapping for open positions
#use this to feed into the somatic genomes to
#get started in generating valid positions via wcu distributions
#gs is the sequence and length map
def open_pos(mnv_vcam, sv_vcam, gs):
    mnv_pos_map = vcam_to_pos_map(mnv_vcam)
    sv_pos_map  = vcam_to_pos_map(sv_vcam)
    open_mnv = {k:[[0,gs[k]]] for k in gs}
    open_sv  = {k:[[0,gs[k]]] for k in gs}
    for k in gs:
        if k in mnv_pos_map:
            open_mnv[k] = utils.LRF_1D(open_mnv[k],mnv_pos_map[k])[2]
        if k in sv_pos_map:
            open_sv[k]  = utils.LRF_1D(open_sv[k],sv_pos_map[k])[2]
    return open_mnv,open_sv

#given a wcu and a list of open positions
#merge and reweight the wcu such that the closed
#positions have a 0.0 probabilty wieghting once scaled
def update_wcu(wcu,open_pos):
    return True

#check the class proportion or saturation
#given a wcu that already has the classes embeded and a mut_pos_map M
#check the number of alterations in M that hit members of each class
def prop_class_mut_pos(M,wcu,gene_map):
    gen_wcu,wcu_i,c,B = {}, {}, get_wcu_class(wcu),{}  # get first class
    if len(c)>0:
        c = c[0]
        for k in M:
            gen_wcu[k],wcu_i[k],B[k] = {},{},{}
            for t in M[k]:
                gen_wcu[k][t] = [[M[k][t][i][0], M[k][t][i][1], 0.0, {t: set(['%s_%s' % (t, i + 1)])}] for i in range(len(M[k][t]))]
                wcu_i[k][t] = {}  # wcu can be genes and gens can overlap so you may get more than one index back
                if k in wcu:
                    for i in range(len(wcu[k])):
                        idx = wcu[k][i][3]
                        for g in idx[c]:
                            if g in wcu_i[k][t]: wcu_i[k][t][g] += [i]
                            else:                wcu_i[k][t][g]  = [i]
        I = {}
        for k in gen_wcu:
            for t in gen_wcu[k]:
                if t in I: I[t][k] = gen_wcu[k][t]
                else:      I[t] = {k:gen_wcu[k][t]}
        I.update({c:wcu})
        merge = merge_wcu(I)  # this does all k in wcu, M
        for k in M:  # only update those events that intersect with the class to update them
            G = weight_graph(merge[k])
            P = scan_graph(G)
            for t in M[k]:
                C,T,D,N = {},{},{},[]
                for p in P:
                    if c in p[3] and t in p[3]:  # classes intersect
                        for i in p[3][c]:  # for each gene in the class
                            for j in p[3][t]:  # map each event to it
                                if i in C: C[i].add(j)
                                else:      C[i] = set([j])
                                if j in T: T[j].add(i)
                                else:      T[j] = set([i])
                for j in T:
                    if len(T[j]) > 1:
                        T[j] = set([T[j].pop()])
                for i in C:
                    if len(C[i]) > 1:  # grab the first event only
                        for e in C[i]:
                            if i in T[e]:
                                D[i] = set([C[i].pop()])
                                break
                    else: D[i] = C[i]
                if k in wcu: B[k][t] = [len(D),len(wcu[k]),len(M[k][t]),1.0]
                else:        B[k][t] = [0,0,len(M[k][t]),1.0]
            if k in gene_map['seq']:
                B[k]['genes'] = len(gene_map['seq'][k])
    else:
        for k in M:
            B[k] = {}
            if k in gene_map['seq']:
                B[k]['genes'] = len(gene_map['seq'][k])
    return B

#returns the classes of the wcu
def get_wcu_class(wcu):
    c = set([])
    for k in wcu:
        for i in wcu[k]:
            for j in i[3]:
                c.add(j)
    return sorted(list(c))

#given a one-class wcu, check the mut_pos -> M[k][t]
#and extend/merge areas to help simulate full sections of gain/loss
def extend_class_mut_pos(M,t,wcu):
    M2,gen_wcu,wcu_i,c = copy.deepcopy(M),{},{},get_wcu_class(wcu) #get first class
    if len(c)>0:
        c = c[0]
        for k in M:
            if t in M[k]:
                gen_wcu[k] = [[M[k][t][i][0],M[k][t][i][1],0.0,{t: set(['%s_%s'%(t,i+1)])}] for i in range(len(M[k][t]))]
                wcu_i[k] = {} #wcu can be genes and gens can overlap so you may get more than one index back
                for i in range(len(wcu[k])):
                    idx = wcu[k][i][3]
                    for g in idx[c]:
                        if g in wcu_i[k]: wcu_i[k][g] += [i]
                        else:             wcu_i[k][g]  = [i]
        intersect = merge_wcu({c:wcu,t:gen_wcu}) #this does all k in wcu, M
        for k in M: #only update those events that intersect with the class to update them
            C,T,D,N,I = {},{},{},[],[]
            G = weight_graph(intersect[k])
            P = scan_graph(G)
            for p in P:
                if c in p[3] and t in p[3]: #classes intersect
                    for i in p[3][c]:         #for each gene in the class
                        for j in p[3][t]:           #map each event to it
                            if i in C: C[i].add(j)
                            else:      C[i] = set([j])
                            if j in T: T[j].add(i)
                            else:      T[j] = set([i])
            for j in T:
                if len(T[j])>1:
                    T[j] = set([T[j].pop()])
            for i in C:
                if len(C[i])>1: #grab the first event only
                    for e in C[i]:
                        if i in T[e]:
                            D[i] = set([C[i].pop()])
                            break
                else: D[i] = C[i]
            for d in D:
                e = list(D[d])[0]
                N += [[wcu[k][i][0],wcu[k][i][1]] for i in wcu_i[k][d]]
            for p in P:
                if len(p[3])<=1 and t in p[3]:
                    I += [int(x.split('_')[-1])-1 for x in p[3][t]]
            for i in I: N += [M[k][t][i]]
            M2[k][t] = sorted(N,key=lambda x: x[0])
        return M2
    else:
        return M

#given a set of chrom and pos for each for TRA distributed them into a tra_map
#using ranomd uniform samping across all chrom routing pairs
#tra_dict = {'chr11':[[0,100],...],'chr22':[]}
def gen_tra_map(tra_dict):
    T = {}
    while len(list(tra_dict.keys()))>=2:
        i = np.random.choice(range(len([0 for a,b in it.permutations(list(tra_dict.keys()),2)])))
        j,k = [(a,b) for a,b in it.permutations(list(tra_dict.keys()),2)][i]
        if len(tra_dict[j])>=1 and len(tra_dict[k])>=1:
            if j in T:
                route = {'SCHR':k,'SPOS':[],'DPOS':[]}
                route['DPOS'] = tra_dict[j].pop(np.random.choice(range(len(tra_dict[j]))))
                route['SPOS'] = tra_dict[k].pop(np.random.choice(range(len(tra_dict[k]))))
                T[j] += [route]
            else:
                route = {'SCHR':k,'SPOS':[],'DPOS':[]}
                route['DPOS'] = tra_dict[j].pop(np.random.choice(range(len(tra_dict[j]))))
                route['SPOS'] = tra_dict[k].pop(np.random.choice(range(len(tra_dict[k]))))
                T[j] = [route]
        else:
            if len(tra_dict[j])<1: tra_dict.pop(j)
            if len(tra_dict[k])<1: tra_dict.pop(k)
    #sort each by the DPOS of the key
    for k in T:
        T[k] = sorted(T[k],key=lambda x: x['DPOS'])
    return T

def get_num_mut(full_len, mut_len, mut_rate):
    return int(full_len*mut_rate/mut_len)

#given a
def gen_inner_var_calls():
    return []

def gen_genotypes(ploidy,pos,hetro=0.15):
    if ploidy==1:
        G = [[1] for i in range(len(pos))]
    elif ploidy==2:
        G = [[[0,1],[1,0],[1,1]][i] for i in np.random.choice(range(3),len(pos),p=[hetro/2.0,hetro/2.0,1.0-hetro])]
        if hetro<=0.0: #ensure no rounding issues on random selection
            G = [[1,1] for i in range(len(pos))]
    elif ploidy==3:
        a = [[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]]
        p = [1,1,2,1,2,2,3] #have to compute for tri-ploidy
        G = [a[i] for i in np.random.choice(range(7),len(pos))]
        if hetro<=0.0: #ensure no rounding issues on random selection
            G = [[1,1,1] for i in range(len(pos))]
    else:
        print('ploidy counts over three are not currently supported')
        raise AttributeError
    return G

#update each frmat field to have n identical genotypes
def update_germline_genotypes(vcam,n,gen_delim='\t'):
    for k in vcam:
        for i in range(len(vcam[k])):
            g = [x for x in vcam[k][i].frmat.split('|')]
            vcam[k][i].frmat = gen_delim.join(['|'.join(g) for x in range(n)])
    return vcam

#given a list of positions and the genotypes G,
#and a list of inner positions I, associate and pass the genotype on
def get_inner_genotypes(mut_pos,genotypes,inner_pos):
    G = []
    if len(genotypes)>0:
        G = [[0 for i in range(len(genotypes[0]))] for j in range(len(inner_pos))]
        for i in range(len(inner_pos)):
            for j in range(len(mut_pos)):
                if intersect(inner_pos[i],mut_pos[j]):
                    G[i] = genotypes[j]
    return G

#using a simulated CloneTree as input T
#can make use of the nodes decendant function
#to set up partitions and alleles across one vcam pos mapping
#part_map = {p1:[1,5,6,8,9,10], p2:[2,3,4,7], ...} have to use index_only=True with partition_pos
#AM = {node_id:set([1,2,3,4,5,6,7,...]), node_id:set([1,2,3,4]), ...}
def gen_clone_allele_map(part_map,CT):
    A,AM,clones = {},{},sorted(list(CT.nodes.keys()),key=lambda x: int(x.split('_')[-1]))
    for clone in clones: #get the decendants of each clone
        A[clone]  = sorted(CT.ancestors(clone),key=lambda x: int(x.split('_')[-1]))
        AM[clone] = []
    for i in range(len(clones)):
        if i in part_map:
            if clones[i] in AM: AM[clones[i]] += copy.deepcopy(part_map[i])
            else:               AM[clones[i]]  = copy.deepcopy(part_map[i])
    for clone in clones: #if there are not enough calls to share for partitions
        for k in A[clone]: AM[clone] += AM[k]
    for clone in clones:   AM[clone] = set(AM[clone])
    return AM

#[2] for each position in the uniform distribution apply the edit
#    and update the current position (DEL = -L, INS = +L, SUB = 0, INV = 0, DUP = L*CN)
#    the final result should throw and error or not tabulate the edit
#    if the simple merging step consumes preselected positions
#take in the hetro value and use it to add the GT info in the FORMAT VCA field
def gen_var_calls(refdict,chrom,mut_type,mut_pos,genotype,
                  tra_map=None,allele_map=None,insert_seqs=None,
                  small=int(1E2),verbose=False,ref_path=None,gen_delim='\t'):
    full_ref = refdict[chrom]
    vca,info,i = [],'',0     #init position list pointer
    if len(genotype)>0:
        ploidy = len(genotype[0])
        if allele_map is not None:
            empty = '|'.join([str(0) for i in range(ploidy)])
            clone_order = sorted(allele_map,key=lambda x: int(x.split('_')[-1]))
    while i < len(mut_pos): #now check positions for 'N' masked content....
        pos = mut_pos[i]
        ref = full_ref[pos[0]:pos[1]]
        geno = '|'.join([str(g) for g in genotype[i]])
        if allele_map is not None: # 0|1, 0|1, 0|0, ...
            clone_geno = []
            for clone in clone_order:
                if i in allele_map[clone]: clone_geno += [geno]
                else:                      clone_geno += [empty]
            geno = gen_delim.join(clone_geno)
        #check the proportion of masked content > 0.5
        if ref.count('N')>0.5*len(ref) and len(ref)>small: #this gets checked by mut_pos too...
            if verbose: print('ref region is >0.5 proportion masked at seq=%s,pos=%s-%s'%(chrom,pos[0],pos[1]))
            i += 1
        else:
            if type(mut_type) is dict and 'DUP' in mut_type: #only apply dups at pos or prior to pos for distal...
                dup_type = np.random.choice(list(mut_type['DUP']['TYPE']))
                d_cn     = np.random.choice(list(mut_type['DUP']['CN'])) #randomly CN is the number of duplications CN2 is a ICN==4
                alt      = '<DUP>'
                if dup_type=='DISPERSED': #injected backward for clarity and correctness
                    if i+d_cn<len(mut_pos):
                        d_pos   = mut_pos[i:i+d_cn]
                        pos     = d_pos[-1]      #right most is the dup area on reference
                        ref     = full_ref[pos[0]:pos[1]]
                        d_pos   = d_pos[:-1] #injected to these areas maintaining disjoint adherence
                        dup_pos = ','.join([str(j[0]) for j in d_pos])
                        info = 'SVTYPE=DUP;END=%s;SVLEN=%s;DUP_POS=%s;DUP_TYPE=%s;ICN=%s'%\
                               (pos[1],pos[1]-pos[0],dup_pos,dup_type,2*d_cn)
                        vca += [gen_var_call(chrom,pos[0],ref,alt,info,geno)]
                    else: #dupication boundry would overwrite the current buffer
                        if verbose: print('not able to generate DUP: %s-%s'%(pos[0],info)) #will only print once
                    i += d_cn #consume the dup and all DIS dups
                else: #TANDEM or DEL
                    dup_pos = pos[0]
                    info = 'SVTYPE=DUP;END=%s;SVLEN=%s;DUP_POS=%s;DUP_TYPE=%s;ICN=%s'%\
                           (pos[1],pos[1]-pos[0],dup_pos,dup_type,2*d_cn)
                    vca += [gen_var_call(chrom,pos[0],ref,alt,info,geno)]
                    i   += 1 #consume just the TANDEM DUP
            elif type(mut_type) is dict and 'TRA' in mut_type:   #SHOULD ENABLE DEL OR INS TYPES
                tra_type = np.random.choice(mut_type['TRA']['TYPE'])
                if tra_type=='DEL':
                    a,b = pos[0],pos[1] #svlen = length of missing CHR1 section
                else: #is an INS
                    a,b = pos[0],pos[0] #svlen = 0 => orginial is stil there and CHR2 was inserted
                    ref = 'N'
                if ref_path is not None: #now you have the option to load only the base chrom and fetch TRA inserts
                    alt    = ru.read_fasta_chrom_pos(ref_path,tra_map[i]['SCHR'],
                                                     tra_map[i]['SPOS'][0],tra_map[i]['SPOS'][1])
                else:
                    alt    = refdict[tra_map[i]['SCHR']][tra_map[i]['SPOS'][0]:tra_map[i]['SPOS'][1]]
                info = 'SVTYPE=TRA;END=%s;SVLEN=%s;TRA_TYPE=%s;CHR2=%s;POS2=%s;END2=%s'%\
                       (b,b-a,tra_type,
                        tra_map[i]['SCHR'],tra_map[i]['SPOS'][0],tra_map[i]['SPOS'][1])
                vca  += [gen_var_call(chrom,a,ref,alt,info,geno)]
                i    += 1
            elif mut_type=='SUB':
                alt  = gen_alt_seq(ref)
                info = 'SVTYPE=SUB;END=%s;SVLEN=%s;'%(pos[1],pos[1]-pos[0])
                vca += [gen_var_call(chrom,pos[0],ref,alt,info,geno)]
                i   += 1
            elif mut_type=='INS':
                if insert_seqs is not None:
                    alt = insert_seqs[i] #should be a map from the vcam index to the string
                else:
                    alt  = gen_alt_seq(ref)
                ref  = 'N'
                info = 'SVTYPE=INS;END=%s;SVLEN=%s;'%(pos[0],pos[1]-pos[0])
                vca += [gen_var_call(chrom,pos[0],ref,alt,info,geno)]
                i   += 1
            elif type(mut_type) is dict and 'INV' in mut_type:
                alt  = utils.get_reverse_complement(ref)
                info = 'SVTYPE=INV;END=%s;SVLEN=%s;INV_TYPE=%s'%\
                       (pos[1],pos[1]-pos[0],np.random.choice(mut_type['INV']['TYPE']))
                vca += [gen_var_call(chrom,pos[0],ref,alt,info,geno)]
                i   += 1
            elif mut_type=='DEL':
                alt  = '<DEL>'
                info = 'SVTYPE=DEL;END=%s;SVLEN=%s;'%(pos[1],pos[1]-pos[0])
                vca += [gen_var_call(chrom,pos[0],ref,alt,info,geno)]
                i   += 1
    return vca

def vca_to_pos_list(vca):
    pos_list = []
    for vc in vca:
        pos_list += [[vc.pos,get_info_end(vc.info)]]
    return pos_list

def vcam_to_pos_map(vcam):
    pos_map = {}
    for k in vcam:
        pos_map[k] = []
        for vc in vcam[k]:
            pos_map[k] += [[vc.pos,get_info_end(vc.info)]]
    return pos_map

#given a possible empty string = ref and a length = default is 0
#either produce a new string that has a different random nuclideotide
#for use as a substitution generator or just a full random insertion
def gen_alt_seq(ref,alphabet=['A','C','G','T'],length=0):
    alt = ''
    l = len(ref)
    sigma = set(alphabet) #do set difference for INS/SUB
    for i in range(0,l):  #have to do this for each position
        cand = np.sort(list(sigma-set([ref[i]]))) #candidates to make a sigle mut
        alt += np.random.choice(cand)
    return alt

#one sigma look-back
def gen_alt_seq2(ref,length=0):
    alt = ''
    l = len(ref)
    dna = set(['A','C','G','T']) #do set difference for INS/SUB
    if l>0: alt += np.random.choice(list(dna-set([ref[0]])))[0] #get one new one
    for i in range(1,l):         #have to do this for each position
        cand = np.sort(list(dna-set([ref[i]]))) #candidates to make a snp
        prob = [0.475,0.475,0.475]     #one will be set to 0.2 to penalize similiar
        for j in range(0,3):     #slight random selection
            if alt[i-1]==cand[j]: prob[j]=0.05  #adjustment via single lookback chain
        alt += np.random.choice(cand) #get one new one
    if l<=0: alt = np.random.choice(list(dna),size=length).tostring()
    return alt

#[3] Keep track pf the applyed edits by making a VariantCall
#    alt, chrom, filter, format, id, info, pos, qual, ref, samples
def gen_var_call(chrom,pos,ref,alt,info,frmt):
    vc = VariantCall(chrom,pos,'.',ref,alt,None,None,info,frmt)
    return vc

def get_genotype_vca(vca,g=0,index=0):
    g_vca = []
    for i in range(len(vca)):
        if get_genotype(vca[i].frmat,index)[g]==1:
            g_vca += [vca[i]]
    return g_vca

def set_genotype(vca,geno):
    g_vca = copy.deepcopy(vca)
    for i in range(len(vca)):
        g_vca[i].frmat = '|'.join(geno)
    return g_vca

def dup(s,i):
    return ''.join([s for j in range(i)])

def is_dup(info):
    if info.rsplit('SVTYPE=')[-1].rsplit(';')[0]=='DUP':
        return True
    else:
        return False

def get_info_type(info):
    sv_type = info.rsplit('SVTYPE=')[-1].rsplit(';')[0]
    return sv_type

def get_info_end(info):
    end = int(info.rsplit('END=')[-1].rsplit(';')[0])
    return end

#works on vc.info specifically
def set_info_end(info,end):
    info_s = info.split('END=')
    return info_s[0]+'END=%s;'%end+';'.join(info_s[-1].rsplit(';')[1:-1])

def set_info_len(info,svlen):
    info_s = info.split('SVLEN=')
    return info_s[0]+'SVLEN=%s;'%svlen+';'.join(info_s[-1].rsplit(';')[1:-1])

def get_info_len(info):
    svlen = abs(int(info.rsplit('SVLEN=')[-1].rsplit(';')[0]))
    return svlen

def get_genotype(frmat,index=0):
    form = frmat.split('\t')[index]
    if form.find(r'|')>0:    geno = [int(i) for i in form.split(r'|')]
    elif form.find(r'/')>0:  geno = [int(i) for i in form.split(r'/')]
    else:                    geno = [int(form.split(r'|')[0].split(r'/')[0])]
    return geno

#use the alt and info fields to return an actionable data structure
def get_dup(info): #we assume ICN / len(DUP_POS) >= 3...
    icn    = info.rsplit('ICN=')[-1].rsplit(';')[0]
    d_pos  = info.rsplit('DUP_POS=')[-1].rsplit(';')[0].rsplit(',')
    d_type = info.rsplit('DUP_TYPE=')[-1].rsplit(';')[0]
    return {'TYPE':d_type,'CN':int(icn),'POS':[int(d) for d in d_pos]}

def get_inv(info):
    return info.rsplit('INV_TYPE=')[-1].rsplit(';')[0]

def get_tra(info):
    tra_type = info.rsplit('TRA_TYPE=')[-1].rsplit(';')[0]
    tra_chrom = info.rsplit('CHR2=')[-1].rsplit(';')[0]
    tra_pos   = int(info.rsplit('POS2=')[-1].rsplit(';')[0])
    tra_end   = int(info.rsplit('END2=')[-1].rsplit(';')[0])
    return {'TYPE':tra_type,'CHR2':tra_chrom,'POS2':tra_pos,'END2':tra_end}

#COMPLEX ranges are used for additive SV inside the boundries
def get_complex_inv_pos(vca):
    complex_list = []
    for i in range(len(vca)):
        if get_info_type(vca[i].info)=='INV' and get_inv(vca[i].info)=='COMPLEX':
            complex_list += [vca[i]]
    return vca_to_pos_list(complex_list)

#we assume the disjoint vca is applied and then adjust the pos_list
def update_pos_list(old_list,vca,g=0,index=0): #g is the genotype position either 0 or 1
    pos_list = copy.deepcopy(old_list)
    for i in range(len(vca)):
        sv_type,sv_len = get_info_type(vca[i].info),get_info_len(vca[i].info)
        c1 = [vca[i].pos,get_info_end(vca[i].info)]
        if get_genotype(vca[i].frmat,index)[g]==1: #allele is present
            if sv_type=='INS': #add to pos_list where needed
                for pos in pos_list: #ins to left will update pos_list to right
                    if pos[0]>c1[0]:
                        pos[0] += sv_len
                        pos[1] += sv_len
            elif sv_type=='DUP': #{'CN':int(icn),'POS':[int(d) for d in d_pos]}
                dup = get_dup(vca[i].info)
                if len(dup['POS'])>=1 and dup['POS'][0]!=vca[i].pos:#this one is dispersed---
                    for pos in pos_list:
                        for d in dup['POS']:
                            if pos[0]>d:     #doesn't split boundries
                                pos[0] += sv_len #assume one disjoint
                                pos[1] += sv_len #assume one disjoint
                elif len(dup['POS'])>=1: #this one is the tandem---------------------------
                    for pos in pos_list:
                        for d in dup['POS']:
                            if pos[0]>d:  #doesn't split boundries
                                pos[0] += sv_len*(dup['CN']/2-1) #assume all tandem CN
                                pos[1] += sv_len*(dup['CN']/2-1) #assume all tandem CN
            elif sv_type=='TRA':
                tra = get_tra(vca[i].info)
                tra_len = tra['END2']-tra['POS2']
                for pos in pos_list:
                    if pos[0]>c1[0]:
                        pos[0] += tra_len-sv_len #could be different sizes...
                        pos[1] += tra_len-sv_len #could be different sizes...
            elif sv_type=='DEL':                #sub from pos_list where needed
                for pos in pos_list:
                    if pos[0]>c1[0]:           #doesn't split boundries
                        pos[0] -= sv_len
                        pos[1] -= sv_len
    return pos_list

#update vca_a positions to the new coordinate space
#made by the application of the vc in vca_b
#similiar in principle to the update_pos_list function above
def update_vca_pos(vca_a,vca_b,g=0,index=0):
    vca = copy.deepcopy(vca_a)
    for i in range(len(vca_b)):
        sv_type, sv_len = get_info_type(vca_b[i].info), get_info_len(vca_b[i].info)
        c1 = [vca_b[i].pos,get_info_end(vca_b[i].info)]
        if get_genotype(vca_b[i].frmat, index)[g] == 1:              # allele is present
            if sv_type == 'INS':                          # add to pos_list where needed
                for j in range(len(vca)):  # ins to left will update pos_list to right
                    if vca[j].pos > c1[0]:
                        vca[j].pos += sv_len
                        vca[j].info = set_info_end(vca[j].info,
                                                   get_info_end(vca[j].info) + sv_len)
            elif sv_type == 'DUP':  # {'CN':int(icn),'POS':[int(d) for d in d_pos]}
                dup = get_dup(vca_b[i].info)
                if len(dup['POS']) >= 1 and dup['POS'][0] != vca_b[i].pos:  # this one is dispersed---
                    for j in range(len(vca)):
                        for d in dup['POS']:
                            if vca[j].pos > d:        # doesn't split boundries
                                vca[j].pos += sv_len      # assume one disjoint
                                vca[j].info = set_info_end(vca[j].info,
                                                           get_info_end(vca[j].info) + sv_len)
                elif len(dup['POS']) >= 1:  # this one is the tandem---------------------------
                    for j in range(len(vca)):
                        for d in dup['POS']:
                            if vca[j].pos > d:                     # doesn't split boundries
                                vca[j].pos += sv_len*(dup['CN']/2-1)  # assume all tandem CN
                                vca[j].info = set_info_end(vca[j].info,
                                                           get_info_end(vca[j].info)+sv_len*(dup['CN']/2-1))
            elif sv_type == 'TRA':
                tra = get_tra(vca_b[i].info)
                tra_len = tra['END2']-tra['POS2']
                for j in range(len(vca)):
                    if vca[j].pos > c1[0]:
                        vca[j].pos += tra_len - sv_len  # could be different sizes...
                        vca[j].info = set_info_end(vca[j].info,
                                                   get_info_end(vca[j].info) + (tra_len - sv_len))
            elif sv_type == 'DEL':  #    sub from pos_list where needed
                for j in range(len(vca)):
                    if vca[j].pos > c1[0]:  # doesn't split boundries
                        vca[j].pos -= sv_len
                        vca[j].info = set_info_end(vca[j].info,
                                                   get_info_end(vca[j].info) - sv_len)
    return vca

#given the full refrence string = full_ref and the
#variant call array = vca produce the mutated string = mut
#mut_type = set(['INS','SUB',{'DEL':0},{'DUP':{'CN':3,'POS'=-100]}
def apply_var_calls(full_ref,vca,g=0,index=0):
    mut,r_pos,offset = '',0,0 #mut is the full mutation sequence, pos is old until the end of the if...
    for i in range(len(vca)): #where to stat = pos, what the ref looks like, what the alternate is
        ref,alt,sv_type = vca[i].ref,vca[i].alt,get_info_type(vca[i].info)
        geno = get_genotype(vca[i].frmat,index)
        if len(geno)>g and geno[g]==1:#allele present
            if sv_type=='INS':                                        #INS
                mut    += full_ref[r_pos:vca[i].pos]+alt
                offset += len(alt)
                r_pos   = vca[i].pos            #no skiping for ins in ref
            elif sv_type=='DEL':
                mut    += full_ref[r_pos:vca[i].pos]+''               #DEL
                offset -= len(ref)
                r_pos   = vca[i].pos+len(ref) #skip the del section in ref
            elif sv_type=='INV' or sv_type=='SUB':
                mut    += full_ref[r_pos:vca[i].pos]+alt                #INV
                r_pos   = vca[i].pos+len(alt)     #move past the INV section
            elif sv_type=='TRA':
                mut    += full_ref[r_pos:vca[i].pos]+alt
                offset += len(alt)-len(ref)
                r_pos   = vca[i].pos+len(ref)
            #-------------------------------------------------------DUP has <DUP> in alt VC field
            elif sv_type=='DUP':#............................................................
                alt  = vca[i].ref
                mut += full_ref[r_pos:vca[i].pos] #grab the section up to the dup pos in the ref space
                dp   = get_dup(vca[i].info)       #get the duplication instructions
                if dp['TYPE']=='DISPERSED' and len(dp['POS'])>=1 and dp['POS'][0]!=vca[i].pos:#DUP->DISPERSED---
                    d_pos = []
                    for d in dp['POS']:
                        d_pos  += [d+offset]
                        offset += len(alt)
                    mut = alt.join(split_pos(mut,dp['POS']))+alt
                    offset += len(alt)      #normal section here
                #DUP->DISPERSED----------------------------------------------------------------------------------
                elif dp['TYPE']=='INV':
                   mut    += ''.join([utils.get_reverse_complement(vca[i].ref) for j in range(dp['CN']/2)])
                   offset += len(vca[i].ref)*(dp['CN']/2-1)
                elif dp['TYPE']=='TANDEM' and len(dp['POS'])>=1:#DUP->TANDEM:::::::::::::::::::::::::::::::::::::
                    mut    += ''.join([vca[i].ref for j in range(dp['CN']/2)]) #the normal ref section
                    offset += len(vca[i].ref)*(dp['CN']/2-1)
                #DUP->TANDEM:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                r_pos = vca[i].pos+len(ref) #move past the DUP section in the ref space
            #---------------------------------------------------
            else: raise KeyError
    mut += full_ref[r_pos:] # left over sections on RHS
    #need to propogate the positions in any unapplied vcas...
    return mut

#given a string S and pos_list pos, split S at each location in pos
def split_pos(S,pos):
    T = []
    if len(pos)>1:
        pos = [0]+pos
        T   = [S[pos[i]:pos[i+1]] for i in range(len(pos)-1)]+[S[pos[-1]:]]
        pos = pos[1:]
    elif len(pos)>=1:
        T   = [S[:pos[0]],S[pos[0]:]]
    return T

#given the string S and query q, return all positions
#in S that match q exactly (no edit distances or approximations)
def get_seq_pos(S,q):
    pos = [r.start() for r in re.finditer(q,S)]
    return pos

def merge_filter_sort_vcam(vcam,aneuploidy,small_cut=50,large_cut=int(260E6)):
    final,dki = [],{}
    for k in vcam:
        if k in aneuploidy:
            print('there is chrom %s aneuploidy'%k)
            for j in range(len(aneuploidy[k])):
                an_svtype = get_info_type(aneuploidy[k][j].info)
                if an_svtype=='DEL':
                    print('loss of chrom %s has been detected'%k)
                    an_geno = [[int(e) for e in g.rsplit('|')] \
                               for g in aneuploidy[k][j].frmat.rsplit('\t')]
                    print('aneuploidy genotype is %s'%an_geno)
                    for i in range(len(vcam[k])):
                        vcam[k][i].frmat = vcam[k][i].frmat.replace('/','|')
                        try:
                            sv_geno = [[int(e) for e in g.rsplit('|')] \
                                       for g in vcam[k][i].frmat.rsplit('\t')]
                            for c in range(len(sv_geno)): #reset all the DEL unsafe genos
                                for p in range(1,len(sv_geno[c])+1):
                                    if an_geno[c][p-1]==1 and sv_geno[c][p-1]==1:
                                        # print('clearing calls for ploidy %s'%p)
                                        sv_geno[c][p-1]=0
                            vcam[k][i].frmat='\t'.join(['|'.join([str(e) for e in c]) for c in sv_geno])
                        except Exception as E:
                            if k in dki: dki[k] += [(j,i,an_geno,sv_geno)]
                            else:        dki[k]  = [(j,i,an_geno,sv_geno)]
                            pass
                        # print('resulting genotype is %s'%sv_geno)
                        if sum([sum(c) for c in sv_geno])>0:
                            sv_len = get_info_end(vcam[k][i].info)-vcam[k][i].pos+1
                            if sv_len >= small_cut and sv_len <= large_cut:
                                final += [vcam[k][i]]
                        # else:
                        #     print('events were lost to aneuploidy: %s : %s : %s : %s'%\
                        #           (vcam[k][i].chrom,vcam[k][i].pos,vcam[k][i].info,vcam[k][i].frmat))
                else: #all vcam[k] elements are DEL safe
                    for i in range(len(vcam[k])):
                        sv_len = get_info_end(vcam[k][i].info)-vcam[k][i].pos+1
                        if sv_len >= small_cut and sv_len <= large_cut:
                            final += [vcam[k][i]]
                final += [aneuploidy[k][j]] #now put in the aneuploidy call

        else:
            final += vcam[k]
    if len(dki)>0:
        print('atempting to dump simulation state...')
        dump_tag = 'ano_issue_' + hashlib.md5(str(len(dki))).hexdigest()[0:5].upper()
        ru.dump_state({'vcam':vcam,'aneuploidy': aneuploidy, 'dki':dki},
                      dump_tag, os.path.expanduser('~/'))
        raise AttributeError
    #:::TODO:::
    return sorted(final,key=lambda x:(x.chrom.zfill(100),x.pos))

#write a VCF file to the outpath, given a vca
#that contains SNVs,SVs, auto builds the header
#using the seqs dict mapping and a data default header template
#usually located in the data directory that needs to have seq key
#injected nefore the first row of data is written to the file
#::TO DO:: this is an easy place to correct the ICN field->/2 and then multiple by the genotype
#::TO DO:: this should be done for the DEL to so you get ;ICN=0 and ICN=1, would be nice!!!
def write_vcf(samples,vca,seqs,ref_path,out_path):
    data_path   = os.path.dirname(os.path.abspath(__file__))+'/data/'
    header_path = data_path+'/header_template.vcf'
    header      = []
    with open(header_path,'r') as f:
        for line in f:
            if line.startswith('#'):
                header += [line.replace('\n','')]
    header[1] += datetime.datetime.now().strftime ("%D")
    header[2] += ref_path.split('/')[-1].rsplit('.fa')[0]
    for k in seqs:
        header += ['##contig=<ID=%s,length=%s>'%(k,len(seqs[k]))]
    header += ['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s'%'\t'.join(samples)]
    s,x = '\n'.join(header)+'\n',1
    for vc in sorted(vca,key=lambda x: (x.chrom.zfill(100),x.pos)):
        if len(vc.ref)>=1:
            s += '\t'.join([vc.chrom,str(vc.pos),'gv_%s'%x,vc.ref[0],vc.alt,'.','PASS',vc.info,'GT',vc.frmat])+'\n'
        x += 1
    with open(out_path,'w') as f:
        f.write(s)
        return True
    return False

#takes in the CT and adds pedegree tags
def write_somatic_vcf(CT,vca,seqs,ref_path,out_path):
    samples = sorted(list(CT.nodes.keys()),key=lambda x: int(x.rsplit('_')[-1]))
    tree = copy.deepcopy(CT.tree)
    data_path   = os.path.dirname(os.path.abspath(__file__))+'/data/'
    header_path = data_path+'/header_template.vcf'
    header      = []
    with open(header_path,'r') as f:
        for line in f:
            if line.startswith('#'):
                header += [line.replace('\n','')]
    header[1] += datetime.datetime.now().strftime ("%D")
    header[2] += ref_path.split('/')[-1].rsplit('.fa')[0]
    for k in seqs:
        header += ['##contig=<ID=%s,length=%s>'%(k,len(seqs[k]))]

    pedigree = ['##PEDIGREE=<Derived=%s,Original=%s>'%(samples[0],samples[0].rsplit('_')[0])] #germline is the base
    for clone in samples:
        for child in tree[clone]:
            pedigree += ['##PEDIGREE=<Derived=%s,Original=%s>'%(child,clone)]
    header += pedigree

    header += ['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s'%'\t'.join(samples)]
    s,x = '\n'.join(header)+'\n',1
    for vc in sorted(vca,key=lambda x: (x.chrom.zfill(100),x.pos)):
        if len(vc.ref)>=1:
            s += '\t'.join([vc.chrom,str(vc.pos),'gv_%s'%x,vc.ref[0],vc.alt,'.','PASS',vc.info,'GT',vc.frmat])+'\n'
        x += 1
    with open(out_path,'w') as f:
        f.write(s)
        return True
    return False

def write_merged_somatic_vcf(normal_vcf,somatic_vcf,out_path):
    t_col,t_head,t_ped,n_col,n_data,n_sample,t_data = [],[],[],[],[],'',[]
    print('reading normal vcf file %s ...'%normal_vcf)
    with open(normal_vcf,'r') as f:
        n_raw = [x.replace('\n','').rsplit('\t') for x in f.readlines()]
        for n in n_raw:
            if n[0].startswith('#CHROM'):     n_col   =  n
            elif not n[0].startswith('#'):    n_data += [n]
        if len(n_col)>0:
            n_sample = n_col[-1] #should be last column
    print('reading somatic vcf file %s ...'%somatic_vcf)
    with open(somatic_vcf,'r') as f:
        t_raw = [x.replace('\n','').rsplit('\t') for x in f.readlines()]
        for t in t_raw:
            if t[0].startswith('#'):          t_head += [t]
            if t[0].startswith('##PEDIGREE'): t_ped  +=  t
            elif t[0].startswith('#CHROM'):   t_col   = t[0]
            elif not t[0].startswith('#'):    t_data += [t]
    print('processing pedigree inforation for genotype inference')
    CT = vcf_pedigree_to_CT(t_ped)
    cs = count_pedigree_clones(t_ped)
    if n_sample in CT: #hash the coordinate type tuple and insert new sorted list, ig will be renumbered gv_x, ...
        print('merging unique VCF entries from normal and somatic...')
        M,i,s = {},1,''
        for h in t_head:
            s += '\t'.join(h)+'\n'
        for n in n_data: #(chrom,start,end,type,geno)
            idx = (n[0],int(n[1]),get_info_end(n[7]),get_info_type(n[7]),tuple(get_genotype(n[9],0)))
            if idx not in M: M[idx] = n[:9]+[n[9] for c in range(cs+1)]
        for t in t_data:
            idx = (t[0],int(t[1]),get_info_end(t[7]),get_info_type(t[7]),tuple(get_genotype(t[9],0)))
            if idx not in M: M[idx] = t
        print('sorting by chrom and start coordinate and enumerating VC id fields')
        for m in sorted(list(M.keys()),key=lambda x: (x[0].zfill(100),x[1])):
            s += '\t'.join(M[m])+'\n'
        print('writing %s results to merged vcf file'%len(M))
        with open(out_path,'w') as f:
            f.write(s)
            return True
        return False
    else:
        print('normal file %s and somatic file %s samples do not match'%(normal_vcf,somatic_vcf))
        return False

def write_split_somatic_vcf(normal_vcf,somatic_vcf,out_path):
    t_col,t_head,t_ped,n_col,n_data,n_sample,t_data = [],[],[],[],[],'',[]
    print('reading normal vcf file %s ...'%normal_vcf)
    with open(normal_vcf,'r') as f:
        n_raw = [x.replace('\n','').rsplit('\t') for x in f.readlines()]
        for n in n_raw:
            if n[0].startswith('#CHROM'):     n_col   =  n
            elif not n[0].startswith('#'):    n_data += [n]
        if len(n_col)>0:
            n_sample = n_col[-1] #should be last column
    print('reading somatic vcf file %s ...'%somatic_vcf)
    with open(somatic_vcf,'r') as f:
        t_raw = [x.replace('\n','').rsplit('\t') for x in f.readlines()]
        for t in t_raw:
            if t[0].startswith('#'):          t_head += [t]
            if t[0].startswith('##PEDIGREE'): t_ped  +=  t
            elif t[0].startswith('#CHROM'):   t_col   = t[0]
            elif not t[0].startswith('#'):    t_data += [t]
    print('processing pedigree inforation for genotype inference')
    CT = vcf_pedigree_to_CT(t_ped)
    cs = count_pedigree_clones(t_ped)
    if n_sample in CT: #hash the coordinate type tuple and insert new sorted list, ig will be renumbered gv_x, ...
        print('removing normal VCF entries from somatic...')
        N,M,i,s = {},{},1,''
        for h in t_head:
            s += '\t'.join(h) + '\n'
        for n in n_data:  # (chrom,start,end,type,geno)
            idx = (n[0],int(n[1]),get_info_end(n[7]),get_info_type(n[7]),tuple(get_genotype(n[9],0)))
            if idx not in N: N[idx] = n
        for t in t_data:
            idx = (t[0],int(t[1]),get_info_end(t[7]),get_info_type(t[7]),tuple(get_genotype(t[9],0)))
            if idx not in N and idx not in M: M[idx] = t
        print('sorting by chrom and start coordinate and enumerating VC id fields')
        for m in sorted(list(M.keys()), key=lambda x: (x[0].zfill(100), x[1])):
            s += '\t'.join(M[m]) + '\n'
        print('writing %s results to merged vcf file' % len(M))
        with open(out_path, 'w') as f:
            f.write(s)
            return True
        return False
    else:
        print('normal file %s and somatic file %s samples do not match'%(normal_vcf,somatic_vcf))
        return False

# list of vcf pedigree strings:
# '##PEDIGREE=<Derived=EDB09649ED_0,Original=EDB09649ED>\n'
def vcf_pedigree_to_CT(ped_ls):
    CT = {}
    for p in ped_ls:
        original = p.rsplit('Original=')[-1].rsplit('>')[0]
        derived  = p.rsplit('Derived=')[-1].rsplit(',')[0]
        if original in CT: CT[original] += [derived]
        else:              CT[original]  = [derived]
    return CT

def count_pedigree_clones(ped_ls):
    cs,x = set([]),0
    for p in ped_ls:
        original = p.rsplit('Original=')[-1].rsplit('>')[0]
        derived = p.rsplit('Derived=')[-1].rsplit(',')[0]
        if original.find('_')>0: cs.add(int(original.rsplit('_')[-1]))
        if derived.find('_') > 0: cs.add(int(derived.rsplit('_')[-1]))
    if len(cs)>0:
        x = sorted(list(cs))[-1]
    return x

#will read the PEDIGREE and FORMAT fields to split out calls
def write_clone_vcf(somatic_vcf,out_path):
    raw,data,header,pedigree,fields,samples,i_k,k_i,outs = [],[],[],[],[],[],{},{},{}
    with open(somatic_vcf,'r') as f:
        for line in f:
            if line.startswith('##PEDIGREE'):
                pedigree += [line.replace('\n','')]
            elif line.startswith('##'):
                header += [line.replace('\n','')]
            elif line.startswith('#'):
                fields += [line.replace('\n','').split('\t')]
                samples += line.replace('\n','').split('\t')[9:]
                for sample in sorted(samples,key=lambda x: int(x.rsplit('_')[-1])):
                    i = int(sample.rsplit('_')[-1])
                    i_k[i] = sample
                    k_i[sample] = i
            else:
                data += [line.replace('\n','').split('\t')]
    #now parse out header for PEDIGREE and last line for FORMAT
    for row in data:
        for i in range(0,len(row[9:])):
            if row[9:][i].find('1')>=0:
                if i_k[i] in outs: outs[i_k[i]] += [row[:9]+[row[9+k_i[i_k[i]]]]]
                else:              outs[i_k[i]]  = [row[:9]+[row[9+k_i[i_k[i]]]]]
    for sample in outs:
        sample_vcf = out_path+'/%s.somatic_S0.vcf'%sample
        with open(sample_vcf,'w') as f:
            s = '\n'.join(header+['\t'.join(fields[0][:9]+[fields[0][9+k_i[sample]]])])+'\n'
            f.writelines(s+'\n'.join(['\t'.join(row) for row in outs[sample]])+'\n')
            print('finished writing %s variation events to %s VCF file'%(len(outs[sample]),sample_vcf))

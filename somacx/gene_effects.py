#!/usr/bin/env python
import os
import glob
import argparse
import variant_utils as vu
import utils as utils

#def samplot_gene_rename(samplot_dir)

def vcf_replace_ids(vcf_in,vcf_out,prefix='vcf'):
    header,data,raw = [],[],[]
    with open(vcf_in,'r') as f:
        raw = [row.replace('\n','').rsplit('\t') for row in f.readlines()]
    for row in raw:
        if row[0].startswith('#'): header += [row]
        else:                      data   += [row]
    for i in range(len(data)):
        data[i][2] = '%s_%s'%(prefix,i+1)
    s = '\n'.join(['\t'.join(h) for h in header+data])+'\n'
    with open(vcf_out,'w') as f:
        f.write(s)
        return True
    return False

def count_validated(wcu,svtypes=['DEL','DUP','INV','TRA']):
    V = {t:{} for t in svtypes}
    for k in wcu:
        for l in wcu[k]:
            for i in range(len(wcu[k][l])):
                t = wcu[k][l][i][3].keys()[0]
                if t not in V:      V[t]     = {k:1}
                elif k not in V[t]: V[t][k]  = 1
                else:               V[t][k] += 1
    return V

def count_totals(vcf_path,svtypes=['DEL','DUP','INV'],label='vcf',filter_size=[50,int(1E9)]):
    W = vu.vcf_to_wcu(vcf_path,w=1.0,remove_dup=True,label=label,filter_size=filter_size)
    V = count_validated(W,svtypes)
    return {t:sum([V[t][k] for k in V[t]]) for t in svtypes}

def partition(wcu):
    W = {}
    for k in wcu:
        for i in range(len(wcu[k])):
            t = wcu[k][i][3].keys()[0]
            if t in W:
                if k in W[t]: W[t][k] += [wcu[k][i]]
                else:         W[t][k]  = [wcu[k][i]]
            else:             W[t]     = {k:[wcu[k][i]]}
    return W

def coordinate_overlap(a,b):
    l = (a[1]-a[0]+1)+(b[1]+b[0]+1)
    u = float(min(l,max(a[1],b[1])-min(a[0],b[0])+1))
    i = 1.0*abs(a[0]-b[0])+abs(a[1]-b[1])
    x = max(0.0,u-i)/u
    return x

def score(wcu1,wcu2,r=0.5):
    S = {t: {'prec':0.0,'rec':0.0,'f1':0.0,'n':0,'m':0,'j':0.0} for t in set(wcu1.keys()).union(set(wcu2.keys()))}
    for t in wcu1:
        if t in wcu2:
            ns = 0
            for k in wcu1[t]:
                if k in wcu2[t]:
                    for i in range(len(wcu1[t][k])):
                        for j in range(len(wcu2[t][k])):
                            y = coordinate_overlap(wcu1[t][k][i],wcu2[t][k][j])
                            if y >= r:
                                ns += 1
                                break
                else: wcu2[t][k] = []
                S[t]['n'] += len(wcu1[t][k])
            if S[t]['n']>0: S[t]['prec'] = 1.0*ns/S[t]['n']
        else:
            wcu2[t] = {k:[]}
            S[t]['n'] += sum([len(wcu1[t][k]) for k in wcu1[t]])
    for t in wcu2:
        if t not in wcu1: wcu1[t] = {k:[]}
        ms = 0
        for k in wcu2[t]:
            if k not in wcu1[t]: wcu1[t][k] = []
            for i in range(len(wcu2[t][k])):
                for j in range(len(wcu1[t][k])):
                    y = coordinate_overlap(wcu2[t][k][i],wcu1[t][k][j])
                    if y >= r:
                        ms += 1
                        break
            S[t]['m'] += len(wcu2[t][k])
        if S[t]['m']>0: S[t]['rec'] = 1.0*ms/S[t]['m']
    for t in wcu1:
        for k in wcu1[t]:
            if k not in wcu2[t]: wcu2[t][k] = []
            F = utils.LRF_1D(wcu1[t][k],wcu2[t][k])
            I,U = sum([abs(f[1]-f[0]) for f in F[0]]),sum([abs(f[1]-f[0]) for f in F[1]])
            if U>0: S[t]['j'] = (1.0*I)/(1.0*U)
            if S[t]['rec']>0.0: S[t]['f1'] = 2.0*(S[t]['prec']*S[t]['rec'])/(S[t]['prec']+S[t]['rec'])
    return S

#given two vcf compute the prec, rec, f1, j scores.
def metrics(vcf1,vcf2,filter_size=[50,int(1E9)],merge=True,r=0.5):
    wcu1 = vu.vcf_to_wcu(vcf1,remove_dup=True,w=1.0,label='vcf1',filter_size=filter_size)      #read and wcu
    wcu1 = partition({k:wcu1[k]['vcf1'] for k in wcu1})                                        #by sv types
    if merge: wcu1 = {t:{k:utils.merge_1D(sorted(wcu1[t][k])) for k in wcu1[t]} for t in wcu1} #merge if needed
    else:     wcu1 = {t:{k:sorted(wcu1[t][k]) for k in wcu1[t]} for t in wcu1}                 #many callers do
    wcu2 = vu.vcf_to_wcu(vcf2,remove_dup=True,w=1.0,label='vcf2',filter_size=filter_size)      # read and wcu
    wcu2 = partition({k:wcu2[k]['vcf2'] for k in wcu2})                                        # by sv types
    if merge: wcu2 = {t:{k:utils.merge_1D(sorted(wcu2[t][k])) for k in wcu2[t]} for t in wcu2} # merge if needed
    else:     wcu2 = {t:{k:sorted(wcu2[t][k]) for k in wcu2[t]} for t in wcu2}                 # many callers do
    S = score(wcu1,wcu2,r=r)
    return S

def g1kp3_nn(vcf1,g1kp3_dir,svtypes=['DEL','DUP','INV'],filter_size=[50,int(1E9)],merge=True,r=0.5,metric='f1'):
    NN,n = {},1.0/(len(svtypes)*1.0)
    g1kp3_vcfs = glob.glob(g1kp3_dir+'/*/*_S0.vcf')
    for g1kp3_vcf in g1kp3_vcfs:
        sample = g1kp3_vcf.rsplit('/')[-1].rsplit('_S0.vcf')[0]
        print('calculating sample= %s'%sample)
        S = metrics(vcf1,g1kp3_vcf,filter_size=filter_size,merge=merge,r=r)
        NN[sample] = sum([n*S[t][metric] for t in S])
    return NN

#given vcf1 and vcf2 return those vcf data rows of vcf1 that have less than
#r=0.5 overlap by svtype with calls that appear in vcf2 and write results to vcf_out
def vcf_difference(vcf1,vcf2,filter_size=[50,int(1E9)],r=0.5,vcf_out=None):
    #load as wcu assuming you have vcf_ids you will be able to pull them out later
    wcu1 = vu.vcf_to_wcu(vcf1,remove_dup=True,w=1.0,label='vcf1',filter_size=filter_size)      #read and wcu
    wcu1 = partition({k:wcu1[k]['vcf1'] for k in wcu1})                                        #by sv types
    wcu1 = {t:{k:sorted(wcu1[t][k]) for k in wcu1[t]} for t in wcu1}                           #many callers do
    wcu2 = vu.vcf_to_wcu(vcf2,remove_dup=True,w=1.0,label='vcf2',filter_size=filter_size)      # read and wcu
    wcu2 = partition({k:wcu2[k]['vcf2'] for k in wcu2})                                        # by sv types
    wcu2 = {t:{k:sorted(wcu2[t][k]) for k in wcu2[t]} for t in wcu2}                           # many callers do
    #find calls in wcu1 that haev < r=0.5 overlap and pull the ids which is more robust
    S = {t:set([]) for t in set(wcu1.keys()).union(set(wcu2.keys()))}
    for t in wcu1:
        if t in wcu2:
            for k in wcu1[t]:
                if k in wcu2[t]:
                    for i in range(len(wcu1[t][k])):
                        for j in range(len(wcu2[t][k])):
                            y = coordinate_overlap(wcu1[t][k][i],wcu2[t][k][j])
                            if y >= r:
                                ids = sorted(list(wcu1[t][k][i][3][t]))
                                for idx in ids: S[t].add(idx)
                                break
    #retrieve data rows from the difference of ids
    T = {t:set([]) for t in set(wcu1.keys()).union(set(wcu2.keys()))}
    header,data,raw = [],[],[]
    with open(vcf1, 'r') as f:
        raw = [row.replace('\n','').rsplit('\t') for row in f.readlines()]
    for row in raw:
        if row[0].startswith('#'): header += [row]
        else:                      data   += [row]
    for t in wcu1:
        for k in wcu1[t]:
            for i in range(len(wcu1[t][k])):
                ids = sorted(list(wcu1[t][k][i][3][t]))
                for idx in ids: T[t].add(idx)
    Q = {t:set([]) for t in set(wcu1.keys()).union(set(wcu2.keys()))}
    D = []
    for t in T:
        Q[t] = T[t].difference(S[t])
        for i in range(len(data)):
            if data[i][2] in Q[t]: D += [data[i]]
    D = sorted(D,key=lambda x: (x[0].zfill(100),int(x[1])))
    for t in Q: print('%s : %s r<0.5 unique calls in vcf1=%s'%(t,len(Q[t]),vcf1))
    if vcf_out is not None:
        s = '\n'.join(['\t'.join(h) for h in header+D])+'\n'
        with open(vcf_out,'w') as f:
            f.write(s)
            return True
        return False

def vcf_intersect(vcf1,vcf2,filter_size=[50,int(1E9)],r=0.5,vcf_out=None):
    #load as wcu assuming you have vcf_ids you will be able to pull them out later
    wcu1 = vu.vcf_to_wcu(vcf1,remove_dup=True,w=1.0,label='vcf1',filter_size=filter_size)      #read and wcu
    wcu1 = partition({k:wcu1[k]['vcf1'] for k in wcu1})                                        #by sv types
    wcu1 = {t:{k:sorted(wcu1[t][k]) for k in wcu1[t]} for t in wcu1}                           #many callers do
    wcu2 = vu.vcf_to_wcu(vcf2,remove_dup=True,w=1.0,label='vcf2',filter_size=filter_size)      # read and wcu
    wcu2 = partition({k:wcu2[k]['vcf2'] for k in wcu2})                                        # by sv types
    wcu2 = {t:{k:sorted(wcu2[t][k]) for k in wcu2[t]} for t in wcu2}                           # many callers do
    #find calls in wcu1 that haev < r=0.5 overlap and pull the ids which is more robust
    S = {t:set([]) for t in set(wcu1.keys()).union(set(wcu2.keys()))}
    for t in wcu1:
        if t in wcu2:
            for k in wcu1[t]:
                if k in wcu2[t]:
                    for i in range(len(wcu1[t][k])):
                        for j in range(len(wcu2[t][k])):
                            y = coordinate_overlap(wcu1[t][k][i],wcu2[t][k][j])
                            if y >= r:
                                ids = sorted(list(wcu1[t][k][i][3][t]))
                                for idx in ids: S[t].add(idx)
                                break
    header,data,raw = [],[],[]
    with open(vcf1, 'r') as f:
        raw = [row.replace('\n','').rsplit('\t') for row in f.readlines()]
    for row in raw:
        if row[0].startswith('#'): header += [row]
        else:                      data   += [row]
    D = []
    for t in S:
        for i in range(len(data)):
            if data[i][2] in S[t]: D += [data[i]]
    D = sorted(D,key=lambda x: (x[0].zfill(100),int(x[1])))
    for t in S: print('%s : %s r>=0.5 similiar calls in vcf1=%s'%(t,len(S[t]),vcf1))
    if vcf_out is not None:
        s = '\n'.join(['\t'.join(h) for h in header+D])+'\n'
        with open(vcf_out,'w') as f:
            f.write(s)
            return True
        return False

def vcf_to_gene_intersect(vcf1,gene_list,gene_map,filter_size=[50,int(1E9)],label='',vcf_out=None):
    wcu1 = vu.vcf_to_wcu(vcf1,remove_dup=True,w=1.0,label='vcf1',filter_size=filter_size)      #read and wcu
    wcu1 = partition({k:wcu1[k]['vcf1'] for k in wcu1})                                        #by sv types
    wcu1 = {t:{k:sorted(wcu1[t][k]) for k in wcu1[t]} for t in wcu1}                           #many callers do
    wcu2 = vu.gene_list_to_wcu(gene_list,[1.0 for i in range(len(gene_list))],gene_map,label=label)
    #find calls in wcu1 that haev < r=0.5 overlap and pull the ids which is more robust
    S = {t:set([]) for t in wcu1}
    for t in wcu1:
        for k in wcu1[t]:
            if k in wcu2:
                for i in range(len(wcu1[t][k])):
                    for j in range(len(wcu2[k])):
                        y = coordinate_overlap(wcu1[t][k][i],wcu2[k][j])
                        if y > 0.0:
                            ids = sorted(list(wcu1[t][k][i][3][t]))
                            for idx in ids: S[t].add(idx)
                            break
    header,data,raw = [],[],[]
    with open(vcf1, 'r') as f:
        raw = [row.replace('\n','').rsplit('\t') for row in f.readlines()]
    for row in raw:
        if row[0].startswith('#'): header += [row]
        else:                      data   += [row]
    D = []
    for t in S:
        for i in range(len(data)):
            if data[i][2] in S[t]: D += [data[i]]
    D = sorted(D,key=lambda x: (x[0].zfill(100),int(x[1])))
    for t in S: print('%s : %s r>=0.5 gene_list hits in vcf1=%s'%(t,len(S[t]),vcf1))
    if vcf_out is not None:
        s = '\n'.join(['\t'.join(h) for h in header+D])+'\n'
        with open(vcf_out,'w') as f:
            f.write(s)
            return True
        return False

des = """
Variant Gene Effects - Compute the effected gene lists for normal versus tumor samples
Timothy James Becker, PhD candidate, UCONN 06/19/2018-05/01/2019\n"""
parser = argparse.ArgumentParser(description=des,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-i', '--in_dir',type=str, help='vcf input folder\t[None]')
parser.add_argument('-o', '--out_dir',type=str, help='analysis output folder\t[None]')
args = parser.parse_args()

in_dir  = args.in_dir
out_dir = args.out_dir
if not os.path.exists(out_dir): os.mkdir(out_dir)
if not os.path.exists(out_dir+'/meta/'): os.mkdir(out_dir+'/meta/')

if not os.path.exists(out_dir+'/meta/gene_map.json'):
    print('didn\'t find a pre-processed gene map, building a new one...')
    gene_map = vu.build_ucsc_gene_exon_map(vu.get_local_path('refGene.txt.gz'))
    vu.write_json_gene_map(out_dir+'/meta/gene_map.json',gene_map)
    print('completed building a python serialized pickle and a json gene map')
gene_map = vu.read_json_gene_map(out_dir+'/meta/gene_map.json')

onco_gl  = vu.read_gene_list(vu.get_local_path('onco_bushman_gene_list.txt'))
mmej_gl  = vu.read_gene_list(vu.get_local_path('mmej_sharma_gene_list.txt'))
nhej_gl  = vu.read_gene_list(vu.get_local_path('nhej_davis_gene_list.txt'))
apot_gl  = vu.read_gene_list(vu.get_local_path('apatosis_thermofisher_gene_list.txt'))
mitcp_gl = vu.read_gene_list(vu.get_local_path('mitcp_giam_gene_list.txt'))
g1kp3_gl = vu.read_gene_list(vu.get_local_path('g1kp3_filtered_gene_list.txt'))

onco_wcu  = vu.gene_list_to_wcu(onco_gl,[1.0 for g in onco_gl],gene_map,'onco')
mmej_wcu  = vu.gene_list_to_wcu(mmej_gl,[1.0 for g in mmej_gl],gene_map,'mmej')
nhej_wcu  = vu.gene_list_to_wcu(nhej_gl,[1.0 for g in nhej_gl],gene_map,'nhej')
apot_wcu  = vu.gene_list_to_wcu(apot_gl,[1.0 for g in apot_gl],gene_map,'apot')
mitcp_wcu = vu.gene_list_to_wcu(mitcp_gl,[1.0 for g in mitcp_gl],gene_map,'mitcp')
g1kp3_wcu = vu.gene_list_to_wcu(g1kp3_gl,[1.0 for g in g1kp3_gl],gene_map,'g1kp3')


GS = {}
vcf_glob = glob.glob('/media/data/tcrboa_validated/SMX*S_S-1.vcf')
# for vcf_path in sorted(vcf_glob,key=lambda x: int(x.rsplit('_F')[-1].rsplit('n')[0])):
for vcf_path in vcf_glob:
    print('starting analysis of %s'%vcf_path)
    vcf_wcu  = vu.vcf_to_wcu(vcf_path,w=1.0,label=None)
    gain,loss,total = vu.g1kp3_to_gene_list(vcf_path,gene_map,filter=[0.0,1.0,True])
    G,L,T = {'onco':[],'mmej':[],'nhej':[],'apot':[],'mitcp':[],'g1kp3':[]},\
            {'onco':[],'mmej':[],'nhej':[],'apot':[],'mitcp':[],'g1kp3':[]},\
            {'onco':[],'mmej':[],'nhej':[],'apot':[],'mitcp':[],'g1kp3':[]}
    A = {'onco':onco_gl,'mmej':mmej_gl,'nhej':nhej_gl,'apot':apot_gl,'mitcp':[],'g1kp3':g1kp3_gl}
    for a in A:
        for i in gain:
            if i in A[a]: G[a] += [i]
        for i in loss:
            if i in A[a]: L[a] += [i]
        for i in total:
            if i in A[a]: T[a] += [i]
    print('TOTAL : DEL,DUP,INV genes-----------------------------------------------------------')
    for t in T:
        prop = 0.0
        if len(A[t])>0: prop = 1.0*len(T[t])/len(A[t])
        print('%s/%s genes in %s (%s)'%(len(T[t]),len(A[t]),t,prop))
    prop = 0.0
    if len(T['g1kp3'])>0: prop = 1.0*len(T['onco'])/len(T['g1kp3'])
    print('onco prop to g1kp3 prop = %s'%prop)
    print('FROM %s VCF calls-------------------------------------------------------------------'%\
          sum([len(vcf_wcu[k]) for k in vcf_wcu]))
    GS[vcf_path.rsplit('/')[-1].rsplit('.vcf')[0]] = T['onco']





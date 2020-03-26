#Timothy James Becker, PhD candidate, UCONN 01/10/2017-03/20/2020
import copy
import json
import hashlib
import gc
import glob
import time
import socket
import sys
import os

import numpy as np

import somacx.read_utils as ru
import utils
import somacx.variant_utils as vu

class Clone:
    def __init__(self,node_id,model=0.0,branch=2,decay=0.1,depth=None):
        self.node_id = node_id
        self.model   = model
        self.branch  = branch
        self.decay   = decay
        self.depth   = depth

    def __enter__(self):
        return self

    def __del__(self):
        return 0

    def __exit__(self, type, value, traceback):
        return 0

class CloneTree:
    def __init__(self,ref_path=None,out_dir=None,
                 root_id=None,model=0.0,branch=2,decay=0.1,n=10,
                 sim_cycle=None,json_path=None,gz=True):
        if gz: ext_pat,idx_pat = '.fa.gz','i'
        else:  ext_pat,idx_pat = '.fa','.fai'
        if json_path is None:
            #get the germline genome root_id in the out_dir
            self.ref_path    = ref_path
            self.gz          = gz
            self.out_dir     = out_dir
            self.genome_list = {i.rsplit('/')[-1].replace(ext_pat,''):i for i in glob.glob(out_dir+'/*'+ext_pat)}
            self.vcf_list    = {i.rsplit('/')[-1].replace('_S0.vcf',''):i for i in glob.glob(out_dir+'/*_S0.vcf')}
            self.root_id     = root_id+'_0' #now they can be sorted by casting the '_' split
            if root_id not in self.genome_list or root_id not in self.vcf_list:
                print('root_id was not found')
                raise AttributeError
            #load the chrom names and lengths from the ref_path/germline
            rs = ru.get_fasta_seq_names_lens(ref_path)                       #all possible
            ks = ru.get_fasta_seq_names_lens(self.genome_list[root_id]) #from generated
            self.gs = {}
            for k in ks:
                q = ''.join(k.rsplit('_')[:-1])
                if q in rs: self.gs[q] = rs[q]
            self.cycle,self.m,self.n,self.d = 0,0,n,1
            self.new,self.past  = [],[] #nodes that are born and die and are added here
            self.tree  = {self.root_id:[]}   #tree structure
            self.alive = {self.root_id:True} #life status
            self.nodes = {self.root_id:Clone(root_id,model,branch,decay,depth=0)} #data access
            self.freq  = {self.root_id:1}
            self.depth = {self.root_id:0}
            if sim_cycle is not None:
                self.simulate(sim_cycle)
        else:
            print('using clone_tree.json constructor form')
            #load and unpack the json object once then propigate
            clone_tree = self.read_json_clone_tree(json_path)
            print(clone_tree)
            self.json_to_clone_tree(clone_tree)
            self.ref_path = ref_path
            self.out_dir  = out_dir
            # [A] ref_path must be provided by command line
            # [B] out_dir must be provided by command line
            print('completed json CloneTree loading')

    def nodes_to_json(self):
        N = {}
        for node in self.nodes:
            N[node] = {'node_id': self.nodes[node].node_id,
                       'model':   self.nodes[node].model,
                       'branch':  self.nodes[node].branch,
                       'decay':   self.nodes[node].decay,
                       'depth':   self.nodes[node].depth}
        return N

    def clone_tree_to_json(self):
        clone_tree = {'ref_path':    self.ref_path,
                      'out_dir':     self.out_dir,
                      'genome_list': self.genome_list,
                      'vcf_list':    self.vcf_list,
                      'root_id':     self.root_id,
                      'gs':          self.gs,
                      'cycle':       self.cycle,
                      'm':           self.m,
                      'n':           self.n,
                      'd':           self.d,
                      'new':         self.new,
                      'past':        self.past,
                      'tree':        self.tree,
                      'alive':       self.alive,
                      'nodes':       self.nodes_to_json(),
                      'freq':        self.freq,
                      'depth':       self.depth,
                      'gz':          self.gz}
        return clone_tree

    #needed for conversion JSON to Python class
    def json_to_nodes(self,N):
        nodes = {}
        for node in N:
            nodes[node] = Clone(node_id=N[node]['node_id'],
                                model=N[node]['model'],
                                branch=N[node]['branch'],
                                decay=N[node]['decay'],
                                depth=N[node]['depth'])
        return nodes

    #used to convert a clone_tree json file into an actual CloneTree instance
    def json_to_clone_tree(self,clone_tree):
        self.ref_path     = clone_tree['ref_path']    #for forward testing
        self.out_dir      = clone_tree['out_dir']     #for forward testing
        self.genome_list  = clone_tree['genome_list']
        self.vcf_list     = clone_tree['vcf_list']
        self.root_id      = clone_tree['root_id']
        self.gs           = clone_tree['gs']
        self.cycle        = clone_tree['cycle']
        self.m            = clone_tree['m']
        self.n            = clone_tree['n']
        self.d            = clone_tree['d']
        self.new          = clone_tree['new']
        self.past         = clone_tree['past']
        self.tree         = clone_tree['tree']
        self.alive        = clone_tree['alive']
        self.nodes        = self.json_to_nodes(clone_tree['nodes'])
        self.freq         = clone_tree['freq']
        self.depth        = clone_tree['depth']
        self.gz           = clone_tree['gz']

    def write_json_clone_tree(self,json_path):
        with open(json_path,'w') as f:
            f.write(json.dumps(self.clone_tree_to_json()))
            return True
        return False

    def read_json_clone_tree(self,json_path):
        clone_tree = {}
        with open(json_path,'r') as f:
            clone_tree = json.load(f,object_pairs_hook=vu.str_hook)
        return clone_tree

    def new_node(self,parent_id):
        self.m += 1  # update the new clone enumaration
        node_id = '%s_%s' % (self.root_id.split('_')[0], self.m)  # get the str id
        self.tree[parent_id] += [node_id]
        self.tree[node_id]    = []
        self.alive[node_id]   = True
        model  = self.nodes[parent_id].model
        branch = self.nodes[parent_id].branch * self.nodes[parent_id].model
        decay  = self.nodes[parent_id].decay
        depth  = self.nodes[parent_id].depth+1
        decay  = min(1.0,max(1E-9,decay*((2.0*model)-decay+1E-9)))
        self.nodes[node_id]   = Clone(node_id,model,branch,decay,depth)
        self.depth[node_id] = depth

    def del_node(self,node_id):
        alive = self.living()
        if self.living()>(self.nodes[self.root_id].decay)*self.n:
            self.alive[node_id] = False
        return True

    def dft(self,node_id):
        if node_id != [] and node_id in self.nodes:
            if self.alive[node_id]: #use the clones params to generate possible subclones
                if np.random.choice([True,False],p=[min(1.0,max(0.0,self.nodes[node_id].branch)),1.0-min(1.0,max(0.0,self.nodes[node_id].branch))]):
                    self.new += [node_id]
            if np.random.choice([True,False],p=[min(1.0,max(0.0,self.nodes[node_id].decay)),1.0-min(1.0,max(0.0,self.nodes[node_id].decay))]):
                self.past += [node_id]
            for k in self.tree[node_id]: self.dft(k)

    def generate_cycle(self):
        self.dft(self.root_id)
        for node_id in self.new:
            if self.living()<self.n:
                self.new_node(node_id)
        for node_id in self.past:
            self.del_node(node_id)
        self.past, self.new = [], []
        self.cycle += 1

    def simulate(self,cycles=1E3):
        i = 0
        while self.living()<self.n and i<cycles:
            self.generate_cycle()
            i += 1
        F = {}
        self.allele_freq(self.root_id,F)
        self.freq = F

    #for each living node id compute the frequency of its alleles
    def allele_freq(self,node_id,F={}):
        F[node_id] = 1
        if self.tree[node_id]!=[]:
            for k in self.tree[node_id]:
                self.allele_freq(k,F)
                F[node_id] += F[k]

    def max_depth(self):
        return max([self.nodes[i].depth for i in self.nodes])

    def living(self):
        return sum([1 if self.alive[k] else 0 for k in self.alive])

    #returns the list of decendant node_ids
    def decendants(self,node_id):
        D = []
        if self.tree[node_id]!=[]:
            D = []
            for k in self.tree[node_id]: D += [k] + self.decendants(k)
        return D

    def ancestors(self,node_id):
        A = []
        for k in self.tree:
            if node_id in self.tree[k]:
                A += [k] + self.ancestors(k)
        return A

    #:::TO DO::: add aditional distributions :::TO DO:::
    def mut_prop(self,mode='lin-desc'):
        if mode=='lin-desc':
            P = [1.0/(x+1) for x in range(len(self.nodes))][::-1]
        return P

    def __enter__(self):
        return self

    def __del__(self):
        return 0

    def __exit__(self, type, value, traceback):
        return 0

# spooling ensures you are only loading one seq into memory at a time
# in combination with single seq loading and fecthing for TRA
# this should lower memory usaged from the germline workflow
# by ~ 6GB for human genomes due to ref=3GB and mut2-3GB
# som_fa_map = {'NA12892_0':'/data/NA12892_0'}
# clone 0 chr1 ploidy 1 fasta will be: /data/NA12892_0.1.chr1.fa
# clone 0 chr1 ploidy 2 fasta will be: /data/NA12892_0.2.chr1.fa
# clone 1 chr1 ploidy 1 fasta will be: /data/NA12892_1.1.chr1.fa
# clone 1 chr1 ploidy 2 fasta will be: /data/NA12892_1.2.chr1.fa
#:::TO DO::: May have to handle the part_map being starved due to lack of vc...
def spool_fasta_clones(som_fa_map,chrom,ploidy,vca,
                       full_ref=None,complex=None,anueuploidy=None,
                       faidx=True,gz=True):
    if gz: ext_pat,idx_pat = '.fa.gz','i'
    else:  ext_pat,idx_pat = '.fa','.fai'
    clones = sorted(som_fa_map,key=lambda x: int(x.split('_')[-1]))
    if complex is None:
        for i in range(len(clones)): #clones are enumerated and sorted started from root=0
            for p in range(1,ploidy+1): #use a dot notation for seq duplicates
                fasta_spool = som_fa_map[clones[i]]+'.%s.%s'%(p,chrom)+ext_pat
                clone_id   = clones[i].rsplit('_')[-1]
                if os.path.exists(fasta_spool):
                    try:
                        mut = ru.read_fasta_chrom(fasta_spool,chrom+'.%s.%s'%(p,clone_id))
                        mut = vu.apply_var_calls(mut,vca,g=p-1,index=i)
                        utils.write_fasta({chrom+'.%s.%s'%(p,clone_id):mut},fasta_spool,index=faidx,gz=gz)
                    except Exception as E:
                        print(E.message)
                        print('fasta_spool=%s chrom=%s p=%s i=%s clone=%s,complex=%s'%(fasta_spool,chrom,p,i,clone_id,complex))
                        raise AttributeError
                else: #have to assume a germline inputs to split out
                    mut = copy.deepcopy(full_ref)
                    mut = vu.apply_var_calls(mut,vca,g=p-1,index=i)
                    utils.write_fasta({chrom+'.%s.%s'%(p,clone_id):mut},fasta_spool,index=faidx,gz=gz)
        print('applied %s shared SNV|MNV to seq %s'%(len(vca),chrom))
    else:
        for i in range(len(clones)): #clones are enumerated and sorted started from root=0
            for p in range(1,ploidy+1): #use a dot notation for seq duplicates
                fasta_spool = som_fa_map[clones[i]]+'.%s.%s'%(p,chrom)+ext_pat
                clone_id = clones[i].rsplit('_')[-1]
                if os.path.exists(fasta_spool):
                    try:
                        mut = ru.read_fasta_chrom(fasta_spool,chrom+'.%s.%s'%(p,clone_id))
                        mut = vu.apply_var_calls(mut,vu.update_vca_pos(vca,complex,g=p-1,index=i),g=p-1,index=i)
                        utils.write_fasta({chrom+'.%s.%s'%(p,clone_id):mut},fasta_spool,index=faidx,gz=gz)
                    except Exception as E:
                        print(E.message)
                        print('fasta_spool=%s chrom=%s p=%s i=%s clone=%s,complex=%s'%(fasta_spool,chrom,p,i,clone_id,complex))
                        raise AttributeError
                else: #have to assume a germline inputs to split out
                    mut = copy.deepcopy(full_ref)
                    #update the vca_pos---------------------------
                    mut = vu.apply_var_calls(mut,vu.update_vca_pos(vca,complex,g=p-1,index=i),g=p-1,index=i)
                    utils.write_fasta({chrom+'.%s.%s'%(p,clone_id):mut},fasta_spool,index=faidx,gz=gz)
        print('applied %s shared SNV|MNV to seq %s'%(len(vca),chrom))
    if anueuploidy is not None: #assume the other chrom have been written
        for i in range(len(clones)):
            clone_id = int(clones[i].rsplit('_')[-1]) #clone id and ploidy 1 to ploidy+1
            print('clone_id is %s'%clone_id)
            if chrom in anueuploidy:
                for vca in anueuploidy[chrom]:
                    geno = vu.get_genotype(vca.frmat,clone_id)
                    print('geno is %s'%(geno,))
                    for j in range(len(geno)):
                        if geno[j]==1:
                            print('anueuploidy vca being applied')
                            sv_len = vu.get_info_len(vca.info)
                            sv_type = vu.get_info_type(vca.info)
                            fasta_spool = som_fa_map[clones[i]]+'.%s.%s'%(j+1,chrom)+ext_pat
                            print('clone=%s,ploidy=%s,chrom=%s,len=%s,type=%s'%\
                                  (clone_id,len(geno),chrom,sv_len,sv_type))
                            print(fasta_spool)
                            if sv_type=='DEL':    #delete the chrom
                                for fa in glob.glob(fasta_spool+'*'): os.remove(fa)
                            elif sv_type=='DUP': #add another chrom ploidy=len(geno)+1
                                sv_dup = vu.get_dup(vca.info)['CN']
                                spool_pat = fasta_spool.rsplit('/')[-1].rsplit('_')[0]+\
                                            '_%s.*.%s'%(clone_id,chrom)+ext_pat
                                high = max([int(x.rsplit('/')[-1].rsplit('.')[1]) \
                                            for x in glob.glob('/'.join(fasta_spool.rsplit('/')[0:-1])+spool_pat)])
                                mut = ru.read_fasta_chrom(fasta_spool,chrom+'.%s.%s'%(j+1,clone_id))
                                for d in range(sv_dup):
                                    new_spool = '_'.join(fasta_spool.rsplit('_')[:-1])+\
                                                '_%s.%s.%s'%(clone_id,high+d,chrom)+ext_pat
                                    utils.write_fasta({chrom+'.%s.%s'%(high+d,clone_id):mut},new_spool,index=faidx,gz=gz)
    return True

#merge the input fasta files (optionally apply aneuploidy probabilties to each clone)
#living will not add those sequences from the dead clonal tissues and
#clean parameter will delete all the input files to save disk space as they are
#loading into memory one at a time and appended into the destination file
#this can help with large somatic clone nodes and provides provisions for managing disk
def merge_fasta_clones(som_fa_map,fa_out_path,aneuploidy=None,living=None,clean=False,gz=True):
    fastas,chroms = [],{}
    if gz: ext_pat,idx_pat = '.fa.gz','i'
    else:  ext_pat,idx_pat = '.fa','.fai'
    for clone in som_fa_map:
        fastas += glob.glob(som_fa_map[clone]+'*.*.*'+ext_pat) #get just chrom.ploidy.clone.fa...
    fastas = list(set(fastas))
    for i in range(len(fastas)):
        seqs = ru.get_fasta_seq_names_lens(fastas[i])
        for k in seqs:
            if k in chroms: chroms[k] += [fastas[i]]
            else:           chroms[k]  = [fastas[i]]
    # print(chroms)
    for i in range(len(fastas)):
        clone = fastas[i].rsplit('/')[-1].rsplit('.')[0]
        clone_num = clone.rsplit('_')[-1]
        seqs = ru.get_fasta_seq_names_lens(fastas[i])
        for k in seqs:
            if k==list(seqs.keys())[-1] and i==(len(fastas)-1):
                print('\nlast seq, now indexing...\n')
                utils.write_fasta({k:ru.read_fasta_chrom(fastas[i],k)},fa_out_path,mode='a',index=True,gz=gz)
            else:
                utils.write_fasta({k:ru.read_fasta_chrom(fastas[i],k)},fa_out_path,mode='a',index=False,gz=gz)
        if clean:
            for fa in glob.glob(fastas[i]+'*'): os.remove(fa)
    return True

#will merge all the individual seqs for each sample into single-cell type fasta files
#downstream sequencing simulation than can be done with premixed or with single-cell
def merge_each_fasta_clone(som_fa_map,out_dir,aneuploidy=None,living=None,clean=False,gz=True):
    fastas = []
    if gz: ext_pat,idx_pat = '.fa.gz','i'
    else:  ext_pat,idx_pat = '.fa','.fai'
    for clone in som_fa_map:
        fastas += glob.glob(som_fa_map[clone]+'*.*.*'+ext_pat)
    fastas = list(set(fastas))
    for i in range(len(fastas)):
        clone = fastas[i].rsplit('/')[-1].rsplit('.')[0]
        clone_num = clone.rsplit('_')[-1]
        seqs = ru.get_fasta_seq_names_lens(fastas[i])
        for k in seqs:
            if k==list(seqs.keys())[-1] and i==(len(fastas)-1): #index after the last one has been added
                utils.write_fasta({k:ru.read_fasta_chrom(fastas[i],k)},out_dir+'/%s'%clone+ext_pat,mode='a',index=True,gz=gz)
            else:
                utils.write_fasta({k:ru.read_fasta_chrom(fastas[i],k)},out_dir+'/%s'%clone+ext_pat,mode='a',index=False,gz=gz)
        if clean:
            for fa in glob.glob(fastas[i]+'*'): os.remove(fa)
    return True

#for each type display, rate, loss and gain
def display_rlg(ks,rate,loss,gain):
    for k in ks:
        if k in rate:
            for t in rate[k]:
                print('rate contig=%s, type=%s:'%(k,t))
                for s in sorted(rate[k][t]):
                    print('%s:%s'%(s,rate[k][t][s]))
        #this is k and then c for class
        if k in loss:
            for c in loss[k]:
                print('loss contig=%s, class=%s:'%(k,c))
                print(loss[k][c])
        if k in gain:
            for c in gain[k]:
                print('gain contig=%s, class=%s:'%(k,c))
                print(gain[k][c])

def wcu_enrichment(M, loss_wcu, gain_wcu, gene_map):
    gain, loss = {}, {}
    for k in M:
        loss[k], gain[k] = {}, {}
        for c in loss_wcu[k]:
            loss[k][c] = vu.prop_class_mut_pos({k: M[k]}, {k: loss_wcu[k][c]}, gene_map)
        for c in gain_wcu[k]:
            gain[k][c] = vu.prop_class_mut_pos({k: M[k]}, {k: gain_wcu[k][c]}, gene_map)
    return loss, gain

def write_genome_from_vcam(ref_path,vcam,sample,out_dir,rs,gene_map,
                           write_snv_indel=True,small_cut=50,gz=True):
    if gz: ext_pat,idx_pat = '.fa.gz','i'
    else:  ext_pat,idx_pat = '.fa','.fai'
    fa_out_path    = out_dir+'/'+sample+ext_pat
    vcf_out_path   = out_dir+'/'+sample+'_S0.vcf'
    loss,gain,rate = {},{},{} #loads these from vcam => wcu frameworks
    loss_wcu,gain_wcu = {},{} #loads these from analysis of VCF file
    ks = set([])
    for l in vcam:
        for k in vcam[l]:
            if k not in ks: ks.add(k)
    ks = sorted(list(ks))
    rs = sorted(list(ks))
    f_start = time.time()
    for k in rs: #vcams already have viable alt and ref strings
        if k in ks:
            loss_wcu[k],gain_wcu[k] = [],[] #do analysis---------------------------------------
            seqs = {k:ru.read_fasta_chrom(ref_path,k)}
            ploidy = 2
            if 'Y' in ks and 'X' in ks and (k == 'Y' or k == 'X'): ploidy = 1
            print('processing seq=%s with ploidy=%s----------------' % (k, ploidy))
            mut1 = copy.deepcopy(seqs[k])
            if ploidy > 1:
                mut2 = copy.deepcopy(seqs[k])
            #-----------------------------------------------------------------------------------
            l = 1
            if l in vcam and k in vcam[l]: #variation layer 1:linked MNVs
                print('L1:linked MNV editing mut1, mut2 strings for chrom %s' % k)
                mut1 = vu.apply_var_calls(mut1,vcam[l][k],g=0)
                if ploidy > 1:
                    if k == 'Y'  and len(vcam[l][k])>0 and vcam[l][k][0].chrom == 'Y' or \
                       k == 'MT' and len(vcam[l][k])>0 and vcam[l][k][0].chrom == 'MT':
                        print('key=%s, chrom=%s' % (k,vcam[l][k][0].chrom))
                    mut2 = vu.apply_var_calls(mut2,vcam[l][k],g=1)
                print('L1:applied %s linked SNV, MNV to seq %s' % (len(vcam[l][k]),k))
            l = 2
            if l in vcam and k in vcam[l]: #variation layer 2: SVs,cxSVs
                print('calculating expected l1:l2 differences')
                x,F = [0, 0],{'DEL':[],'DUP':[],'INS':[],'INV':[],'TRA':[]}
                for vc in vcam[l][k]:
                    svtype = vu.get_info_type(vc.info)
                    svlen  = vu.get_info_len(vc.info)
                    svgen  = vu.get_genotype(vc.frmat,0)
                    if svtype=='DEL':
                        if svgen[0]==1:                  x[0] -= svlen-1
                        if len(svgen)>1 and svgen[1]==1: x[1] -= svlen-1
                    elif svtype=='DUP':
                        if svgen[0]==1:                  x[0] += (vu.get_dup(vc.info)['CN']/4)*svlen-1
                        if len(svgen)>1 and svgen[1]==1: x[1] += (vu.get_dup(vc.info)['CN']/4)*svlen-1
                    elif svtype=='INS':
                        if svgen[0]==1:                  x[0] += len(vc.alt)
                        if len(svgen)>1 and svgen[1]==1: x[1] += len(vc.alt)
                    F[svtype] += [svlen]
                for f in F: print('L2:applied %s %s of mean-len=%s for SVs to seq %s'%\
                                  (len(F[f]),f, (0 if len(F[f])<1 else np.mean(F[f])), k))
                y = [0,0]
                print('L2:SV editing mut1, mut2 strings for chrom %s'%k)
                sv1 = vu.apply_var_calls(mut1,vcam[l][k],g=0)
                y[0] += len(sv1)-len(mut1)
                if ploidy > 1:
                    sv2 = vu.apply_var_calls(mut2,vcam[l][k],g=1)
                    y[1] += len(sv2)-len(mut2)
                if abs(x[0]-y[0])>0 or abs(x[1]-y[1])>0:
                    print('@@VCF@@-------------representation violation = %s : %s --------------------@@VCF@@'%(x,y))
                    ru.dump_state(vcam[l][k],'vcf_%s_violation_error'%k,out_dir+'/')
                else:
                    print('actual l1:l2 differences match expected')
                print('L2: applied %s SVs to seq %s' % (len(vcam[l][k]), k))
            l = 3
            if l in vcam and k in vcam[l]:  # variation layer 1:linked MNVs
                print('L3:unlinked MNV editing mut1, mut2 strings for chrom %s' % k)
                mut1 = vu.apply_var_calls(sv1,vcam[l][k],g=0)
                if ploidy > 1:
                    if k == 'Y'  and len(vcam[l][k])>0 and vcam[l][k][0].chrom == 'Y' or \
                       k == 'MT' and len(vcam[l][k])>0 and vcam[l][k][0].chrom == 'MT':
                        print('key=%s, chrom=%s' % (k,vcam[l][k][0].chrom))
                    mut2 = vu.apply_var_calls(sv2,vcam[l][k],g=1)
                print('L3: applied %s unlinked SNV, MNV to seq %s' % (len(vcam[l][k]), k))
            #------------------------------------------------------------------------------------
            print('writing the results to a fasta')
            if ploidy > 1:
                utils.write_fasta({k+'.1':mut1,k+'.2':mut2},fa_out_path,mode='a',index=True,gz=gz)
            else:
                utils.write_fasta({k+'.1':mut1},fa_out_path,mode='a',index=True,gz=gz)
            print('%s--------------------------------\n'%(''.join(['-' for i in range(len(k))])))
    seqs = {k:ru.read_fasta_chrom(ref_path,k) for k in rs} #now load up all seqs for vcf
    # write a g1k style VCF 4.2 file for FusorSV S0 inputs for big DEL,DUP,INV call
    if not write_snv_indel:
        print('saving the large SVs to VCF')
        final_sv_vca = vu.merge_filter_sort_vcam(vcam[2], {}, small_cut=small_cut)
        vu.write_vcf([sample], final_sv_vca,seqs,ref_path,vcf_out_path)
        vu.wcu_to_bedgraph(loss_wcu, out_dir + '/%s_loss.bedgraph' % sample, sign=-1.0)
        vu.wcu_to_bedgraph(gain_wcu, out_dir + '/%s_gain.bedgraph' % sample, sign=1.0)
    else:
        print('saving SNV/MNV,INDEL and large SVs to VCF')
        snv_indel_vcam = vu.vcam_union(vcam[1], vcam[3])
        final_vcam = vu.vcam_union(vcam[2], snv_indel_vcam)
        final_vca = vu.merge_filter_sort_vcam(final_vcam, {}, small_cut=0)
        vu.write_vcf([sample], final_vca, seqs, ref_path, vcf_out_path)
        vu.wcu_to_bedgraph(loss_wcu, out_dir + '/%s_loss.bedgraph' % sample, sign=-1.0)
        vu.wcu_to_bedgraph(gain_wcu, out_dir + '/%s_gain.bedgraph' % sample, sign=1.0)
    f_stop = time.time()
    print('genome generated in %s sec' % round(f_stop - f_start, 2))
    return sample,vcam,loss,gain,rate

# rs is the chrom regions that will be written to VCF/fasta
# ks is the chroms that will get SVs
# seqs is the reference dict
def germline_genome(ref_path,out_dir,rs,ks,mut_p,loss_wcu,gain_wcu,gene_map,
                    gen_method='fast',gz=True,write_snv_indel=True,small_cut=50):
    if gz: ext_pat,idx_pat = '.fa.gz','i'
    else:  ext_pat,idx_pat = '.fa','.fai'
    if sys.version_info.major < 3:
        sample = hashlib.md5(socket.gethostname()+str(time.time())).hexdigest()[:10].upper()
    else:
        sample = hashlib.md5(str(socket.gethostname()+str(time.time())).encode('utf_8')).hexdigest()[:10].upper()
    fa_out_path        = out_dir+'/'+sample+ext_pat
    vcf_out_path       = out_dir+'/'+sample+'_S0.vcf'
    print('loading all ref seq into memory for TRA generation')
    seqs = {k: ru.read_fasta_chrom(ref_path, k) for k in rs}
    vcam = {k:{} for k in mut_p} #get the layers used for generation: 1,2,3
    # --------------------------------------------------------------------------------------------
    M,G,T,loss,gain,rate = {},{},{},{},{},{}  # generate non-overlapping regions for SVs in layer 2 of mut_p
    for k in rs:
        if k in ks:  # make this mut_l and mut_n a distribution [1E3,1E4,1E5,1E6],[30,20,10,5],ect...
            print('building initial SV distributions for seq %s' % k)
            l_pos,l_w = vu.wcu_to_pos_w(loss_wcu[k],len(seqs[k])) #option to set background and scaling
            g_pos,g_w = vu.wcu_to_pos_w(gain_wcu[k],len(seqs[k])) #option to set background and scaling
            class_p = {'DEL':[l_pos,l_w],'INV':[l_pos,l_w],
                       'INS':[l_pos,l_w],'TRA':[l_pos,l_w],'DUP':[g_pos,g_w]}
            if gen_method=='fast':
                M[k] = vu.gen_class_mut_pos_map(seqs[k],copy.deepcopy(class_p),mut_p,l=2)
            else:
                M[k] = vu.gen_class_mut_pos_map_slow(seqs[k],copy.deepcopy(class_p),mut_p,l=2)
            loss[k],gain[k] = {},{}
            for c in loss_wcu[k]:
                loss[k][c] = vu.prop_class_mut_pos({k:M[k]},{k:loss_wcu[k][c]},gene_map)
            for c in gain_wcu[k]:
                gain[k][c] = vu.prop_class_mut_pos({k:M[k]},{k:gain_wcu[k][c]},gene_map)
            rate[k] = vu.mut_p_rate(mut_p,M[k],len(seqs[k]),l=2)
    # set genotypes on the corrected events in M------------------------------
    for k in M:
        G[k] = {}
        for t in M[k]:
            ploidy = 2
            if 'Y' in ks and 'X' in ks and (k == 'Y' or k == 'X'): ploidy = 1
            G[k][t] = vu.gen_genotypes(ploidy, M[k][t], hetro=mut_p[2][t]['h'])
    print('generating TRA mapping from available seqs\n')
    TK = {}
    for k in M:
        if 'TRA' in M[k]: TK[k] = copy.deepcopy(M[k]['TRA'])
    T = vu.gen_tra_map(TK)

    f_start = time.time()
    for k in rs:
        vcam[1][k],vcam[2][k],vcam[3][k],vca = [],[],[],[]

        ploidy = 2
        if 'Y' in ks and 'X' in ks and (k == 'Y' or k == 'X'):ploidy = 1
        if 'chrY' in ks and 'chrX' in ks and (k == 'chrY' or k == 'chrX'): ploidy = 1
        print('processing seq=%s with ploidy=%s----------------' % (k, ploidy))

        # start with a some SNV before SV ----------------------------------------------
        pos,w = vu.wcu_to_pos_w(loss_wcu[k],len(seqs[k]))
        class_p = {'SUB':[pos,w]}
        mut_pos = vu.gen_class_mut_pos_map(seqs[k],class_p,mut_p,l=1)
        if 'SUB' in mut_pos:
            mut_pos = mut_pos['SUB']
            mut_gen = vu.gen_genotypes(ploidy, mut_pos, hetro=mut_p[1]['SUB']['h'])
            vca = vu.gen_var_calls(seqs, k, 'SUB', mut_pos, mut_gen)  # gen_var_calls should have a heterogenity prob
            print('L1: generated %s SNV, MNVs for seq %s' % (len(mut_pos), k))
            mut1 = vu.apply_var_calls(seqs[k], vca, g=0)
            if ploidy > 1:
                mut2 = vu.apply_var_calls(seqs[k], vca, g=1)
            print('L1: applied %s SNV, MNV to seq %s' % (len(vca), k))
            vcam[1][k] += vca
            vcam[1][k] = sorted(vcam[1][k],key=lambda x:(x.chrom.zfill(100),x.pos))
        else:
            mut1 = copy.deepcopy(seqs[k])
            if ploidy > 1:
                mut2 = copy.deepcopy(seqs[k])
        # now set total SVs and partition into DEL,DUP,INV,TRA
        ins_vca,del_vca,dup_vca,tra_vca,inv_vca = [],[],[],[],[]
        inside_del_vca,inside_ins_vca,comp_del_vca,comp_ins_vca,flanking_dup_vca = [],[],[],[],[]
        if k in ks:  # the target chrom to test against, only these get the layer 2 SVs
            print('L2:generating disjoint pos list for DEL, DUP, INV SVs')
            # (1) apply  INS::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            if 'INS' in M[k]:
                ins_vca = vu.gen_var_calls(seqs,k,'INS',M[k]['INS'],G[k]['INS'])  # del_pos
                print('L2:applied %s INS mean-len=%s SVs to seq %s' % \
                      (len(ins_vca), vu.get_mean_size(M[k]['INS']), k))
            # (2) apply  DEL::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            if 'DEL' in M[k]:
                del_vca = vu.gen_var_calls(seqs, k, 'DEL', M[k]['DEL'], G[k]['DEL'])  # del_pos
                print('L2:applied %s DEL mean-len=%s SVs to seq %s' % \
                      (len(del_vca), vu.get_mean_size(M[k]['DEL']), k))
            # (3) apply large DUP:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            if 'DUP' in M[k]:
                dup_params = {'DUP': {'TYPE': mut_p[2]['DUP']['TYPE'],
                                      'TP':   mut_p[2]['DUP']['TP'],
                                      'CN':   mut_p[2]['DUP']['CN'],
                                      'CNP':  mut_p[2]['DUP']['CNP']}}
                dup_vca = vu.gen_var_calls(seqs, k, dup_params, M[k]['DUP'], G[k]['DUP'])  # dup_pos
                print('L2:applied %s DUP mean-len=%s SVs to seq %s' % \
                      (len(dup_vca), vu.get_mean_size(M[k]['DUP']), k))
            # (4) apply  INV->(inner DEL INS and SUB):::::::::::::::::::::::::::::::::::::
            if 'INV' in M[k]:
                complex_pos = []
                inv_params = {'INV': {'TYPE': mut_p[2]['INV']['TYPE'],
                                      'TP':   mut_p[2]['INV']['TP'],
                                      's:p':  mut_p[2]['INV']['s:p']}}
                inv_vca = vu.gen_var_calls(seqs, k, inv_params, M[k]['INV'], G[k]['INV'])  # inv pos
                print('L2:applied %s INV mean-len=%s SVs to seq %s' % \
                      (len(inv_vca), vu.get_mean_size(M[k]['INV']), k))
                complex_pos = vu.get_complex_inv_pos(inv_vca)
                dup_params = {'DUP': {'CN': [2, 3], 'TYPE': ['TANDEM']}}
                sm_l = max(list(mut_p[2]['INV']['s:p'].keys()))
                n_l,flank_l,flank_r = mut_p[2]['INV']['s:p'][sm_l]
                dup_flank_pos = vu.gen_mut_pos_flanking(complex_pos,int(sm_l),flank_l,flank_r)
                dup_flank_gen = vu.gen_genotypes(ploidy, dup_flank_pos, hetro=mut_p[2]['INV']['h'])  # same rate as INV
                if len(dup_flank_pos) > 0:
                    flanking_dup_vca = vu.gen_var_calls(seqs, k, dup_params, dup_flank_pos, dup_flank_gen)
                    if len(flanking_dup_vca) > 0:
                        print('L2:applied %s flanking DUP mean-len=%s SVs to seq %s' % \
                              (len(flanking_dup_vca), vu.get_mean_size(dup_flank_pos), k))
            # (5) apply TRA in either INS or DEL form
            if 'TRA' in M[k] and k in T:
                print('starting TRA calls for destination chrom %s' % k)
                tra_pos = [i['DPOS'] for i in T[k]]
                tra_params = {'TRA': {'TYPE': mut_p[2]['TRA']['TYPE'],
                                      'TP': mut_p[2]['TRA']['TP']}}
                tra_vca = vu.gen_var_calls(seqs, k, tra_params, tra_pos, G[k]['TRA'], tra_map=T[k])
                print('L2:applied %s TRA mean-len=%s SVs to seq %s' % \
                      (len(tra_vca), vu.get_mean_size(M[k]['TRA']), k))
            vca = sorted(ins_vca + del_vca + dup_vca + inv_vca + flanking_dup_vca + tra_vca,
                         key=lambda x: (x.chrom.zfill(100), x.pos))
            print('L2:editing mut1, mut2 strings for chrom %s' % k)
            mut1 = vu.apply_var_calls(mut1, vca, g=0)
            if ploidy > 1:
                mut2 = vu.apply_var_calls(mut2, vca, g=1)
            vcam[2][k] += vca
            vcam[2][k] = sorted(vcam[2][k], key=lambda x: (x.chrom.zfill(100), x.pos))
            # now work on the inner variants that are attached to complex SVs
            if 'INV' in M[k]:  # check for complex inv now
                inside_del_vca1,inside_del_vca2,inside_ins_vca1,inside_ins_vca2 = [],[],[],[]
                sm_l = min(mut_p[2]['INV']['s:p'])
                n_l,flank_l,flank_r = mut_p[2]['INV']['s:p'][sm_l]
                print('L2:Inner:acquiring inner positions for complex INV calls on chrom %s' % k)
                inner_pos = vu.gen_mut_pos_inside(complex_pos,mut_l=int(sm_l),mut_n=int(n_l),
                                                  break_ends=flank_l,low_bp=25)
                print('L2:Inner:partitioning for inner DEL and INS calls on chrom %s' % k)
                part_pos = vu.partition_pos(inner_pos,p=2,prop=[0.4, 0.6])  # split among INS,DEL
                if 0 in part_pos:  # del string modification and vcf generation
                    del_pos = part_pos[0]
                    del_gen = vu.get_inner_genotypes(complex_pos, G[k]['INV'],del_pos)
                    comp_del_pos = vu.reverse_complement_pos(complex_pos,del_pos)

                    del_vca = vu.gen_var_calls(seqs,k,'DEL',del_pos,del_gen)
                    comp_del_vca = vu.gen_var_calls(seqs, k, 'DEL', comp_del_pos, del_gen)
                    inside_del_vca1 = vu.update_vca_pos(del_vca,vcam[2][k],g=0)

                    if ploidy > 1:
                        inside_del_vca2 = vu.update_vca_pos(del_vca,vcam[2][k],g=1)

                        if len(inside_del_vca1) + len(inside_del_vca2) > 0:
                            print('L2:Inner:appplied %s small inner DEL to INV on seq %s' % \
                                  (len(inside_del_vca1) + len(inside_del_vca2), k))
                    else:
                        print('L2:Inner:appplied %s small inner DEL to INV on seq %s' % \
                              (len(inside_del_vca1), k))
                if 1 in part_pos:  # ins string modification and vcf generation
                    ins_pos = part_pos[1]
                    comp_ins_pos = vu.reverse_complement_pos(complex_pos, ins_pos)
                    ins_gen = vu.get_inner_genotypes(complex_pos, G[k]['INV'], ins_pos)

                    ins_vca = vu.gen_var_calls(seqs,k,'INS',ins_pos,ins_gen)
                    ins_str = {i:ins_vca[i].alt[::-1] for i in range(len(ins_vca))} #reversed insertions
                    comp_ins_vca = vu.gen_var_calls(seqs, k,'INS',comp_ins_pos,ins_gen,
                                                    insert_seqs=ins_str)
                    inside_ins_vca1 = vu.update_vca_pos(ins_vca,vcam[2][k],g=0)
                    if ploidy > 1:
                        inside_ins_vca2 = vu.update_vca_pos(ins_vca,vcam[2][k],g=1)
                        if len(inside_ins_vca1) + len(inside_ins_vca2) > 0:
                            print('L2:Inner:applied %s small inner INS to INV on seq %s' % \
                                  (len(inside_ins_vca1) + len(inside_ins_vca2), k))
                    else:
                        print('L2:Inner:applied %s small inner INS to INV on seq %s' % \
                              (len(inside_ins_vca1), k))
                # change this for each genotype to keep pos in mut1 and mut2 in tact
                vca1 = sorted(inside_del_vca1 + inside_ins_vca1, key=lambda x: x.pos)
                if ploidy > 1:
                    vca2 = sorted(inside_del_vca2 + inside_ins_vca2, key=lambda x: x.pos)
                print('L2:Inner:editing mut1, mut2 strings for chrom %s' % k)
                mut1 = vu.apply_var_calls(mut1, vca1, g=0)
                if ploidy > 1:
                    mut2 = vu.apply_var_calls(mut2, vca2, g=1)

                # correction of coordinates from each genotype
                vcam[2][k] = sorted(vcam[2][k] + comp_del_vca + comp_ins_vca,
                                    key=lambda x: (x.chrom.zfill(100), x.pos))
        # second round of SNV that are post SV
        pos,w = vu.wcu_to_pos_w(loss_wcu[k], len(seqs[k]))
        class_p,vca = {'SUB': [pos,w]},[]
        mut_pos = vu.gen_class_mut_pos_map(seqs[k], class_p,mut_p,l=3)
        if 'SUB' in mut_pos:
            mut_pos = mut_pos['SUB']
            mut_gen = vu.gen_genotypes(ploidy, mut_pos, hetro=mut_p[3]['SUB']['h'])
            vca = vu.gen_var_calls(seqs, k, 'SUB', mut_pos, mut_gen)
            vcam[3][k] += vca
            vcam[3][k] = sorted(vcam[3][k],key=lambda x:(x.chrom.zfill(100),x.pos))
            print('L3:applied %s SNVs to seq %s' % (len(mut_pos), k))
            mut1 = vu.apply_var_calls(mut1,vca,g=0)  # all ref_chroms get some SNV
            if ploidy > 1:
                mut2 = vu.apply_var_calls(mut2,vca,g=1)  # all ref_chroms get some SNV
        # test a modification that does chr1_DNA1, chr1_DNA2, ect
        print('writing the results to a fasta')
        if ploidy > 1:
            utils.write_fasta({k+'.1':mut1,k+'.2':mut2},fa_out_path,mode='a',index=True,gz=gz)
        else:
            utils.write_fasta({k+'.1':mut1},fa_out_path,mode='a',index=True,gz=gz)
        print('%s--------------------------------\n' % (''.join(['-' for i in range(len(k))])))
    # write a g1k style VCF 4.2 file for FusorSV S0 inputs for big DEL,DUP,INV call
    if not write_snv_indel:
        print('saving the large SVs to VCF')
        final_sv_vca = vu.merge_filter_sort_vcam(vcam[2],{},small_cut=small_cut)
        vu.write_vcf([sample],final_sv_vca,seqs,ref_path,vcf_out_path)
        vu.wcu_to_bedgraph(loss_wcu,out_dir+'/%s_loss.bedgraph'%sample,sign=-1.0)
        vu.wcu_to_bedgraph(gain_wcu, out_dir+'/%s_gain.bedgraph'%sample,sign=1.0)
    else:
        print('saving SNV/MNV,INDEL and large SVs to VCF')
        snv_indel_vcam = vu.vcam_union(vcam[1],vcam[3])
        final_vcam     = vu.vcam_union(vcam[2],snv_indel_vcam)
        final_vca      = vu.merge_filter_sort_vcam(final_vcam,{},small_cut=0)
        vu.write_vcf([sample],final_vca,seqs,ref_path,vcf_out_path)
        vu.wcu_to_bedgraph(loss_wcu,out_dir+'/%s_loss.bedgraph'%sample,sign=-1.0)
        vu.wcu_to_bedgraph(gain_wcu,out_dir+'/%s_gain.bedgraph'%sample,sign=1.0)
    f_stop = time.time()
    print('genome generated in %s sec' % round(f_stop-f_start, 2))
    return sample,vcam,loss,gain,rate

#somatic workflow that generates heterogenious data through a CSC or sub clonal model
def somatic_genomes(ref_path,out_dir,sample,gs,vcam,mut_p,loss_wcu,gain_wcu,gene_map,
                    gen_method='fast',gz=True,clean=True,gen_center=False,write_snv_indel=True,
                    small_cut=50,clone_tree_path=None,clone_tree_params=None,aneuploidy=None):
    if gz: ext_pat,idx_pat = '.fa.gz','i'
    else:  ext_pat,idx_pat = '.fa','.fai'
    #[1] generate the clone tree using model=0.0->CSC or model=1.0->sub clonal for hetergeneity
    if clone_tree_path is None and clone_tree_params is None:
        print('need to pass a clonal_tree json file path or use parameters to build a new one...')
        raise AttributeError
    if clone_tree_path is not None:
        CT = CloneTree(ref_path=ref_path,out_dir=out_dir,json_path=clone_tree_path)
    if clone_tree_params is not None: #command lines will overide the json file
        CT = CloneTree(ref_path,out_dir,sample,
                       clone_tree_params['model'],
                       clone_tree_params['branch'],
                       clone_tree_params['decay'],
                       clone_tree_params['cov']/2,
                       clone_tree_params['cov']/2,
                       gz=gz)
    print('%s clones are active'%len(CT.nodes))
    ks = list(set(list(vcam[1].keys())+list(vcam[2].keys())+list(vcam[3].keys())))
    s_vcam    = {k:{} for k in vcam}
    s_vcam[1] = vu.update_germline_genotypes(vcam[1],len(CT.nodes))
    s_vcam[2] = vu.update_germline_genotypes(vcam[2],len(CT.nodes))
    s_vcam[3] = vu.update_germline_genotypes(vcam[3],len(CT.nodes))
    germ_fa      = out_dir+'/'+sample+ext_pat
    fa_out_path  = out_dir+'/'+sample+'.somatic'+ext_pat #final merged fasta name
    vcf_out_path = out_dir+'/'+sample+'.somatic_S0.vcf'     #final somatic VCF only has novel somatic calls?
    som_fa_map   = {clone:out_dir+'/'+clone for clone in CT.nodes}
    # --------------------------------------------------------------------------------------------

    M,G,T,rate = {},{},{} ,{} # generate non-overlapping regions for SVs in layer 2 of mut_p
    for k in gs:
        print('building initial SV distributions for seq %s' % k)
        seqs = {k:ru.read_fasta_chrom(ref_path,k)} #one at a time now
        l_pos,l_w = vu.wcu_to_pos_w(loss_wcu[k],len(seqs[k])) #option to set background and scaling
        g_pos,g_w = vu.wcu_to_pos_w(gain_wcu[k],len(seqs[k])) #option to set background and scaling
        class_p = {'DEL':copy.deepcopy([l_pos,l_w]),'INV':copy.deepcopy([l_pos,l_w]),
                   'INS':copy.deepcopy([l_pos,l_w]),'TRA':copy.deepcopy([l_pos,l_w]),
                   'DUP':copy.deepcopy([g_pos,g_w])}
        if gen_method=='fast':
            M[k] = vu.gen_class_mut_pos_map(seqs[k],copy.deepcopy(class_p),mut_p,
                                            l=2,size_prop=2.0,center=gen_center,germ=False)
        else:
            M[k] = vu.gen_class_mut_pos_map_slow(seqs[k],copy.deepcopy(class_p),mut_p,
                                                 l=2,size_prop=2.0,center=gen_center)
        rate[k] = vu.mut_p_rate(mut_p,M[k],gs[k],l=2)
    l_pos,l_w,g_pos,g_w = [],[],[],[]

    #check overlapping events and filter as needed
    # old_pos_map,over = vu.vcam_to_pos_map(vcam[2]),{}
    # for k in old_pos_map:
    #     over[k] = []
    #     if M.has_key(k):
    #         N = copy.deepcopy(M[k])
    #         N['germ'] = []
    #         for pos in old_pos_map[k]: N['germ']  += [pos]
    #         N['germ'] = sorted(N['germ'],key=lambda x: x[0])
    #         over[k] = vu.mut_overlap(N)
    #         #if [i[3].keys() for i in over]
    #         if not vu.mut_over_detected(over[k]):
    #             print('no somatic to germline overlap detected for seq=%s'%k)
    #         else:
    #             print('somatic to germline overlap detected for seq=%s'%k)
    #             print(over[k])
    #return seqs[k], copy.deepcopy(class_p), mut_p
    old_pos_map,over = {},{}

    # correct areas of gain for onco genes------------------------------------
    print('correcting SV events to boundaries for seq %s' % k)
    extend_wcu = {k:gain_wcu[k]['onco'] if 'onco' in gain_wcu[k] else [] for k in gain_wcu}
    M = vu.extend_class_mut_pos(M,'DUP',copy.deepcopy(extend_wcu))
    loss,gain = wcu_enrichment(M,loss_wcu,gain_wcu,gene_map)
    extend_wcu = {}

    # set genotypes on the corrected events in M------------------------------
    # for somatic will have 'h' = 1.0 due to DNA repair and replication events
    for k in M:
        G[k] = {}
        for t in M[k]:
            ploidy = 2
            if 'Y' in ks and 'X' in ks and (k == 'Y' or k == 'X'): ploidy = 1
            G[k][t] = vu.gen_genotypes(ploidy, M[k][t], hetro=mut_p[2][t]['h'])
    print('generating TRA mapping from available seqs\n')
    TK = {}
    for k in M:
        if 'TRA' in M[k]: TK[k] = copy.deepcopy(M[k]['TRA'])
    T = vu.gen_tra_map(TK)

    s_A = None
    if aneuploidy is not None:
        print('aneuploidy is being modeled')
        ploidy_map = {k:aneuploidy[k]['n'] for k in aneuploidy}
        print('base aneuploidy is %s'%aneuploidy)
        mod_aneuploidy = vu.adjust_aneuploidy_effect(aneuploidy,ploidy_map,loss,dist='uniform')
        print('modified aneuploidy is %s'%mod_aneuploidy)
        s_A = vu.gen_aneuploidy({k:ru.read_fasta_chrom(ref_path,k) for k in ks},ploidy_map,CT,mod_aneuploidy)
        print('aneuploidy has been generated for chrom %s'%list(s_A.keys()))
        print(s_A)

    f_start = time.time()

    for k in gs: #reuse the old mnv_vcam and sv_vcam
        for l in s_vcam:
            if not k in s_vcam[l]: s_vcam[l][k] = []
        seqs = {k: ru.read_fasta_chrom(ref_path,k)}
        ploidy = 2
        if 'Y' in ks and 'X' in ks and (k == 'Y' or k == 'X'): ploidy = 1
        print('processing seq=%s with ploidy=%s----------------' % (k, ploidy))
        # start with a some SNV before SV ----------------------------------------------
        pos,w   = vu.wcu_to_pos_w(loss_wcu[k],len(seqs[k]))
        class_p = {'SUB':[pos,w]}
        mut_pos = vu.gen_class_mut_pos_map(seqs[k],class_p,mut_p,l=1)
        if 'SUB' in mut_pos: mut_pos = mut_pos['SUB']
        else:                mut_pos = []
        pos,w = [],[]
        mut_gen    = vu.gen_genotypes(ploidy, mut_pos, hetro=mut_p[1]['SUB']['h'])
        part_map   = vu.partition_pos(mut_pos,len(CT.nodes),prop=CT.mut_prop(),index_only=True)
        allele_map = vu.gen_clone_allele_map(part_map,CT)
        vca = vu.gen_var_calls(seqs,k,'SUB',mut_pos,mut_gen,
                               allele_map=allele_map,ref_path=ref_path)
        mut_pos = []
        #spool the first layer---------------------------------------------------
        s_vcam[1][k] += vca
        s_vcam[1][k] = sorted(s_vcam[1][k],key=lambda x:(x.chrom.zfill(100),x.pos))

        spool_fasta_clones(som_fa_map,k,ploidy,s_vcam[1][k],seqs[k],gz=gz)
        #spool the first layer---------------------------------------------------
        # now set total SVs and partition into DEL,DUP,INV,TRA
        ins_vca,del_vca,dup_vca,tra_vca,inv_vca = [],[],[],[],[]
        inside_del_vca,inside_ins_vca,comp_del_vca,comp_ins_vca,flanking_dup_vca = [],[],[],[],[]
        if k in gs:  # the target chrom to test against, only these get the layer 2 SVs
            print('L2:generating disjoint pos list for DEL, DUP, INV SVs')
            # (1) apply  INS::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            if 'INS' in M[k]:
                part_map   = vu.partition_pos(M[k]['INS'],len(CT.nodes),prop=CT.mut_prop(),index_only=True)
                allele_map = vu.gen_clone_allele_map(part_map,CT)
                ins_vca = vu.gen_var_calls(seqs,k,'INS',M[k]['INS'],G[k]['INS'],
                                           allele_map=allele_map,ref_path=ref_path)  # del_pos
                print('L2:applied %s INS mean-len=%s SVs to seq %s' % \
                      (len(ins_vca), vu.get_mean_size(M[k]['INS']), k))
            # (2) apply  DEL::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            if 'DEL' in M[k]:
                part_map   = vu.partition_pos(M[k]['DEL'],len(CT.nodes),prop=CT.mut_prop(),index_only=True)
                allele_map = vu.gen_clone_allele_map(part_map,CT)
                del_vca = vu.gen_var_calls(seqs,k,'DEL',M[k]['DEL'],G[k]['DEL'],
                                           allele_map=allele_map,ref_path=ref_path)  # del_pos
                print('L2:applied %s DEL mean-len=%s SVs to seq %s' % \
                      (len(del_vca), vu.get_mean_size(M[k]['DEL']), k))
            # (3) apply large DUP:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            if 'DUP' in M[k]:
                part_map   = vu.partition_pos(M[k]['DUP'],len(CT.nodes),prop=CT.mut_prop(),index_only=True)
                allele_map = vu.gen_clone_allele_map(part_map,CT)
                dup_params = {'DUP': {'TYPE': mut_p[2]['DUP']['TYPE'],
                                      'TP':   mut_p[2]['DUP']['TP'],
                                      'CN':   mut_p[2]['DUP']['CN'],
                                      'CNP':  mut_p[2]['DUP']['CNP']}}
                dup_vca = vu.gen_var_calls(seqs,k,dup_params,M[k]['DUP'],G[k]['DUP'],
                                           allele_map=allele_map,ref_path=ref_path)  # dup_pos
                print('L2:applied %s DUP mean-len=%s SVs to seq %s' % \
                      (len(dup_vca), vu.get_mean_size(M[k]['DUP']), k))
            # (4) apply  INV->(inner DEL INS and SUB):::::::::::::::::::::::::::::::::::::
            if 'INV' in M[k]:
                part_map    = vu.partition_pos(M[k]['INV'],len(CT.nodes),prop=CT.mut_prop(),index_only=True)
                allele_map  = vu.gen_clone_allele_map(part_map,CT)
                complex_pos = []
                inv_params  = {'INV': {'TYPE': mut_p[2]['INV']['TYPE'],
                                       'TP': mut_p[2]['INV']['TP'],
                                       's:p': mut_p[2]['INV']['s:p']}}
                inv_vca = vu.gen_var_calls(seqs,k,inv_params,M[k]['INV'],G[k]['INV'],
                                           allele_map=allele_map,ref_path=ref_path)  # inv pos
                print('L2:applied %s INV mean-len=%s SVs to seq %s' % \
                      (len(inv_vca), vu.get_mean_size(M[k]['INV']), k))
                complex_pos = vu.get_complex_inv_pos(inv_vca)
                if len(complex_pos)>0:
                    dup_params = {'DUP': {'CN': [2, 3], 'TYPE': ['TANDEM']}}
                    sm_l = max(list(mut_p[2]['INV']['s:p'].keys()))
                    n_l,flank_l,flank_r = mut_p[2]['INV']['s:p'][sm_l]
                    dup_flank_pos = vu.gen_mut_pos_flanking(complex_pos,int(sm_l),flank_l,flank_r)
                    dup_flank_gen = vu.gen_genotypes(ploidy, dup_flank_pos, hetro=mut_p[2]['INV']['h'])  # same rate as INV
                    if len(dup_flank_pos) > 0:
                        part_map = vu.partition_pos(dup_flank_pos,len(CT.nodes),
                                                    prop=[1.0/(x+1) for x in range(len(CT.nodes))][::-1],
                                                    index_only=True)
                        allele_map = vu.gen_clone_allele_map(part_map,CT)
                        flanking_dup_vca = vu.gen_var_calls(seqs,k,dup_params,dup_flank_pos,dup_flank_gen,
                                                            allele_map=allele_map,ref_path=ref_path)
                        if len(flanking_dup_vca) > 0:
                            print('L2:applied %s flanking DUP mean-len=%s SVs to seq %s' % \
                                  (len(flanking_dup_vca), vu.get_mean_size(dup_flank_pos), k))
            # (5) apply TRA in either INS or DEL form
            if 'TRA' in M[k] and k in T:
                print('L2:starting TRA calls for destination chrom %s' % k)
                tra_pos  = [i['DPOS'] for i in T[k]]
                part_map = vu.partition_pos(tra_pos,len(CT.nodes),prop=CT.mut_prop(),index_only=True)
                allele_map = vu.gen_clone_allele_map(part_map,CT)
                tra_params = {'TRA': {'TYPE': mut_p[2]['TRA']['TYPE'],
                                      'TP':   mut_p[2]['TRA']['TP']}}
                tra_vca = vu.gen_var_calls(seqs,k,tra_params,tra_pos,G[k]['TRA'],
                                           tra_map=T[k],allele_map=allele_map,ref_path=ref_path)
                print('L2:applied %s TRA mean-len=%s SVs to seq %s' % \
                      (len(tra_vca), vu.get_mean_size(M[k]['TRA']), k))
            vca = sorted(ins_vca +del_vca+dup_vca+inv_vca+flanking_dup_vca+tra_vca,
                         key=lambda x: (x.chrom.zfill(100), x.pos))
            #spool the second layer---------------------------------------------------
            s_vcam[2][k] += vca
            s_vcam[2][k] = sorted(s_vcam[2][k], key=lambda x: (x.chrom.zfill(100), x.pos))
            spool_fasta_clones(som_fa_map,k,ploidy,s_vcam[2][k],seqs[k],gz=gz)
            #spool the second layer---------------------------------------------------
            # now work on the inner variants that are attached to complex SVs
            if 'INV' in M[k]:  # check for complex inv now
                inv_del_vca,comp_del_vca,inv_ins_vca,comp_ins_vca = [],[],[],[]
                sm_l = min(mut_p[2]['INV']['s:p'])
                n_l,flank_l,flank_r = mut_p[2]['INV']['s:p'][sm_l]
                print('L2:Inner:acquiring inner positions for complex INV calls on chrom %s' % k)
                inner_pos = vu.gen_mut_pos_inside(complex_pos,mut_l=int(sm_l),mut_n=int(n_l),
                                                  break_ends=flank_l,low_bp=25)
                print('L2:Inner:partitioning for inner DEL and INS calls on chrom %s' % k)
                part_pos = vu.partition_pos(inner_pos,p=2,prop=[0.4, 0.6])  # split among INS,DEL
                if 0 in part_pos:  # del string modification and vcf generation
                    inv_del_pos  = part_pos[0]
                    inv_del_gen  = vu.get_inner_genotypes(complex_pos, G[k]['INV'],inv_del_pos)
                    comp_del_pos = vu.reverse_complement_pos(complex_pos,inv_del_pos)
                    part_map     = vu.partition_pos(inv_del_pos,len(CT.nodes),prop=CT.mut_prop(),index_only=True)
                    del_al_map   = vu.gen_clone_allele_map(part_map,CT)
                    inv_del_vca  = vu.gen_var_calls(seqs,k,'DEL',inv_del_pos,inv_del_gen,
                                                    allele_map=del_al_map,ref_path=ref_path)
                    comp_del_vca = vu.gen_var_calls(seqs,k,'DEL',comp_del_pos,inv_del_gen,
                                                    allele_map=del_al_map,ref_path=ref_path)
                    print('L2:Inner:appplied %s small inner DEL to INV on seq %s'%(len(inv_del_vca),k))
                if 1 in part_pos:  # ins string modification and vcf generation
                    inv_ins_pos  = part_pos[1]
                    inv_ins_gen  = vu.get_inner_genotypes(complex_pos,G[k]['INV'],inv_ins_pos)
                    comp_ins_pos = vu.reverse_complement_pos(complex_pos,inv_ins_pos)
                    part_map     = vu.partition_pos(inv_ins_pos,len(CT.nodes),prop=CT.mut_prop(),index_only=True)
                    ins_al_map   = vu.gen_clone_allele_map(part_map,CT)
                    inv_ins_vca  = vu.gen_var_calls(seqs,k,'INS',inv_ins_pos,inv_ins_gen,
                                                    allele_map=ins_al_map,ref_path=ref_path)
                    comp_ins_vca = vu.gen_var_calls(seqs,k,'INS',comp_ins_pos,inv_ins_gen,
                                                    allele_map=ins_al_map,ref_path=ref_path)
                    print('L2:Inner:appplied %s small inner INS to INV on seq %s'%(len(ins_vca),k))
                #now write out the complex INV inner events
                vca          = sorted(del_vca+ins_vca,key=lambda x:(x.chrom.zfill(100),x.pos))
                r_vca        = sorted(comp_del_vca+comp_ins_vca,key=lambda x:(x.chrom.zfill(100),x.pos))
                s_vcam[2][k] = sorted(s_vcam[2][k]+r_vca,key=lambda x: (x.chrom.zfill(100),x.pos))

                #check for nhej and mmej here-----------------------------------------
                #s_vcam[2][k] = vu.adjust_del_effect()
                #s_vcam[2][k] = vu.adjust_ins_effect()

                spool_fasta_clones(som_fa_map,k,ploidy,vca,seqs[k],complex=s_vcam[2][k],gz=gz)
        # second round of SNV that are post SV
        pos,w = vu.wcu_to_pos_w(loss_wcu[k], len(seqs[k]))
        class_p = {'SUB': [pos,w]}
        mut_pos = vu.gen_class_mut_pos_map(seqs[k], class_p,mut_p,l=3)
        if 'SUB' in mut_pos: mut_pos = mut_pos['SUB']
        else:                mut_pos = []
        pos,w = [],[]
        mut_gen    = vu.gen_genotypes(ploidy, mut_pos, hetro=mut_p[3]['SUB']['h'])
        part_map   = vu.partition_pos(mut_pos,len(CT.nodes),prop=CT.mut_prop(),index_only=True)
        allele_map = vu.gen_clone_allele_map(part_map,CT)
        vca = vu.gen_var_calls(seqs,k,'SUB',mut_pos,mut_gen,
                               allele_map=allele_map,ref_path=ref_path)
        mut_pos = []
        s_vcam[3][k] += vca
        s_vcam[3][k] = sorted(s_vcam[3][k],key=lambda x:(x.chrom.zfill(100),x.pos))
        spool_fasta_clones(som_fa_map,k,ploidy,s_vcam[3][k],seqs[k],anueuploidy=s_A,gz=gz)

    # test a modification that does chr1_DNA1, chr1_DNA2, ect
    m_start = time.time()
    print('merging all spooled seqs to the somatic fasta')
    merge_fasta_clones(som_fa_map,fa_out_path,aneuploidy=None,living=None,clean=clean,gz=gz)
    m_stop = time.time()
    print('spooled seq were merged in %s sec'%round(m_stop-m_start,2))
    if not clean:
        c_start = time.time()
        print('merging indivdual clone sequences to clone fastas')
        merge_each_fasta_clone(som_fa_map,out_dir,aneuploidy=None,living=None,clean=clean,gz=gz)
        c_stop = time.time()
        print('individual seq were merged in %s sec' % round(c_stop-c_start, 2))
    # write a g1k style VCF 4.2 file for FusorSV S0 inputs for big DEL,DUP,INV call
    seqs = {k:ru.read_fasta_chrom(ref_path,k) for k in gs}
    if not write_snv_indel:
        print('saving the large SVs to VCF')
        final_vca = vu.merge_filter_sort_vcam(s_vcam[2],s_A,small_cut=small_cut)
        vu.write_somatic_vcf(CT,final_vca,seqs,ref_path,vcf_out_path)  #has all events not in germline
        vu.write_clone_vcf(vcf_out_path,out_dir)  #has events only for each clone
        vu.wcu_to_bedgraph(loss_wcu,out_dir+'/%s_somatic_loss.bedgraph'%sample,sign=-1.0)
        vu.wcu_to_bedgraph(gain_wcu,out_dir+'/%s_somatic_gain.bedgraph'%sample,sign=1.0)
    else:
        print('saving SNV/MNV,INDEL and large SVs to VCF')
        final_vcam     = vu.vcam_union(vu.vcam_union(s_vcam[1],s_vcam[3]),s_vcam[2])
        final_vca      = vu.merge_filter_sort_vcam(final_vcam,s_A,small_cut=0)
        vu.write_somatic_vcf(CT,final_vca,seqs,ref_path,vcf_out_path)  #has all events not in germline
        if not clean:
            vu.write_clone_vcf(vcf_out_path,out_dir)  #has events only for each clone
        vu.wcu_to_bedgraph(loss_wcu,out_dir+'/%s_somatic_loss.bedgraph'%sample,sign=-1.0)
        vu.wcu_to_bedgraph(gain_wcu,out_dir+'/%s_somatic_gain.bedgraph'%sample,sign=1.0)
    f_stop = time.time()
    print('genome generated in %s sec'%round(f_stop-f_start,2))
    return CT,M,T,s_vcam,loss,gain,rate,over

#!/usr/bin/env python
#Timothy James Becker, PhD candidate, UCONN 01/10/2017-03/20/2020
#soMaCX: Somatic Complex SV Diploid Genome Generator, produces disjoint SVs on demand and then
#allows for secondary SV layers to be introduced using several parameters
#supported types are INS, DEL, SUB, TANDEM DUP, DISPERSED DUP, INV DUP, INV,COMPLEX INV, TRA (chroma->chromB)
#each of these can have inner layers of INS, DEL, DUP, ect using a distribution and
#a separate breakpoint probability that allows for easily setting DEL, INS directly on
#DUP or INV break ends, thereby providing a realiable, distribution-based test platform for simulation
#as FASTA and VCF formated files provide easy benchmarking

import argparse
import os
import sys
import glob
import gzip
import json
import time
import numpy as np
import somacx.read_utils as ru
import somacx.variant_utils as vu
import somacx.sim as sim

if sys.version_info.major<3:
    def str_hook(obj):
        return {k.encode('utf-8') if isinstance(k,unicode) else k:
                v.encode('utf-8') if isinstance(v, unicode) else v for k,v in obj}
else:
    def str_hook(obj):
        return {k:v for k,v in obj}

#util functions for unpackaging and dispersing json data
def disperse_complex_generator_json(json_path,out_dir):
    if json_path.endswith('.gz'):
        with gzip.open(json_path,'rb') as f:
            raw = json.load(f,object_pairs_hook=str_hook)
    else:
        with open(json_path,'r') as f:
            raw = json.load(f,object_pairs_hook=str_hook)
    if not os.path.exists(out_dir+'/meta/'): os.makedirs(out_dir+'/meta/')
    for k in raw:
        with open('%s/%s.json'%(out_dir+'/meta/',k),'w') as f:
            f.write(json.dumps(raw[k]))
            print('copying %s to %s'%(k,'%s/%s.json'%(out_dir+'/meta/',k)))
    return True

#util functions for packaging json data into a single full json file
def write_complex_generator_json(json_path,in_dir,gene_map,g_var_map,g_loss_wild,g_gain_wild,
                                 s_var_map,s_loss_wild,s_gain_wild,s_anueuploidy,clone_tree,gz=True):

    C = {'gene_map':in_dir+gene_map,'g_var_map':in_dir+g_var_map,'s_var_map':in_dir+s_var_map,
         'clone_tree':in_dir+clone_tree,'somatic_aneueploidy':in_dir+s_anueuploidy}
    W = [in_dir+g_loss_wild,in_dir+g_gain_wild,in_dir+s_loss_wild,in_dir+s_gain_wild]
    for c in C:
        print('loading %s'%c)
        if C[c].endswith('.gz'):
            with gzip.open(C[c],'rb') as f:
                C[c.rsplit('/')[-1].rsplit('.json.gz')[0]] = json.load(f,object_pairs_hook=str_hook)
        else:
            with open(C[c],'r') as f:
                C[c.rsplit('/')[-1].rsplit('.json')[0]] = json.load(f,object_pairs_hook=str_hook)
    for w in W:
        for k in glob.glob(w):
            print('loading %s'%k)
            if k.endswith('.gz'):
                with gzip.open(k,'rb') as f:
                    C[k.rsplit('/')[-1].rsplit('.json.gz')[0]] = json.load(f,object_pairs_hook=str_hook)
            else:
                with open(k,'r') as f:
                    C[k.rsplit('/')[-1].rsplit('.json.gz')[0]] = json.load(f,object_pairs_hook=str_hook)
    if gz:
        if not json_path.endswith('.gz'): json_path += '.gz'
        if sys.version_info.major<3:
            with gzip.open(in_dir+json_path,'w') as f:
                f.write(json.dumps(C))
                return True
        else:
            with gzip.open(in_dir+json_path,'wb') as f:
                f.write(bytes(json.dumps(C),'utf_8'))
                return True
        return False
    else:
        with open(in_dir+json_path,'w') as f:
            f.write(json.dumps(C))
            return True
        return False

des = """soMaCX: Generator v0.1.0, 01/01/2017-02/02/2020 Timothy James Becker"""
parser = argparse.ArgumentParser(description=des)
parser.add_argument('-r','--ref_path',type=str,help='reference fasta input file\t[None]')
parser.add_argument('-o','--out_dir',type=str,help='output directory\t[None]')
parser.add_argument('-c','--chroms',type=str,help='comma seperate chrom list to create large DEL,DUP,INV on\t[11,20]')
parser.add_argument('-C','--ref_chroms',type=str,help='comma seperate chrom list to create fasta with\t[1-22,X,Y,MT]')
parser.add_argument('-j','--complex_generator_json',type=str,help='full somacx complex generator json file.\t[None]')
parser.add_argument('--vcf',type=str,help='provide a FusorSV format VCF4 file for SVs generation at germline step\t[None]')
parser.add_argument('--method',type=str,help='fast or slow method for generation\t[fast]')
parser.add_argument('--uncomp',action='store_true',help='don\'t use htslib/bgzip .gz compression for read/write\t[False]')
parser.add_argument('--single',action='store_true',help='don\'t delete single clone genomes\t[False]')
parser.add_argument('--center',type=bool,help='for somatic position generation, center the ranges\t[False]')
parser.add_argument('--model',type=str,help='[0.0,1.0] value where 0.0=>CSC and 1.0=>sub clonal,csvs will indicate a random range\t[0.5]')
parser.add_argument('--branch',type=str,help='[0.0,1.0] branching factor or chance to branch each cycle, where model=1.0 and >= 1.0 => binary tree,csvs will indicate a random range\t[0.5]')
parser.add_argument('--decay',type=str,help='[0.0,1.0] decay factor or chance that a cell node dies each cycle, where 1.0 is certain cell death,csvs will indicate a random range\t[0.001]')
parser.add_argument('--cov',type=str,help='>2 total coverage/ploidy that will be simulated, used as a lower bound so data is not clipped past noise level,csvs will indicate a random range\t[10]')
parser.add_argument('--small_cut',type=int,help='filter out variants that are less than this size in bp\t[1]')
parser.add_argument('--verbose',action='store_true',help='output more result mertics to stdout\t[False]')
args = parser.parse_args()

print('starting: \n'+des+'\n')
if args.ref_path is not None:
    ref_path = args.ref_path
    fasta_path = ref_path
    print('reading reference meta data')
    S = ru.get_fasta_seq_names_lens(ref_path)
    ratios = vu.get_ratios(S) #ratios of the length in scale to largest AKA chr1
    print('fasta=%s meta data read'%fasta_path.rsplit('/')[-1])
else:
    print('reference fasta was not found')
    raise IOError
if args.vcf is not None:
    prior_vcf = args.vcf
    compatible = True
    vcf_seqs = vu.vcf_chroms(prior_vcf) #
    for v in vcf_seqs:
        if v not in S: compatible = False
    if not compatible:
        print('VCF input file is not compatible with the reference fasta!')
        print('ref seqs: %s'%(S))
        print('vcf seqs: %s'%(vcf_seqs.keys()))
        raise AttributeError
else:
    prior_vcf = None
if args.out_dir is not None:
    out_dir = args.out_dir
    if not os.path.exists(out_dir): os.makedirs(out_dir)
else:
    print('output directory not specified')
    raise IOError
if args.chroms is not None:
    ks = args.chroms.split(',')
else:
    ks = [str(i) for i in range(1,23)]+['X','Y','MT']
if args.ref_chroms is not None:
    rs = args.ref_chroms.split(',')
    kk = list(S.keys())
    for k in kk:
        if k not in rs: S.pop(k)
    kk = sorted(S,key=lambda x: x.zfill(max([len(k) for k in S])))
else:#1-22,X,Y,MT and chr checking prefix 
    rs,ss = [],[str(i) for i in range(1,23)]+['X','Y','MT']
    kk = sorted(S,key=lambda x: x.zfill(max([len(k) for k in S])))
    for k in kk:
        for s in ss:
            if k.upper()==s.upper():           rs += [k]
            elif k.upper()==('chr'+s).upper(): rs += [k]
    kk = list(S.keys())
    for k in kk:
        if k not in rs: S.pop(k)
    kk = sorted(S,key=lambda x: x.zfill(max([len(k) for k in S])))
print('starting SV generation with seqs=%s'%kk)
if args.complex_generator_json == 'DEFAULT':
    disperse_complex_generator_json(vu.get_local_path()+'full.json',out_dir)
    args.complex_generator_json = vu.get_local_path()+'full.json'
elif args.complex_generator_json is not None:
    disperse_complex_generator_json(args.complex_generator_json,out_dir)
else:
    full_json = None
if args.method is not None:          method = args.method
else:                                method = 'fast'
if args.uncomp:                      gz     = False
else:                                gz     = True
if args.single:                      clean  = False
else:                                clean  = True
if args.center is not None:          center = args.center
else:                                center = False
if args.model is not None:           model  = sorted([float(x) for x in args.model.rsplit(',')])
else:                                model  = [0.25,0.75]
if len(model)>0 and len(model)<=1:   model  = model[0]
else:                                model  = np.random.choice(np.arange(model[0],model[1],(model[1]-model[0])/10.0),1)[0]
if args.branch is not None:          branch = sorted([float(x) for x in args.branch.rsplit(',')])
else:                                branch = [0.5,1.0]
if len(branch)>0 and len(branch)<=1: branch = branch[0]
else:                                branch = np.random.choice(np.arange(branch[0],branch[1],(branch[1]-branch[0])/10.0),1)[0]
if args.decay is not None:           decay  = sorted([float(x) for x in args.decay.rsplit(',')])
else:                                decay  = [0.001,0.0001]
if len(decay)>0 and len(decay)<=1:   decay  = decay[0]
else:                                decay  = np.random.choice(np.arange(decay[0],decay[1],(decay[1]-decay[0])/10.0),1)[0]
if args.cov is not None:             cov    = sorted([float(x) for x in args.cov.rsplit(',')])
else:                                cov    = [15,35]
if len(cov)>0 and len(cov)<=1:       cov    = cov[0]
else:                                cov    = np.random.choice(np.arange(cov[0],cov[1],(cov[1]-cov[0])/10.0),1)[0]
if args.small_cut is not None: small_cut = args.small_cut
else:                          small_cut = 1
#if a user puts in any clone tree params, this will overide the full.json file
clone_tree_params = None
if args.model is not None or args.branch is not None or \
    args.decay is not None or  args.cov is not None:
    clone_tree_params = {'model':model,'branch':branch,'decay':decay,'cov':cov}
    clone_tree_path   = None
else: clone_tree_path = out_dir+'/meta/clone_tree.json'
start = time.time()
print('reading gene map...')
if not os.path.exists(out_dir+'/meta/gene_map.json'):
    print('didn\'t find a pre-processed gene map, building a new one...')
    gene_map = vu.build_ucsc_gene_exon_map(vu.get_local_path('refGene.txt.gz'))
    vu.write_json_gene_map(out_dir+'/meta/gene_map.json',gene_map)
    print('completed building a json gene map')
gene_map    = vu.read_json_gene_map(out_dir+'/meta/gene_map.json')

#::UPDATING DISTRIBUTIONS::######################################################################################
#read a genelist and store as a wcu json file for building the full json
mitcp_genes = vu.read_gene_list(vu.get_local_path('mitcp_giam_gene_list.txt'))
mitcp_wcu   = vu.gene_list_to_wcu(mitcp_genes,[0.0 for i in range(len(mitcp_genes))],gene_map,label='mitcp')
vu.write_json_wcu(out_dir+'/meta/'+'germline_mitcp_gain_wcu.json',mitcp_wcu)
mitcp_wcu   = vu.gene_list_to_wcu(mitcp_genes,[0.0 for i in range(len(mitcp_genes))],gene_map,label='mitcp')
vu.write_json_wcu(out_dir+'/meta/'+'germline_mitcp_loss_wcu.json',mitcp_wcu)
mitcp_wcu   = vu.gene_list_to_wcu(mitcp_genes,[0.0 for i in range(len(mitcp_genes))],gene_map,label='mitcp')
vu.write_json_wcu(out_dir+'/meta/'+'somatic_mitcp_gain_wcu.json',mitcp_wcu)
mitcp_wcu   = vu.gene_list_to_wcu(mitcp_genes,[1E2 for i in range(len(mitcp_genes))],gene_map,label='mitcp')
vu.write_json_wcu(out_dir+'/meta/'+'somatic_mitcp_loss_wcu.json',mitcp_wcu)
# #write out a full complex_generator_json from several seperate ones

print('updating the full param file')
write_complex_generator_json(in_dir=out_dir+'/meta/',
                             json_path='full.json',
                             gene_map='gene_map.json',
                             g_var_map='g_var_map.json',
                             g_loss_wild='germline_*_loss_wcu.json',
                             g_gain_wild='germline_*_gain_wcu.json',
                             s_var_map='s_var_map.json',
                             s_loss_wild='somatic_*_loss_wcu.json',
                             s_gain_wild='somatic_*_gain_wcu.json',
                             s_anueuploidy='somatic_anueuploidy.json',
                             clone_tree='clone_tree.json')
print('finished with full param file')
disperse_complex_generator_json(args.complex_generator_json,out_dir) #copies all params to meta folder
#::UPDATING DISTRIBUTIONS::######################################################################################

#load the first germline mutation map derived from g1kp3 distributions
germline_var_map = vu.read_json_mut_map(out_dir+'/meta/g_var_map.json')
#alu_wcu = vu.vcf_to_wcu(vu.get_local_path('melt2_prior_mei.tar.gz') #have to use tarfile lib

#default weighted class units for germline regions that are amenable to variations
g_loss_raw = {g.rsplit('germline_')[-1].rsplit('_loss_wcu.json')[0]:vu.read_json_wcu(g) \
              for g in glob.glob(out_dir+'/meta/germline_*_loss_wcu.json')}
loss_wcu = vu.merge_wcu(g_loss_raw)
g_gain_raw = {g.rsplit('germline_')[-1].rsplit('_gain_wcu.json')[0]:vu.read_json_wcu(g) \
              for g in glob.glob(out_dir+'/meta/germline_*_gain_wcu.json')}
gain_wcu = vu.merge_wcu(g_gain_raw)

#add empty keys for chroms that are not affected by this
for k in rs:
    if k not in loss_wcu: loss_wcu[k] = {}
    if k not in gain_wcu: gain_wcu[k] = {}

if prior_vcf is not None: #via FusorSV_VCF file
    print('loading prior vcf file: %s'%prior_vcf)
    sample = prior_vcf.split('/')[-1].split('_')[0]
    vcf_vcam = vu.vcf_to_vcam(prior_vcf,ref_path,ks) #replace with VCF SVs
    if 'Y' not in vcf_vcam and 'Y' in ks: ks.remove('Y') #update ks
    if 'Y' not in vcf_vcam and 'Y' in rs: rs.remove('Y') #update rs
    if 'Y' not in vcf_vcam and 'Y' in S: S.pop('Y')
    print('checking VCF call conflicts')
    vcf_vcam = vu.vcam_remove_conflicts(vcf_vcam)
    print('variantion is present on: %s'%list(vcf_vcam.keys()))

print('ks = %s rs = %s at germline_genome start'%(ks,rs))
#not fully expressed, but getting closer to germline goal once transposons are in place
sample,vcam,g_loss,g_gain,g_rate = sim.germline_genome(ref_path,out_dir,rs,ks,germline_var_map,loss_wcu,gain_wcu,
                                                       gene_map,gen_method=method,gz=gz,small_cut=small_cut)

if prior_vcf is not None: #via FusorSV_VCF file
    vcf_sample = prior_vcf.split('/')[-1].split('_')[0]
    print('hybrid %s-L1:MNV + %s-L2:SV + %s-L3:MNV'%(sample,vcf_sample,sample))
    if 2 in vcam: vcam.pop(2)
    vcam[2] = vcf_vcam
    vcam = vu.alter_lvcam_genotypes(vcam,h=0.25) #for now assumes no ploidy 1

    # sample = vcf_sample    #
    # write_snv_indel = True #
    # small_cut = 50         #
    # gz = True              #
    sample,vcam,g_loss,g_gain,g_rate = sim.write_genome_from_vcam(ref_path,vcam,vcf_sample,out_dir,rs,gene_map)
#need to write function to apply fsvcam to ref_path+out_dir

#default weighted class units for somatic regions that include onco genes and nhej pathways...
somatic_var_map = vu.read_json_mut_map(out_dir+'/meta/s_var_map.json')
mnv_wcu = vu.vcam_to_wcu(vcam[1],0.0,'germ')
sv_wcu  = vu.vcam_to_wcu(vcam[2],0.0,'germ')
s_loss_raw = {g.rsplit('somatic_')[-1].rsplit('_loss_wcu.json')[0]:vu.read_json_wcu(g) \
              for g in glob.glob(out_dir + '/meta/somatic_*_loss_wcu.json')}

s_loss_raw['germ'] = sv_wcu
loss_wcu = vu.merge_wcu(s_loss_raw)
#functional hits can increase mut_p rates, aneuploidy

s_gain_raw = {g.rsplit('somatic_')[-1].rsplit('_gain_wcu.json')[0]:vu.read_json_wcu(g) \
              for g in glob.glob(out_dir + '/meta/somatic_*_gain_wcu.json')}
gain_wcu['germ'] = sv_wcu
gain_wcu = vu.merge_wcu(s_gain_raw)

for k in rs:
    if k not in loss_wcu: loss_wcu[k] = {}
    if k not in gain_wcu: gain_wcu[k] = {}

#sim.somatic_genomes test case-----------
anueuploidy = vu.read_json_anueuploidy(out_dir+ '/meta/somatic_anueuploidy.json')
mut_p,gs = somatic_var_map,S
print('starting initialization of somatic clonal tree using params:')
print('sample=%s clone_tree_params=%s clone_tree_path=%s'%(sample,clone_tree_params,clone_tree_path))
CT,M,T,s_vcam,s_loss,s_gain,s_rate,s_over = sim.somatic_genomes(ref_path,out_dir,sample,gs,vcam,
                                                                somatic_var_map,loss_wcu,gain_wcu,gene_map,
                                                                gen_method=method,gz=gz,clean=clean,
                                                                gen_center=center,
                                                                clone_tree_path=clone_tree_path,
                                                                clone_tree_params=clone_tree_params,
                                                                anueuploidy=anueuploidy,
                                                                small_cut=small_cut)
# print('germline generation enrichment targets:')
# sim.display_rlg(ks,g_rate,g_loss,g_gain)
# print('somatic generation enrichment targets:')
# sim.display_rlg(ks,s_rate,s_loss,s_gain)
"""
#loss
driver,loss_types = 'mitcp',['DEL','TRA','INS','INV']
g_total = {}
for k in g_loss:
    if g_loss[k].has_key(driver):
        v = g_loss[k][driver][k]
        for t in loss_types:
            if v.has_key(t) and v[t][0]>0:
                if g_total.has_key(k): g_total[k] += 1
                else:                  g_total[k]  = 1

s_total = {}
for k in s_loss:
    if s_loss[k].has_key(driver):
        v = s_loss[k][driver][k]
        for t in loss_types:
            if v.has_key(t) and v[t][0]>0:
                if s_total.has_key(k): s_total[k] += 1
                else:                  s_total[k]  = 1

print('somatic driver %s loss to germline loss is %s:%s'%\
      (driver,sum([s_total[k] for k in s_total]),sum([g_total[k] for k in g_total])))
"""

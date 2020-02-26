#!/usr/bin/env python
#Timothy James Becker, PhD candidate, UCONN 02/02/2020
#read simulation alignment automation for aligning reads using minimap2 and alignment using art, pbsim, lrsim
import argparse
import os
import glob
try:
    import subprocess32 as subprocess
except Exception as E:
    import subprocess
import gzip
import json
import time

des = """soMaCX: FASTQ Alignment Tool v0.1.0, 02/20/2020 Timothy James Becker"""
parser = argparse.ArgumentParser(description=des)
parser.add_argument('-f','--fastq',type=str,help='fastq input files\t[None]')
parser.add_argument('-r','--ref',type=str,help='reference fasta input file\t[None]')
parser.add_argument('-o','--out',type=str,help='reference fasta input file\t[fastq/out]')
parser.add_argument('--tools',type=str,help='minimap2,samtools,sambamba tools path\t[None]')
parser.add_argument('--threads',type=int,help='reference fasta input file\t[4]')
parser.add_argument('--platform',type=str,help='illumina,pacbio\t[illumina]')
parser.add_argument('--clean',type=bool,help='delete fastq files as you proceed to save disk\t[False]')
args = parser.parse_args()

if args.fastq is not None:
    fastq_path = args.fastq+'/'
    files =  glob.glob(fastq_path+'*.fq')+glob.glob(fastq_path+'*.fq.gz')+\
             glob.glob(fastq_path+'*.fastq')+glob.glob(fastq_path+'*.fasq.gz')
    if len(files)<1:
        print('fastq files were not found!')
        raise IOError
else:
    print('fastq file input was not specified')
    raise IOError
if args.ref is not None:
    ref_path = args.ref
    if not os.path.exists(ref_path):
        print('ref %s was not found'%ref_path)
        raise IOError
    # else:
    #     if os.path.exists(ref_path.rsplit('.fa')[-1]+'.mmi'):
    #         ref_path = ref_path.rsplit('.fa')[-1]+'.mmi' #use mmi
else:
    print('reference path was not specified!')
    raise IOError
if args.platform is not None:
    platforms = ['illumina','pacbio']
    platform = args.platform
    if platform not in platforms:
        print('%s platform is not valid, please choose on of: %s'%(platform,','.join(platforms)))
else: platform = 'illumina'
if args.tools is not None:
    tools = args.tools
    if not os.path.exists(tools+'samtools'):
        print('%s was not a valid tools path to minimap2,samtools,sambamba'%tools)
else: tools = ''
if args.out is not None:
    out_path = args.out
    if not os.path.exists(out_path): os.mkdir(out_path)
else:
    out_path = fastq_path+'/out/'
    if not os.path.exists(out_path): os.mkdir(out_path)
if args.threads is not None:  threads = args.threads
else:                         threads = 4

if __name__ == '__main__':
    if platform=='illumina':
        print('illumina platform selected...')
        #[1] gather fastq input read patterns---------------------------------------------
        fqs      = glob.glob(fastq_path+'*1.fq')+glob.glob(fastq_path+'*2.fq')
        fqgzs    = glob.glob(fastq_path+'*1.fq.gz')+glob.glob(fastq_path+'*2.fq.gz')
        fastqs   = glob.glob(fastq_path+'*1.fastq')+glob.glob(fastq_path+'*2.fastq')
        fastqgzs = glob.glob(fastq_path+'*1.fastq.gz')+glob.glob(fastq_path+'*2.fastq.gz')
        if len(fqs)>0:      ext = 'fq'
        if len(fqgzs)>0:    ext = 'fq.gz'
        if len(fastqs)>0:   ext = 'fastq'
        if len(fastqgzs)>0: ext = 'fastq.gz'
        F,fq1,fq2 = [],glob.glob(fastq_path+'*1.%s'%ext),glob.glob(fastq_path+'*2.%s'%ext)
        print('ext was: %s'%ext)
        for f1 in fq1:
            prefix1 = f1.rsplit('1.%s'%ext)[0]
            for f2 in fq2:
                prefix2 = f2.rsplit('2.%s'%ext)[0]
                if prefix1==prefix2: F += [[f1,f2]]
        print('found the following fastq pairs:\n%s'%'\n'.join([','.join(f) for f in F]))
        #[2] for each read pattern, (a) get SM (b) make RG (c) minimap2, sort, delete temp+unsorted
        x,B = 1,[]
        for f in sorted(F): #will call the first lane1 etc...
            sm = f[0].rsplit('/')[-1].rsplit('.')[0].rsplit('_')[0] #prefixes match so take pair1
            rg = r'"@RG\tID:%s'%sm+"_lane%s"%x+r"\tLB:%s"%sm+r"_lane%s"%x+r"\tPL:"+\
                 platform.upper()+r"\tPU:%s"%sm+r"_lane%s"%x+r'\tSM:%s"'%sm

            # $BIO/minimap2 -ax sr -Y -t $TH -R $RG1 $REF $DIR/lane1.1.fq.gz $DIR/lane1.2.fq.gz | $BIO/samtools view - -Shb > $DIR/$SM.lane1_N.bam
            command = [tools+'minimap2','-ax sr','-Y','-t %s'%threads,'-R %s'%rg, ref_path,
                       f[0],f[1],'|',tools+'samtools view','-','-Shb','>','%s/%s.lane%s.bam'%(out_path,sm,x)]
            try: out = subprocess.check_output(' '.join(command),shell=True)
            except Exception as E: print(command,E)

            # $BIO/samtools sort -l 9 -@ $TH -T _sort -o $DIR/$SM.lane1_N.sorted.bam $DIR/$SM.lane1_N.bam
            command = [tools+'samtools sort','-l 9','-@ %s'%threads, '-T %s'%out_path+'/_sort','-o',
                       '%s/%s.lane%s.sorted.bam'%(out_path,sm,x),'%s/%s.lane%s.bam'%(out_path,sm,x)]
            try: out = subprocess.check_output(' '.join(command),shell=True)
            except Exception as E: print(command,E)

            #rm $DIR/$SM.lane1_N.bam
            command = ['rm','%s/%s.lane%s.bam'%(out_path,sm,x)]
            try: out = subprocess.check_output(' '.join(command),shell=True)
            except Exception as E: print(command,E)
            B += ['%s/%s.lane%s.sorted.bam'%(out_path,sm,x)]
            x += 1

        #[3] merge all sorted bams, delete each bam
        # $BIO/samtools merge -l 9 -@ $TH -O BAM $DIR/$SM.N.merged.bam $DIR/$SM.lane1_N.sorted.bam $DIR/$SM.lane2_N.sorted.bam
        command = [tools+'samtools merge','-l 9','-@ %s'%threads, '-O  BAM',
                   '%s/%s.merged.bam'%(out_path,sm)]+sorted(B)
        try: out = subprocess.check_output(' '.join(command),shell=True)
        except Exception as E: print(command,E)

        #rm $DIR/$SM.lane1_N.sorted.bam $DIR/$SM.lane2_N.sorted.bam
        command = ['rm']+sorted(B)
        try: out = subprocess.check_output(' '.join(command),shell=True)
        except Exception as E: print(command,E)

        #[4] markdups on merged-sorted bam, delete merged bam
        #$BIO/sambamba markdup -l 9 -t $TH --tmpdir=_sort --sort-buffer-size=8096 --overflow-list-size=2000000 $DIR/$SM.N.merged.bam $DIR/$SM.N.final.bam
        command = [tools+'sambamba markdup','-l 9','-t %s'%threads,'--tmpdir=%s'%out_path+'/_sort',
                   '--sort-buffer-size=8096 --overflow-list-size=2000000',
                   '%s/%s.merged.bam'%(out_path,sm),'%s/%s.final.bam'%(out_path,sm)]
        try: out = subprocess.check_output(' '.join(command),shell=True)
        except Exception as E: print(command,E)
        command = ['rm','%s/%s.merged.bam'%(out_path,sm)]
        if os.path.exists(out_path+'/_sort'): command += [out_path+'/_sort']
        try: out = subprocess.check_output(' '.join(command),shell=True)
        except Exception as E: print(command,E)
    if platform=='pacbio':
        print('pacbio platform selected...')
        #[1] get the fastq reads
        fqs      = glob.glob(fastq_path+'*.fq')
        fqgzs    = glob.glob(fastq_path+'*.fq.gz')
        fastqs   = glob.glob(fastq_path+'*.fastq')
        fastqgzs = glob.glob(fastq_path+'*.fastq.gz')
        if len(fqs)>0:      ext = 'fq'
        if len(fqgzs)>0:    ext = 'fq.gz'
        if len(fastqs)>0:   ext = 'fastq'
        if len(fastqgzs)>0: ext = 'fastq.gz'
        F = glob.glob(fastq_path+'*%s'%ext)
        print(F)
        print('found the following fastq files:\n%s'+'\n'.join(F))
        x,B = 1,[]
        for f in sorted(F): #will call the first lane1 etc...
            sm = f.rsplit('/')[-1].rsplit('.')[0].rsplit('_')[0] #prefixes match so take pair1
            rg = r'"@RG\tID:%s'%sm+"_lane%s"%x+r"\tLB:%s"%sm+r"_lane%s"%x+r"\tPL:"+\
                 platform.upper()+r"\tPU:%s"%sm+r"_lane%s"%x+r'\tSM:%s"'%sm

            # $BIO/minimap2 -ax sr -Y -t $TH -R $RG1 $REF $DIR/lane1.1.fq.gz $DIR/lane1.2.fq.gz | $BIO/samtools view - -Shb > $DIR/$SM.lane1_N.bam
            command = [tools+'minimap2','-ax map-pb','-Y','-t %s'%threads,'-R %s'%rg, ref_path,
                       f,'|',tools+'samtools view','-','-Shb','>','%s/%s.lane%s.bam'%(out_path,sm,x)]
            try: out = subprocess.check_output(' '.join(command),shell=True)
            except Exception as E: print(command,E)

            # $BIO/samtools sort -l 9 -@ $TH -T _sort -o $DIR/$SM.lane1_N.sorted.bam $DIR/$SM.lane1_N.bam
            command = [tools+'samtools sort','-l 9','-@ %s'%threads, '-T %s'%out_path+'/_sort','-o',
                       '%s/%s.lane%s.sorted.bam'%(out_path,sm,x),'%s/%s.lane%s.bam'%(out_path,sm,x)]
            try: out = subprocess.check_output(' '.join(command),shell=True)
            except Exception as E: print(command,E)

            #rm $DIR/$SM.lane1_N.bam
            command = ['rm','%s/%s.lane%s.bam'%(out_path,sm,x)]
            try: out = subprocess.check_output(' '.join(command),shell=True)
            except Exception as E: print(command,E)
            B += ['%s/%s.lane%s.sorted.bam'%(out_path,sm,x)]
            x += 1
        #[3] merge all sorted bams, delete each bam
        # $BIO/samtools merge -l 9 -@ $TH -O BAM $DIR/$SM.N.merged.bam $DIR/$SM.lane1_N.sorted.bam $DIR/$SM.lane2_N.sorted.bam
        command = [tools+'samtools merge','-l 9','-@ %s'%threads, '-O BAM',
                   '%s/%s.merged.bam'%(out_path,sm)]+sorted(B)
        try: out = subprocess.check_output(' '.join(command),shell=True)
        except Exception as E: print(command,E)

        #rm $DIR/$SM.lane1_N.sorted.bam $DIR/$SM.lane2_N.sorted.bam
        command = ['rm']+sorted(B)
        try: out = subprocess.check_output(' '.join(command),shell=True)
        except Exception as E: print(command,E)

        #[4] markdups on merged-sorted bam, delete merged bam
        #$BIO/sambamba markdup -l 9 -t $TH --tmpdir=_sort --sort-buffer-size=8096 --overflow-list-size=2000000 $DIR/$SM.N.merged.bam $DIR/$SM.N.final.bam
        command = [tools+'sambamba markdup','-l 9','-t %s'%threads,'--tmpdir=%s'%out_path+'/_sort',
                   '--sort-buffer-size=8096 --overflow-list-size=2000000',
                   '%s/%s.merged.bam'%(out_path,sm),'%s/%s.final.bam'%(out_path,sm)]
        try: out = subprocess.check_output(' '.join(command),shell=True)
        except Exception as E: print(command,E)
        command = ['rm -rf','%s/%s.merged.bam'%(out_path,sm)]
        if os.path.exists(out_path+'/_sort'): command += [out_path+'/_sort']
        try: out = subprocess.check_output(' '.join(command),shell=True)
        except Exception as E: print(command,E)



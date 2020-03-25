#Timothy James Becker, PhD candidate, UCONN 01/10/2017-03/20/2020
#cython edit distance algorithms, string processing utils, weighted_random distribution
#c imports
cimport cython
cimport numpy as np
#regular imports
import math
import random
import itertools as it
import numpy as np
import pysam
__version__ = '0.1.1'
#uses the boundry in A to rescale the values in R
#destructively edits values in R as a result
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def rescale_scan(double [:,:] A, double [:] R): 
    cdef int j1,j2,n1,n2
    j1,j2,n1,n2 = 0,0,len(A),len(R)
    while j1<n1 and j2<n2:
        if R[j2]>=A[j1][0] and R[j2]<A[j1][1]:
            R[j2] = (R[j2]-A[j1][0])/(A[j1][1]-A[j1][0])*(A[j1][3]-A[j1][2])+A[j1][2]
            j2 += 1
        if R[j2]>=A[j1][1]:
            j1 += 1
        if j1>=n1: j1,j2 = n1,j2+1
        if j2>=n2: j2,j1 = n2,j1+1
        
#weighted algorithm for optimal calculation speed
#[1] from left to right multiple the range of each pos by its weight and add making projected line
#[2] sample n values from the projected line (much faster)
#[3] scale each of the n values back using the mapping in A
#w_pos = [[pos_0,pos_1],[],[],[]] #sorted by pos[i][0], len(pos) == len(w), pos[i][1]>=pos[i][0] for all i
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def weighted_random(long [:,:] pos, double [:] w, long n,
                    bint normalize=False,bint sort=False): #get a view of numpy w_pos
    cdef int i,j,k
    cdef double a,b
    cdef double [:,:] A = np.zeros((len(pos),4), dtype=np.double)
    cdef np.ndarray[double, ndim=1] R = np.zeros([n,],dtype=np.double)
    a,b,k = 0.0,0.0,len(pos)
    if len(pos)==len(w) and k>0 and len(pos[0])>0:
        if normalize:
            a = 0.0 #set weights to sum to 1.0
            for i in range(k): a += w[i]
            for i in range(k): w[i] /= a
        a,b = pos[0][0],pos[0][0]
        for i in range(k):
            A[i][0] = b
            A[i][2], A[i][3] = <double>pos[i][0],<double>pos[i][1]            
            A[i][1] = b+(A[i][3]-A[i][2])*w[i]
            b = A[i][1]
        R[:] = np.random.uniform(low=a,high=b,size=n)
        if sort:
            R = np.sort(R)
            rescale_scan(A,R)
            np.random.shuffle(R)
        else:
            for j in range(n):
                for i in range(k):
                    if R[j]>=A[i][0] and R[j]<A[i][1]:
                        R[j] = (R[j]-A[i][0])/(A[i][1]-A[i][0])*(A[i][3]-A[i][2])+A[i][2]
                        break
    return np.array(R,dtype=long)

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def weighted_random_alt(long [:,:] pos, double [:] w, long n,
                        bint normalize=True, inner_dist='uniform'):
    cdef long i,j
    cdef double a,x
    cdef np.ndarray[long, ndim=1] A = np.zeros((n,),dtype=long)
    cdef np.ndarray[long, ndim=1] R = np.zeros((n,),  dtype=long)
    k = len(pos)
    if len(pos)==len(w) and k>0 and len(pos[0])>0:
        if normalize:
            a = 0.0 #set weights to sum to 1.0
            for i in range(k): a += w[i]
            for i in range(k): w[i] /= a
        A[:] = np.random.choice(range(k),n,p=w)
        if inner_dist=='uniform':
            for i in range(n):
                R[i] = <long>np.random.uniform(pos[A[i]][0],pos[A[i]][1],1)
        elif inner_dist=='triangular':
            for i in range(n):
                R[i] = <long>np.random.triangular(pos[A[i]][0],pos[A[i]][1],1)
        #can do more distribution here too...
    return R
    
#this one is without a GIL and assumes that |c1|>=|c2|
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def edit_dist(const unsigned char[::1] c1, const unsigned char[::1] c2,
              unsigned int w_mat, unsigned int w_ins, unsigned int w_del, unsigned int w_sub):
    cdef unsigned int i,j,k,u,v,x,y,z
    u,v,k = len(c1),len(c2),2
    cdef np.ndarray[unsigned int, ndim=2] D = np.zeros([u+1,k], dtype=np.uint32)
    for i in range(u+1):   D[i][0] = i
    for j in range(v%k+1): D[0][j] = j
    for j in range(1,v+1):
        for i in range(1,u+1):
            if c1[i-1] == c2[j-1]:    
                D[i][j%k] = D[i-1][(j-1)%k]+w_mat  #matching
            else:                                  #mismatch, del, ins, sub
                x,y,z = D[i-1][j%k]+w_del,D[i][(j-1)%k]+w_ins,D[i-1][(j-1)%k]+w_sub
                if x<=y and x<=z:   D[i][j%k] = x
                elif y<=x and y<=z: D[i][j%k] = y
                else:               D[i][j%k] = z
    return D[u][v%k]

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def edit_dist_str(str s1, str s2, list w):
    cdef unsigned int i,j,k,u,v
    u,v,k = len(s1),len(s2),2
    if u<v: u,v,s1,s2 = v,u,s2,s1 #flip to longest of the two
    cdef np.ndarray[unsigned int, ndim=2] D = np.zeros([u+1,k], dtype=np.uint32)
    for i in range(u+1):   D[i][0] = i
    for j in range(v%k+1): D[0][j] = j
    for j in range(1,v+1):
        for i in range(1,u+1):
            if s1[i-1] == s2[j-1]:    
                D[i][j%k] = D[i-1][(j-1)%k] #matching
            else:                           #mismatch, del, ins, sub 
                D[i][j%k] = min(D[i-1][j%k]+w[0],D[i][(j-1)%k]+w[1],D[i-1][(j-1)%k]+w[2])
    return D[u][v%k]

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def reverse_complement(unsigned char[::1] s1):
    cdef unsigned int i,j,k,n
    #[1] build an byte map out of a C array:
    #[2] fill the complement values for lower and upper case conversion of a,A,c,C,g,G,t,T,N
    return 0

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def write_fasta(dict seqs, str fasta_path, str mode='w',bint index=False, bint gz=True, int linewidth=80):
    cdef int i,y
    cdef str s,k
    if not gz:
        with open(fasta_path,mode) as fasta:
            for k in seqs:
                s = '\n'.join(['>%s'%k]+[seqs[k][i:(i+linewidth)] for i in range(0,len(seqs[k]),linewidth)])+'\n'
                fasta.write(s)
        if index: pysam.faidx(fasta_path) #reindex
    else:
        if mode=='w' or mode=='a': mode = mode+'b'
        else: return False
        if index:
            with pysam.BGZFile(fasta_path,mode=mode,index=fasta_path+'i') as fasta:
                for k in seqs:
                    s = '\n'.join(['>%s'%k]+[seqs[k][i:(i+linewidth)] for i in range(0,len(seqs[k]),linewidth)])+'\n'
                    fasta.write(s)
            pysam.faidx(fasta_path)
        else:
            with pysam.BGZFile(fasta_path,mode=mode) as fasta:
                for k in seqs:
                    s = '\n'.join(['>%s'%k]+[seqs[k][i:(i+linewidth)] for i in range(0,len(seqs[k]),linewidth)])+'\n'
                    fasta.write(s)
    return True 

#given a fastq with a certain number of reads, subsample a certain number of them
#paire=True implies your fastq_in_path will have the 1 or 2 file name ending
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def sample_fastq(str fastq_in_path, str fastq_out_path, double p, bint paired=False):
    cdef int x,y
    x,y = 0,0
    if not paired:
        with open(fastq_in_path,'r') as q_in_1:
            with open(fastq_out_path,'w') as q_out_1:
                for r1l1,r1l2,r1l3,r1l4 in it.izip_longest(*[q_in_1]*4):
                    y += 1
                    if random.random()<=p:
                        x += 1
                        q_out_1.write(''.join([r1l1,r1l2,r1l3,r1l4]))
    else:
        with open(fastq_in_path+'1.fq','r') as q_in_1:
            with open(fastq_in_path+'2.fq','r') as q_in_2:
                with open(fastq_out_path+'1.fq','w') as q_out_1:
                    with open(fastq_out_path+'2.fq','w') as q_out_2:
                        for r1l1,r2l1,r1l2,r2l2,r1l3,r2l3,r1l4,r2l4 in it.izip(*[q_in_1,q_in_2]*4):
                            y += 1
                            if random.random()<=p:
                                x += 1
                                q_out_1.write(''.join([r1l1,r1l2,r1l3,r1l4]))
                                q_out_2.write(''.join([r2l1,r2l2,r2l3,r2l4]))
    return x,y

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def get_reverse_complement(str seq):
    cdef dict m = {'A':'T','a':'t','T':'A','t':'a','C':'G','c':'g','G':'C','g':'c','N':'N'}
    return ''.join([m[k] for k in seq[::-1]])

@cython.boundscheck(False)
@cython.nonecheck(False)
def merge_1D(list C):
    cdef long i,j,b,n
    cdef list M = []
    n = len(C)
    if n > 0:
        i,W = 0,[]
        while i < n-1:
            j = i+1
            b = C[i][1]
            while b+1 >= C[j][0] and j < n-1:
                if b < C[j][1]: b = C[j][1]
                j += 1
            M += [[C[i][0],b]]
            i = j                              #----------span of i to j here-------------
        if len(M)>0 and M[-1][1]+1>=C[i][0]:   #----------potential span of i-1 to i here-
            if M[-1][1]<C[i][1]:
                M[-1][1] = C[i][1]
        else:                                  #------------only i is here----------------
            M += [[C[i][0],C[i][1]]]
    return M

@cython.boundscheck(False)
@cython.nonecheck(False)
def LRF_1D(list C1, list C2):
    cdef long n1,n2,j1,j2,upper
    cdef list I,U,D1,D2
    j1,j2,upper = 0,0,0         #initializations and padding
    n1,n2 = len(C1)+1,len(C2)+1 #boundries here
    I,U,  = [[-2,-2]],[[-2,-2]]
    D1,D2 = [[-2,-2]],[[-2,-2]]
    if n1 > 1 and n2 > 1:
        upper = max(C1[-1][1],C2[-1][1])
        C1 += [[upper+2,upper+2],[upper+4,upper+4]] #pad out the end of C1
        C2 += [[upper+2,upper+2],[upper+4,upper+4]] #pad out the end of C2
        while j1+j2 < n1+n2:  #pivioting dual ordinal indecies scan left to right on C1, C2
            a = C1[j1][0]-C2[j2][0]
            b = C1[j1][0]-C2[j2][1]
            c = C1[j1][1]-C2[j2][0]
            d = C1[j1][1]-C2[j2][1]
            if    C1[j1][0:2]==C2[j2][0:2]: #[7] c1 and c2 are equal on x
                orientation_7(C1,j1,C2,j2,I,U,D1,D2)
                j1 += 1
                j2 += 1 
            elif  c<0:                      #[1] c1 disjoint of left of c2
                orientation_1(C1,j1,C2,j2,I,U,D1,D2)                   
                j1 += 1    
            elif  b>0:                      #[6] c1 disjoint right of c2
                orientation_6(C1,j1,C2,j2,I,U,D1,D2)             
                j2 += 1 
            elif  a<0 and d<0:              #[2] c1 right overlap to c2 left no envelopment
                orientation_2(C1,j1,C2,j2,I,U,D1,D2)           
                j1 += 1 
            elif  a>0 and d>0:              #[4] c1 left overlap of c2 right no envelopment
                orientation_4(C1,j1,C2,j2,I,U,D1,D2) 
                j2 += 1 
            elif  a<=0 and d>=0:            #[3] c1 envelopment of c2
                orientation_3(C1,j1,C2,j2,I,U,D1,D2)
                j2 += 1 
            elif  a>=0 and d<=0:            #[5] c1 enveloped by c2
                orientation_5(C1,j1,C2,j2,I,U,D1,D2)
                j1 += 1 
            if j1>=n1: j1,j2 = n1,j2+1 #sticky indecies wait for eachother
            if j2>=n2: j2,j1 = n2,j1+1 #sticky indecies wait for eachother
        #pop off extras for each features (at most two at the end)
        while len(C1) > 0 and C1[-1][0] > upper:  C1.pop()    
        while len(C2) > 0 and C2[-1][0] > upper:  C2.pop()
        while len(I)  > 0 and I[-1][0]>upper:      I.pop()
        if len(I) > 0 and I[-1][1]>upper:         I[-1][1] = upper
        while len(U) > 0 and U[-1][0]>upper:      U.pop()
        if len(U) > 0 and U[-1][1]>upper:         U[-1][1] = upper
        while len(D1) > 0 and D1[-1][0]>upper:    D1.pop()
        if len(D1) > 0 and D1[-1][1]>upper:       D1[-1][1] = min(C2[-1][0]-1,C1[-1][1])
        while len(D2) > 0 and D2[-1][0]>upper:    D2.pop()
        if len(D2) > 0 and D2[-1][1]>upper:       D2[-1][1] = min(C1[-1][0]-1,C2[-1][1])
    else:
        if   n1==1:  
            if n2>1: U,D2 = U+C2,D2+C2
        elif n2==1:
            if n1>1: U,D1 = U+C1,D1+C1 
    return I[1:],U[1:],D1[1:],D2[1:]     
    
#[1] c1 disjoint of left of c2
@cython.boundscheck(False)
@cython.nonecheck(False)
cdef void orientation_1(list C1,long j1,list C2,long j2,list I,list U,list D1,list D2):
    #[i]----------------------[i] #no intersection
    #[u]----------------------[u]
    if U[-1][1]+1 >= C1[j1][0]:
        U[-1][1] = C1[j1][1]
    else:                         
        U += [[C1[j1][0],C1[j1][1]]]
    #[d1]--------------------[d1]
    if D1[-1][1]+1 >= C1[j1][0]:  #extend segment
        if D1[-1][1]+1!=C2[j2-1][0]:
            D1[-1][1] = C1[j1][1]
    else:                         #new segment                        
        D1 += [[C1[j1][0],C1[j1][1]]]            
    #[d2]--------------------[d2] #no set two difference

#[2] c1 right overlap to c2 left no envelopment
@cython.boundscheck(False)
@cython.nonecheck(False)
cdef void orientation_2(list C1,long j1,list C2,long j2,list I,list U,list D1,list D2):
    #[i]----------------------[i]
    if I[-1][1]+1 >= C2[j2][0]: 
        I[-1][1] = C1[j1][1] #was C2[j2][1]
    else:                       
        I += [[C2[j2][0],C1[j1][1]]]                 
    #[u]----------------------[u]
    if U[-1][1]+1 >= C1[j1][0]: 
        U[-1][1] = C2[j2][1]
    else:                       
        U += [[C1[j1][0],C2[j2][1]]]
    #[d1]--------------------[d1]
    if D1[-1][1]+1 >= C1[j1][0]: 
        D1[-1][1] = C1[j1][1]
        if D1[-1][1] > C2[j2][0]-1:
            D1[-1][1] = C2[j2][0]-1
            if D1[-1][1] < D1[-1][0]: D1.pop()
    else:
        D1 += [[C1[j1][0],C2[j2][0]-1]]            
    #[d2]--------------------[d2] 
    D2 += [[C1[j1][1]+1,C2[j2][1]]]

#[3] c1 envelopment of c2
@cython.boundscheck(False)
@cython.nonecheck(False)
cdef void orientation_3(list C1,long j1,list C2,long j2,list I,list U,list D1,list D2):
    #[i]----------------------[i]
    if I[-1][1]+1 >= C2[j2][0]: 
        I[-1][1] = C2[j2][1]          
    else:                       
        I += [[C2[j2][0],C2[j2][1]]]                        
    #[u]----------------------[u]
    if U[-1][1]+1 >= C1[j1][0]: 
        U[-1][1] = C1[j1][1]
    else:                       
        U += [[C1[j1][0],C1[j1][1]]]
    #[d1]--------------------[d1] 
    if D1[-1][1]+1 >= C1[j1][0]:
        D1[-1][1] = C2[j2][0]-1
        if C2[j2][1] < C1[j1][1]:
            D1 += [[C2[j2][1]+1,C1[j1][1]]]
    elif D1[-1][1] >= C2[j2][0]:
        D1[-1][1] = C2[j2][0]-1
        if C2[j2][1] < C1[j1][1]:  #has a right side
            D1 += [[C2[j2][1]+1,C1[j1][1]]]
    else:
        if C1[j1][0] < C2[j2][0]:  #has a left side
            D1 += [[C1[j1][0],C2[j2][0]-1]]
        if C2[j2][1] < C1[j1][1]:  #has a right side
            D1 += [[C2[j2][1]+1,C1[j1][1]]]


#[4] c1 left overlap of c2 right no envelopment            
@cython.boundscheck(False)
@cython.nonecheck(False)
cdef void orientation_4(list C1,long j1,list C2,long j2,list I,list U,list D1,list D2):
    #[i]----------------------[i]
    if I[-1][1]+1 >= C1[j1][0]: 
        I[-1][1] = C2[j2][1]        
    else:                       
        I += [[C1[j1][0],C2[j2][1]]]           
    #[u]----------------------[u]
    if U[-1][1]+1 >= C2[j2][0]: 
        U[-1][1] = C1[j1][1]
    else:                       
        U += [[C2[j2][0],C1[j1][1]]]       
    #[d1]--------------------[d1]
    D1 += [[C2[j2][1]+1,C1[j1][1]]]            
    #[d2]--------------------[d2] 
    if D2[-1][1]+1 >= C2[j2][0]: 
        D2[-1][1] = C2[j2][1]
        if D2[-1][1] > C1[j1][0]-1:
            D2[-1][1] = C1[j1][0]-1
            if D2[-1][1] < D2[-1][0]: D2.pop()
    else:                        
        D2 += [[C2[j2][0],C1[j1][0]-1]]


#[5] c1 enveloped by c2
@cython.boundscheck(False)
@cython.nonecheck(False)
cdef void orientation_5(list C1,long j1,list C2,long j2,list I,list U,list D1,list D2):
    #[i]----------------------[i]
    if I[-1][1]+1 >= C1[j1][0]: 
        I[-1][1] = C1[j1][1]          
    else:                       
        I += [[C1[j1][0],C1[j1][1]]]
                                      
    #[u]----------------------[u]
    if U[-1][1]+1 >= C2[j2][0]: 
        U[-1][1] = C2[j2][1]
    else:                       
        U += [[C2[j2][0],C2[j2][1]]]
    #[d1]--------------------[d1] #no set one difference
    #[d2]--------------------[d2]
    if D2[-1][1]+1 >= C2[j2][0]: 
        D2[-1][1] = C1[j1][0]-1
        if C1[j1][1] < C2[j2][1]:
            D2 += [[C1[j1][1]+1,C2[j2][1]]]
    elif D2[-1][1] >= C1[j1][0]:
        D2[-1][1] = C1[j1][0]-1
        if C1[j1][1] < C2[j2][1]:  #has a right side
            D2 += [[C1[j1][1]+1,C2[j2][1]]]
    else:
        if C2[j2][0] < C1[j1][0]:  #has a left side
            D2 += [[C2[j2][0],C1[j1][0]-1]]
        if C1[j1][1] < C2[j2][1]:  #has a right side
            D2 += [[C1[j1][1]+1,C2[j2][1]]]

#[6] c1 disjoint right of c2
@cython.boundscheck(False)
@cython.nonecheck(False)
cdef void orientation_6(list C1,long j1,list C2,long j2,list I,list U,list D1,list D2):
    #[i]----------------------[i] #no instersection
    if U[-1][1]+1 >= C2[j2][0]:   
        U[-1][1] = C2[j2][1]
    else:                         
        U += [[C2[j2][0],C2[j2][1]]]
    #[d1]--------------------[d1] #no set one difference
    #[d2]--------------------[d2] 
    if D2[-1][1]+1 >= C2[j2][0]:
        if D2[-1][1]+1!=C1[j1-1][0]:
            D2[-1][1] = C2[j2][1]
    else:                        
        D2 += [[C2[j2][0],C2[j2][1]]] 

#[7] c1 and c2 are equal on x
@cython.boundscheck(False)
@cython.nonecheck(False)
cdef void orientation_7(list C1,long j1,list C2,long j2,list I,list U,list D1,list D2):
    #[i]----------------------[i]
    if I[-1][1]+1 >= C1[j1][0]:
        I[-1][1] = C1[j1][1]
    else:
        I += [[C1[j1][0],C1[j1][1]]]
               
    #[u]----------------------[u]
    if U[-1][1]+1 >= C1[j1][0]:
        U[-1][1] = C1[j1][1]
    else:
        U += [[C1[j1][0],C1[j1][1]]]
    #[d1]----------------------[d1]
    #[d2]----------------------[d2]    
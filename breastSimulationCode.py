#! /usr/bin/python

######################################################################################
## A Python script to simulate 3D peripherally dominated tumor growth and 			##		
## multi-region sequencing data via an agent-based model. Deme subdivision is       ##
## assumed in order to model cell mixing and spatial contraint. 20 samples across   ##
## all octants of the deme are recorded to profile simulated spatial tumor          ##
## heterogeneity. This version assumes at most two drivers occur on the same        ##
## lineage; the fitness of Tier1 and Tier2 driver lineages is 1+s and (1+s)^2,      ##
## respectively.                 													##
## 																					##
## Author for breast-specific modifications: Katherine McNamara (Curtis Lab at Stanford)##                                           ##
## Original author: Zheng Hu (Curtis Lab at Stanford)                               ##                                                      ##
######################################################################################


import sys,os,math,random
import numpy as np
from collections import Counter
import sets

class deme():
    def __init__(self):
        self.present= 0         ## whether the deme is empty or occupied: 0-empty;1-occupied
        self.neutral = []    ## the neutral founder lineage after tumor tranformation
        self.advant1 = []       ## the advantageous cells having one driver mutation
        self.advant2 = []       ## the advantageous cells having two driver mutations

def createLattice(d):
    """
    Create a 3D cubic lattice with side length of 2d+1 where each site contains a empty deme.
    """
    lattice = {}
    for x in range(0,2*d+1):
        for y in range(0,2*d+1):
            for z in range(0,2*d+1):
                lattice[(x,y,z)] = deme()
    return lattice


def neighbor26((a,b,c)):
    """
    Moore neighbourhood: 26 neighbour sites of (a,b,c)
    """
    neighbor = [(a+i,b+j,c+k)
                for i in [-1,0,1]
                for j in [-1,0,1]
                for k in [-1,0,1]
                if not (i==0 and j==0 and k==0)]
    
    return neighbor


#def neighbor26((a,b,c)):
#    """
#    Moore neighbourhood: 26 neighbour sites of (a,b,c).
#    """
#    neighbor = [(a-1, b-1, c-1),(a-1, b-1, c),(a-1, b-1, c+1),(a-1, b, c-1),(a-1, b, c),(a-1, b, c+1),(a-1, b+1, c-1),(a-1, b+1, c),(a-1, b+1, c+1),(a, b-1, c-1),(a, b-1, c),(a, b-1, c+1),(a, b, c-1),(a, b, c+1),(a, b+1, c-1),(a, b+1, c),(a, b+1, c+1),(a+1, b-1, c-1),(a+1, b-1, c),(a+1, b-1, c+1),(a+1, b, c-1),(a+1, b, c),(a+1, b, c+1),(a+1, b+1, c-1),(a+1, b+1, c),(a+1, b+1, c+1)]
#    return neighbor


def neighbor6((a,b,c)):
    """
    von Neumann neighbourhood: 6 neighbour sites of (a,b,c).
    """
    neighbor = [(a-1, b, c),(a+1, b, c),(a, b-1, c),(a, b+1, c),(a, b, c-1),(a, b, c+1)]
    return neighbor


def localNeighbor((a,b,c),r):
    """
    A function to search the local neighbour sites of (a,b,c) within an area of radius r in the 3D cubic lattice.
    """
    neighbor = []
    for x in range(-r,r+1):
        for y in range(-r,r+1):
            for z in range(-r,r+1):
                if pow(x,2)+pow(y,2)+pow(z,2) < pow(r+1,2):
                    neighbor += [(a+x,b+y,c+z)]
    return neighbor


def traceLineage(mlineage,mutid):
    """
    A function to obtain the mutational lineage of a cell from the mutation id of the most recently occurred mutation in the cell. 
    For example, the input ID (most recently occurred mutation) of target cell is "100" and the output is "1-12-35-56-100", which is the mutation lineage of the cell
    
    mlineage - the list that could be used to recover the mutational lineage given the most recent mutation id of a lineage
    mutid - the mutation ID of the most recently occurred mutation in the cell
    """
    recent_muts = mutid.split(',')  # it is possible that multiple mutations occur during in a cell division. For instance, the mutation id of most recently occurred mutations is "100,101"
    recent_muts = [int(t) for t in recent_muts]
    first_mut = recent_muts[0]      # the first mutation id in a multi-mutation event
    trace = []
    while first_mut > 0:
        trace += recent_muts
        recent_muts = mlineage[first_mut].split(',')
        recent_muts = [int(t) for t in recent_muts]
        first_mut = recent_muts[0]
    return trace

    
def lowerORupper(value):
    """
    A function to choose the upper or lower integral value given a non-integral number
    """
    lower_int = int(value)
    upper_int = lower_int+1
    if random.random() < value-lower_int:
        return upper_int
    else:
        return lower_int


def initiateFirstDeme_v1(maxsize,lineage,current_id,sfit):
    """
    The growth of the initial deme from a single transformed tumor cell via a random discrete-time birth-death process
    v1 - one-tier driver model
    maxsize - size limit of a deme
    lineage - a list that stores the lineage information of mutations
    current_id - the starting mutation ID
    sfit - selection fitness of advantageous mutations
    """
    neu_list = [str(current_id)]
    adv_list = []
    current_deme_size = 1
    while current_deme_size < maxsize:
        n1,n2 = len(neu_list),len(adv_list)                         #n1 and n2 are the current number of neutral founder cells and advantageous cells, respectively
        neu_divcells =  int(n1*birth_rate+1)                        #number of dividing cells of neutral lineage in this generation. The other cells will die in the next generation
        neu_list = random.sample(neu_list,neu_divcells)*2
        if n2 > 0:
            adv_divcells = lowerORupper(n2*birth_rate*(1+sfit))   #number of dividing cells of advantageous lineage in this generation        
            adv_list = random.sample(adv_list,adv_divcells)*2
        n1,n2 = len(neu_list),len(adv_list)
        current_deme_size = n1+n2
        if n1 > 0:
            new_mut1 = np.random.poisson(mut_rate*n1)               # the total number of mutations occurring in a generation follows Poission distribution with lambda=u*n
            mut_assig1 = Counter(np.random.choice(n1,new_mut1))
            for x1 in mut_assig1.keys():
                nmut = mut_assig1[x1]
                new_mut1 = range(current_id+1,current_id+1+nmut)
                mut_str = ",".join(map(str,new_mut1))
                #if nmut > 1:
                #    for t in new_mut1:
                #        multi_events[str(t)] = mut_str
                for xn in range(0,nmut):
                    current_id += 1
                    lineage += [neu_list[x1]]
                neu_list[x1] = mut_str
        if n2 > 0:
            new_mut2 = np.random.poisson(mut_rate*n2)
            mut_assig2 = Counter(np.random.choice(n2,new_mut2))
            for x2 in mut_assig2.keys():
                nmut = mut_assig2[x2]
                new_mut2 = range(current_id+1,current_id+1+nmut)
                mut_str = ",".join(map(str,new_mut2))
                #if nmut > 1:
                #    for t in new_mut2:
                #        multi_events[str(t)] = mut_str
                for xn in range(0,nmut):
                    current_id += 1
                    lineage += [adv_list[x2]]
                adv_list[x2] = mut_str
        
        if random.random() < adv_rate*n1:                           # occurence of advantageous mutation on the neutral lineage
            current_id += 1
            current_n1 = len(neu_list)
            lineage += [str(neu_list[current_n1-1])]
            adv_list += [str(current_id)]
            neu_list = neu_list[0:current_n1-1]
    
    return neu_list,adv_list,current_id,lineage


def initiateFirstDeme_v2(maxsize,lineage,current_id,sfit):
    """
    The growth of the initial deme from a single transformed tumor cell via a random discrete-time birth-death process
    v2 - two-tier driver model
    maxsize - size limit of a deme
    lineage - a list that stores the lineage information of mutations
    current_id - the starting mutation ID
    sfit - selection fitness of advantageous mutations
    """
    neu_list = [str(current_id)]
    adv1_list = []
    adv2_list = []
    current_deme_size = 1
    while current_deme_size < maxsize:
        n1,n2,n3 = len(neu_list),len(adv1_list),len(adv2_list)
        neu_divcells =  int(n1*birth_rate+1)                            #number of dividing cells in this generation
        neu_list = random.sample(neu_list,neu_divcells)*2
        if n2 > 0:
            adv1_divcells = lowerORupper(n2*birth_rate*(1+sfit))        #number of dividing cells in this generation        
            adv1_list = random.sample(adv1_list,adv1_divcells)*2
        if n3 > 0:
            adv2_divcells = lowerORupper(n3*birth_rate*pow(1+sfit,2))   #number of dividing cells in this generation        
            adv2_list = random.sample(adv2_list,adv2_divcells)*2
        n1,n2,n3 = len(neu_list),len(adv1_list),len(adv2_list)
        current_deme_size = n1+n2+n3
        if n1 > 0:
            new_mut1 = np.random.poisson(mut_rate*n1)
            mut_assig1 = Counter(np.random.choice(n1,new_mut1))
            for x1 in mut_assig1.keys():
                nmut = mut_assig1[x1]
                new_mut1 = range(current_id+1,current_id+1+nmut)
                mut_str = ",".join(map(str,new_mut1))
                #if nmut > 1:
                #    for t in new_mut1:
                #        multi_events[str(t)] = mut_str
                for xn in range(0,nmut):
                    current_id += 1
                    lineage += [neu_list[x1]]
                neu_list[x1] = mut_str
        if n2 > 0:
            new_mut2 = np.random.poisson(mut_rate*n2)
            mut_assig2 = Counter(np.random.choice(n2,new_mut2))
            for x2 in mut_assig2.keys():
                nmut = mut_assig2[x2]
                new_mut2 = range(current_id+1,current_id+1+nmut)
                mut_str = ",".join(map(str,new_mut2))
                #if nmut > 1:
                #    for t in new_mut2:
                #        multi_events[str(t)] = mut_str
                for xn in range(0,nmut):
                    current_id += 1
                    lineage += [adv1_list[x2]]
                adv1_list[x2] = mut_str
        if n3 > 0:
            new_mut3 = np.random.poisson(mut_rate*n3)
            mut_assig3 = Counter(np.random.choice(n3,new_mut3))
            for x3 in mut_assig3.keys():
                nmut = mut_assig3[x3]
                new_mut3 = range(current_id+1,current_id+1+nmut)
                mut_str = ",".join(map(str,new_mut3))
                #if nmut > 1:
                #    for t in new_mut2:
                #        multi_events[str(t)] = mut_str
                for xn in range(0,nmut):
                    current_id += 1
                    lineage += [adv2_list[x3]]
                adv2_list[x3] = mut_str
        
        if random.random() < adv_rate*n1:
            current_id += 1
            current_n1 = len(neu_list)
            lineage += [str(neu_list[current_n1-1])]
            adv1_list += [str(current_id)]
            neu_list = neu_list[0:current_n1-1]
        
        if random.random() < adv_rate*n2:
            current_id += 1
            current_n2 = len(adv1_list)
            lineage += [str(adv1_list[current_n2-1])]
            adv2_list += [str(current_id)]
            adv1_list = adv1_list[0:current_n2-1]
    
    
    return neu_list,adv1_list,adv2_list,current_id,lineage


def demeGrowthFission_v1(neu_list,adv_list,lineage,current_id,current_deme_number,sfit):
    """
    A function to simulate deme expansion and fission and keep track of the mutational lineages
    v1 - one-tier driver model
    """
    current_deme_size = len(neu_list)+len(adv_list)
    while current_deme_size < 2*deme_size:                          #when the deme size doubles, it will split into two offspring demes
        n1,n2 = len(neu_list),len(adv_list)
        neu_divcells =  lowerORupper(n1*birth_rate)                 #number of dividing cells in this generation
        neu_list = random.sample(neu_list,neu_divcells)*2
        if n2 > 0:
            adv_divcells =  lowerORupper(n2*birth_rate*(1+sfit))  #number of dividing cells in this generation
            adv_list = random.sample(adv_list,adv_divcells)*2
        n1,n2 = len(neu_list),len(adv_list)
        current_deme_size = n1+n2
        if current_deme_number < 5*pow(10,7)/deme_size:             #stop mutation occurring when the tumor size is larger than 5*10^7 cells. The reason is that late occuring mutations have very small chance to present at detectable frequency even under selection.
            if n1 > 0:
                new_mut1 = np.random.poisson(mut_rate*n1)
                mut_assig1 = Counter(np.random.choice(n1,new_mut1))
                for x1 in mut_assig1.keys():
                    nmut = mut_assig1[x1]
                    new_mut1 = range(current_id+1,current_id+1+nmut)
                    mut_str = ",".join(map(str,new_mut1))
                    #if nmut > 1:
                    #    for t in new_mut1:
                    #        multi_events[str(t)] = mut_str
                    for xn in range(0,nmut):
                        current_id += 1
                        lineage += [neu_list[x1]]
                    neu_list[x1] = mut_str
            if n2 > 0:
                new_mut2 = np.random.poisson(mut_rate*n2)
                mut_assig2 = Counter(np.random.choice(n2,new_mut2))
                for x2 in mut_assig2.keys():
                    nmut = mut_assig2[x2]
                    new_mut2 = range(current_id+1,current_id+1+nmut)
                    mut_str = ",".join(map(str,new_mut2))
                    #if nmut > 1:
                    #    for t in new_mut2:
                    #        multi_events[str(t)] = mut_str
                    for xn in range(0,nmut):
                        current_id += 1
                        lineage += [adv_list[x2]]
                    adv_list[x2] = mut_str
            
            if random.random() < adv_rate*n1:
                current_id += 1
                current_n1 = len(neu_list)
                lineage += [str(neu_list[current_n1-1])]
                adv_list += [str(current_id)]
                neu_list = neu_list[0:current_n1-1]
            #n1,n2 = len(neu_list),len(adv_list)
    random.shuffle(neu_list)
    if len(neu_list) > 0:
        offspring_neu = np.random.binomial(len(neu_list),0.5)       # the offpring deme size is determined by a Binomial distribution B(n,0.5)
    else:
        offspring_neu = 0
    neu_list1=neu_list[0:offspring_neu]
    neu_list2=neu_list[offspring_neu:len(neu_list)]
    random.shuffle(adv_list)
    if len(adv_list) > 0:
        offspring_adv = np.random.binomial(len(adv_list),0.5)
    else:
        offspring_adv = 0
    adv_list1=adv_list[0:offspring_adv]
    adv_list2=adv_list[offspring_adv:len(adv_list)]
    
    return neu_list1,neu_list2,adv_list1,adv_list2,current_id,lineage


def demeGrowthFission_v2(neu_list,adv1_list,adv2_list,lineage,current_id,current_deme_number,sfit):
    """
    A function to simulate deme growth and fission and keep track of the mutational lineages
    v2 - two-tier driver model
    """
    current_deme_size = len(neu_list)+len(adv1_list)+len(adv2_list)
    while current_deme_size < 2*deme_size:
        n1,n2,n3 = len(neu_list),len(adv1_list),len(adv2_list)
        if n1 > 0:
            neu_divcells =  lowerORupper(n1*birth_rate)             #number of dividing cells in this generation
            neu_list = random.sample(neu_list,neu_divcells)*2
        if n2 > 0:
            adv1_divcells =  lowerORupper(n2*birth_rate*(1+sfit))   #number of dividing cells in this generation
            adv1_list = random.sample(adv1_list,adv1_divcells)*2
        if n3 > 0:
            adv2_divcells =  lowerORupper(n3*birth_rate*pow(1+sfit,2)) #number of dividing cells in this generation
            adv2_list = random.sample(adv2_list,adv2_divcells)*2
        
        n1,n2,n3 = len(neu_list),len(adv1_list),len(adv2_list)
        current_deme_size = n1+n2+n3
        if current_deme_number < 5*pow(10,7)/deme_size: # stop mutation occurence when the tumor size is larger than 10^4*5000 = 5*10^7
            if n1 > 0:
                new_mut1 = np.random.poisson(mut_rate*n1)
                mut_assig1 = Counter(np.random.choice(n1,new_mut1))
                for x1 in mut_assig1.keys():
                    nmut = mut_assig1[x1]
                    new_mut1 = range(current_id+1,current_id+1+nmut)
                    mut_str = ",".join(map(str,new_mut1))
                    #if nmut > 1:
                    #    for t in new_mut1:
                    #        multi_events[str(t)] = mut_str
                    for xn in range(0,nmut):
                        current_id += 1
                        lineage += [neu_list[x1]]
                    neu_list[x1] = mut_str
            if n2 > 0:
                new_mut2 = np.random.poisson(mut_rate*n2)
                mut_assig2 = Counter(np.random.choice(n2,new_mut2))
                for x2 in mut_assig2.keys():
                    nmut = mut_assig2[x2]
                    new_mut2 = range(current_id+1,current_id+1+nmut)
                    mut_str = ",".join(map(str,new_mut2))
                    #if nmut > 1:
                    #    for t in new_mut2:
                    #        multi_events[str(t)] = mut_str
                    for xn in range(0,nmut):
                        current_id += 1
                        lineage += [adv1_list[x2]]
                    adv1_list[x2] = mut_str
            if n3 > 0:
                new_mut3 = np.random.poisson(mut_rate*n3)
                mut_assig3 = Counter(np.random.choice(n3,new_mut3))
                for x3 in mut_assig3.keys():
                    nmut = mut_assig3[x3]
                    new_mut3 = range(current_id+1,current_id+1+nmut)
                    mut_str = ",".join(map(str,new_mut3))
                    #if nmut > 1:
                    #    for t in new_mut2:
                    #        multi_events[str(t)] = mut_str
                    for xn in range(0,nmut):
                        current_id += 1
                        lineage += [adv2_list[x3]]
                    adv2_list[x3] = mut_str
            
            if random.random() < adv_rate*n1:
                current_id += 1
                current_n1 = len(neu_list)
                lineage += [str(neu_list[current_n1-1])]
                adv1_list += [str(current_id)]
                neu_list = neu_list[0:current_n1-1]
            
            if random.random() < adv_rate*n2:
                current_id += 1
                current_n2 = len(adv1_list)
                lineage += [str(adv1_list[current_n2-1])]
                adv2_list += [str(current_id)]
                adv1_list = adv1_list[0:current_n2-1]
    
    random.shuffle(neu_list)
    if len(neu_list) > 0:
        offspring_neu = np.random.binomial(len(neu_list),0.5)
    else:
        offspring_neu = 0
    neu_list1=neu_list[0:offspring_neu]
    neu_list2=neu_list[offspring_neu:len(neu_list)]
    
    random.shuffle(adv1_list)
    if len(adv1_list) > 0:
        offspring_adv1 = np.random.binomial(len(adv1_list),0.5)
    else:
        offspring_adv1 = 0
    adv1_list1=adv1_list[0:offspring_adv1]
    adv1_list2=adv1_list[offspring_adv1:len(adv1_list)]
    
    random.shuffle(adv2_list)
    if len(adv2_list) > 0:
        offspring_adv2 = np.random.binomial(len(adv2_list),0.5)
    else:
        offspring_adv2 = 0
    adv2_list1=adv2_list[0:offspring_adv2]
    adv2_list2=adv2_list[offspring_adv2:len(adv2_list)]
    
    return neu_list1,neu_list2,adv1_list1,adv1_list2,adv2_list1,adv2_list2,current_id,lineage


def seqProcessing(sp,sample_keys,mlineage,size_par,mean_depth,purity):
    """
    Model the random sampling process in NGS and report the sequencing allele frequencies in a sample of cells
    
    sp- the lattice space
    sample_keys- the locations for the demes in a bulk sample
    size_par- variance parameter for negative-binomial distribution
    mean_depth- the mean depth of the sequencing
    purity- tumor purity
    """
    all_cur_id = []                                     # all most recently occurred mutations
    all_mut_id = []                                     # all mutations in the sampled cells
    for key in sample_keys:
        smuts = list(sp[key].neutral + sp[key].advant1 + sp[key].advant2)
        all_cur_id += smuts
    sample_size = 10000                                 # the number of cells for sequencing analysis
    sample_id = random.sample(all_cur_id,sample_size)
    id_count = Counter(sample_id)
    for x in id_count.keys():
        xlineage = traceLineage(mlineage,x)
        all_mut_id += xlineage*id_count[x]
    mut_count = Counter(all_mut_id)
    prob_par=size_par*1.0/(size_par+mean_depth)
    sampleAF = {}                                       # a dictionary storing the mutation IDs and corresponding depth and allele frequency the seq data
    for x in mut_count.keys():
        true_af = mut_count[x]*0.5*purity/sample_size   # the true allele frequency in the sample
        if true_af > 0.001:                             # filter mutations with very low frequency that is not detectable by ~100X sequencing depth
            site_depth = np.random.negative_binomial(size_par,prob_par)
            if site_depth >= 15:                        # seq depth cutoff for "calling" a mutation
                var_reads = np.random.binomial(site_depth,true_af)
                seq_af = var_reads*1.0/site_depth
                if var_reads >= 4:                      # variant reads cutof for "calling" a mutation
                    sampleAF[str(x)] = (site_depth,seq_af)
    return sampleAF


def highMuts(sp,position,mlineage,cutoff):
    """
    Obtain the high-frequency mutations (vaf>cutoff) in a particular deme
    
    sp - the lattice space
    position - the location of the deme
    mlineage - mutation lineage dictionary
    cutoff - the VAF cutoff for a "high-frequency" mutation, e.g. 0.4
    """
    all_cur_id = sp[position].neutral + sp[position].advant1 + sp[position].advant2
    all_mut_id = []
    sample_size = 100
    sample_id = random.sample(all_cur_id,sample_size)
    id_count = Counter(sample_id)
    for y in id_count.keys():
        xlineage = traceLineage(mlineage,y)
        all_mut_id += xlineage*id_count[y]
    mut_count = Counter(all_mut_id)
    highAF_muts = []
    for x in mut_count.keys():
        allele_freq = mut_count[x]*1.0/sample_size
        if allele_freq > cutoff:
            highAF_muts += [int(x)]
    
    return highAF_muts


def pubMutGenerator(n,size_par,mean_depth,purity):
    """
    A function to generate the public clonal mutations occured during the multi-step tumorigenesis before transformation.
    
    n- number of clonal mutations
    size_par- variation parameter in the negative binomial distribution
    mean_death- mean seq depth
    """
    prob_par=size_par*1.0/(size_par+mean_depth)
    mean_af = 0.5*purity
    depth_pub = []
    maf_pub = []
    for k in range(0,n):
        correct = 0
        while correct == 0:
            site_depth = np.random.negative_binomial(size_par,prob_par)
            if site_depth >= 15:
                correct =1
        var_reads = np.random.binomial(site_depth,mean_af)
        site_maf = var_reads*1.0/site_depth
        depth_pub += [site_depth]
        maf_pub += [site_maf]
    return depth_pub,maf_pub


def localSampling(region,sample_number,cutoff):
    """
    A function to sampling the locations of multiple bulk samples in a local region.
    """
    success = 0
    while success == 0:
        locations = random.sample(region,sample_number)
        repeat = sample_number*(sample_number-1)
        minall = 999
        for x in range(0,repeat):
            rs = random.sample(locations,2)
            min_distance = min([abs(rs[0][0]-rs[1][0]),abs(rs[0][1]-rs[1][1]),abs(rs[0][2]-rs[1][2])])
            if min_distance < minall:
                minall = min_distance
        if min_distance > 2*cutoff:
            success = 1
    return locations


def bulkTissueSampling(sp,location,radius):
    """
    A function to sampling a bulk sample in a local region.
    """
    local_region = localNeighbor(location,radius)
    bulk_tissue = []
    for x in local_region:
        if sp[x].present == 1:
            bulk_tissue += [x]
    return bulk_tissue


def lineageDashLink(mlist):
    """
    Transform the mutation lineage from list (e.g [1,3,10,20]) to dash-linked string (e.g. 1-3-10-20)
    """
    if len(mlist) > 0:
        dstring = str(mlist[0])
        for x in mlist[1:len(mlist)]:
            dstring += "-"
            dstring += str(x)
        return dstring
    else:
        return "0"


def missingDepth(vafdata,absent_muts,mean_depth):
    """
    Randomly generate the sequencing depth for the mutation-absent sites across samples
    """
    for x in absent_muts:
        done = 0
        while done == 0:
            missing_depth = np.random.negative_binomial(2,2.0/(2+mean_depth))
            if missing_depth >= 15:
                done = 1
        vafdata[str(x)] = (missing_depth,0)
    
    return vafdata



#############main script to simulate a tumor and multi-region sequencing data#########
###parameter intiation###
s_coef = float(sys.argv[1])         # the selection coefficient, ranges from 0 to 0.5
repl = int(sys.argv[2])             # replication of simulation

rd = 60                             # the side length of the 3D space, may need to be increased for very small demes
deme_size = 5000                    # the deme size, ranges from 1K to 50K
final_tumor_size = pow(10,9)        # the number of cells in the final tumor (1 billion cells)
final_deme_number = final_tumor_size/deme_size    # the final number of demes in the tumor
birth_rate = 0.55                   # the birth probability at each cell generation during tumor growth

mut_rate = 0.6                      # the neutral mutation rate at whole exonic region
adv_rate = pow(10,-5) if s_coef > 0 else 0
percentage = int(s_coef*100)        # the percentage form of the selection

npub=50                             # the number of public mutation to be generated, based on our breast cohort
seq_depth=100                       # the average sequencing depth, based on utilized mutaiton from our breast cohort

mut_id = 0
mutlineage = ['0']                  # the lineage tracer
######################################################################################

first_neu,first_adv1,first_adv2,mut_id,mutlineage = initiateFirstDeme_v2(deme_size,mutlineage,mut_id,s_coef)  #the growth of the fisrt deme from single transformed cell

print "No. of neutral, tie1 and tie2 advantageous mutations in the first deme:",len(first_neu),len(first_adv1),len(first_adv2)

space = createLattice(rd)
space[(rd,rd,rd)].present = 1
space[(rd,rd,rd)].neutral = list(first_neu)
space[(rd,rd,rd)].advant1 = list(first_adv1)
space[(rd,rd,rd)].advant2 = list(first_adv2)
current_keys = [(rd,rd,rd)]
current_deme_number = 1                                 #current deme number
surface_keys = [(rd,rd,rd)]
surface_deme_number = 1
deme_time_generation = 0

while current_deme_number < final_deme_number:
    new_keys = []
    for w in range(0,surface_deme_number):             # deme expansion occurs in the surface of a tumor
        ckey = random.choice(surface_keys)
        if space[ckey].present == 1:
            #rx,ry,rz = ckey[0],ckey[1],ckey[2]
            nei_sites = neighbor26(ckey)  # neighbor sites of (rx,ry,rz)
            empty_sites = [key for key in nei_sites if space[key].present == 0]                    # the empty neighbor sites
            empty_site_number = len(empty_sites)
            if empty_site_number > 0:
                rand_prob = random.random()
                if rand_prob < 1-math.exp(-empty_site_number*0.25): # the probability for a deme to grow and divide is proportional to the # of empty neighbor sites
                    pre_neu = list(space[ckey].neutral)
                    pre_adv1 = list(space[ckey].advant1)
                    pre_adv2 = list(space[ckey].advant2)
                    post_neu_l1,post_neu_l2,post_adv1_l1,post_adv1_l2,post_adv2_l1,post_adv2_l2,mut_id,mutlineage = demeGrowthFission_v2(pre_neu,pre_adv1,pre_adv2,mutlineage,mut_id,current_deme_number,s_coef)
                    space[ckey].neutral = list(post_neu_l1)
                    space[ckey].advant1 = list(post_adv1_l1)
                    space[ckey].advant2 = list(post_adv2_l1)
                    
                    nkey = random.choice(empty_sites)
                    space[nkey].neutral = list(post_neu_l2)
                    space[nkey].advant1 = list(post_adv1_l2)
                    space[nkey].advant2 = list(post_adv2_l2)
                    
                    space[nkey].present = 1
                    current_keys += [nkey]
                    current_deme_number += 1
                    new_keys += [nkey]
        else:
            print "something is wrong!"     
    
    ###update surface###
    surface_update = list(surface_keys+new_keys)
    surface_keys = []
    for fkey in surface_update:
        neisites = neighbor26(fkey)
        random.shuffle(neisites)
        for key in neisites:
            if space[key].present == 0:
                surface_keys += [fkey]
                break
    
    surface_deme_number = len(surface_keys)
    current_deme_number = len(current_keys)
    deme_time_generation += 1

    print "generation, no. of demes in surface and whole tumor=",deme_time_generation,surface_deme_number,current_deme_number
    
# sample from the periphery to get demes that are farther apart ...  using surface keys

quadrant1,quadrant2,quadrant3,quadrant4,quadrant5,quadrant6,quadrant7,quadrant8 = [],[],[],[],[],[],[],[] #surface demes in the eight quadrants
for pky in surface_keys:
    if pky[0] > rd and pky[1] > rd and pky[2] > rd:
        quadrant1 += [pky]
    if pky[0] < rd and pky[1] < rd and pky[2] < rd:
        quadrant2 += [pky]
    if pky[0] < rd and pky[1] > rd and pky[2] > rd:
        quadrant3 += [pky]
    if pky[0] > rd and pky[1] < rd and pky[2] < rd:
        quadrant4 += [pky]
    
    if pky[0] > rd and pky[1] > rd and pky[2] < rd:
        quadrant5 += [pky]
    if pky[0] < rd and pky[1] < rd and pky[2] > rd:
        quadrant6 += [pky]
    if pky[0] > rd and pky[1] < rd and pky[2] > rd:
        quadrant7 += [pky]
    if pky[0] < rd and pky[1] > rd and pky[2] < rd:
        quadrant8 += [pky]


# sample 3 times each from the first 4 demes and 2x each from the 2nd 4 demes for a total of 20 samples
locat1_1 = random.choice(quadrant1) # location of bulk tissue1
locat1_2 = random.choice(quadrant1) # location of bulk tissue1
locat1_3 = random.choice(quadrant1) # location of bulk tissue1
locat2_1 = random.choice(quadrant2)
locat2_2 = random.choice(quadrant2)
locat2_3 = random.choice(quadrant2)
locat3_1 = random.choice(quadrant3)
locat3_2 = random.choice(quadrant3)
locat3_3 = random.choice(quadrant3)
locat4_1 = random.choice(quadrant4)
locat4_2 = random.choice(quadrant4)
locat4_3 = random.choice(quadrant4)

locat5_1 = random.choice(quadrant5)
locat5_2 = random.choice(quadrant5)
locat6_1 = random.choice(quadrant6)
locat6_2 = random.choice(quadrant6)
locat7_1 = random.choice(quadrant7)
locat7_2 = random.choice(quadrant7)
locat8_1 = random.choice(quadrant8)
locat8_2 = random.choice(quadrant8)
sample20 = [locat1_1,locat1_2,locat1_3,locat2_1,locat2_2,locat2_3,locat3_1,locat3_2,locat3_3,locat4_1,locat4_2,locat4_3,locat5_1,locat5_2,locat6_1,locat6_2,locat7_1,locat7_2,locat8_1,locat8_2]

# sample a radius of 3 demes around the center of the sample, the radius can be adjusted to generate larger or smaller samples
tissue1_1 = bulkTissueSampling(space,sample20[0],3)
tissue1_2 = bulkTissueSampling(space,sample20[1],3)
tissue1_3 = bulkTissueSampling(space,sample20[2],3)
tissue2_1 = bulkTissueSampling(space,sample20[3],3)
tissue2_2 = bulkTissueSampling(space,sample20[4],3)
tissue2_3 = bulkTissueSampling(space,sample20[5],3)
tissue3_1 = bulkTissueSampling(space,sample20[6],3)
tissue3_2 = bulkTissueSampling(space,sample20[7],3)
tissue3_3 = bulkTissueSampling(space,sample20[8],3)
tissue4_1 = bulkTissueSampling(space,sample20[9],3)
tissue4_2 = bulkTissueSampling(space,sample20[10],3)
tissue4_3 = bulkTissueSampling(space,sample20[11],3)
tissue5_1 = bulkTissueSampling(space,sample20[12],3)
tissue5_2 = bulkTissueSampling(space,sample20[13],3)
tissue6_1 = bulkTissueSampling(space,sample20[14],3)
tissue6_2 = bulkTissueSampling(space,sample20[15],3)
tissue7_1 = bulkTissueSampling(space,sample20[16],3)
tissue7_2 = bulkTissueSampling(space,sample20[17],3)
tissue8_1 = bulkTissueSampling(space,sample20[18],3)
tissue8_2 = bulkTissueSampling(space,sample20[19],3)

maf1_1 = seqProcessing(space,tissue1_1,mutlineage,2,seq_depth,1)
maf1_2 = seqProcessing(space,tissue1_2,mutlineage,2,seq_depth,1)
maf1_3 = seqProcessing(space,tissue1_3,mutlineage,2,seq_depth,1)
maf2_1 = seqProcessing(space,tissue2_1,mutlineage,2,seq_depth,1)
maf2_2 = seqProcessing(space,tissue2_2,mutlineage,2,seq_depth,1)
maf2_3 = seqProcessing(space,tissue2_3,mutlineage,2,seq_depth,1)
maf3_1 = seqProcessing(space,tissue3_1,mutlineage,2,seq_depth,1)
maf3_2 = seqProcessing(space,tissue3_2,mutlineage,2,seq_depth,1)
maf3_3 = seqProcessing(space,tissue3_3,mutlineage,2,seq_depth,1)
maf4_1 = seqProcessing(space,tissue4_1,mutlineage,2,seq_depth,1)
maf4_2 = seqProcessing(space,tissue4_2,mutlineage,2,seq_depth,1)
maf4_3 = seqProcessing(space,tissue4_3,mutlineage,2,seq_depth,1)

maf5_1 = seqProcessing(space,tissue5_1,mutlineage,2,seq_depth,1)
maf5_2 = seqProcessing(space,tissue5_2,mutlineage,2,seq_depth,1)
maf6_1 = seqProcessing(space,tissue6_1,mutlineage,2,seq_depth,1)
maf6_2 = seqProcessing(space,tissue6_2,mutlineage,2,seq_depth,1)
maf7_1 = seqProcessing(space,tissue7_1,mutlineage,2,seq_depth,1)
maf7_2 = seqProcessing(space,tissue7_2,mutlineage,2,seq_depth,1)
maf8_1 = seqProcessing(space,tissue8_1,mutlineage,2,seq_depth,1)
maf8_2 = seqProcessing(space,tissue8_2,mutlineage,2,seq_depth,1)

MAF_file = open("tumor"+str(repl)+"_simVAF_deme"+str(deme_size)+"_s"+str(percentage)+"percent.txt","w")
MAF_file.write("mut_id"+" "+"public"+" "+"depth1_1"+" "+"maf1_1"+" "+"depth1_2"+" "+"maf1_2"+" "+"depth1_3"+" "+"maf1_3"+" "+"depth2_1"+" "+"maf2_1"+" "+"depth2_2"+" "+"maf2_2"+" "+"depth2_3"+" "+"maf2_3"+" "+"depth3_1"+" "+"maf3_1"+" "+"depth3_2"+" "+"maf3_2"+" "+"depth3_3"+" "+"maf3_3"+" "+"depth4_1"+" "+"maf4_1"+"depth4_2"+" "+"maf4_2"+"depth4_3"+" "+"maf4_3"+"depth5_1"+" "+"maf5_1"+"depth5_2"+" "+"maf5_2"+"depth6_1"+" "+"maf6_1"+"depth6_2"+" "+"maf6_2"+"depth7_1"+" "+"maf7_1"+"depth7_2"+" "+"maf7_2"+"depth8_1"+" "+"maf8_1"+"depth8_2"+" "+"maf8_2")
MAF_file.write("\n")

for k in range(0,npub):
    pdepth,pmaf = pubMutGenerator(20,2,seq_depth,1)
    MAF_file.write("c"+str(k+1)+" "+"1"+" "+str(pdepth[0])+" "+str(pmaf[0])+" "+str(pdepth[1])+" "+str(pmaf[1])+" "+str(pdepth[2])+" "+str(pmaf[2])+" "+str(pdepth[3])+" "+str(pmaf[3])+" "+str(pdepth[4])+" "+str(pmaf[4])+" "+str(pdepth[5])+" "+str(pmaf[5])+" "+str(pdepth[6])+" "+str(pmaf[6])+" "+str(pdepth[7])+" "+str(pmaf[7])+" "+str(pdepth[8])+" "+str(pmaf[8])+" "+str(pdepth[9])+" "+str(pmaf[9])+" "+str(pdepth[10])+" "+str(pmaf[10])+" "+str(pdepth[11])+" "+str(pmaf[11])+" "+str(pdepth[12])+" "+str(pmaf[12])+" "+str(pdepth[13])+" "+str(pmaf[13])+" "+str(pdepth[14])+" "+str(pmaf[14])+" "+str(pdepth[15])+" "+str(pmaf[15])+" "+str(pdepth[16])+" "+str(pmaf[16])+" "+str(pdepth[17])+" "+str(pmaf[17])+" "+str(pdepth[18])+" "+str(pmaf[18])+" "+str(pdepth[19])+" "+str(pmaf[19]))
    MAF_file.write("\n")

muts_all = sets.Set(maf1_1.keys()) | sets.Set(maf1_2.keys()) | sets.Set(maf1_3.keys()) | sets.Set(maf2_1.keys()) | sets.Set(maf2_2.keys()) | sets.Set(maf2_3.keys()) | sets.Set(maf3_1.keys()) | sets.Set(maf3_2.keys()) | sets.Set(maf3_3.keys()) | sets.Set(maf4_1.keys()) | sets.Set(maf4_2.keys()) | sets.Set(maf4_3.keys()) | sets.Set(maf5_1.keys()) | sets.Set(maf5_2.keys()) | sets.Set(maf6_1.keys()) | sets.Set(maf6_2.keys()) | sets.Set(maf7_1.keys()) | sets.Set(maf7_2.keys()) | sets.Set(maf8_1.keys()) | sets.Set(maf8_2.keys()) 

absent1_1 = muts_all-sets.Set(maf1_1.keys())
absent1_2 = muts_all-sets.Set(maf1_2.keys())
absent1_3 = muts_all-sets.Set(maf1_3.keys())
absent2_1 = muts_all-sets.Set(maf2_1.keys())
absent2_2 = muts_all-sets.Set(maf2_2.keys())
absent2_3 = muts_all-sets.Set(maf2_3.keys())
absent3_1 = muts_all-sets.Set(maf3_1.keys())
absent3_2 = muts_all-sets.Set(maf3_2.keys())
absent3_3 = muts_all-sets.Set(maf3_3.keys())
absent4_1 = muts_all-sets.Set(maf4_1.keys())
absent4_2 = muts_all-sets.Set(maf4_2.keys())
absent4_3 = muts_all-sets.Set(maf4_3.keys())
absent5_1 = muts_all-sets.Set(maf5_1.keys())
absent5_2 = muts_all-sets.Set(maf5_2.keys())
absent6_1 = muts_all-sets.Set(maf6_1.keys())
absent6_2 = muts_all-sets.Set(maf6_2.keys())
absent7_1 = muts_all-sets.Set(maf7_1.keys())
absent7_2 = muts_all-sets.Set(maf7_2.keys())
absent8_1 = muts_all-sets.Set(maf8_1.keys())
absent8_2 = muts_all-sets.Set(maf8_2.keys())

maf1_1 = missingDepth(maf1_1,absent1_1,seq_depth)
maf1_2 = missingDepth(maf1_2,absent1_2,seq_depth)
maf1_3 = missingDepth(maf1_3,absent1_3,seq_depth)
maf2_1 = missingDepth(maf2_1,absent2_1,seq_depth)
maf2_2 = missingDepth(maf2_2,absent2_2,seq_depth)
maf2_3 = missingDepth(maf2_3,absent2_3,seq_depth)
maf3_1 = missingDepth(maf3_1,absent3_1,seq_depth)
maf3_2 = missingDepth(maf3_2,absent3_2,seq_depth)
maf3_3 = missingDepth(maf3_3,absent3_3,seq_depth)
maf4_1 = missingDepth(maf4_1,absent4_1,seq_depth)
maf4_2 = missingDepth(maf4_2,absent4_2,seq_depth)
maf4_3 = missingDepth(maf4_3,absent4_3,seq_depth)
maf5_1 = missingDepth(maf5_1,absent5_1,seq_depth)
maf5_2 = missingDepth(maf5_2,absent5_2,seq_depth)
maf6_1 = missingDepth(maf6_1,absent6_1,seq_depth)
maf6_2 = missingDepth(maf6_2,absent6_2,seq_depth)
maf7_1 = missingDepth(maf7_1,absent7_1,seq_depth)
maf7_2 = missingDepth(maf7_2,absent7_2,seq_depth)
maf8_1 = missingDepth(maf8_1,absent8_1,seq_depth)
maf8_2 = missingDepth(maf8_2,absent8_2,seq_depth)


# output the mafs for the 20 samples, the first number indicates which octant of the tumor the sample is take from
for mt in sorted(muts_all):
    MAF_file.write(str(mt)+" "+"0"+" "+str(maf1_1[mt][0])+" "+str(maf1_1[mt][1])+" "+str(maf1_2[mt][0])+" "+str(maf1_2[mt][1])+" "+str(maf1_3[mt][0])+" "+str(maf1_3[mt][1])+" "+str(maf2_1[mt][0])+" "+str(maf2_1[mt][1])+" "+str(maf2_2[mt][0])+" "+str(maf2_2[mt][1])+" "+str(maf2_3[mt][0])+" "+str(maf2_3[mt][1])+" "+str(maf3_1[mt][0])+" "+str(maf3_1[mt][1])+" "+str(maf3_2[mt][0])+" "+str(maf3_2[mt][1])+" "+str(maf3_3[mt][0])+" "+str(maf3_3[mt][1])+" "+str(maf4_1[mt][0])+" "+str(maf4_1[mt][1])+" "+str(maf4_2[mt][0])+" "+str(maf4_2[mt][1])+" "+str(maf4_3[mt][0])+" "+str(maf4_3[mt][1])+" "+str(maf5_1[mt][0])+" "+str(maf5_1[mt][1])+" "+str(maf5_2[mt][0])+" "+str(maf5_2[mt][1])+" "+str(maf6_1[mt][0])+" "+str(maf6_1[mt][1])+" "+str(maf6_2[mt][0])+" "+str(maf6_2[mt][1])+" "+str(maf7_1[mt][0])+" "+str(maf7_1[mt][1])+" "+str(maf7_2[mt][0])+" "+str(maf7_2[mt][1])+" "+str(maf8_1[mt][0])+" "+str(maf8_1[mt][1])+" "+str(maf8_2[mt][0])+" "+str(maf8_2[mt][1])+" ")
    MAF_file.write("\n")


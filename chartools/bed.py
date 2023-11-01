import subprocess
import numpy as np


def make_tx_to_ix(tx_dict, genes_dict, genes_order):
    for k, v in genes_dict.items():
        tx_dict[k]=v
        
    gene_to_ix={}

    for i, k in enumerate(genes_order):
        gene_to_ix[k]=i
        
    tx_to_ix={}
    for k, v in tx_dict.items():
        if v[0] in gene_to_ix:
            tx_to_ix[k]=gene_to_ix[v[0]]

    return tx_to_ix

def sites_table(pileup, tx_to_ix, site_to_ix, cmd_filter=None, nmax=-1):
    # make_tx_to_ix(tx_dict, genes_dict, genes_order)
    cmd="gunzip -c "+pileup
    if not cmd_filter is None:
        cmd+=" | "+cmd_filter
    print(cmd)
    nok=0

    nsites=max((v for _, v in site_to_ix.items()))+1
    ngenes=max((v for _, v in tx_to_ix.items()))+1
    

    M=np.zeros((ngenes, nsites))

    
    
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
    nok=0
    for _, line in enumerate(p.stdout):
        if nmax>0 and nok>nmax:
            break
        read_data=line.strip().split("\t")
#         print(line)
        siteName=read_data[3]
        geneName=read_data[5]
        i=tx_to_ix.get(geneName, -1)
        j=site_to_ix.get(siteName, -1)
#         L=int(read_data[-2])
        if (i>-1) and (j>-1):
            M[i,j]+=1
#         print("%g\t%g"%(int(read_data[3]), int(read_data[-1])))
        nok+=1
    p.terminate()
    del p
    
    return M

def sites_dist(pileup, site_to_ix, tx_to_ix=None, cmd_filter=None, nmax=-1):
    cmd="gunzip -c "+pileup
    if not cmd_filter is None:
        cmd+=" | "+cmd_filter
    print(cmd)

    nok=0
    nsites=max((v for _, v in site_to_ix.items()))+1
    sites_vector=np.zeros(nsites)

    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
    nok=0
    if not tx_to_ix is None:
        for _, line in enumerate(p.stdout):
            if nmax>0 and nok>nmax:
                break
            read_data=line.strip().split("\t")
            if read_data[5] in tx_to_ix:
                siteName=read_data[3]
                j=site_to_ix.get(siteName, -1)
                if j>-1:
                    sites_vector[j]+=1
            nok+=1
        p.terminate()
        del p
    else:
        for _, line in enumerate(p.stdout):
            if nmax>0 and nok>nmax:
                break
            read_data=line.strip().split("\t")
            siteName=read_data[3]
            j=site_to_ix.get(siteName, -1) 
            if j>-1:
                sites_vector[j]+=1
            nok+=1
        p.terminate()
        del p

    return sites_vector


def filter_mad(v, site_to_ix=None, madx_thresh=10):
    me=np.median(v)
    mad=np.median(np.abs(v-me))
    isok=np.where(np.abs(v-me)<madx_thresh*mad)[0]

    print("Med=%g, MAD=%g, Found %g sites over %g x MAD [%g total]"%(me,mad,len(v)-1-len(isok),madx_thresh, len(v)-1))

    v_clean=v[isok]

    # if not site_to_ix is None:
    #     site_to_ix2={k:v for k, v in site_to_ix.items() if v in isok}

    
    #     return v_clean, site_to_ix2
    # else:
    #isok are indeces
    return v_clean, isok

def genes_dist(pileup, tx_to_ix, site_to_ix=None, cmd_filter=None, nmax=-1):
    cmd="gunzip -c "+pileup
    if not cmd_filter is None:
        cmd+=" | "+cmd_filter
    print(cmd)

    nok=0
    ngenes=max((v for _, v in tx_to_ix.items()))+1
    genes_vector=np.zeros(ngenes)

    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
    nok=0
    if not site_to_ix is None:
        for _, line in enumerate(p.stdout):
            if nmax>0 and nok>nmax:
                break
            read_data=line.strip().split("\t")
            if read_data[3] in site_to_ix:
                geneName=read_data[5]
                j=tx_to_ix.get(geneName, -1)
                if j>-1:
                    genes_vector[j]+=1
            nok+=1
        p.terminate()
        del p
    else:
        for _, line in enumerate(p.stdout):
            if nmax>0 and nok>nmax:
                break
            read_data=line.strip().split("\t")
            geneName=read_data[5]
            j=tx_to_ix.get(geneName, 0)
            if j>-1:
                genes_vector[j]+=1
            nok+=1
        p.terminate()
        del p

    return genes_vector




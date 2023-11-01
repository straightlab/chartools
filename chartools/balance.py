import pandas as pd
import subprocess
import numpy as np
from scipy import linalg
# def dpnII_compensate(O, dpnII_vector, inplace=True):
#     Y=O
#     m=dpnII_vector.mean()
#     if not inplace:
#         Y=np.array(O)
#     Y=Y/dpnII_vector[np.newaxis,...]
#     return Y





def _balance_trans_global(O, rows_chr_ids, col_chr_blocks, bkg_compensation=None, inplace=True, expect_only=False): 
    # dispatches the counts over all trans chromosomes

    #O=oberved RPKM MATRIX
    #row_chr_blocks [0,10,10,...]--> block definition [chr1, chr2, ....]
    #col_chr_blocks 
    ngenes, nbins=O.shape
    nchr=len(col_chr_blocks)-1
    Y=O


    if bkg_compensation in None:
        bkg_compensation=np.ones(nchr) #compensate 

    if not inplace:
        Y=np.array(O)
    
    
    if expect_only:
        Y=np.zeros(O.shape)
        for i in range(ngenes):
            s1=col_chr_blocks[rows_chr_ids[i]]
            s2=col_chr_blocks[rows_chr_ids[i]+1]
            e=(O[i,0:s1].sum()+O[i,s2:-1].sum())/(nbins-s2+s1)*dpnII_compensation[rows_chr_ids[i]] #average over all other chr
            Y[i,0:s1]=e
            Y[i,s2:-1]=e
    else:
        for i in range(ngenes):
            s1=col_chr_blocks[rows_chr_ids[i]]
            s2=col_chr_blocks[rows_chr_ids[i]+1]
            e=(O[i,0:s1].sum()+O[i,s2:-1].sum())/(nbins-s2+s1)*dpnII_compensation[rows_chr_ids[i]]
            Y[i,0:s1]=O[i,0:s1]/e
            Y[i,s2:-1]=O[i,s2:-1]/e
    return Y
    
def _balance_trans_chr(O, rows_chr_ids, col_chr_blocks, trans_proba=None, inplace=True, expect_only=False): #O=oberved RPKM MATRIX
#row_chr_blocks [0,10,10,...]--> block definition [chr1, chr2, ....]
    #col_chr_blocks 
    ngenes, nbins=O.shape
    nchr=len(col_chr_blocks)-1
    Y=O

    if trans_proba is None:
        trans_proba=np.ones((nchr, nchr))
    if not inplace:
        Y=np.array(O)
    
    if expect_only:
        Y=np.ones(O.shape)
        for i in range(ngenes):
            cischr=rows_chr_ids[i]
            s1=col_chr_blocks[cischr]
            s2=col_chr_blocks[cischr+1]
#             for j in range(s1):                
            e_src=np.mean(O[i,s1:s2]) #cis rate
            for j in range(nchr):
                if not (j==rows_chr_ids[i]):
                
                    e=e_src*trans_proba[cischr,j]
                    Y[i,col_chr_blocks[j]:(col_chr_blocks[j+1])]=e
                else:
                    Y[i,s1:s2]=O[i,s1:s2]
                   
    else:
        for i in range(ngenes):
            cischr=rows_chr_ids[i]
            s1=col_chr_blocks[cischr]
            s2=col_chr_blocks[cischr+1]
#             for j in range(s1):                
            e_src=np.mean(O[i,s1:s2]) #cis rate
            for j in range(nchr):
                if not (j==rows_chr_ids[i]):
                
                    e=e_src*trans_proba[cischr,j]
                    Y[i,col_chr_blocks[j]:(col_chr_blocks[j+1])]=O[i,col_chr_blocks[j]:(col_chr_blocks[j+1])]/e
                else:
                    Y[i,s1:s2]=1.0
    return Y

# def _balance_cis_chr(O, rows_chr_ids, col_chr_blocks, rows_binoforigin, decay_model=None, inplace=True, expect_only=False):
    
#     ngenes, nbins=O.shape
#     nchr=len(col_chr_blocks)-1
#     Y=O

#     ix0=int((len(decay_model)-1)/2)+1
    

#     if not inplace:
#         Y=np.array(O)
    
#     if expect_only:
#         Y=np.zeros(O.shape)
#         for i in range(ngenes):
#             j=rows_chr_ids[i] #only modify these
#             O_cis=O[i,col_chr_blocks[j]:(col_chr_blocks[j]+1)]
#             # e=O[i,col_chr_blocks[j]:(col_chr_blocks[j]+1)].sum() #take the little slicey
#             delta=ix0
#             # O[i,col_chr_blocks[j]:(col_chr_blocks[j]+1)]
#             ixboi=rows_binoforigin[i]
#             Y[i,col_chr_blocks[j]:(col_chr_blocks[j]+1)]=O_cis.sum()*decay_model[(ix0-ixboi):(ix0-ixboi+nbins)]
                   
#     else:
#         for i in range(ngenes):
#            O_cis=O[i,col_chr_blocks[j]:(col_chr_blocks[j]+1)]
#             # e=O[i,col_chr_blocks[j]:(col_chr_blocks[j]+1)].sum() #take the little slicey
#             delta=ix0
#             # O[i,col_chr_blocks[j]:(col_chr_blocks[j]+1)]
#             ixboi=rows_binoforigin[i]
#             Y[i,col_chr_blocks[j]:(col_chr_blocks[j]+1)]=Y[i,col_chr_blocks[j]:(col_chr_blocks[j]+1)]/(O_cis.sum()*decay_model[(ix0-ixboi):(ix0-ixboi+nbins)])
#     return Y



def get_sampling_distribution(O, rows_chr_ids, col_chr_blocks, bychr=True):
    # returns Nj / N in block j
    nchr=len(col_chr_blocks)-1
    b=O.sum(axis=0)
    # ngenes, nbins=cmap.shape
    # b=np.zeros(nbins)
    Ntrans=np.zeros(nchr)
    for i in range(nchr):
#         this_chr_rng=
        b[col_chr_blocks[i]:(col_chr_blocks[i+1])]-=O[rows_chr_ids==i,col_chr_blocks[i]:col_chr_blocks[i+1]].sum(axis=0)
        Ntrans[i]=b[col_chr_blocks[i]:(col_chr_blocks[i+1])].sum()
        if bychr:
            b[col_chr_blocks[i]:(col_chr_blocks[i+1])]=b[col_chr_blocks[i]:(col_chr_blocks[i+1])]/Ntrans[i]

    if not bychr:
        b=b/Ntrans.sum()

    return b, Ntrans

def get_trans_matrix(O, rows_chr_ids, col_chr_blocks=None):
    # returns Nj / N in block j
    # if col_chr_blocks is None:
    #     col_chr_blocks=np.arange(())
    # rows_chr_ids=genes_df['chr_id'].loc[cmap.index].values.astype(int)
    nchr=len(col_chr_blocks)-1
    # b=np.array(cmap.sum(axis=0).values)
    ngenes, nbins=O.shape
    # b=np.zeros(nbins)
    Mtrans=np.zeros((nchr,nchr))
    
    for i in range(ngenes):
        for j in range(nchr):
            Mtrans[rows_chr_ids[i],j]+=O[i, col_chr_blocks[j]:(col_chr_blocks[j+1])].sum()

    return Mtrans



# def _get_expected(O, rows_chr_ids, col_chr_blocks, sampling_dist, trans_proba=None, inplace=True, expect_only=False, bychr=True, src_relative=True): #O=oberved RPKM MATRIX
# #row_chr_blocks [0,10,10,...]--> block definition [chr1, chr2, ....]
#     #col_chr_blocks 
#     ngenes, nbins=O.shape
#     nchr=len(col_chr_blocks)-1
#     Y=O

#     if trans_proba is None:
#         trans_proba=np.ones((nchr, nchr))
#     if not inplace:
#         Y=np.array(O)
    
#     if expect_only:
#         Y=np.ones(O.shape)
#         for i in range(ngenes):
#             cischr=rows_chr_ids[i]
#             s1=col_chr_blocks[cischr]
#             s2=col_chr_blocks[cischr+1]
# #             for j in range(s1):                
#             e_src=np.sum(O[i,s1:s2]) #cis contacts rate
#             for j in range(nchr):
#                 if not (j==rows_chr_ids[i]):
#                     if src_relative:
#                         e=e_src*trans_proba[cischr,j]
#                     else:
#                         e=np.sum(Y[i,col_chr_blocks[j]:(col_chr_blocks[j+1])])
#                     Y[i,col_chr_blocks[j]:(col_chr_blocks[j+1])]=e*sampling_dist[col_chr_blocks[j]:(col_chr_blocks[j+1])]
#                 else:
#                     Y[i,s1:s2]=O[i,s1:s2]
                   
#     else:
#         for i in range(ngenes):
#             cischr=rows_chr_ids[i]
#             s1=col_chr_blocks[cischr]
#             s2=col_chr_blocks[cischr+1]
# #             for j in range(s1):                
#             e_src=np.sum(O[i,s1:s2]) #cis contacts rate
#             for j in range(nchr):
#                 if not (j==rows_chr_ids[i]):
#                     if src_relative:
#                         e=e_src*trans_proba[cischr,j]
#                     else:
#                         e=np.sum(Y[i,col_chr_blocks[j]:(col_chr_blocks[j+1])])
#                     Y[i,col_chr_blocks[j]:(col_chr_blocks[j+1])]=O[i,col_chr_blocks[j]:(col_chr_blocks[j+1])]/(e*sampling_dist[col_chr_blocks[j]:(col_chr_blocks[j+1])])
#                 else:
#                     Y[i,s1:s2]=1.0
#     return Y

def _get_expected(O, rows_chr_ids, col_chr_blocks, sampling_dist, trans_proba=None, inplace=True, expect_only=False, bychr=True, src_relative=True, delta=0.0): #O=oberved RPKM MATRIX
#row_chr_blocks [0,10,10,...]--> block definition [chr1, chr2, ....]
    #col_chr_blocks 
    ngenes, nbins=O.shape
    nchr=len(col_chr_blocks)-1
    Y=O

    if trans_proba is None:
        trans_proba=np.ones((nchr, nchr))
    if not inplace:
        Y=np.array(O)
    
    if expect_only:
        Y=np.ones(O.shape)
        for i in range(ngenes):
            cischr=rows_chr_ids[i]
            s1=col_chr_blocks[cischr]
            s2=col_chr_blocks[cischr+1]
#             for j in range(s1):                
            e_src=np.sum(O[i,s1:s2]) #cis contacts rate
            for j in range(nchr):
                if not (j==rows_chr_ids[i]):
                    if src_relative:
                        e=e_src*trans_proba[cischr,j]
                    else:
                        e=np.sum(O[i,col_chr_blocks[j]:(col_chr_blocks[j+1])])
                    Y[i,col_chr_blocks[j]:(col_chr_blocks[j+1])]=e*sampling_dist[col_chr_blocks[j]:(col_chr_blocks[j+1])]
                else:
                    Y[i,s1:s2]=O[i,s1:s2]
                   
    else:
        for i in range(ngenes):
            cischr=rows_chr_ids[i]
            s1=col_chr_blocks[cischr]
            s2=col_chr_blocks[cischr+1]
#             for j in range(s1):                
            e_src=np.sum(O[i,s1:s2]) #cis contacts rate
            for j in range(nchr):
                if not (j==rows_chr_ids[i]):
                    if src_relative:
                        e=e_src*trans_proba[cischr,j]
                    else:
                        e=np.sum(O[i,col_chr_blocks[j]:(col_chr_blocks[j+1])])
                    Y[i,col_chr_blocks[j]:(col_chr_blocks[j+1])]=(delta+O[i,col_chr_blocks[j]:(col_chr_blocks[j+1])])/(delta+(e*sampling_dist[col_chr_blocks[j]:(col_chr_blocks[j+1])]))
                else:
                    Y[i,s1:s2]=1.0
    return Y


def _get_expected_total(O, expected_cis, rows_chr_ids, col_chr_blocks, sampling_dist, trans_proba=None, inplace=True, expect_only=False, bychr=True, src_relative=True, delta=0.0): #O=oberved RPKM MATRIX
#row_chr_blocks [0,10,10,...]--> block definition [chr1, chr2, ....]
    #col_chr_blocks 
    ngenes, nbins=O.shape
    nchr=len(col_chr_blocks)-1
    Y=O

    if trans_proba is None:
        trans_proba=np.ones((nchr, nchr))
    if not inplace:
        Y=np.array(O)
    
    gene_blocks=np.hstack([0,np.cumsum(np.array([k.shape[0] for k in expected_cis]))])
    if expect_only:
        Y=np.ones(O.shape)
        for i in range(ngenes):
            cischr=rows_chr_ids[i]
            s1=col_chr_blocks[cischr]
            s2=col_chr_blocks[cischr+1]
#             for j in range(s1):                
            e_src=np.sum(O[i,s1:s2]) #cis contacts rate
            for j in range(nchr):
                if not (j==rows_chr_ids[i]):
                    if src_relative:
                        e=e_src*trans_proba[cischr,j]
                    else:
                        e=np.sum(O[i,col_chr_blocks[j]:(col_chr_blocks[j+1])])
                    Y[i,col_chr_blocks[j]:(col_chr_blocks[j+1])]=e*sampling_dist[col_chr_blocks[j]:(col_chr_blocks[j+1])]
                else:
                    Y[i,s1:s2]=expected_cis[j][i-gene_blocks[j],:] #.iloc[i,:].values
                    Y[i,s1:s2]=Y[i, s1:s2]*(np.sum(O[i, s1:s2])/np.sum(Y[i, s1:s2]))
                   
        # for j in range(nchr):
        #     Y[gene_blocks[j]:gene_blocks[j+1],col_chr_blocks[j]:(col_chr_blocks[j+1])]=expected_cis[j]
    else:
        for i in range(ngenes):
            cischr=rows_chr_ids[i]
            s1=col_chr_blocks[cischr]
            s2=col_chr_blocks[cischr+1]
#             for j in range(s1):                
            e_src=np.sum(O[i,s1:s2]) #cis contacts rate
            for j in range(nchr):
                if not (j==rows_chr_ids[i]):
                    if src_relative:
                        e=e_src*trans_proba[cischr,j]
                    else:
                        e=np.sum(O[i,col_chr_blocks[j]:(col_chr_blocks[j+1])])
                    Y[i,col_chr_blocks[j]:(col_chr_blocks[j+1])]=(delta+O[i,col_chr_blocks[j]:(col_chr_blocks[j+1])])/(delta+(e*sampling_dist[col_chr_blocks[j]:(col_chr_blocks[j+1])]))
                else:
                    Y[i,s1:s2]=expected_cis[j][i-gene_blocks[j],:] #.iloc[i,:].values
                    Y[i,s1:s2]=O[i,s1:s2]/(Y[i, s1:s2]*(np.sum(O[i, s1:s2])/np.sum(Y[i, s1:s2])))
        # for j in range(nchr):
        #     Y[gene_blocks[j]:gene_blocks[j+1],col_chr_blocks[j]:(col_chr_blocks[j+1])]=O[gene_blocks[j]:gene_blocks[j+1],col_chr_blocks[j]:(col_chr_blocks[j+1])]/expected_cis[j]
    return Y

def fly_distribution_global(pairs, chr_dict, bins=np.linspace(-1000000,1000000,201), pairs_filter=None, nmax=0): #dict of vectors, which are targeted chr, first numer is # match,
# second is # of dna match

    if pairs_filter is None:
        cmd = "cat "+pairs+" | awk '($25==1 && $9>15){print $4, substr($2,3), $5, $3}'"
    else:
        cmd = "cat "+pairs+" | awk '("+pairs_filter+"){print $4, substr($2,3), $5, $3}'"

    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
    nok=0

    nchr=len(chr_dict)
    nbins=len(bins)+1
    

    out=np.zeros(nbins,int)
    
    
        # read all reads in the file unitl reach end of file
    for _, line in enumerate(p.stdout):
        if nmax>0 and nok>nmax:
            break
        read_data=line.strip().split(" ")
        tgt_chr=chr_dict.get(read_data[0],-1)
        src_chr=chr_dict.get(read_data[1],-1)

        if tgt_chr>-1 and (src_chr==tgt_chr): #only if this is good stuff
            
            nok+=1
            
            fly=int(read_data[2])-int(read_data[3])
                # src_pos=int(read_data[3])
            out[np.searchsorted(bins,fly)]+=1
    
            
    # dw=np.hstack([1,np.diff(bins),1])
    # out=out*(1000000.0/dw)*(1000000.0/nok)

    p.terminate()
    del p

    return out, nok

def fly_distribution_egdgeCorrected(pairs, chr_dict, chr_len_vec=None, bins=np.linspace(-1000000,1000000,201), pairs_filter=None, nmax=0): #dict of vectors, which are targeted chr, first numer is # match,
# second is # of dna match

    if pairs_filter is None:
        cmd = "cat "+pairs+" | awk '($25==1 && $9>15){print $4, substr($2,3), $5, $3}'"
    else:
        cmd = "cat "+pairs+" | awk '("+pairs_filter+"){print $4, substr($2,3), $5, $3}'"

    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
    nok=0

    nchr=len(chr_dict)
    nbins=len(bins)+1
    

    out=np.zeros(nbins,int)
    possible_contacts=np.zeros(nbins,int)
    
        # read all reads in the file unitl reach end of file
    for _, line in enumerate(p.stdout):
        if nmax>0 and nok>nmax:
            break
        read_data=line.strip().split(" ")
        tgt_chr=chr_dict.get(read_data[0],-1)
        src_chr=chr_dict.get(read_data[1],-1)

        if tgt_chr>-1 and (src_chr==tgt_chr): #only if this is good stuff
            
            nok+=1
            
            fly=int(read_data[2])-int(read_data[3])
            src_pos=int(read_data[3])
            out[np.searchsorted(bins,fly)]+=1
            maxflight=chr_len_vec[tgt_chr]-src_pos
            minflight=-src_pos
            possible_contacts[np.searchsorted(bins,minflight):(np.searchsorted(bins,maxflight)+1)]+=1
            
    # dw=np.hstack([1,np.diff(bins),1])
    # out=out*(1000000.0/dw)*(1000000.0/nok)

    p.terminate()
    del p

    out=out/np.maximum(1,possible_contacts)
    out=out/np.sum(out)
    return out, nok, possible_contacts

# def finite_chr_correction(genes_df, binsize, col_chr_blocks):

def fly_distribution_bychr(pairs, chr_dict, chr_len_vec=None, bins=np.linspace(-1000000,1000000,201), pairs_filter=None, nmax=0): #dict of vectors, which are targeted chr, first numer is # match,
# second is # of dna match

    if pairs_filter is None:
        cmd = "cat "+pairs+" | awk '($25==1 && $9>15){print $4, substr($2,3), $5, $3}'"
    else:
        cmd = "cat "+pairs+" | awk '("+pairs_filter+"){print $4, substr($2,3), $5, $3}'"

    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
    nok=0

    nchr=len(chr_dict)
    nbins=len(bins)+1
    

    out=np.zeros((nchr,nbins),int)
    possible_contacts=np.zeros((nchr,nbins),int)
    NN=np.zeros(nchr,int)
    
        # read all reads in the file unitl reach end of file
    for _, line in enumerate(p.stdout):
        if nmax>0 and nok>nmax:
            break
        read_data=line.strip().split(" ")
        tgt_chr=chr_dict.get(read_data[0],-1)
        src_chr=chr_dict.get(read_data[1],-1)

        if tgt_chr>-1 and (src_chr==tgt_chr): #only if this is good stuff
            NN[src_chr]+=1
            nok+=1
            
            fly=int(read_data[2])-int(read_data[3])
            src_pos=int(read_data[3])
            out[src_chr,np.searchsorted(bins,fly)]+=1
            maxflight=chr_len_vec[tgt_chr]-src_pos
            minflight=-src_pos
            possible_contacts[src_chr, np.searchsorted(bins,minflight):(np.searchsorted(bins,maxflight)+1)]+=1
            
    # dw=np.hstack([1,np.diff(bins),1])
    # out=out*(1000000.0/dw)*(1000000.0/nok)

    p.terminate()
    del p

    out=out/np.maximum(1,possible_contacts)
    out=out/(np.sum(out, axis=1)[...,np.newaxis])

    return out, possible_contacts, NN


def make_expected_cis_(pairs, chr_len_dict, chr_to_ID, gene_to_ID, tx_dict, binsize=1000000, bygene=True, modelP=None, modelM=None, qmin_rna=0, qmin_dna=15, noambiv=True, nmax=0): #dict of vectors, which are targeted chr, first numer is # match,
# second is # of dna match'
    if noambiv:
        cmd="cat "+pairs+" | awk '($25==1 && $8>"+str(qmin_rna)+" && $9>"+str(qmin_dna)+"){print $18, $4, substr($2,3), $5, $3}'"
    else:
        cmd = "cat "+pairs+" | awk '($8>"+str(qmin_rna)+" && $9>"+str(qmin_dna)+"){print $18, $4, substr($2,3), $5, $3}'"
    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
    nok=0

    ngenes=max([v[0] for _, v in gene_to_ID.items()])+1
    
    # mP=interpolate.interp1d(model_x, model_yP, kind='linear')
    # mM=interpolate.interp1d(model_x, model_yM, kind='linear')

    totalchrlen=np.sum([v for _, v in chr_len_dict.items()])
    nchr=max([v for _, v in chr_to_ID.items()])+1
    ID_to_chr={v:k for k, v in chr_to_ID.items()}
    chr_len_vec=np.zeros(nchr)

    NN=np.zeros(nchr)
    
    bins=np.arange(0,250000000, binsize)
    
    for i in range(nchr):
        chr_len_vec[i]=chr_len_dict.get(ID_to_chr[i],0)
    chr_offset_vec=np.hstack([0,np.cumsum(chr_len_vec)[0:-1]])
    
    ngenes_bychr=np.zeros(nchr, int)
    for _, v in gene_to_ID.items():
        ngenes_bychr[v[0]]=max(ngenes_bychr[v[0]],v[1]+1)
    
    nbins=[int(np.floor(((clen)*1.0)/binsize)) for clen in chr_len_vec]
#     print(nbins)
    out=[np.zeros((ngenes_bychr[k],nbins[k])) for k in range(nchr)]
    ntotal=0
    for _, line in enumerate(p.stdout):
        ntotal+=1
        if nmax>0 and nok>nmax:
            break
        if (ntotal%100000)==0:
            print("Million reads processed: %g"%(ntotal/1000000))
        read_data=line.strip().split(" ")
        T_id=read_data[0]
        T_annot=tx_dict.get(T_id,[])
        # tgt_chr_offset=chr_off_dict.get(read_data[1],-1)
        # src_chr_offset=chr_off_dict.get(read_data[2],-1)
        tgt_chr=chr_to_ID.get(read_data[1],-1)
        src_chr=chr_to_ID.get(read_data[2],-1)
        if bygene:
            Tid=T_annot[0] if len(T_annot)>0 else "*"
        else:
            Tid="*" if (skip_unknown_tx and len(T_annot)<2) else T_id

        gene_IDs=gene_to_ID.get(Tid,[-1,-1])
        gene_ID=gene_IDs[1]
        gchr=gene_IDs[0]
        # print(chr_offset_vec)

        if (gene_ID>-1) and (tgt_chr>-1) and (src_chr==tgt_chr) and (src_chr==gchr): #only if this is good stuff
            
            NN[src_chr]+=1
            cout=out[src_chr]
            # tgt_binid=int(np.floor((read_data[3]-1)/binsize))
            src_pos=int(read_data[4])-1
            
            src_binid=int(np.floor(src_pos/binsize))
            if src_binid<nbins[src_chr]:
                nok+=1
                # xedges=bins[(src_binid+1):]-src_pos
                valsP=modelP(bins[(src_binid+1):]-src_pos)
                valsM=modelM(-bins[0:(src_binid+1)]+src_pos)
                cout[gene_ID,src_binid]+=(valsP[0]+valsM[-1])
                #nbins_src_chr-src_bin_id-1

                # xx=np.diff(valsP[0:(nbins[src_chr]-src_binid)])
    #             print("%g_%g_%g_%g_%g"%(src_chr,src_binid,cout.shape[1],nbins[src_chr],len(xx)))
                cout[gene_ID,(src_binid+1):]=cout[gene_ID,(src_binid+1):]+valsP[1:(nbins[src_chr]-src_binid)]-valsP[0:(nbins[src_chr]-src_binid-1)]#np.diff(valsP[0:(nbins[src_chr]-src_binid)])
                cout[gene_ID,0:src_binid]=cout[gene_ID,0:src_binid]+valsM[(-src_binid-1):-1]-valsM[(-src_binid):]
#             cout[gene_ID,0:src_binid+1]+=np.diff(valsP[0:(nbins[src_chr]-src_binid)])
            # tgt_binid=int(np.floor((read_data[3]-1)/binsize))

            # out[gene_ID,binid]+=1
    
    p.terminate()
    del p


    for k, cout in enumerate(out):
        cout=cout*(1000000.0/binsize)*(1000000.0/NN[k])

    return out, NN

def make_expected_cis_by_chr(pairs, chr_len_dict, chr_to_ID, gene_to_ID, tx_dict, binsize=1000000, bygene=True, modelP=None, modelM=None, qmin_rna=0, qmin_dna=15, noambiv=True, nmax=0): #dict of vectors, which are targeted chr, first numer is # match,
# second is # of dna match'
    if noambiv:
        cmd="cat "+pairs+" | awk '($25==1 && $8>"+str(qmin_rna)+" && $9>"+str(qmin_dna)+"){print $18, $4, substr($2,3), $5, $3}'"
    else:
        cmd = "cat "+pairs+" | awk '($8>"+str(qmin_rna)+" && $9>"+str(qmin_dna)+"){print $18, $4, substr($2,3), $5, $3}'"
    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
    nok=0

    ngenes=max([v[0] for _, v in gene_to_ID.items()])+1
    
    # mP=interpolate.interp1d(model_x, model_yP, kind='linear')
    # mM=interpolate.interp1d(model_x, model_yM, kind='linear')

    totalchrlen=np.sum([v for _, v in chr_len_dict.items()])
    nchr=max([v for _, v in chr_to_ID.items()])+1
    ID_to_chr={v:k for k, v in chr_to_ID.items()}
    chr_len_vec=np.zeros(nchr)

    NN=np.zeros(nchr)
    
    bins=np.arange(0,250000000, binsize)
    
    for i in range(nchr):
        chr_len_vec[i]=chr_len_dict.get(ID_to_chr[i],0)
    chr_offset_vec=np.hstack([0,np.cumsum(chr_len_vec)[0:-1]])
    
    ngenes_bychr=np.zeros(nchr, int)
    for _, v in gene_to_ID.items():
        ngenes_bychr[v[0]]=max(ngenes_bychr[v[0]],v[1]+1)
    
    nbins=[int(np.floor(((clen)*1.0)/binsize)) for clen in chr_len_vec]
#     print(nbins)
    out=[np.zeros((ngenes_bychr[k],nbins[k])) for k in range(nchr)]
    ntotal=0
    for _, line in enumerate(p.stdout):
        ntotal+=1
        if nmax>0 and nok>nmax:
            break
        if (ntotal%100000)==0:
            print("Million reads processed: %g"%(ntotal/1000000))
        read_data=line.strip().split(" ")
        T_id=read_data[0]
        T_annot=tx_dict.get(T_id,[])
        # tgt_chr_offset=chr_off_dict.get(read_data[1],-1)
        # src_chr_offset=chr_off_dict.get(read_data[2],-1)
        tgt_chr=chr_to_ID.get(read_data[1],-1)
        src_chr=chr_to_ID.get(read_data[2],-1)
        if bygene:
            Tid=T_annot[0] if len(T_annot)>0 else "*"
        else:
            Tid="*" if (skip_unknown_tx and len(T_annot)<2) else T_id

        gene_IDs=gene_to_ID.get(Tid,[-1,-1])
        gene_ID=gene_IDs[1]
        gchr=gene_IDs[0]
        # print(chr_offset_vec)

        if (gene_ID>-1) and (tgt_chr>-1) and (src_chr==tgt_chr) and (src_chr==gchr): #only if this is good stuff
            
            NN[src_chr]+=1
            cout=out[src_chr]
            # tgt_binid=int(np.floor((read_data[3]-1)/binsize))
            src_pos=int(read_data[4])-1
            
            src_binid=int(np.floor(src_pos/binsize))
            if src_binid<nbins[src_chr]:
                nok+=1
                # xedges=bins[(src_binid+1):]-src_pos
                valsP=modelP[src_chr](bins[(src_binid+1):]-src_pos)
                valsM=modelM[src_chr](-bins[0:(src_binid+1)]+src_pos)
                cout[gene_ID,src_binid]+=(valsP[0]+valsM[-1])
                #nbins_src_chr-src_bin_id-1

                # xx=np.diff(valsP[0:(nbins[src_chr]-src_binid)])
    #             print("%g_%g_%g_%g_%g"%(src_chr,src_binid,cout.shape[1],nbins[src_chr],len(xx)))
                cout[gene_ID,(src_binid+1):]=cout[gene_ID,(src_binid+1):]+valsP[1:(nbins[src_chr]-src_binid)]-valsP[0:(nbins[src_chr]-src_binid-1)]#np.diff(valsP[0:(nbins[src_chr]-src_binid)])
                cout[gene_ID,0:src_binid]=cout[gene_ID,0:src_binid]+valsM[(-src_binid-1):-1]-valsM[(-src_binid):]
#             cout[gene_ID,0:src_binid+1]+=np.diff(valsP[0:(nbins[src_chr]-src_binid)])
            # tgt_binid=int(np.floor((read_data[3]-1)/binsize))

            # out[gene_ID,binid]+=1
    
    p.terminate()
    del p


    for k, cout in enumerate(out):
        cout=cout*(1000000.0/binsize)*(1000000.0/NN[k])

    return out, NN


def make_heatmap(pairs, chr_len_dict, chr_to_ID, gene_to_ID, tx_dict, bygene=True, binsize=1000000, skip_unknown_tx=True, qmin_rna=0, qmin_dna=10, nmax=0, noambiv=False): #dict of vectors, which are targeted chr, first numer is # match,
# second is # of dna match'
    if noambiv:
        cmd="cat "+pairs+" | awk '($25==1 && $8>"+str(qmin_rna)+" && $9>"+str(qmin_dna)+"){print $18, $4, substr($2,3), $5, $3}'"
    else:
        cmd = "cat "+pairs+" | awk '($8>"+str(qmin_rna)+" && $9>"+str(qmin_dna)+"){print $18, $4, substr($2,3), $5, $3}'"
    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
    nok=0
    nokok=0
    ngenes=max([v for _, v in gene_to_ID.items()])+1
    

    totalchrlen=np.sum([v for _, v in chr_len_dict.items()])
    nchr=max([v for _, v in chr_to_ID.items()])+1
    ID_to_chr={v:k for k, v in chr_to_ID.items()}
    chr_len_vec=np.zeros(nchr)
    for i in range(nchr):
        chr_len_vec[i]=chr_len_dict.get(ID_to_chr[i],0)
    
    chr_offset_vec=np.hstack([0,np.cumsum(chr_len_vec)[0:-1]])
    
    nbins_bychr=[int(np.floor(((clen)*1.0)/binsize)) for clen in chr_len_vec]
    nbins=np.sum(nbins_bychr) #int((totalchrlen*1.0)/binsize)+1
    bins_offset=np.hstack([0,np.cumsum(nbins_bychr)[0:-1]])

    out=np.zeros((ngenes,nbins))
    out_chr=np.zeros((ngenes,nchr)) #Ninbounds, Noffchrtgt
    

    for _, line in enumerate(p.stdout):
        if nmax>0 and nok>nmax:
            break
        read_data=line.strip().split(" ")
        T_id=read_data[0]
        T_annot=tx_dict.get(T_id,[])
        # tgt_chr_offset=chr_off_dict.get(read_data[1],-1)
        # src_chr_offset=chr_off_dict.get(read_data[2],-1)
        tgt_chr=chr_to_ID.get(read_data[1],-1)
        src_chr=chr_to_ID.get(read_data[2],-1)
        if bygene:
            Tid=T_annot[0] if len(T_annot)>0 else "*"
        else:
            Tid="*" if (skip_unknown_tx and len(T_annot)<2) else T_id

        gene_ID=gene_to_ID.get(Tid,-1)
        # print(chr_offset_vec)
        if tgt_chr>-1 and src_chr>-1:
            nok+=1
        if gene_ID>-1 and tgt_chr>-1 and src_chr>-1: #only if this is good stuff
            nokok+=1
           
            # print(tgt_chr)
            out_chr[gene_ID,tgt_chr]+=1
            
            # binid=int((int(read_data[3])+chr_offset_vec[tgt_chr]-1)/binsize)
            binid_relative=int(np.floor((int(read_data[3])*1.0)/binsize))
            if binid_relative<nbins_bychr[tgt_chr]:
                nokok+=1
                binid=binid_relative+bins_offset[tgt_chr]
            
                out[gene_ID,binid]+=1


    # we want all the counts in counts per million bp per million input read
    out=out*(1000000.0/binsize)*(1000000.0/nokok)
    
    p.terminate()
    del p

    if not chr_len_vec is None:
        out_chr=out_chr/(chr_len_vec/1000000.0)/(nok/1000000.0)
    
        # out_chr_df=pd.DataFrame(out_chr*(1000000.0/nok),index=classes_lut,columns=)
        # same_others_df=pd.DataFrame(out_nsame_others*(1000000.0/nok),index=classes_lut,columns=(['N_same_chr','N_other_chr']))
        # # so far all offchr_df are "per mb", normalize per million input contacts

    # offchr_df.sort_index(axis=1, inplace=True)
    return out, out_chr, nok, nokok




def blocks_to_contact_map(blocks, genes_df, g2keep=None, renormalize=None, sampling_dist_allchr=None, col_chr_blocks=None, merge=True):
    if renormalize is None:
        renormalize=1.0
    if isinstance(renormalize,list):
        renormalize=np.array(renormalize)
    if not isinstance(renormalize,np.ndarray):
        renormalize=renormalize*np.ones(len(blocks))

    # print(renormalize)

    # sampling_dist=[1.0]*len(blocks) #*np.ones(len(blocks))
    # if not sampling_dist_allchr is None:
        
    #     for k in range(len(blocks)):
    #         sampling_dist[k]=sampling_dist_allchr[col_chr_blocks[k]:col_chr_blocks[k+1]].reshape(1, blocks[k].shape[1])
    #         sampling_dist[k]=sampling_dist[k]/np.sum(sampling_dist[k])
    # else:
    #     sampling_dist=[1.0]*len(blocks)
   
    if sampling_dist_allchr is None:
        if g2keep is None:
            blocks_df=[pd.DataFrame((blocks[k]*renormalize[k]), index=genes_df.loc[genes_df['chr_id']==k].index) for k in range(len(blocks))]
            if merge:
                contact_map=pd.DataFrame(linalg.block_diag(*[e.values for e in blocks_df]), index=genes_df.index)
            else:
                contact_map=[e.values for e in blocks_df]
        else:
            gdf=genes_df.loc[g2keep]
            blocks_df=[select_genes(pd.DataFrame((blocks[k]*renormalize[k]), index=genes_df.loc[genes_df['chr_id']==k].index), gdf.loc[genes_df['chr_id']==k].index) for k in range(len(blocks))]
            if merge:
                contact_map=pd.DataFrame(linalg.block_diag(*[e.values for e in blocks_df]), index=gdf.index)
            else:
                contact_map=[e.values for e in blocks_df]
        return contact_map
    else:
        sampling_dist=[1.0]*len(blocks)
        for k in range(len(blocks)):
            sampling_dist[k]=sampling_dist_allchr[col_chr_blocks[k]:col_chr_blocks[k+1]].reshape(1, blocks[k].shape[1])
            sampling_dist[k]=sampling_dist[k]/np.sum(sampling_dist[k])
        if g2keep is None:
            blocks_df=[[]]*len(blocks)
            for k in range(len(blocks)):
                Nk=np.sum(blocks[k], axis=1)*renormalize[k]
                vk=blocks[k]*(renormalize[k]*sampling_dist[k])
                normfactor=(Nk/(np.sum(vk, axis=1))).reshape((blocks[k].shape[0],1))
       
                blocks_df[k]=pd.DataFrame(vk*normfactor, index=genes_df.loc[genes_df['chr_id']==k].index)
            # blocks_df=[pd.DataFrame((blocks[k]*(renormalize[k]*sampling_dist[k])/, index=genes_df.loc[genes_df['chr_id']==k].index) for k in range(len(blocks))]
            if merge:
                contact_map=pd.DataFrame(linalg.block_diag(*[e.values for e in blocks_df]), index=genes_df.index)
            else:
                contact_map=[e.values for e in blocks_df]
        else:
            gdf=genes_df.loc[g2keep]
            blocks_df=[[]]*len(blocks)
            for k in range(len(blocks)):
                Nk=np.sum(blocks[k], axis=1)*renormalize[k]
                vk=blocks[k]*(renormalize[k]*sampling_dist[k])
                normfactor=(Nk/(np.sum(vk, axis=1))).reshape((blocks[k].shape[0],1))
                
                blocks_df[k]=select_genes(pd.DataFrame(vk*normfactor, index=genes_df.loc[genes_df['chr_id']==k].index),gdf.loc[genes_df['chr_id']==k].index)

            # blocks_df=[select_genes(pd.DataFrame((blocks[k]*renormalize[k])*sampling_dist[k], index=genes_df.loc[genes_df['chr_id']==k].index), gdf.loc[genes_df['chr_id']==k].index) for k in range(len(blocks))]
            if merge:
                contact_map=pd.DataFrame(linalg.block_diag(*[e.values for e in blocks_df]), index=gdf.index)
            else:
                contact_map=[e.values for e in blocks_df]
        return contact_map
    
def rescale(m, rng=(1,1e6)): #this allows us to avoid zeros and such
    minv=np.min(m.values)
    maxv=np.max(m.values)
    mr=(m-minv)/(maxv-minv)*(rng[1]-rng[0])+rng[0]
    return mr


def to_cpm(m):
#     e=m.sum(axis=1)
#     N=e.sum()
    N=np.sum(m.values)
    mc=m/N*1000000
    return mc

# def normalize_bygenes_locus(m, genes_df, minval=0, rng=1e6):
#     gene_src=genes_df['source_pos'].loc[m.index].values
#     gene_src_val=m.values[np.arange(len(m)),gene_src]

#     i2keep=gene_src_val>minval
#     norm=gene_src_val[i2keep]
#     mnorm=m.loc[i2keep].div(norm,axis=0)*rng

# #     minv=np.min(mnorm.values)
# #     maxv=np.max(mnorm.values)
# #     mnorm=(mnorm-minv)/(maxv-minv)*(rng[1]-rng[0])+rng[0]

#     return mnorm


def normalize_bygenes_locus(m, genes_df, minval=0, rng=1e6):
    gene_src=genes_df['source_pos'].loc[m.index].values
    gene_src_val=m.values[np.arange(len(m)),gene_src]

    i2keep=gene_src_val>minval
    norm=gene_src_val[i2keep]
    mnorm=m.loc[i2keep].div(norm,axis=0)*rng

#     minv=np.min(mnorm.values)
#     maxv=np.max(mnorm.values)
#     mnorm=(mnorm-minv)/(maxv-minv)*(rng[1]-rng[0])+rng[0]

    return mnorm

def normalize_bygenes(m, norm_df):
    ix=pd.concat([m.iloc[:,0], norm_df], join='inner').index

    norm_val=norm_df.loc[ix].values

    mnorm=m.loc[ix].div(norm_val, axis=1)
    return mnorm

def select_genes(m, g2keep, keeporder=True):
    if keeporder:
        msub=m.loc[m.index.isin(g2keep)]
    else:
        msub=m.loc[g2keep]
    return msub

def plot_map(m, ax=None, annotations=None, plotsrc=False, label_Y=0, label_chr=False, bins_ix_chrdict=None, genes_df=None, cmap='YlGnBu',vmin=2, vmax=6, chrRNA=None, chrDNA=None):
    if not chrRNA is None:
        if type(chrRNA) is str:
            m=select_genes(m, genes_df.loc[genes_df['chr']==chrRNA].index)
        else:
            m=select_genes(m, genes_df.loc[genes_df['chr'].isin(chrRNA)].index)
            
    if not chrDNA is None:
        if type(chrDNA) is str:
            m=m.iloc[:,bins_ix_chrdict[chrDNA][0]:bins_ix_chrdict[chrDNA][1]].copy()
        else:
            m=m.iloc[:,bins_ix_chrdict[chrDNA[0]][0]:bins_ix_chrdict[chrDNA[-1]][1]].copy()
    #left, right, bottom, top
    xmin=m.columns[0]-0.5
    xmax=m.columns[0]+m.shape[1]-0.5
    ymin=-0.5
    ymax=m.shape[0]-0.5
    im = ax.imshow(np.log10(m.values), cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto', extent=(m.columns[0]-0.5, m.columns[0]+m.shape[1]-0.5, m.shape[0]-0.5, -0.5))
    if plotsrc:
        ax.plot(genes_df['source_pos'].loc[m.index].values,np.arange(len(m)), linestyle=':', color='k', linewidth=1)
    
    if label_chr:
        bins_ix_low=np.array([v[0] for v in bins_ix_chrdict.values()])-0.5
        bins_ix_low.sort()
#         bins_ix_high=np.array([v[1] for v in bins_ix_chrdict.values()])
#         bins_ix_high.sort()
#         xtk=np.hstack bins_ix_low
#         px=(bins_ix_low+bins_ix_high)/2
        tks=bins_ix_low[(bins_ix_low>=xmin) & (bins_ix_low<=xmax)]
        ax.set_xticks(tks)
        ax.set_xticklabels([])
        for k, v in bins_ix_chrdict.items():
            x_label=np.mean(v)
            if ((x_label>=xmin) and (x_label<=xmax)):
                ax.annotate(k,xy=(x_label+1.5,ymax+1), horizontalalignment='center', verticalalignment='top', rotation=90, fontsize=8)
        ax.set_ylim([ymax+1, ymin])
    
                
    if label_Y>0:
        
        pos=pd.DataFrame(np.array(range(len(m))), index=m.index, columns=['pos'])
        
        mygenes=pd.concat([pos,genes_df[['chr']]], axis=1, join='inner').groupby('chr').apply(lambda x: x.sort_values('pos').head(1)).sort_values('pos')
#         print(mygenes)
        mygenes_chr=mygenes['chr'].values
        pos_chr=mygenes['pos'].values
        pos_chr2=np.hstack([pos_chr,m.shape[0]])
#         print(pos_chr)
        pos_label=(pos_chr2[1:]+pos_chr2[0:-1])/2
        
#         for ii, c in enumerate(mygenes_chr):
#             ax.annotate(c, xy=(xmax+1,pos_label[ii]), horizontalalignment='left', verticalalignment='center')
    
        for ii, c in enumerate(pos_chr):
            ax.annotate(mygenes_chr[ii], xy=(xmax+1,c), horizontalalignment='left', verticalalignment='top', fontsize=8)
            if label_Y==2:
                ax.plot([xmin, xmax],[c-0.5,c-0.5],color='k',linestyle=':')
        
    ax.set_xlim([xmin, xmax+1])

    if not annotations is None:
        pos=pd.DataFrame(np.array(range(len(m))), index=m.index, columns=['pos'])
        pos_genes=pos['pos'].loc[pos.index.isin(annotations.index)]

        mygenes=genes_df[['name','type']].loc[pos_genes.index].copy()
        labs=[r['name']+", "+r['type'] for _, r in mygenes.iterrows()]

        
        ax.set_yticks(pos_genes.values)
        ax.set_yticklabels(labs)
#     im = ax.imshow(np.log10(m.values), cmap='YlGnBu',vmin=2, vmax=6, aspect='auto')


    if label_chr>1:
        for kkk in tks:
            ax.plot([kkk,kkk], ax.get_ylim(), color='k',linestyle='-.')

def submap(m, genes_df, bins_ix_chrdict, chrRNA=None, chrDNA=None):
    m2=m
    if not chrRNA is None:
        if type(chrRNA) is str:
            m2=select_genes(m, genes_df.loc[genes_df['chr']==chrRNA].index)
        else:
            m2=select_genes(m, genes_df.loc[genes_df['chr'].isin(chrRNA)].index)
            
    if not chrDNA is None:
        if type(chrDNA) is str:
            m2=m2.iloc[:,bins_ix_chrdict[chrDNA][0]:bins_ix_chrdict[chrDNA][1]].copy()
        else:
            m2=m2.iloc[:,bins_ix_chrdict[chrDNA[0]][0]:bins_ix_chrdict[chrDNA[-1]][1]].copy()
    return m2
            
    # if not chrDNA is None:
    #     if type(chrDNA) is str:
    #         m2=m.loc[:,bins_ix_chrdict[chrDNA][0]:bins_ix_chrdict[chrDNA][1]]
    #     else:
    #         m2=m.loc[:,bins_ix_chrdict[chrDNA[0]][0]:bins_ix_chrdict[chrDNA[-1]][1]]
    # else:
    #     m2=m

        
    # if not chrRNA is None:
    #     if type(chrRNA) is str:
    #         m2=select_genes(m2, genes_df.loc[genes_df['chr']==chrRNA].index)
    #     else:
    #         m2=select_genes(m2, genes_df.loc[genes_df['chr'].isin(chrRNA)].index)
            
    # return m2

def block_summary(m, genes_df, bins_ix_chrdict):
    
    bins_ix_low=np.array([v[0] for v in bins_ix_chrdict.values()])
    bins_ix_low.sort()
    cols_chr_id=[bisect.bisect_right(bins_ix_low,v)-1 for v in m.columns]
    cols_chr_df=pd.DataFrame(cols_chr_id, columns=['chrID'])
    cols_chr_df['iiloc']=np.arange(len(cols_chr_df))
    cols_chr_df2=cols_chr_df.groupby('chrID').aggregate(['first','last'])
#     cols_chr_df2.columns=cols_chr_df2.columns.droplevel(0)
#     cols_chr_df2['n']=cols_chr_df['last']-cols_chr_df['first']+1
    cols_chr_M=cols_chr_df2.values
    
    
    chr_to_id={k:i for i,k in enumerate(list(bins_ix_chrdict.keys()))}
    
    genes_src=genes_df['source_pos'].loc[m.index].values
    genes_chr=genes_df['chr'].loc[m.index].values
    genes_chr2=[chr_to_id[v] for v in genes_chr]
    
    x=np.zeros((len(m),6)) #total, total src bin, total src chr, total other, Nrc, N other
    
    i=-1
    for ix, r in m.iterrows():
        i+=1
        cchr=genes_chr2[i]
        
        x[i,0]=r.sum()
        x[i,1]=m.at[ix,genes_src[i]]
        x[i,2]=np.sum(m.values[i,cols_chr_M[cchr,0]:(cols_chr_M[cchr,1]+1)])
        x[i,3]=x[i,0]-x[i,2]
        x[i,4]=(cols_chr_M[cchr,1]+1)-cols_chr_M[cchr,0]
        x[i,5]=len(r)-x[i,4]
    
    bs=pd.DataFrame(x,index=m.index, columns=['e','e_locus','e_src','e_other','N_src','N_other'])
    bs['chr']=genes_chr
    return bs
    
    
def sort_by_chr(bs, by_value, ascending=False):
    ixs=bs.reset_index().groupby('chr').apply(lambda x: x.sort_values(by_value, ascending=ascending))['index'].values
    return ixs

def summarize_chrmatrix_(m, genes_df, chr_len_vec): #assumes the genes are in the same order
    
    S=np.zeros((len(genes_df),4)) #e src, e others, e total, t-score, entropy total, entropy others
    m_sum=m.sum(axis=1)
    chr_norm=chr_len_vec.sum()-chr_len_vec
    for k in range(len(genes_df)):
        chrid=int(genes_df['chr_id'].iloc[k])
        v=m[k,:]*chr_len_vec/1000000.0
        S[k, 2]=v.sum()
        S[k,0]=v[chrid]
        S[k,1]=S[k, 2]-S[k, 0]
        S[k,3]=(S[k, 1]/chr_norm[chrid])/(S[k, 0]/chr_len_vec[chrid])
    S1=pd.DataFrame(S, index=genes_df.index)
    S1.columns=['e_cis','e_trans','e_total','t_score']
    Sdf=pd.concat([genes_df[['type','name','chr_id']],S1], axis=1, join='outer').loc[genes_df.index].copy()
#     Sdf.col
#         S[k,4]=entropy(m[k,:])
#         S[k,5]=
#         entro[k,0]=entropy(m.values[k,:])
#         entro[k,1]=entropy(m.values[k,bins_ix_chrdict[chrs[k]][0]:bins_ix_chrdict[chrs[k]][1]])
#     s=pd.DataFrame(entro, index=ix, columns=['entro_all','entro_src'])
#     s['e']=m.sum(axis=1)
#     s['n']=m.std(axis=1)/s['e']
    
#     tot_e=s['e'].sum(axis=0)
#     s['e']=s['e']/tot_e*1000000
    return Sdf

def read_csv_data(filename, binsize, colnum, col_chr_blocks, chr_to_ID, chr_len_vec):
    # totalchrlen=np.sum([v for _, v in chr_len_dict.items()])
    # nchr=max([v for _, v in chr_to_ID.items()])+1
    # ID_to_chr={v:k for k, v in chr_to_ID.items()}
    # chr_len_vec=np.zeros(nchr)
    # for i in range(nchr):
    #     chr_len_vec[i]=chr_len_dict.get(ID_to_chr[i],0)
    
    # chr_offset_vec=np.hstack([0,np.cumsum(chr_len_vec)[0:-1]])
    
    nbins_bychr=[int(np.floor(((clen)*1.0)/binsize)) for clen in chr_len_vec]
    nbins=np.sum(nbins_bychr) #int((totalchrlen*1.0)/binsize)+1
    # # bins_offset=np.hstack([0,np.cumsum(nbins_bychr)[0:-1]])

    # nbins=col_chr_blocks[-1]+1

    out=np.zeros(nbins, float)
    out_n=np.zeros(nbins, int)
    with open(filename, 'r') as fh:
        for line in fh:
            read_data=line.strip().split("\t")
            
            tgt_chr=chr_to_ID.get(read_data[0], -1)
            binid_relative=int(np.floor((int(read_data[1])*1.0)/binsize))
            if (binid_relative<nbins_bychr[tgt_chr]) and (tgt_chr>-1):

                binid=binid_relative+col_chr_blocks[tgt_chr]
            
                if not(read_data[colnum]=='NA'):
                    out[binid]+=float(read_data[colnum])
                    out_n[binid]+=1
    out=out/np.maximum(out_n,1)
    return out


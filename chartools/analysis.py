import subprocess
import numpy as np
import pandas as pd
import csv
from scipy import stats
def make_tx_dict_fromCSV(annot_file):
    annot={}
    with open(annot_file,'r') as f:
        myreader=csv.reader(f, delimiter='\t')
        for row in myreader:
            annot[row[0]]=row[1:]
    return annot

def make_chr_length_dict_fromCSV(chr_file, delim='\t'):
    chrdict = {}
    with open(chr_file, 'r') as f:
        myreader = csv.reader(f, delimiter=delim)
        for _, row in enumerate(myreader):
            chrdict[row[0]] = int(row[1])
    return chrdict

def load_index_fromCSV(file):
    index={}
    i=0
    with open(file,'r') as f:
        myreader=csv.reader(f, delimiter='\t')
        for row in myreader:
            index[row[0]]=i
            i+=1
    return index

def tx_binding_table(pairs, chr_dict, tx_dict, bygene=True, skip_unknown_tx=True, chr_len_vec=None, qmin_rna=0, qmin_dna=10, nmax=0, noambiv=False): #dict of vectors, which are targeted chr, first number is # match,
# second is # of dna match, also returns another dataframe with mean distace for given RNA, std

        #tx_dict: tx_id (key), gene_id, promoter_pos, tx_total_lenght

        #NEW OUTPUT 0=readID, 1=chrR, 2=posR, 3=chrD, 4=posD, 5=strandR, 6=strandD, 
        # 7=QR, 8=QD, 9=flagR, 10=flagD, 11=gapR, 12=gapD, 13=query alignment length RNA, 14=query aligment length DNA, 
        # 15=nR, 16=nD, 17=annot_ref(ENST), 18=annot_pos, 19=annot_name(ex: gapdh), 20=annot_type, 21=gS (#of genomic aligmnet with compatible annotation RNA side), 22=aS (#of annotations), 23=ai, 24=aI, 25=ENSG
 
    # annot id, d pos, tgt chr, useless, src chr, readid --> annotid, tgt chr, src chr, tgp pos, src pos
        if noambiv:
            cmd = "cat "+pairs+" | awk '($25==1 && $8>"+str(qmin_rna)+" && $9>"+str(qmin_dna)+"){print $18, $4, substr($2,3), $5, $3}'"
        else:
            cmd = "cat "+pairs+" | awk '($8>"+str(qmin_rna)+" && $9>"+str(qmin_dna)+"){print $18, $4, substr($2,3), $5, $3}'"
        m={}
        nchr=len(chr_dict)
        print(cmd)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
        nok=0
        nsame=0
        # read all reads in the file unitl reach end of file
        for _, line in enumerate(p.stdout):
            if nmax>0 and nok>nmax:
                break
            read_data=line.strip().split(" ")
            T_id=read_data[0]
            T_annot=tx_dict.get(T_id,[])
            tgt_chr=chr_dict.get(read_data[1],-1)
            src_chr=chr_dict.get(read_data[2],-1)
            

            if bygene:
                Tid=T_annot[0] if len(T_annot)>0 else "*"
            else:
                Tid="*" if (skip_unknown_tx and len(T_annot)<2) else T_id

            

            if len(Tid)>1 and tgt_chr>-1 and src_chr>-1 and (len(T_annot)>0): #only if this is good stuff
                nok+=1
                if Tid in m:
                    x=m[Tid]
                    x[tgt_chr+1]+=1
                    if tgt_chr==src_chr:
                        nsame+=1
                        flight=(int(read_data[3])-int(read_data[4]))/10000.0
                        flight2=(int(read_data[3])-int(T_annot[2]))/10000.0
                        x[nchr+1]+=flight
                        x[nchr+2]+=abs(flight)
                        x[nchr+3]+=flight*flight
                        x[nchr+4]+=flight2
                        x[nchr+5]+=abs(flight2)
                        x[nchr+6]+=flight2*flight2
                else:
                    x=np.zeros(nchr+1+6)
                    x[tgt_chr+1]=1
                    x[0]=src_chr
                    if tgt_chr==src_chr:
                        nsame+=1
                        flight=(int(read_data[3])-int(read_data[4]))/10000.0
                        flight2=(int(read_data[3])-int(T_annot[2]))/10000.0
                    else:
                        flight=0
                        flight2=0

                    x[nchr+1]=flight
                    x[nchr+2]=abs(flight)
                    x[nchr+3]=flight*flight
                    x[nchr+4]=flight2
                    x[nchr+5]=abs(flight2)
                    x[nchr+6]=flight2*flight2
                    m[Tid]=x

        p.terminate()
        del p
        
        
        inv_chr = {v: k for k, v in chr_dict.items()} 
        df=pd.DataFrame.from_dict(m,orient='index')     
        df.columns=['src']+[inv_chr[i] for i in range(len(inv_chr))]+['D_sum','Dabs_sum','D2_sum','P_sum','Pabs_sum','P2_sum']
        
        if not chr_len_vec is None:
            norm_vec=1000000.0/chr_len_vec #counts per mbp
            df.values[:,1:-6]=df.values[:,1:-6]*norm_vec[None,:]

        x=df.values
        # x[:,1:(nchr+1)].sum(axis=1)-
        n_samechr=np.maximum(1,x[np.arange(x.shape[0]),x[:,0].astype(int)+1])
        # print(n_samechr.shape)
        
        x[:,(nchr+1):]=(x[:,(nchr+1):]/n_samechr[:,np.newaxis])
        x[:,nchr+3]=np.sqrt(x[:,nchr+3]-x[:,nchr+1]*x[:,nchr+1])
        x[:,nchr+6]=np.sqrt(x[:,nchr+6]-x[:,nchr+4]*x[:,nchr+4])
        x[:,(nchr+1):]=x[:,(nchr+1):]*10
        
        if bygene:
            df.index.names=['ENSG']
        else:
            df.index.names=['ENST']
        return df

def _binding_statistics_aux(v, chr_len_vec=None):
    #on self, mean on other, std on other, entropy total, entropy other
    if not chr_len_vec is None:
        norm_vec=1000000.0/chr_len_vec #counts per mbp
        # df.values[:,1:-6]=df.values[:,1:-6]*norm_vec[None,:]
        norm_total=1000000.0/chr_len_vec.sum()
    else:
        norm_total=1

    nr, nc = v.shape
    s_all=np.log(nc-1)
    s_others=np.log(nc-2)
    curr_src=0
    mask_others=np.ones(nc-1,dtype=bool)

    sc=np.zeros((nr,6))
    for i in range(nr):
        e_chr= v[i,1:]#signal by chromosome
        mask_others[curr_src]=True
        
        curr_src=int(v[i,0])

        mask_others[curr_src]=False
        sc[i,0]=np.sum(e_chr)*norm_total #total counts of this tx per M targets
        
        if not chr_len_vec is None:
            e_chr=e_chr*norm_vec

        if sc[i,0]>0:    
            sc[i,1]=e_chr[curr_src]/sc[i,0] #counts on chr i /M target over counts all chr per M target = density on self over total density
            sc[i,2]=np.mean(e_chr[mask_others])/sc[i,0]
            sc[i,3]=np.std(e_chr[mask_others])/sc[i,0]
            sc[i,5]=stats.entropy(e_chr[mask_others])/s_others
            sc[i,4]=stats.entropy(e_chr)/s_all

    s=sc[:,0].sum() #total count per M target
    sc[:,0]=sc[:,0]*1000000.0/s #total counts per M targets per M reads
    return sc

    #e_sother_mean * e_total = RPK other / RPK avg * (RPK avg/(total RPK avg)*1e6)
    #RPK other / total RPK avg *1e6
    # density on other / total densit 


def binding_statistics(data, chr_len_vec=None):
    s=pd.DataFrame(_binding_statistics_aux(data.values[:,0:-6], chr_len_vec=chr_len_vec),columns=['e_total','e_src','e_others_mean','e_others_std','S_total','S_other'])
    s.index=data.index
    s['noise_ratio']=s['e_others_mean']/s['e_src']
    # for ix, col_ix in enumerate(columns_ixs):
    #     fun2apply=lambda x: "|".join(sorted(set([annot_dict.get(geneid,["*"])[col_ix] for geneid in x.split("|")])))
    #     s.insert(ix,columns[ix],list(map(fun2apply,s.index)))
    return s


def intersect_bed(pairs, bed, tx_dict=None, chr_index=None, genes_classes_dict=None, bygene=True, normalize=True, chr_len_vec=None, qmin_rna=0, qmin_dna=10, nmax=0, noambiv=False, skip_unknown_tx=True, maskout=0):
    index, chr_index = bed_to_index(bed, chr_index=chr_index, ignore_other_chr=True)
    
    offsets=np.hstack([0,1,1+np.cumsum([int(len(ix)/2)+1 for ix in index])])
    ws=None
    if normalize:
        if chr_len_vec is None:
            chr_len_vec=[0]*len(index)

        ws=[np.zeros(0) for _ in range(len(index))] #only for the bins not for the pos
        for i, ix in enumerate(index):
            if len(ix)>0:
                ws[i]=np.diff(ix)[0::2]
                tot=ws[i].sum()
                ws[i]=np.hstack([chr_len_vec[i]-tot,ws[i]])
            else:
                ws[i]=np.array([chr_len_vec[i]])
        wsall=np.hstack(ws)
    nbins=offsets[-1]
    print(offsets)
    print(nbins)
    print(len(offsets))

    nok=0
    if noambiv:
        cmd="cat "+pairs+" | awk '($25==1 && $8>"+str(qmin_rna)+" && $9>"+str(qmin_dna)+"){print $18, $4, substr($2,3), $5, $3}'"
    else:
        cmd = "cat "+pairs+" | awk '($8>"+str(qmin_rna)+" && $9>"+str(qmin_dna)+"){print $18, $4, substr($2,3), $5, $3}'"
    
    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
    nok=0
    nok2=0
    nok3=0
    
    out={}
    
        # read all reads in the file until reach end of file
    for _, line in enumerate(p.stdout):
        if nmax>0 and nok>nmax:
            break
        read_data=line.strip().split(" ")
        T_id=read_data[0]
        T_annot=tx_dict.get(T_id,[])
        tgt_chr=chr_index.get(read_data[1],-1)
        src_chr=chr_index.get(read_data[2],-1)
        if bygene:
            Tid=T_annot[0] if len(T_annot)>0 else "*"
        else:
            Tid="*" if (skip_unknown_tx and len(T_annot)<2) else T_id

        


        # if not genes_classes_dict is None:
        #     Tid=genes_classes_dict.get(Tid,"*")

        nok+=1
        # print(T_id)
        # print(Tid)
        if not(Tid=="*") and tgt_chr>-1: #only if this is good stuff
            nok2+=1
            distrib=out.get(Tid,None)
            new=False
            if distrib is None:
                new=True
                distrib=np.zeros(nbins)
         
                if (len(T_annot)>2) and (src_chr>-1):
                    promoter_pos=int(T_annot[2])
                    bin_src=np.searchsorted(index[src_chr],promoter_pos)
                    if bin_src%2==0:
                        print('yep')
                        distrib[0]=offsets[src_chr+1]
                    else:
                        
                        distrib[0]=int((bin_src+1)/2)+offsets[src_chr+1]
                else:
                    distrib[0]=-1

            pos=np.searchsorted(index[tgt_chr],int(read_data[3]))
     
            ##
            mask=False
            if maskout>0 and (tgt_chr==src_chr):
                flight=(int(read_data[3])-int(read_data[4]))
                mask=True if (abs(flight)<maskout) else False
            

            if not mask:
                nok3+=1
                if pos%2==0:
                    distrib[offsets[tgt_chr+1]]+=1
                else:
                    distrib[int((pos+1)/2)+offsets[tgt_chr+1]]+=1

                if new:
                    out[Tid]=distrib

    if normalize:
        for k, v in out.items():
            out[k]=np.hstack([v[0],v[1:]*(1000000.0/wsall)*(1000000.0/nok)])

    # we want all the counts in counts per million bp per million input read
   
    p.terminate()
    del p


    return out, nok, nok2, nok3, ws, offsets

def summarize_bed_intersection(out, offsets):
    summary_mean={}
    summary_std={}
    summary_entropy={}
    summary_src={}
    # summary_esrc={}
    for k, v in out.items():
        summary=np.zeros((3,len(offsets)-2))
        for i in range(len(offsets)-2):
            # print(offsets[i+2])
            # print(v[(offsets[i+1]+1):offsets[i+2]])
            summary[0,i]=np.mean(v[(offsets[i+1]+1):offsets[i+2]])
            summary[1,i]=np.std(v[(offsets[i+1]+1):offsets[i+2]])/np.sqrt((offsets[i+2]-offsets[i+1]+1)) #stderr
            nelem=np.sum(v[(offsets[i+1]+1):offsets[i+2]])
            summary[2,i]=stats.entropy(v[(offsets[i+1]+1):offsets[i+2]])/np.log2(nelem)
        src_chr=np.searchsorted(offsets,v[0]+1)-2
        summary_src[k]=np.array([summary[0,src_chr],summary[1,src_chr],summary[2,src_chr]])
        summary[0,src_chr]=np.nan
        summary[1,src_chr]=np.nan
        summary[2,src_chr]=np.nan

        summary_mean[k]=summary[0,:]
        summary_std[k]=summary[1,:]
        summary_entropy[k]=summary[2,:]
        
    summary_mean_df=pd.DataFrame.from_dict(summary_mean, orient='index')
    summary_mean_df.columns=["chr%g_mean"%(i+1) for i in range(len(offsets)-2)]
    summary_std_df=pd.DataFrame.from_dict(summary_std, orient='index')
    summary_std_df.columns=["chr%g_se"%(i+1) for i in range(len(offsets)-2)]
    summary_entropy_df=pd.DataFrame.from_dict(summary_std, orient='index')
    summary_entropy_df.columns=["chr%g_S"%(i+1) for i in range(len(offsets)-2)]
    summary_src_df=pd.DataFrame.from_dict(summary_src, orient='index')
    summary_src_df.columns=["mean","se","S"]
    summary_df=pd.concat([summary_mean_df,summary_entropy_df,summary_std_df,summary_src_df],axis=1,keys=['mean','S','se','src'])
    
    # xx.columns=["chr%g_mean"%(i+1) for i in range(23)]
    return summary_df

def bed_to_index(bed, nmax=-1, chr_index=None, ignore_other_chr=True):
    reader=open(bed,'r')
    
    if chr_index is None:
        chr_index={}
        ignore_other_chr=False
        nchr=-1
        index=[]
    else:
        chr_index2={k:v for k, v in chr_index.items()}
        chr_index=chr_index2
        nchr=max((v for _, v in chr_index.items()))
        index=[[] for _ in range(nchr+1)]
        
    curr_chr=""
    curr_chr_ix=-1
    nok=0
    while (nmax<1 or nok<nmax):
        try:
            line=reader.readline()
        except StopIteration:
            reader.close()
            return index, chr_index
        nok+=1
        read_data=line.strip().split("\t")
        if len(read_data)<3:
            reader.close()
            return index, chr_index
        if not (curr_chr==read_data[0]):
            curr_chr=read_data[0]
            curr_chr_ix=chr_index.get(curr_chr,-1)
            if (curr_chr_ix<0):
                if not(ignore_other_chr):
                    nchr+=1
                    chr_index[curr_chr]=nchr
                    curr_chr_ix=nchr
                    index+=[[]]
#         print(curr_chr_ix)
        if curr_chr_ix>-1:
            # print(curr_chr_ix)
            index[curr_chr_ix]+=[int(read_data[1]),int(read_data[2])]
    for ix in index:
        ix.sort()
    reader.close()

    # offsets=np.hstack([0,np.cumsum([int(len(ix)/2)+1 for ix in index])])

    return index, chr_index

def make_heatmap(pairs, chr_len_dict, chr_to_ID, gene_to_ID, tx_dict, bygene=True, binsize=1000000, skip_unknown_tx=True, qmin_rna=0, qmin_dna=10, nmax=0, noambiv=False): #dict of vectors, which are targeted chr, first numer is # match,
# second is # of dna match'
    if noambiv:
        cmd="cat "+pairs+" | awk '($25==1 && $8>"+str(qmin_rna)+" && $9>"+str(qmin_dna)+"){print $18, $4, substr($2,3), $5, $3}'"
    else:
        cmd = "cat "+pairs+" | awk '($8>"+str(qmin_rna)+" && $9>"+str(qmin_dna)+"){print $18, $4, substr($2,3), $5, $3}'"
    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
    nok=0

    ngenes=max([v for _, v in gene_to_ID.items()])+1
    

    totalchrlen=np.sum([v for _, v in chr_len_dict.items()])
    nchr=max([v for _, v in chr_to_ID.items()])+1
    ID_to_chr={v:k for k, v in chr_to_ID.items()}
    chr_len_vec=np.zeros(nchr)
    for i in range(nchr):
        chr_len_vec[i]=chr_len_dict.get(ID_to_chr[i],0)
    chr_offset_vec=np.hstack([0,np.cumsum(chr_len_vec)[0:-1]])
    nbins=int((totalchrlen*1.0)/binsize)+1

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
           
            # print(tgt_chr)
            out_chr[gene_ID,tgt_chr]+=1
            
            binid=int((int(read_data[3])+chr_offset_vec[tgt_chr]-1)/binsize)

            out[gene_ID,binid]+=1


    # we want all the counts in counts per million bp per million input read
    out=out*(1000000.0/binsize)*(1000000.0/nok)
    
    p.terminate()
    del p

    if not chr_len_vec is None:
        out_chr=out_chr/(chr_len_vec/1000000.0)/(nok/1000000.0)
    
        # out_chr_df=pd.DataFrame(out_chr*(1000000.0/nok),index=classes_lut,columns=)
        # same_others_df=pd.DataFrame(out_nsame_others*(1000000.0/nok),index=classes_lut,columns=(['N_same_chr','N_other_chr']))
        # # so far all offchr_df are "per mb", normalize per million input contacts

    # offchr_df.sort_index(axis=1, inplace=True)
    return out, out_chr

#distribution by genename: for each gene in dictionary, target distance for that gene
#bybedfile/source: for each interval in bed file, source-target dist when source is in there, or befile/target when target is in there
#
def fly_statistics(pairs, chr_dict, tx_dict, genes_classes_dict, bygene=True, bins=np.linspace(-1000000,1000000,201), mask_dist=[0, 10000, 30000, 100000, 300000,1000000,3000000, 10000000], from_promoter=False, direction_relative=True, skip_unknown_tx=True, chr_len_vec=None, qmin_rna=0, qmin_dna=10, nmax=0, noambiv=False): #dict of vectors, which are targeted chr, first numer is # match,
# second is # of dna match'
    if noambiv:
        cmd="cat "+pairs+" | awk '($25==1 && $8>"+str(qmin_rna)+" && $9>"+str(qmin_dna)+"){print $18, $4, substr($2,3), $5, $3}'"
    else:
        cmd = "cat "+pairs+" | awk '($8>"+str(qmin_rna)+" && $9>"+str(qmin_dna)+"){print $18, $4, substr($2,3), $5, $3}'"
    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
    nok=0

    print(cmd)
    #create an index to map class to ouput row
    classes_index={v:ix for ix, v in enumerate(set(vv for _,vv in genes_classes_dict.items()))} #protein coding: 0, etcx...
    genes_dict={k:classes_index.get(v,-1) for k, v in genes_classes_dict.items()}

    classes_index_inv={ix:k for k,ix in classes_index.items()}
    classes_lut=[classes_index_inv.get(ix) for ix in range(len(classes_index_inv))]

    nclasses=len(classes_index)
    nbins=len(bins)+1
    
    nmasks=len(mask_dist)
    # w=bins[-1]-bins[0]
    dw=np.hstack([1,np.diff(bins),1])
    if not chr_len_vec is None:
        
        # w_mb=w/1000000.0
        # n2add_in=1/w_mb
        #10kb, 50kb, 100kb, 500kb, 1Mb, 5Mb, 10Mb
        n2add_total=1000000.0/(chr_len_vec.sum())
        n2add_other=1000000.0/(chr_len_vec.sum()-chr_len_vec)
        # n2add_same=[1000000.0/(chr_len_vec-2*wmask) for wmask in mask_dist] 
    else:
        n2add_total=1
        n2add_other=np.ones(len(chr_dict))
        # n2add_same=[np.ones(len(chr_dict)) for _ in range(nmasks)]

    out=np.zeros((nclasses,nbins),int)
    # out_promoter=np.zeros((nclasses,nbins),int)
    out_offchr=np.zeros((nclasses,2+nmasks)) #Ninbounds, Noffchrtgt
    out_nsame_others=np.zeros((nclasses,2), int)
        # read all reads in the file until reach end of file
    for _, line in enumerate(p.stdout):
        if nmax>0 and nok>nmax:
            break
        read_data=line.strip().split(" ")
        T_id=read_data[0]
        T_annot=tx_dict.get(T_id,[])
        tgt_chr=chr_dict.get(read_data[1],-1)
        src_chr=chr_dict.get(read_data[2],-1)

        if bygene:
            Tid=T_annot[0] if len(T_annot)>0 else "*"
        else:
            Tid="*" if (skip_unknown_tx and len(T_annot)<2) else T_id

        T_class=genes_dict.get(Tid,-1)
        if tgt_chr>-1 and src_chr>-1:
            nok+=1
        if T_class>-1 and tgt_chr>-1 and src_chr>-1: #only if this is good stuff
           
            out_offchr[T_class,-1]+=n2add_total
            if tgt_chr==src_chr:
                out_nsame_others[T_class,0]+=1
                if from_promoter:
                    src_pos=int(T_annot[2])
                    fly=int(read_data[3])-int(T_annot[2])
                    # n2add_this=1000000.0/(chr_len_vec-2*wmask)
                else:
                    src_pos=int(read_data[4])
                    fly=int(read_data[3])-int(read_data[4])
                if direction_relative and T_annot[1]=="-":
                    fly=-fly

                out[T_class,np.searchsorted(bins,fly)]+=1
                # out_promoter[T_class,np.searchsorted(bins,fly_promoter)]+=1
                for i, wm  in enumerate(mask_dist):
                    if abs(fly)>= wm: #>bins[-1] or fly<bins[0]:
                        n2add_this=1 if (chr_len_vec is None) else 1000000.0/(max(0.0,chr_len_vec[src_chr]-(src_pos+wm))+max(0.0,src_pos-wm))
                        out_offchr[T_class,i]+=n2add_this #n2add_same[i][src_chr]
                # else:
                    # out_offchr[T_class,0]+=n2add_in
            else:
                out_offchr[T_class,-2]+=n2add_other[src_chr]
                out_nsame_others[T_class,1]+=1

    # we want all the counts in counts per million bp per million input read
    out=out*(1000000.0/dw)*(1000000.0/nok)
    p.terminate()
    del p

    offchr_df=pd.DataFrame(out_offchr*(1000000.0/nok),index=classes_lut,columns=(['N_outof_%ikbp'%(wm/1000) for wm in mask_dist]+['N_other_chr','N']))
    same_others_df=pd.DataFrame(out_nsame_others*(1000000.0/nok),index=classes_lut,columns=(['N_same_chr','N_other_chr']))
    # so far all offchr_df are "per mb", normalize per million input contacts

    # offchr_df.sort_index(axis=1, inplace=True)
    return out, offchr_df, classes_index, nok, same_others_df


def fly_statistics_bychr(pairs, chr_dict, bins=np.linspace(-1000000,1000000,201), chr_len_vec=None, mask_dist=[0, 10000, 30000, 100000, 300000,1000000,3000000, 10000000], qmin_rna=0, qmin_dna=10, nmax=0, noambiv=False): #dict of vectors, which are targeted chr, first numer is # match,
# second is # of dna match
    if noambiv:
        cmd = "cat "+pairs+" | awk '($25==1 && $8>"+str(qmin_rna)+" && $9>"+str(qmin_dna)+"){print $4, substr($2,3), $5, $3}'"
    else:
        cmd = "cat "+pairs+" | awk '($8>"+str(qmin_rna)+" && $9>"+str(qmin_dna)+"){print $4, substr($2,3), $5, $3}'"

    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
    nok=0

    nchr=len(chr_dict)
    nbins=len(bins)+1
    nmasks=len(mask_dist)

    out=np.zeros((nchr,nbins),int)
    out_offchr=np.zeros((nchr,2+nmasks))
    out_nsame_others=np.zeros((nchr,2), int)
    if not chr_len_vec is None:
        n2add_total=1000000.0/(chr_len_vec.sum())
        n2add_other=1000000.0/(chr_len_vec.sum()-chr_len_vec)
    else:
        n2add_total=1
        n2add_other=np.ones(len(chr_dict))

        # read all reads in the file unitl reach end of file
    for _, line in enumerate(p.stdout):
        if nmax>0 and nok>nmax:
            break
        read_data=line.strip().split(" ")
        tgt_chr=chr_dict.get(read_data[0],-1)
        src_chr=chr_dict.get(read_data[1],-1)

        if tgt_chr>-1 and src_chr>-1: #only if this is good stuff
            
            out_offchr[src_chr,-1]+=n2add_total
            nok+=1
            if tgt_chr==src_chr:
                out_nsame_others[src_chr,0]+=1
                fly=int(read_data[2])-int(read_data[3])
                src_pos=int(read_data[3])
                out[src_chr,np.searchsorted(bins,fly)]+=1
    
                for i, wm  in enumerate(mask_dist):
                    if abs(fly)>= wm:
                        n2add_this=1 if (chr_len_vec is None) else 1000000.0/(max(0.0,chr_len_vec[src_chr]-(src_pos+wm))+max(0.0,src_pos-wm))
                        out_offchr[src_chr,i]+=n2add_this #n2add_same[i][src_chr]

            else:
                out_nsame_others[src_chr,1]+=1
                out_offchr[src_chr,-2]+=n2add_other[src_chr]
    dw=np.hstack([1,np.diff(bins),1])
    out=out*(1000000.0/dw)*(1000000.0/nok)

    p.terminate()
    del p

    offchr_df=pd.DataFrame(out_offchr*(1000000.0/nok) ,columns=(['N_outof_%ikbp'%(wm/1000) for wm in mask_dist]+['N_other_chr','N']))
    same_others_df=pd.DataFrame(out_nsame_others*(1000000.0/nok), columns=(['N_same_chr','N_other_chr']))
    return out, offchr_df, nok, same_others_df

def fly_statistics_bychr_coolerdump(cdump, chr_dict, bins=np.linspace(-1000000,1000000,201), chr_len_vec=None, mask_dist=[0, 10000, 30000, 100000, 300000,1000000,3000000, 10000000], nmax=0): #dict of vectors, which are targeted chr, first numer is # match,
# second is # of dna match
    cmd = "cat "+cdump+" | awk '{print $4, $1, $5, $2, $7}'"

    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
    nok=0

    nchr=len(chr_dict)
    nbins=len(bins)+1
    nmasks=len(mask_dist)

    out=np.zeros((nchr,nbins),int)
    out_offchr=np.zeros((nchr,2+nmasks))

    if not chr_len_vec is None:
        n2add_total=1000000.0/(chr_len_vec.sum())
        n2add_other=1000000.0/(chr_len_vec.sum()-chr_len_vec)
    else:
        n2add_total=1
        n2add_other=np.ones(len(chr_dict))

        # read all reads in the file unitl reach end of file
    for _, line in enumerate(p.stdout):
        print(line)
        if nmax>0 and nok>nmax:
            break
        read_data=line.strip().split(" ")
        tgt_chr=chr_dict.get(read_data[0],-1)
        src_chr=chr_dict.get(read_data[1],-1)
        nhits=int(read_data[4])

        if tgt_chr>-1 and src_chr>-1: #only if this is good stuff
            out_offchr[src_chr,-1]+=(n2add_total*nhits)
            nok+=1
            if tgt_chr==src_chr:
                fly=abs(int(read_data[2])-int(read_data[3]))+1 #+1 is so that 0 gives 1bp to fall into first bin, and 500 gives 501 which falls into next bin
                src_pos=int(read_data[3])
                out[src_chr,np.searchsorted(bins,fly)]+=nhits
    
                for i, wm  in enumerate(mask_dist):
                    if abs(fly)>= wm:
                        n2add_this=1 if (chr_len_vec is None) else 1000000.0/(max(0.0,chr_len_vec[src_chr]-(src_pos+wm))+max(0.0,src_pos-wm))
                        out_offchr[src_chr,i]+=(n2add_this*nhits) #n2add_same[i][src_chr]

            else:
                out_offchr[src_chr,-2]+=(n2add_other[src_chr]*nhits)
    dw=np.hstack([1,np.diff(bins),1])
    out=out*(1000000.0/dw)*(1000000.0/nok)

    p.terminate()
    del p

    offchr_df=pd.DataFrame(out_offchr*(1000000.0/nok) ,columns=(['N_outof_%ikbp'%(wm/1000) for wm in mask_dist]+['N_other_chr','N']))

    return out, offchr_df, nok

def fly_statistics_bychr_coolerdump2(cdump, nbins=25000, nmax=0, dw=1000): #dict of vectors, which are targeted chr, first numer is # match,
# second is # of dna match
    cmd = "source activate charseq; cooler dump "+cdump

    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
    nok=0


    out=np.zeros(nbins+2,int)
    # out_offchr=np.zeros((nchr,2+nmasks))



        # read all reads in the file unitl reach end of file
    for _, line in enumerate(p.stdout):
        # print (line)
        if nmax>0 and nok>nmax:
            break
        nok+=1
        read_data=line.strip().split("\t")
        fly=min(abs(int(read_data[1])-int(read_data[0])),nbins)

        out[fly+1]+=int(read_data[2])
    
    out=out*(1000000.0/dw)*(1000000.0/nok)

    p.terminate()
    del p

    return out, nok
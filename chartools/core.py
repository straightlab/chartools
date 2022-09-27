import numpy as np
import csv
from scipy import sparse
import pandas as pd
from scipy import stats as stats


from numpy.random import multinomial, random_sample, randint
#importnp.random.random_sample(10)

class Chartable:
    def __init__(self, file_bed=None, chr_dict=None, annot_chr_dict=None, chr_len_vec=None, bins_name_col=-1):
        # rows are annotations
        
        self.annot_dict={} # name to ix, chrix  #dictionnary
        
        self.annot_LUT=[]
        self.annot_chr=[]
        self.nannot=0

        self.bins_LUT=[]


       
        if not chr_len_vec is None:
            self.chr_len_vec=chr_len_vec.copy()
        else:
             self.chr_len_vec=np.zeros(0, int)

        # columns are bins, blocked by chr
        if not chr_dict is None:
            self.chr_dict={k: v for k,v in chr_dict.items()}
            if file_bed is None:
                self.chr_blocks=[]
                self.nbins=0
            else:
                self.chr_blocks, self.bins_LUT=Chartable.load_chr_blocks(file_bed, chr_dict, bins_name_col=bins_name_col)
                self.nbins=self.chr_blocks[-1] #np.max(self.chr_blocks[:,1])+1

            self.nchr=max((v for v in chr_dict.values()))+1
            
        
            self.chr_LUT=[""]*self.nchr
        
            for k, v in self.chr_dict.items():
                self.chr_LUT[v]=k


            if not annot_chr_dict is None:
                self.load_annot_dict(annot_chr_dict)
            else:
                self.chr_dict={}
                self.chr_blocks=[]
                self.nchr=0
                self.nbins=0
                self.chr_LUT=[]

        # count matrix 
        self.counts=[]

    def copy(self, copydata=True):
        # rows are annotations
        out=Chartable()
        out.annot_dict={k:v for k,v in self.annot_dict.items()} # name to ix, chrix  #dictionnary
        
        out.annot_LUT=self.annot_LUT.copy()
        out.annot_chr=self.annot_chr.copy()
        out.nannot=self.nannot

        # columns are bins, blocked by chr
        out.chr_dict={k:v for k,v in self.chr_dict.items()}
        out.chr_blocks=self.chr_blocks.copy()
        out.chr_len_vec=self.chr_len_vec.copy()
        out.nchr=self.nchr
        out.nbins=self.nbins
        out.bins_LUT=self.bins_LUT.copy()
        
        out.chr_LUT=self.chr_LUT.copy()

        if copydata:
            out.counts=self.counts.copy()
        return out


    # @staticmethod
    # def load_chr_blocks(file_bed, chr_to_chrid):
    #     chr_startstop=-1*np.ones((len(chr_to_chrid),2), int)
        
    #     with open(file_bed) as file:
    #         reader = csv.reader(file, delimiter="\t")
    #         iterator = iter(reader)
    
    #         this_chr=""
    #         this_chr_id=-1
    #         this_row=-1

    #         for row in iterator:
    #             if row:
    #                 this_row+=1
    #                 if not (row[0]==this_chr):
    #                     if this_chr_id>-1:
    #                         chr_startstop[this_chr_id,1]=this_row
    #                         this_chr=row[0]
    #                         this_chr_id=chr_to_chrid.get(this_chr,-1)
    #                         if this_chr_id>-1:
    #                             chr_startstop[this_chr_id,0]=this_row
    #                     else:
    #                         this_chr=row[0]
    #                         this_chr_id=chr_to_chrid.get(this_chr,-1)
    #                         if this_chr_id>-1:
    #                             chr_startstop[this_chr_id,0]=this_row
                                
    #         if this_chr_id>-1:
    #             chr_startstop[this_chr_id,1]=this_row+1
    #     chr_blocks=np.hstack([chr_startstop[:,0].squeeze(), chr_startstop[-1,1]])            
    #     return chr_blocks

    @staticmethod
    def load_chr_blocks(file_bed, chr_to_chrid, bins_name_col=-1):
        nchr=len(chr_to_chrid)
        blocks=np.zeros(nchr+1, int)
        bins_LUT=[]

        with open(file_bed) as file:
            reader = csv.reader(file, delimiter="\t")
            iterator = iter(reader)
    
            this_chr=""
            this_chr_id=-1
            this_row=-1
            # this_chr_id_ok=-1


            for row in iterator:
                if row:
                    this_row+=1
                    if bins_name_col>-1:
                        bins_LUT.append(row[bins_name_col])
                    if not (row[0]==this_chr):
                        # starting a new chr:
                        this_chr=row[0]
                        this_chr_id=chr_to_chrid.get(this_chr,-1)

                        if this_chr_id>-1:
                            # this_chr_id_ok=this_chr_id # last ok chr
                            if blocks[this_chr_id]>0:
                                raise ValueError()
                            else:
                                blocks[this_chr_id]=this_row

#             if blocks[this_chr_id_ok+1]>0:
#                 raise ValueError()
#             else:
#                 blocks[this_chr_id_ok+1]=this_row+1
                
            blocks[0]=0
            blocks[-1]=this_row+1

            for i in range(nchr,0,-1):
                if blocks[i]==0:
                    blocks[i]=blocks[i+1]
        return blocks, bins_LUT


    def load_annot_dict(self, annot_chr_dict):
        self.annot_dict={k: v[0] for k, v in annot_chr_dict.items()}
        self.nannot=max((v[0] for v in annot_chr_dict.values()))+1

        self.annot_LUT=[""]*self.nannot
        self.annot_chr=np.zeros(self.nannot, int)
        
        for k, v in annot_chr_dict.items():
            self.annot_LUT[v[0]]=k
            self.annot_chr[v[0]]=v[1]
            

    # def load_raw_counts(self, file_count, binsize=1, reorganize_genes=True): # counts are sorted by gene #rows are chr, columns are bins
    #     # annot_dict has been loaded
    #     self.binsize=binsize
    #     annot_dict=self.annot_dict


    #     chr_offset_vec=np.hstack([0,np.cumsum(int(np.ceil(self.chr_len_vec/binsize)))])
    #     chr_offset_dict={self.chr_LUT[i]:v for i, v in enumerate(chr_offset_vec[0:-1])}

    #     data=[]
    #     indices=[]
    #     indptr=[] #np.zeros(self.nannot, int)
    #     # counts=sparse.csr_matrix((transfer_data,transfer_indices, transfer_indptr), shape=(self.nbins, new_nbins))

    #     with open(file_count) as file:
    #         reader = csv.reader(file, delimiter="\t")
    #         iterator = iter(reader)
    
    #         this_gene=""
    #         this_gene_ix=-1
    #         this_pos=-1
    #         # this_gene_ndata=0
    #         ix=-1 #track pos of last entered data

    #         ix_gene=-1 #track the local "id" of the last recorded gene, cross referenced to true gene id with genex_ixs
    #         gene_ixs=[]
            
    #         for r in iterator:
    #             if r:
                    
    #                 if (r[2]==this_gene) and (this_gene_ix>0): #same gene as previous gene entered and valid
    #                     new_pos=int(r[1]/binsize)+chr_offset_dict[r[0]]
    #                     if new_pos==this_pos:
    #                         data[ix]+=1
    #                     else:
    #                         data.append(1)
    #                         this_pos=new_pos
    #                         indptr.append(this_pos)
    #                         ix+=1 #last entered data

    #                 else: #new gene

    #                     #reset gene and pos info trackers
    #                     this_gene=r[2]
    #                     this_pos=-1
    #                     this_gene_ix=annot_dict.get(r[2],-1)
                        
    #                     if this_gene_ix>0:
                        
    #                         ix_gene+=1 #increase local gene ix
    #                         gene_ixs.append(this_gene_ix) #cross ref true gene_ix
                             
    #                         this_pos=int(r[1]/binsize)+chr_offset_dict[r[0]]
    #                         data.append(1)
    #                         indptr.append(this_pos)

    #                         ix+=1
    #                         indices.append(ix) #this is where the new gene (this) start
    #     indptr.append(len(data)+1)
        
    #     self.nbins=chr_offset_vec[-1]
    #     counts=sparse.csr_matrix((data,indices, indptr), shape=(ix_gene+1,self.nbins))

    #     print("Found %g genes"%(ix_gene+1))
    #     if reorganize_genes: # if true, we make sure the row id are the same are the annot_dict
    #         remap_indices=np.argsort(gene_ixs)
    #         remap_indptr=np.zeros(self.nannot, int)
    #         remap_indptr[gene_ixs]=1
    #         remap_indptr=np.hstack((0, np.cumsum(remap_indptr)))
    #         remap_data=np.ones(len(gene_ixs))

    #         remap_matrix=sparse.csr_matrix((remap_data,remap_indices, remap_indptr), shape=(self.nannot,ix_gene+1))

    #         counts=remap_matrix*counts
    #         self.counts=counts

    #     else:
    #         self.counts=counts



    #     # now we need to reindex the rows

    

    def load_count_table(self, file_count): # rows are chr, columns are bins
        
        with open(file_count) as file:
            reader = csv.reader(file, delimiter="\t")
            iterator = iter(reader)
    
            # data structures for the coo_matrix
            data = []
            rows = []
            cols = []
    
            for row in iterator:
                if row:
                    ix=self.annot_dict.get(row[0], -1)
                    if ix>-1:
                        rows.append(ix)
                        cols.append(int(row[1])-1) 
                        data.append(int(row[2]))
        
            # create and return the coo_matrix
            data = np.array(data, copy=False)
            rows = np.array(rows, dtype=np.int32, copy=False)
            cols = np.array(cols, dtype=np.int32, copy=False)
            m = sparse.coo_matrix((data, (rows, cols)), shape=(self.nannot, self.nbins))
            self.counts=m.tocsr()
    

    def load_fly_table(self, file_count, bins, reorganize_genes=False, largerorequal=True, gene_strand_dict=None):

        b=bins
        nb=len(b)
        self.nbins=len(b)+1

        annot_dict=self.annot_dict
        chr_len_vec=self.chr_len_vec
        
        # possible=np.zeros((self.nchr,self.nbins),np.int64)

        if largerorequal:
            b=bins-1


        with open(file_count) as file:
            reader = csv.reader(file, delimiter="\t")
            iterator = iter(reader)
    
            # data structures for the coo_matrix
            data=[]
            indices=[]
            indptr=[] 

            
            data_impossible=[]
            this_gene_impossible=np.zeros(nb+1, int)

       
            this_gene=""
            this_gene_ix=-1
            this_pos=-1

            this_chr_ix=-1
            this_chr=""
           
            ix=-1 #track pos of last entered data
           

            ix_gene=-1 #track the local "id" of the last recorded gene, cross referenced to true gene id with genex_ixs
            gene_ixs=[]
            
            directional = False if gene_strand_dict is None else True
            dirflag=1



            for r in iterator:
                if r:
                    if (r[2]==this_gene) and (this_gene_ix>-1): #same gene as previous gene entered and valid
                        new_pos=np.searchsorted(b, dirflag*int(r[1]))
                        if new_pos==this_pos:
                            data[ix]+=1
                        else:
                            
                            data.append(1)
                            this_pos=new_pos
                            indices.append(this_pos)
                            ix+=1 #last entered data



                        max_flight_pos=(chr_len-int(r[3])-1) if dirflag>0 else int(r[3])
                        min_flight_pos=-int(r[3]) if dirflag>0 else (-chr_len-int(r[3])+1) 
                        
                        max_flight=np.searchsorted(b,max_flight_pos)
                        min_flight=np.searchsorted(b,min_flight_pos)

                        # possible[this_chr_ix,min_flight:max_flight]+=1

                        n_impossible_left=min_flight
                        n_impossible_right=nb-max_flight
                    
                        if n_impossible_left>0:
                            this_gene_impossible[0:(min_flight)]+=1
                        if n_impossible_right>0:
                            this_gene_impossible[(max_flight+1):(nb+1)]+=1


                    else: #new gene

                        if this_gene_ix>-1: #previous gene is valid
                            data_impossible+=[sparse.csr_matrix(this_gene_impossible)]
                            
                        #reset gene and pos info trackers
                        this_gene_impossible=np.zeros(nb+1, int)
                        this_gene=r[2]

                        this_pos=-1
                        this_gene_ix=annot_dict.get(r[2],-1)

                        this_chr=r[0]
                        this_chr_ix=self.chr_dict[this_chr]
                        chr_len=chr_len_vec[this_chr_ix]

                        if directional and (gene_strand_dict.get(this_gene, False)): #TRUE means negative strand
                            dirflag=-1
                        else:
                            dirflag=1

                        if this_gene_ix>-1: #valid gene

                            ix_gene+=1 #increase local gene ix
                            gene_ixs.append(this_gene_ix) #cross ref true gene_ix
                             
                            this_pos=np.searchsorted(b, dirflag*int(r[1]))

                            max_flight_pos=(chr_len-int(r[3])-1) if dirflag>0 else int(r[3])
                            min_flight_pos=-int(r[3]) if dirflag>0 else (-chr_len-int(r[3])+1) 
                            
                            max_flight=np.searchsorted(b,max_flight_pos)
                            min_flight=np.searchsorted(b,min_flight_pos)

                            # possible[this_chr_ix,min_flight:max_flight]+=1

                            # n_impossible_left=min_flight-1
                            # n_impossible_right=nb-max_flight

                            # if n_impossible_left>0:
                            #     this_gene_impossible[0:(min_flight-1)]+=1
                            # if n_impossible_right>0:
                            #     this_gene_impossible[max_flight:nb]+=1

                            n_impossible_left=min_flight
                            n_impossible_right=nb-max_flight
                        
                            if n_impossible_left>0:
                                this_gene_impossible[0:(min_flight)]+=1
                            if n_impossible_right>0:
                                this_gene_impossible[(max_flight+1):(nb+1)]+=1


                            data.append(1)
                            indices.append(this_pos)

                            ix+=1
                            indptr.append(ix) #this is where the new gene (this) start
        
        if this_gene_ix>-1:
            data_impossible+=[sparse.csr_matrix(this_gene_impossible)]
        indptr.append(len(data))
        

        
        self.chr_blocks=np.array([0,self.nbins])
        self.bins_LUT=[0 for _ in range(len(b)+1)]
        
        print("Found %g genes"%(ix_gene+1))
        self.counts=sparse.csr_matrix((data,indices, indptr), shape=(ix_gene+1,self.nbins))

        full_annot_LUT=[k for k in self.annot_LUT]
        full_annot_chr=self.annot_chr.copy()
        
        new_annot_LUT=[full_annot_LUT[ix] for ix in (gene_ixs)]
        new_annot_dict={g:i for i, g in enumerate(new_annot_LUT)}
        new_nannot=len(new_annot_LUT)
        new_annot_chr=full_annot_chr[gene_ixs]

        self.annot_LUT=new_annot_LUT
        self.annot_dict=new_annot_dict
        self.nannot=new_nannot
        self.annot_chr=new_annot_chr

        alldata_impossible=sparse.vstack(data_impossible)

        if reorganize_genes: # if true, we make sure the row id are the same are the annot_dict
            self.select_annots(full_annot_LUT, inplace=True)

        # possible=possible[:,1:-1]

        inaccessible_flight=self.copy(copydata=False)
        inaccessible_flight.counts=alldata_impossible

        return inaccessible_flight

    def load_fly_table_density(self, file_count, bins, reorganize_genes=False, largerorequal=True, gene_strand_dict=None):

        b=bins
        nb=len(b)
        self.nbins=len(b)+1

        annot_dict=self.annot_dict
        chr_len_vec=self.chr_len_vec
        
        # possible=np.zeros((self.nchr,self.nbins),np.int64)

        if largerorequal:
            b=bins-1

        deltas=np.diff(bins)

        with open(file_count) as file:
            reader = csv.reader(file, delimiter="\t")
            iterator = iter(reader)
    
            # data structures for the coo_matrix
            data=[]
            indices=[]
            indptr=[] 

            
            data_impossible=[]
            this_gene_impossible=np.zeros(nb+1, int)

       
            this_gene=""
            this_gene_ix=-1
            this_pos=-1

            this_chr_ix=-1
            this_chr=""
           
            ix=-1 #track pos of last entered data
           

            ix_gene=-1 #track the local "id" of the last recorded gene, cross referenced to true gene id with genex_ixs
            gene_ixs=[]
            
            directional = False if gene_strand_dict is None else True
            dirflag=1



            for r in iterator:
                if r:
                    if (r[2]==this_gene) and (this_gene_ix>-1): #same gene as previous gene entered and valid
                        new_pos=np.searchsorted(b, dirflag*int(r[1]))
                        if new_pos==this_pos:
                            data[ix]+=1
                        else:
                            
                            data.append(1)
                            this_pos=new_pos
                            indices.append(this_pos)
                            ix+=1 #last entered data



                        max_flight_pos=(chr_len-int(r[3])-1) if dirflag>0 else int(r[3])
                        min_flight_pos=-int(r[3]) if dirflag>0 else (-chr_len-int(r[3])+1) 
                        
                        max_flight=np.searchsorted(b,max_flight_pos)
                        min_flight=np.searchsorted(b,min_flight_pos)

                        # possible[this_chr_ix,min_flight:max_flight]+=1

                        n_impossible_left=min_flight
                        n_impossible_right=nb-max_flight
                    
                        if n_impossible_left>0:
                            this_gene_impossible[0:(min_flight)]+=1
                        if n_impossible_right>0:
                            this_gene_impossible[(max_flight+1):(nb+1)]+=1


                    else: #new gene

                        if this_gene_ix>-1: #previous gene is valid
                            data_impossible+=[sparse.csr_matrix(this_gene_impossible)]
                            
                        #reset gene and pos info trackers
                        this_gene_impossible=np.zeros(nb+1, int)
                        this_gene=r[2]

                        this_pos=-1
                        this_gene_ix=annot_dict.get(r[2],-1)

                        this_chr=r[0]
                        this_chr_ix=self.chr_dict[this_chr]
                        chr_len=chr_len_vec[this_chr_ix]

                        if directional and (gene_strand_dict.get(this_gene, False)): #TRUE means negative strand
                            dirflag=-1
                        else:
                            dirflag=1

                        if this_gene_ix>-1: #valid gene

                            ix_gene+=1 #increase local gene ix
                            gene_ixs.append(this_gene_ix) #cross ref true gene_ix
                             
                            this_pos=np.searchsorted(b, dirflag*int(r[1]))

                            max_flight_pos=(chr_len-int(r[3])-1) if dirflag>0 else int(r[3])
                            min_flight_pos=-int(r[3]) if dirflag>0 else (-chr_len-int(r[3])+1) 
                            
                            max_flight=np.searchsorted(b,max_flight_pos)
                            min_flight=np.searchsorted(b,min_flight_pos)

                            # possible[this_chr_ix,min_flight:max_flight]+=1

                            # n_impossible_left=min_flight-1
                            # n_impossible_right=nb-max_flight

                            # if n_impossible_left>0:
                            #     this_gene_impossible[0:(min_flight-1)]+=1
                            # if n_impossible_right>0:
                            #     this_gene_impossible[max_flight:nb]+=1

                            n_impossible_left=min_flight
                            n_impossible_right=nb-max_flight
                        
                            if n_impossible_left>0:
                                this_gene_impossible[0:(min_flight)]+=1
                            if n_impossible_right>0:
                                this_gene_impossible[(max_flight+1):(nb+1)]+=1


                            data.append(1)
                            indices.append(this_pos)

                            ix+=1
                            indptr.append(ix) #this is where the new gene (this) start
        
        if this_gene_ix>-1:
            data_impossible+=[sparse.csr_matrix(this_gene_impossible)]
        indptr.append(len(data))
        

        
        self.chr_blocks=np.array([0,self.nbins])
        self.bins_LUT=[0 for _ in range(len(b)+1)]
        
        print("Found %g genes"%(ix_gene+1))
        self.counts=sparse.csr_matrix((data,indices, indptr), shape=(ix_gene+1,self.nbins))

        full_annot_LUT=[k for k in self.annot_LUT]
        full_annot_chr=self.annot_chr.copy()
        
        new_annot_LUT=[full_annot_LUT[ix] for ix in (gene_ixs)]
        new_annot_dict={g:i for i, g in enumerate(new_annot_LUT)}
        new_nannot=len(new_annot_LUT)
        new_annot_chr=full_annot_chr[gene_ixs]

        self.annot_LUT=new_annot_LUT
        self.annot_dict=new_annot_dict
        self.nannot=new_nannot
        self.annot_chr=new_annot_chr

        alldata_impossible=sparse.vstack(data_impossible)

        if reorganize_genes: # if true, we make sure the row id are the same are the annot_dict
            self.select_annots(full_annot_LUT, inplace=True)

        # possible=possible[:,1:-1]

        inaccessible_flight=self.copy(copydata=False)
        inaccessible_flight.counts=alldata_impossible

        return inaccessible_flight

    def reindex_chr(self, chr_LUT, drop_annots=True):

        out=self.copy()
        out.chr_dict={k: i for i, k in enumerate(chr_LUT)}
        new_nchr=len(chr_LUT)
        out.nchr=new_nchr
        out.chr_LUT=chr_LUT.copy()
        
        new_chr_dict=out.chr_dict
        old_chr_LUT=self.chr_LUT
        annot_chr_new=np.array([new_chr_dict.get(old_chr_LUT[k], -1) for k in self.annot_chr])
        
        dropped_annot = []
        if drop_annots:
            keep=(annot_chr_new>-1)
            dropped_annot=[val for is_good, val in zip(keep, self.annot_LUT) if not is_good]

            out.annot_chr=annot_chr_new[keep]
            out.annot_LUT=[val for is_good, val in zip(keep, self.annot_LUT) if is_good]
            out.annot_dict={k: i for i, k in enumerate(out.annot_LUT)}
            out.nannot=len(out.annot_LUT)

        else:
            out.annot_chr=annot_chr_new
            out.annot_LUT=[val for val in self.annot_LUT]
            out.annot_dict={k: i for i, k in enumerate(out.annot_LUT)}
            out.nannot=len(out.annot_LUT)

        chr_new_to_old_LUT= np.zeros(out.nchr, int) #nbin2 in size nbin2[i] = chr
        for i, k in enumerate(chr_LUT): #new ix, name 
            chr_new_to_old_LUT[i]=self.chr_dict[k]
        
        nbins_by_chr=np.diff(self.chr_blocks)
        new_chr_blocks=np.hstack([0,np.cumsum(nbins_by_chr[chr_new_to_old_LUT])])
        
        out.chr_blocks=new_chr_blocks
        
        new_nbins=new_chr_blocks[-1]
        out.nbins=new_nbins
        
        transfer_data=np.ones(new_nbins, int)
        transfer_indptr=np.arange(new_nbins+1)
        transfer_indices=np.zeros(new_nbins, int)
        for i in range(new_nchr):
            this_chr_start=new_chr_blocks[i]
            this_chr_stop=new_chr_blocks[i+1]
            this_chr_nbins=this_chr_stop-this_chr_start
            old_this_chr_start=self.chr_blocks[chr_new_to_old_LUT[i]]
            transfer_indices[this_chr_start:this_chr_stop]=np.arange(old_this_chr_start,old_this_chr_start+this_chr_nbins)
        
        transfer_matrix=sparse.csc_matrix((transfer_data,transfer_indices, transfer_indptr), shape=(self.nbins, new_nbins))

        if drop_annots:
            out.counts=(self.counts[keep,:]*transfer_matrix).tocsr()
        else:
            out.counts=(self.counts*transfer_matrix).tocsr()

        if len(self.bins_LUT)>0:
            new_bins_ixs=(sparse.csc_matrix(np.arange(self.nbins).reshape((1,self.nbins)), shape=(1,self.nbins))*transfer_matrix).toarray().squeeze()
            this_bins_LUT=self.bins_LUT
            new_bins_LUT=[this_bins_LUT[i] for i in new_bins_ixs]
            out.bins_LUT=new_bins_LUT

        chrID_new_to_old = [self.chr_dict[i] for i in chr_LUT]
        new_chr_len_vec = self.chr_len_vec[chrID_new_to_old]
        out.chr_len_vec = new_chr_len_vec
        return out, dropped_annot

    def N_bychr(self):
        # returns vector nrowsxnchr
        ba=self.get_block_summator_bychr()
        N=self.counts*ba
        return N

    def N_cistrans(self, format='sparse', annotate=None):
        
        Nbychr=self.N_bychr()
        Ncis=sparse.csr_matrix([Nbychr[k,c] for k, c in enumerate(self.annot_chr)])
        Ntrans=sparse.csr_matrix(Nbychr.sum(axis=1).reshape((1,self.nannot)))-Ncis # reshape((1,Nbychr.shape[0]))-
        N=sparse.vstack((Ncis, Ntrans))

        if format=='dense':
            N=N.toarray().transpose()
        elif ((format=='pandas') or (not annotate is None)) :
            N=pd.DataFrame(N.toarray().transpose(), columns=['Ncis','Ntrans'], index=self.annot_LUT)
            if not annotate is None:
                N=pd.concat([annotate.reindex(N.index), N], axis=1)
        return N

    def N(self, format='dense', annotate=None):
        N=np.array(self.counts.sum(axis=1)).ravel()
        if ((format=='pandas') or (not annotate is None)) :
            N=pd.DataFrame(N, columns=['Ntotal'], index=self.annot_LUT)
            if not annotate is None:
                N=pd.concat([annotate.reindex(N.index), N], axis=1)
        return N


    def select_annots(self, annots, inplace=False):
        if inplace:
            out=self
        else:
            out=self.copy(copydata=False)
        
        annots_in=[k for k in annots if (k in self.annot_dict)]
        

        annots_ixs=[self.annot_dict[k] for k in annots_in]
        out.counts=self.counts[annots_ixs,:]
        
        out.annot_dict={k:ix for ix, k in enumerate(annots_in)}
        out.nannot=len(annots_in)
        out.annot_LUT=annots_in
        
        out.annot_chr=out.annot_chr[annots_ixs]
        
        
        if not inplace:
            return out


    def select_annots2(self, annots, genes):
        
        out=self.copy(copydata=False)
        annots_ixs=np.array([self.annot_dict.get(k, -1) for k in annots], dtype=int)

        transfer_data=np.ones(len(annots), int)
        transfer_data[annots_ixs<0]=0
        transfer_indices=annots_ixs.copy()
        transfer_indices[annots_ixs<0]=0
        transfer_indptr=np.arange(len(annots)+1)
        annot_select=sparse.csr_matrix((transfer_data,transfer_indices, transfer_indptr), shape=(len(annots), len(self.annot_LUT)))
        out.counts=(annot_select*self.counts)
        out.counts.eliminate_zeros()
        
        out.annot_dict={k:ix for ix, k in enumerate(annots)}
        out.nannot=len(annots)
        out.annot_LUT=annots.copy()
        
        
        out.annot_chr=genes['chr'].loc[annots].copy().map(self.chr_dict).values #annots_chr_id.copy()
        
        return(out)

    # def to_pandas_flat(self, colname = 'N', thr=0):
   
    #     cts=self.counts.tocsr(copy=True)

    #     if thr>0:
    #         cts.data[cts.data<thr] = 0 
    #         cts.eliminate_zeros()

    #     annots_out=[]
    #     bins_out=np.array()
        
        
    #     indices=cts.indices
    #     indptr=cts.indptr
    #     ns=np.diff(indptr)

    #     if len(self.bins_LUT)>0:
    #         bins_LUT = np.array(self.bins_LUT)
    #         for i, a in enumerate(self.annot_LUT):
    #             if ns[i]>0:
    #                 annots_out+=[a]*ns[i]
    #                 bins_out=np.append(bins_out, bins_LUT[indices[indptr[i]:indptr[i+1]]])
    #     else:
    #         for i, a in enumerate(self.annot_LUT):
    #             if ns[i]>0:
    #                 annots_out+=[a]*ns[i]
    #                 bins_out=np.append(bins_out, indices[indptr[i]:indptr[i+1]])

    #     out=pd.DataFrame({'src':annots_out,'tgt':bins_out, colname:cts.data}).set_index(['src', 'tgt'])
    #     return out

    def to_pandas_flat(self, colname = 'N', thr=0):
        print("fast flattening")
        cts=self.counts.tocoo(copy=True)

        if thr>0:
            cts.data[cts.data<thr] = 0 
            cts.eliminate_zeros()
        
        annots_LUT=np.array(self.annot_LUT).copy()
        annots_out = annots_LUT[cts.row]
        if len(self.bins_LUT)>0:
            bins_LUT = np.array(self.bins_LUT).copy()
            bins_out=bins_LUT[cts.col]
        else:
            bins_out = cts.col.copy()
        

        out=pd.DataFrame({'src':annots_out,'tgt':bins_out, colname:cts.data}).set_index(['src', 'tgt'])
        return out



    # def select_bins(self, bins, inplace=False):
    #     if inplace:
    #         out=self
    #     else:
    #         out=self.copy(copydata=False)
        
    #     out.counts=out.counts.tocsc()
    #     if type(bins[0])==str:
    #         old_bins_LUT=self.bins_LUT
    #         if len(old_bins_LUT)==0:
    #             raise ValueError()
            
    #         old_bins_LUD={k: i for i, k in enumerate(old_bins_LUT)}
    #         new_bins_valid=[k for k in bins if k in old_bins_LUD]
    #         bins_ixs=np.array([old_bins_LUD[k] for k in new_bins_valid])
    #         out.bins_LUT=new_bins_valid

    #     else:
    #         bins_ixs=np.array(bins)
        
    #     new_nbins=len(bins_ixs)
    #     out.nbins=new_nbins
        
    #     transfer_data=np.ones(new_nbins, int)
    #     transfer_indptr=np.arange(new_nbins+1)
    #     transfer_indices=bins_ixs
        
        
    #     transfer_matrix=sparse.csc_matrix((transfer_data,transfer_indices, transfer_indptr), shape=(self.nbins, new_nbins))

    #     out.counts=(self.counts*transfer_matrix).tocsr()

    #     old_bin_to_chr_LUT=self.bin_to_chr_LUT()
    #     new_to_old_chr_LUT=(sparse.csc_matrix(np.array(old_bin_to_chr_LUT.reshape((len(old_bin_to_chr_LUT),1))), shape=(len(old_bin_to_chr_LUT),1))*transfer_matrix).toarray().squeeze()

    #     new_chr_blocks=np.array(self.nchr+1)
    #     this_chr=-1
    #     for i in range(len(new_to_old_chr_LUT)):
    #         if new_to_old_chr_LUT[i]!=this_chr:
    #             this_chr=new_to_old_chr_LUT[i]
    #             new_chr_blocks[this_chr]=i
        
    #     new_chr_blocks[-1]=new_nbins

    #     for i in range(self.nchr,0,-1):
    #         if new_chr_blocks[i]==0:
    #             new_chr_blocks[i]=new_chr_blocks[i+1]

    #     out.chr_blocks=new_chr_blocks
        
    #     if not inplace:
    #         return out

    def select_bins(self, bins, inplace=False, reorder=False):
        if inplace:
            out=self
        else:
            out=self.copy(copydata=False)

    #     print(out.counts)
        out.counts=self.counts.tocsc()
        
        if type(bins[0])==str:
            
            old_bins_LUT=self.bins_LUT
            if len(old_bins_LUT)==0:
                raise ValueError()

            if reorder:
                old_bins_LUD={k: i for i, k in enumerate(self.bins_LUT)}
                new_bins_valid=[k for k in bins if k in old_bins_LUD]
                bins_ixs=np.array([old_bins_LUD[k] for k in new_bins_valid])
                out.bins_LUT=new_bins_valid
    #         print(len(bins_LUT))
            else:
                new_bins_LUD={k: i for i, k in enumerate(bins)}
                bins_ixs=np.array([i for i, e in enumerate(self.bins_LUT) if e in new_bins_LUD])
                out.bins_LUT=[k for k in self.bins_LUT if k in new_bins_LUD]
    #         
            

        else:
            bins_ixs=np.array(bins)

        new_nbins=len(bins_ixs)
        #print(new_nbins)
        out.nbins=new_nbins

        transfer_data=np.ones(new_nbins, int)
        transfer_indptr=np.arange(new_nbins+1)
        transfer_indices=bins_ixs


        transfer_matrix=sparse.csc_matrix((transfer_data,transfer_indices, transfer_indptr), shape=(self.nbins, new_nbins))

        out.counts=(self.counts*transfer_matrix).tocsr()
        
        new_chr=np.searchsorted(self.chr_blocks, bins_ixs)
        new_chr_blocks=np.searchsorted(np.sort(new_chr),np.arange(self.nchr))
        out.chr_blocks=new_chr_blocks
    
        if not inplace:
            return out

    def collapse_giset(self, gi_set):
        out=self.copy(copydata=False)
        
        out.counts=(self.counts.tocsc()*gi_set['summing_matrix']).tocsr()
        
        new_chr_blocks=gi_set['summing_matrix'].shape[1]*np.ones(len(self.chr_blocks),int)
        new_chr_blocks[0]=0
        out.chr_blocks=new_chr_blocks
        out.bins_LUT=gi_set['set_names']
        out.nbins=len(out.bins_LUT)

        return out

    def sum(self, axis=0): #sums over annotations
        if axis==0:
             # rows are annotations
            out=Chartable()
            out.annot_dict={'sum':'0'} # name to ix, chrix  #dictionnary
        
            out.annot_LUT=['sum']
            out.annot_chr=np.array([-1], dtype=int)
            out.nannot=1
       

            # columns are bins, blocked by chr
            out.chr_dict={k:v for k,v in self.chr_dict.items()}
            out.chr_blocks=self.chr_blocks.copy()
            out.nchr=self.nchr
            out.nbins=self.nbins
        
            out.chr_LUT=self.chr_LUT.copy()
            out.bins_LUT=self.bins_LUT.copy()

            out.counts=sparse.csr_matrix(self.counts.sum(axis=0))
        else:
            out=Chartable()
            out.annot_dict={k:v for k,v in self.annot_dict.items()} # name to ix, chrix  #dictionnary
            
            out.annot_LUT=self.annot_LUT.copy()
            out.annot_chr=self.annot_chr.copy()
            out.nannot=self.nannot

            # columns are bins, blocked by chr
            out.chr_dict={'all':0}
            out.chr_blocks=np.array([0,1], int)
            out.nchr=1
            out.nbins=1
            out.bins_LUT=['all']
            
            out.chr_LUT=[0]

            out.counts=sparse.csr_matrix(self.counts.sum(axis=1))
        return out

    def to_pandas(self):
        import pandas as pd
        # return a pandas dataframe
        if sparse.issparse(self.counts):
            df=pd.DataFrame(self.counts.toarray(), index=self.annot_LUT)
            if len(self.bins_LUT)>0:
                df.columns=self.bins_LUT.copy()
        else:
            df=pd.DataFrame(self.counts, index=self.annot_LUT)
        return df

    
    # def compute_blocks_ptr(self, blocks=None):
    #     if blocks is None:
    #         blocks=self.chr_blocks

    #     #i is row, goal is to find the first and last index in the data part of the csr matrix where each chromosome starts and stops
    #     blocks_ptr=np.zeros((self.nannot,blocks.shape[0]), int)
    #     for i in range(self.nannot):
    #         blocks_ptr

    def get_block_summator_bychr(self, blocks=None):
        
        if blocks is None:
            data=np.ones(self.nbins, int)
            # for i in range(self.nchr):
            #     data[self.chr_blocks[i,0]:self.chr_blocks[i,1]]=1 #/(self.chr_blocks[i,1]-self.chr_blocks[i,0])
        #         print(data)
            indptr=self.chr_blocks.copy() #np.append(self.chr_blocks[:,0].flatten(),[self.chr_blocks[-1,1]])
            indices=np.arange(self.nbins)
            out=sparse.csr_matrix((data,indices,indptr),shape=(self.nchr,self.nbins)).transpose()
        else:
            ndata=blocks[-1]-blocks[0]
            data=np.ones(ndata, int)
            indptr=blocks.copy()
            indices=np.arange(ndata)
            out=sparse.csr_matrix((data,indices,indptr),shape=(len(blocks)-1,self.nbins)).transpose()
           
        return out


    def get_source_target_chr_data(self): #return the chromosome id of the source
        sourcechr_data=np.zeros(len(self.counts.data), int)
        
        indptr=self.counts.indptr
        annot_chr=self.annot_chr
        for i in range(self.nannot):
            sourcechr_data[indptr[i]:indptr[i+1]]=annot_chr[i]
        
        chr_LUT=self.bin_to_chr_LUT()
        
        targetchr_data=chr_LUT[self.counts.indices]

        return sourcechr_data, targetchr_data

    def cis_trans(self, format='pandas', annotate=None):
        trans=self.copy(copydata=True)
        sourcechr_data, targetchr_data = self.get_source_target_chr_data()
        
        counts=trans.counts.data
        counts[sourcechr_data==targetchr_data]=0

        cis=self.copy(copydata=False)
        
        cis.counts=sparse.csr_matrix((self.counts.data-counts, self.counts.indices.copy(), self.counts.indptr.copy()), shape=(self.nannot, self.nbins))

        cis.counts.eliminate_zeros()
        trans.counts.eliminate_zeros()

        N_cis=np.array(cis.counts.sum(axis=1)).ravel()
        N_trans=np.array(trans.counts.sum(axis=1)).ravel()
        N=np.column_stack([N_cis,N_trans])

        
        if ((format=='pandas') or (not annotate is None)) :
            N=pd.DataFrame(N, columns=['Ncis','Ntrans'], index=self.annot_LUT)
            if not annotate is None:
                N=pd.concat([annotate.reindex(N.index), N], axis=1)

        return cis, trans, N


    def bin_to_chr_LUT(self):
        data=np.ones(self.nbins, int)
        blocks=self.chr_blocks
        for i in range(self.nchr):
            data[blocks[i]:blocks[i+1]]=i
        return data

    

    def normalize(self, bychr=False, inplace=False):
        if inplace:
            out=self
        else:
            out=self.copy(copydata=True)
        if bychr:
            ba=self.get_block_summator_bychr()
            scale=(self.counts*ba).toarray() #.ravel()
            scale[scale==0]=1.0
            bin_to_chr=self.bin_to_chr_LUT()
            c=out.counts.tocoo()
            out.counts.data=self.counts.data/scale[c.row,bin_to_chr[self.counts.indices]]
        else:
            scale=self.counts.sum(axis=1)
            out.counts=sparse.csr_matrix(out.counts/scale)

        if not inplace:
            return out
        
    def normalize_row_bychr(self, normalize=False):
        out=self.copy(copydata=False)
        
        ba=self.get_block_summator_bychr()
        scale=(self.counts*ba).toarray().ravel()
        if normalize:
            scale=scale*np.diff(self.chr_blocks)
        # scale.data=1.0/scale.data


        scale[scale==0]=1.0
        bin_to_chr=self.bin_to_chr_LUT()
        y=bin_to_chr[self.counts.indices]
#             out.counts.data=self.counts.data/scale[y]

        x_change=np.append(np.argwhere(np.diff(y)>0).squeeze(),[len(y)-1]) # 0,0,0,0,1,1,1,1,1,1,1,3,3 [3,10] means y[4]-y[3]>0, y[11]-y[10]>0, so the length of the individual chr marks are 3-0+1=3-(-1), 10-4+1=10-3, final-11+1=final-10. These are for chr y[3], y[10], y[final]

#             print(y.shape)
        ndata_bychr=np.zeros(self.nchr,int)
        ndata_bychr[y[x_change]]=np.diff(np.append([-1],x_change))

        out.counts=sparse.csr_matrix((self.counts.data/scale[y], self.counts.indices.copy(), np.append([0],ndata_bychr.cumsum())), shape=(self.nchr,self.nbins))
#             out.counts.indptr=np.append([0],ndata_bychr.cumsum())
#             out.counts.shape=(self.nchr,self.nbins)

        out.nannot=out.nchr
        out.annot_chr=np.arange(out.nchr)
        out.dict={k:v for k, v in out.chr_dict.items()}
        out.annot_LUT=out.chr_LUT.copy()
       
        return out

#     def block_indices(self, blocks=None):
#         block_summator=self.get_block_summator_bychr(blocks)
#         nblocks=block_summator.get_shape()[1]
#         mask_data=np.ones(len(self.counts.data), int)
#         mask=sparse.csr_matrix((mask_data,self.counts.indices, self.counts.indptr), shape=(self.nannot, self.nbins))

#         out=mask*block_summator
#         return block_indices
        




#     def get_block_ixs(self, ):
#         bin_to_chr=self.bin_to_chr_LUT()
#         y=bin_to_chr[self.counts.indices]
# #             out.counts.data=self.counts.data/scale[y]

#         x_change=np.append(np.argwhere(np.diff(y)>0).squeeze(),[len(y)-1])
    def predict_Ntrans(self, bychr=False, chr_trans_matrix=None):
        if chr_trans_matrix is None:
            if bychr:
                Ntrans_predicted=self.N_bychr()
            else:
                Ntrans_predicted=self.N_cistrans()[1,:].transpose()
        # elif not isinstance(chr_trans_matrix, (list, np.ndarray)): #
        #     Ntrans_predicted=self.N_cistrans().sum(axis=0)*chr_trans_matrix
        # else:
        #     Ntotal=self.N_cistrans().sum(axis=0) #total value
        #     for k in range(self.nchr): #tgt chr
        #         if self.chr_blocks[k,0]>-1:
            #             Ntrans_predicted[:,k]=chr_trans_matrix[self.annot_chr,k]*Ntotal[:,np.newaxis]
            return Ntrans_predicted


    def expected_trans(self, b, bychr=False, chr_trans_matrix=None, exclude_zeros=False): # returns with 
        #3 options: total cis expected (use total cis count as reference), local chr cis expected (use chr specific counts as reference), total expected (uses a chr_trans_matrix or scalar if by chr==False and total observed count
        # )
        # 
        out=self.copy(copydata=False)
        out.counts=np.zeros((self.nannot, self.nbins))

        Ntrans_predicted=self.predict_Ntrans(bychr=bychr, chr_trans_matrix=chr_trans_matrix)

        if bychr:
            bn=b.normalize_row_bychr()
        else:
            bn=b.normalize()

        e=self.copy()
        e.counts=Ntrans_predicted*bn.counts
        return e

    def flat_field(self, b, prj, b_prenormalized=False):
    
        if not b_prenormalized:
            bkg_mat=b.normalize_row_bychr().project(prj).sum(axis=0).counts.toarray().ravel()
        else:
            bkg_mat=b.project(prj).sum(axis=0).counts.toarray().ravel()


            
        regdoms_df2=prj['df'].loc[self.bins_LUT]
        delta=(regdoms_df2['stop']-regdoms_df2['start']).values
        Lvec=self.chr_len_vec[regdoms_df2['chr'].map(self.chr_dict).values]
        delta=delta/Lvec

        dvec=bkg_mat/delta

        alldata=self.counts.data
        indices=self.counts.indices
        indptr=self.counts.indptr

        data_corr = np.zeros(alldata.shape)
        for i in range(len(indptr)-1):
            if indptr[i+1]>indptr[i]:
                data = alldata[indptr[i]:indptr[i+1]]
                Bdata=dvec[indices[indptr[i]:indptr[i+1]]]
                
                s1 = np.sum((data*Bdata))
                if s1>0:
                    data_corr[indptr[i]:indptr[i+1]]=(data*Bdata) * (np.sum(data)/s1)
                else:
                    data_corr[indptr[i]:indptr[i+1]]=0

        cts_flat=self.copy(copydata=True)
        cts_flat.counts.data = data_corr
        return cts_flat, dvec
    

    def replace_cis(self, new_cis, annots=None, genes=None, how='outter',  reorder_genes=True):
        _, out, _ = self.cis_trans()
        new_cis_cis, _, _ = new_cis.cis_trans()

        if (annots is None) and (out.annot_LUT==new_cis_cis.annot_LUT):
            annots_new = None
        else:
            if annots is None:
                if self.annot_LUT == new_cis_cis.annot_LUT:
                    annots_new = None
                else:
                    if how=='inner':
                        annots_new = list((set(self.annots).intersect(set(new_cis_cis.annots))))
                    else:
                        annots_new = list((set(self.annots).union(set(new_cis_cis.annots))))
                    annots_new = genes.loc[genes.index.isin(annots_new)].index.to_list()
            else:
                if (out.annot_LUT==annots) and (new_cis_cis.annot_LUT==annots):
                    annots_new = None
                else:
                    annots_new = annots

        if reorder_genes and (not annots_new is None):
            annots_new = genes.loc[genes.index.isin(annots_new)].index.to_list()

        if not annots_new is None:
            out = out.select_annots2(annots_new, genes=genes)
            new_cis_cis = new_cis_cis.select_annots2(annots_new, genes=genes)

        out.counts = out.counts + new_cis_cis.counts

        return out


    def replace_trans(self, new_trans, annots=None, genes=None, how='outter', reorder_genes=True):
        out, _, _ = self.cis_trans()
        _, new_trans_trans, _ = new_trans.cis_trans()

        if (annots is None) and (out.annot_LUT==new_trans_trans.annot_LUT):
            annots_new = None
        else:
            if annots is None:
                if self.annot_LUT == new_trans_trans.annot_LUT:
                    annots_new = None
                else:
                    if how=='inner':
                        annots_new = list((set(self.annots).intersect(set(new_trans_trans.annots))))
                    else:
                        annots_new = list((set(self.annots).union(set(new_trans_trans.annots))))
                    
            else:
                if (out.annot_LUT==annots) and (new_trans_trans.annot_LUT==annots):
                    annots_new = None
                else:
                    annots_new = annots

        if reorder_genes and (not annots_new is None):
            annots_new = genes.loc[genes.index.isin(annots_new)].index.to_list()

        
        if not annots_new is None:
            out = out.select_annots2(annots_new, genes=genes)
            new_trans_trans = new_trans_trans.select_annots2(annots_new, genes=genes)

        out.counts = out.counts + new_trans_trans.counts

        return out



    def to_bigwig(self, df, filename, normalize=True, normalize_reads=True, global_divide=1.0):
        import pyBigWig as pbw
        bwheader=[(k,self.chr_len_vec[i]) for i,k in enumerate(self.chr_LUT)]
        bw = pbw.open(filename, "w")
        bw.addHeader(bwheader)

        #normfactor=v.sum()*1.0 if normalize else 1.0
        normfactor=1.0
        normfactor_reads=1.0
        if normalize_reads:
            normfactor_reads=(self.counts.sum())/1e6
        for i in range(self.nchr):
            n_on_chr=self.chr_blocks[i+1]-self.chr_blocks[i]
            if n_on_chr>0:
                x_startStop=df.iloc[self.chr_blocks[i]:self.chr_blocks[i+1],1:3].values
                chroms=np.array([self.chr_LUT[i]]*n_on_chr)
                starts=x_startStop[:,0].ravel()
                ends=x_startStop[:,1].ravel()
                if normalize:
                    normfactor=df['L'].iloc[self.chr_blocks[i]:self.chr_blocks[i+1]].values
                values=self.counts[0,self.chr_blocks[i]:self.chr_blocks[i+1]].toarray().ravel().astype(np.float64)/(normfactor*normfactor_reads*global_divide)
                #print(chroms)
                bw.addEntries(chroms,starts, ends=ends, values=values, validate=False)
        bw.close()

    def sparsify(self):
        out = self.copy(copydata = False)
        counts = self.counts
        out_data = np.zeros(counts.data.shape, int)
        for i in range(counts.shape[0]):
            ndata = counts.indptr[i+1]-counts.indptr[i]
            if ndata>0:
                this_data = counts.data[counts.indptr[i]:counts.indptr[i+1]]
                nsample = int(np.ceil(this_data.sum()))
                this_data = this_data/this_data.sum()
                #this_ind = counts.indices[counts.indptr[i]:counts.indptr[i+1]]
                new_data = multinomial(nsample, this_data)
                out_data[counts.indptr[i]:counts.indptr[i+1]]=new_data
                
        out.counts = sparse.csr_matrix((out_data, counts.indices.copy(), counts.indptr.copy()), shape = counts.shape)
        return out

    

    
    # def get_bins(self, binlist):
    #     out=self.copy()
    #     out.counts=self.counts.tocsc()[:,binlist].tocsr()

    #     #update header

    #     # out.chr_dict={k:v for k,v in self.chr_dict.items()}
    #     out.chr_blocks=
    #     # self.chr_blocks.copy()
    #     # out.nchr=self.nchr
    #     out.nbins=len(binlist)
        
    #     out.chr_LUT=self.chr_LUT.copy()



    #     return out


    # def get_chromosomes(self, chrlist, byname=False):
    #     out=self.copy()
    #     out.counts=self.counts.tocsc()[:,binlist].tocsr()
    #     out
    #     self.counts=
    #


class Chartable_binned:
    def __init__(self, chr_dict, chr_len_vec, annot_chr_dict=None, binsize=1):

        self.chr_dict={k: v for k,v in chr_dict.items()}
        self.nchr=max((v for v in chr_dict.values()))+1
        self.chr_LUT=[""]*self.nchr
        for k, v in self.chr_dict.items():
            self.chr_LUT[v]=k

        self.binsize=binsize
        self.chr_len_vec=chr_len_vec.copy()
        self.chr_blocks=np.hstack([0,np.cumsum(np.ceil(self.chr_len_vec/self.binsize).astype(int))])

        self.nbins=self.chr_blocks[-1]


        # loaded with annot_chr_dict in not empty
        self.annot_dict={} # name to ix, chrix  #dictionnary
        self.annot_LUT=[]
        self.annot_chr=[]
        self.nannot=0

        if not annot_chr_dict is None:
            self.load_annot_dict(annot_chr_dict)
        
        # count matrix 
        self.counts=[]

        #self.iscompressed=False
        self.x=np.array(0,int) #compressed position if compressed

    def load_annot_dict(self, annot_chr_dict):
        self.annot_dict={k: v[0] for k, v in annot_chr_dict.items()}
        self.nannot=max((v[0] for v in annot_chr_dict.values()))+1

        self.annot_LUT=[""]*self.nannot
        self.annot_chr=np.zeros(self.nannot, int)
        
        for k, v in annot_chr_dict.items():
            self.annot_LUT[v[0]]=k
            self.annot_chr[v[0]]=v[1]

    def load_annot_from_bed(self, bed_df):
        annot_chr_dict={k: [i, self.chr_dict[bed_df.iat[i,0]]] for i, k in enumerate(bed_df.index.to_list())}
        self.load_annot_dict(annot_chr_dict)


    def copy(self, copydata=True):
        
        out=Chartable_binned(self.chr_dict, self.chr_len_vec, binsize=self.binsize)
        
        out.nbins=self.nbins
        out.chr_len_vec=self.chr_len_vec.copy()
        out.chr_LUT=[k for k in self.chr_LUT]
        out.nchr=self.nchr
        out.chr_dict={k:v for k, v in self.chr_dict.items()}
        out.chr_blocks=self.chr_blocks.copy()


        out.annot_dict={k:v for k,v in self.annot_dict.items()} # name to ix, chrix  #dictionnary
        out.annot_LUT=self.annot_LUT.copy()
        out.annot_chr=self.annot_chr.copy()
        out.nannot=self.nannot
        #out.iscompressed=False
        if copydata:
            out.counts=self.counts.copy()
            out.x=self.x.copy()
            #out.iscompressed=self.iscompressed
        return out

    def load_raw_counts(self, file_count, reorganize_genes=True, compress=True, compress_at=None, idix=2): # counts are sorted by gene #rows are chr, columns are bins
        
        # annot_dict has been loaded
        annot_dict=self.annot_dict
        binsize=self.binsize

        chr_offset_dict={self.chr_LUT[i]:v for i, v in enumerate(self.chr_blocks[0:-1])}

        data=[]
        indices=[]
        indptr=[] 

        with open(file_count) as file:
            reader = csv.reader(file, delimiter="\t")
            iterator = iter(reader)
    
            this_gene=""
            this_gene_ix=-1
            this_pos=-1
            # this_gene_ndata=0
            ix=-1 #track pos of last entered data

            ix_gene=-1 #track the local "id" of the last recorded gene, cross referenced to true gene id with genex_ixs
            gene_ixs=[]
            
            for r in iterator:
                if r:
                    
                    if (r[idix]==this_gene) and (this_gene_ix>-1): #same gene as previous gene entered and valid
                        new_pos=int(int(r[1])/binsize)+chr_offset_dict[r[0]]
                        if new_pos==this_pos:
                            data[ix]+=1
                        else:
                            data.append(1)
                            this_pos=new_pos
                            indices.append(this_pos)
                            ix+=1 #last entered data

                    else: #new gene

                        #reset gene and pos info trackers
                        this_gene=r[idix]
                        # print(this_gene)
                        this_pos=-1
                        this_gene_ix=annot_dict.get(r[idix],-1)
                        # print(this_gene_ix)
                        if this_gene_ix>-1:
                        
                            ix_gene+=1 #increase local gene ix
                            gene_ixs.append(this_gene_ix) #cross ref true gene_ix
                             
                            this_pos=int(int(r[1])/binsize)+chr_offset_dict[r[0]]
                            data.append(1)
                            indices.append(this_pos)

                            ix+=1
                            indptr.append(ix) #this is where the new gene (this) start
        indptr.append(len(data))
        

        print("Found %g genes"%(ix_gene+1))
        self.counts=sparse.csr_matrix((data,indices, indptr), shape=(ix_gene+1,self.nbins))

        comp_counts=None
        if compress:
            comp_counts=self.compress(x=compress_at)

        
        
        full_annot_LUT=[k for k in self.annot_LUT]
        full_annot_chr=self.annot_chr.copy()
        
        new_annot_LUT=[full_annot_LUT[ix] for ix in (gene_ixs)]
        new_annot_dict={g:i for i, g in enumerate(new_annot_LUT)}
        new_nannot=len(new_annot_LUT)
        new_annot_chr=full_annot_chr[gene_ixs]

        self.annot_LUT=new_annot_LUT
        self.annot_dict=new_annot_dict
        self.nannot=new_nannot
        self.annot_chr=new_annot_chr

        if reorganize_genes: # if true, we make sure the row id are the same are the annot_dict
            self.select_annots(full_annot_LUT, inplace=True, keepzeros=True, annots_chr=full_annot_chr)


        return comp_counts

    def loadRNA(self, pairsfile, intervals): 
        import pypairix as px
        binsize = self.binsize
        nbins = self.nbins
        chr_offset_dict = {self.chr_LUT[i]:v for i, v in enumerate(self.chr_blocks[0:-1])}

        data = []
        indices = []

        n_intervals = intervals.shape[0]
        indptr = np.zeros(n_intervals + 1)

        tb=px.open(pairsfile)

        ix = 0

        for i in trange(n_intervals):
            indptr[i]=ix
            querystr = '*|{}:{}-{}'.format(intervals.iat[i,0], intervals.iat[i,1], intervals.iat[i,2])
            #print(querystr)
            it = tb.querys2D(querystr)
            this_rna_pos = -1
            #n_contacts = 0
            for contact in it:
                #print(contact[1][2:])
                rna_pos = int(int(contact[2])/binsize)+chr_offset_dict[contact[1][2:]]
                if rna_pos == this_rna_pos:
                    data[-1] += 1
                else:
                    this_rna_pos = rna_pos
                    data.append(1)
                    indices.append(rna_pos)
                    ix += 1
                
        #print(ix)
        indptr[-1]=ix
        #print(indptr)
        #print(nbins)
        data = np.array(data)
        indices = np.array(indices)
        #print(len(data))
        self.counts=sparse.csr_matrix((data, indices, indptr), shape=(n_intervals,nbins))
        self.compress()
        

    def load_counts_byfile(self, file_count_list, smpl_ids, compress=True, ix_count=3): # counts are sorted by gene #rows are chr, columns are bins
        
        binsize=self.binsize

        chr_offset_dict={self.chr_LUT[i]:v for i, v in enumerate(self.chr_blocks[0:-1])}

        data=[]
        indices=[]
        indptr=[0]

        this_file_ix=-1
        ix=-1

        for file_count in file_count_list:
            this_file_ix+=1
            n=0
            print("Lookig at %s"%(file_count))
            with open(file_count) as file:
                reader = csv.reader(file, delimiter="\t")
                iterator = iter(reader)
                this_pos=-1
               
                for r in iterator:
                    if r:
                        
                        if r[0] in chr_offset_dict:
                            n+=1
                            new_pos=int(int(r[1])/binsize)+chr_offset_dict[r[0]]
                            if new_pos==this_pos:
                                data[ix]+=int(r[ix_count])
                            else:
                                data.append(int(r[ix_count]))
                                this_pos=new_pos
                                indices.append(this_pos)
                                ix+=1 #last entered data
            indptr.append(ix+1)
            

            print("Found %g data in %s"%(n, file_count))
        self.counts=sparse.csr_matrix((data,indices, indptr), shape=(len(smpl_ids),self.nbins))

        self.annot_LUT=smpl_ids.copy()
        self.annot_dict={k: i for i, k in enumerate(smpl_ids)}
        self.nannot=len(smpl_ids)
        self.annot_chr=np.zeros(len(smpl_ids),int)

        comp_counts=None
        if compress:
            comp_counts=self.compress()


        return comp_counts


    def iscompressed(self):
        return ((len(self.x.shape)>0) and (len(self.x)>0))

    def compress(self, x=None, mincount=1):
        
        if self.iscompressed(): #recompress existing (for example after filtering rows)
            # for now ignore the x
            print("Recompressing coordinates")
            idx = np.unique(sparse.find(self.counts)[1], return_counts=True)
            ix=idx[0]
            y=idx[1]
            if mincount>1:
                ix=ix[y>=mincount]
                y=y[y>=mincount]
            self.counts=self.counts[:,ix]
            self.x=self.x[ix]

            chr_ids=np.searchsorted(self.chr_blocks,ix,side='right')-1
            self.chr_blocks=np.searchsorted(chr_ids,np.arange(self.nchr+1))
            self.nbins=self.chr_blocks[-1]
        else:
            #self.iscompressed=True
            if x is None:

                idx = np.unique(sparse.find(self.counts)[1], return_counts=True)
                y=idx[1]
        
                self.counts=self.counts[:,idx[0]]
                self.x=idx[0]
            else:
                self.counts=self.counts[:,x]
                self.x=x
                y=None

            new_chr_blocks=np.searchsorted(self.x,self.chr_blocks)
            self.chr_blocks=new_chr_blocks
            self.nbins=self.chr_blocks[-1]

        return y

    # def evenout_x(self, binsize):
        

    def select_bins_bool(self, ixs_bool):
        
        if self.iscompressed(): #recompress existing (for example after filtering rows)

            print("Selecting bins")
            N0=len(self.x)
            self.counts=self.counts[:,ixs_bool]
            self.x=self.x[ixs_bool]
            self.nbins=len(self.x)
            chr_ids=np.searchsorted(self.chr_blocks,np.arange(N0)[ixs_bool],side='right')-1
            self.chr_blocks=np.searchsorted(chr_ids,np.arange(self.nchr+1))
            
        else:
            print('Not implemented')

    def get_bin_chr(self):
        return (np.searchsorted(self.chr_blocks,np.arange(self.nbins),side='right')-1)

    def select_annots(self, annots, inplace=False, keepzeros=False, annots_chr=None):
        if inplace:
            out=self
        else:
            out=self.copy(copydata=False)
            out.x=self.x.copy()
        
        annots_in=[k for k in annots if (k in self.annot_dict)]

        annots_ixs=np.array([self.annot_dict[k] for k in annots_in])

        out.counts=self.counts[annots_ixs,:]
        
        if keepzeros:
            if out.counts.getformat()!='csr':
                out.counts=out.counts.tocsr()

            new_indptr=np.zeros(len(annots)+1, int)
            new_indptr[np.hstack([np.sort(annots_ixs),len(annots)])]=out.counts.indptr
            # new_indptr[-1]=len(annots_ixs)
            for i in range(len(annots),0,-1):
                if new_indptr[i]==0:
                    new_indptr[i]=new_indptr[i+1]

            

         
            out.counts=sparse.csr_matrix((out.counts.data, out.counts.indices, new_indptr), shape=(len(annots), out.nbins))
            
            out.annot_dict={k:ix for ix, k in enumerate(annots)}
            out.nannot=len(annots)
            out.annot_LUT=[k for k in annots]
            out.annot_chr=annots_chr.copy()
        
        else:

            out.annot_dict={k:ix for ix, k in enumerate(annots_in)}
            out.nannot=len(annots_in)
            out.annot_LUT=annots_in
            
            out.annot_chr=out.annot_chr[annots_ixs]
        

        if not inplace:
            return out


    def select_annots2(self, annots, genes):
        
        out=self.copy(copydata=False)
        out.x = self.x.copy()
        annots_ixs=np.array([self.annot_dict.get(k, -1) for k in annots], dtype=int)

        transfer_data=np.ones(len(annots), int)
        transfer_data[annots_ixs<0]=0
        transfer_indices=annots_ixs.copy()
        transfer_indices[annots_ixs<0]=0
        transfer_indptr=np.arange(len(annots)+1)
        annot_select=sparse.csr_matrix((transfer_data,transfer_indices, transfer_indptr), shape=(len(annots), len(self.annot_LUT)))
        out.counts=(annot_select*self.counts)
        out.counts.eliminate_zeros()
        
        out.annot_dict={k:ix for ix, k in enumerate(annots)}
        out.nannot=len(annots)
        out.annot_LUT=annots.copy()
        
        
        out.annot_chr=genes['chr'].loc[annots].copy().map(self.chr_dict).values #annots_chr_id.copy()
        
        return(out)

    # def select_annots2(self, annots, keepzeros=False, annots_chr=None):

    #     out=self.copy(copydata=False)
    #     out.x=self.x.copy()
        
    #     annots_in=[k for k in annots if (k in self.annot_dict)]

    #     annots_ixs=np.array([self.annot_dict[k] for k in annots_in])

    #     out.counts=self.counts[annots_ixs,:]
        
    #     if keepzeros:
    #         if out.counts.getformat()!='csr':
    #             out.counts=out.counts.tocsr()

    #         transferm=
    #         new_indptr=np.zeros(len(annots)+1, int)
    #         new_indptr[np.hstack([np.sort(annots_ixs),len(annots)])]=out.counts.indptr
    #         # new_indptr[-1]=len(annots_ixs)
    #         for i in range(len(annots),0,-1):
    #             if new_indptr[i]==0:
    #                 new_indptr[i]=new_indptr[i+1]

            

         
    #         out.counts=sparse.csr_matrix((out.counts.data, out.counts.indices, new_indptr), shape=(len(annots), out.nbins))
            
    #         out.annot_dict={k:ix for ix, k in enumerate(annots)}
    #         out.nannot=len(annots)
    #         out.annot_LUT=[k for k in annots]
    #         out.annot_chr=annots_chr.copy()
        
    #     else:

    #         out.annot_dict={k:ix for ix, k in enumerate(annots_in)}
    #         out.nannot=len(annots_in)
    #         out.annot_LUT=annots_in
            
    #         out.annot_chr=out.annot_chr[annots_ixs]
        

    #     if not inplace:
    #         return out

    def bed_to_curtain(self, bedfile, curtain=None, stranded=True, strand_id=-1, feature_id=8, integrate=True, position='left'):

        if curtain is None:
            curtain=np.arange(0,3000).astype(int)
            
        binsize=self.binsize

        if strand_id==-1:
            stranded=False

        if self.iscompressed():
            chr_blocks_refframe=np.hstack([0,np.cumsum(np.ceil(self.chr_len_vec/self.binsize).astype(np.int64))])
            chr_offset_dict={self.chr_LUT[i]:v for i, v in enumerate(chr_blocks_refframe[0:-1])}
            chr_max_pos_dict={self.chr_LUT[i]:v for i, v in enumerate(chr_blocks_refframe[1:])}
        else:
            chr_offset_dict={self.chr_LUT[i]:v for i, v in enumerate(self.chr_blocks[0:-1])}
            chr_max_pos_dict={self.chr_LUT[i]:v for i, v in enumerate(self.chr_blocks[1:])}

        feature_pos=[]
        feature_strand=[]
        feature_name=[]
        feature_chr=[]
        blocks=[]
        ix=0

        this_chr=""
        new_chr_blocks=np.zeros(self.nchr+1)

        data=[]

        pos_offset=1 if integrate else 0

        with open(bedfile) as file:
            reader = csv.reader(file, delimiter="\t")
            iterator = iter(reader)
            
            for r in iterator:
                if r:
                    blocks.append(ix)
                    if r[0]!=this_chr: #we are starting a new chr
                        new_chr_blocks[self.chr_dict[r[0]]]=ix
                        this_chr=r[0]


                    this_pos=int(int(r[1])/binsize)+chr_offset_dict[r[0]]
                    
                    if (stranded and r[strand_id]=="-"):
                        this_pos+=pos_offset
                        cmin=-chr_max_pos_dict[r[0]]+this_pos  
                        cmax=int(int(r[1])/binsize)+1 #=-chr_offset_dict[r[0]]+this_pos+1
                        this_curtain=-curtain[(curtain<cmax) & (curtain>cmin)]+this_pos
                        this_curtain_list=this_curtain.tolist()
                        this_curtain_list.reverse()
                        data+=this_curtain_list
                        feature_strand.append(False)
                        feature_pos.append(this_pos-pos_offset)
                        
                    else:
                        cmin=-int(int(r[1])/binsize)-1 #=chr_offset_dict[r[0]]-this_pos-1
                        cmax=chr_max_pos_dict[r[0]]-this_pos
                        this_curtain=curtain[(curtain<cmax) & (curtain>cmin)]+this_pos
                        
                        
                        data+=this_curtain.tolist()
                        feature_strand.append(True)
                        feature_pos.append(this_pos)
        
                    feature_name+=[r[feature_id]]
                    feature_chr+=[self.chr_dict[r[0]]]
                    ix+=max(0,(len(this_curtain)-pos_offset))
                    # ix+=len(this_curtain)
                    
        
        blocks.append(ix)
        blocks=np.array(blocks)
        
        new_chr_blocks[0]=0
        new_chr_blocks[-1]=ix

        for i in range(self.nchr,0,-1):
            if new_chr_blocks[i]==0:
                new_chr_blocks[i]=new_chr_blocks[i+1]

        if integrate:
            w=np.diff(curtain)*self.binsize
            x=curtain[0:-1]*self.binsize
        else:
            w=self.binsize
            x=curtain*self.binsize

        return {'data':data, 'blocks':blocks, 'feature_name':feature_name, 'feature_chr':feature_chr, 'feature_pos':feature_pos, 'feature_strand':feature_strand, 'chr_blocks':new_chr_blocks, 'x':x, 'w':w, 'integrate':integrate}

    def df_to_curtain(self, df0, curtain=None, stranded=True, feature_info_id=None, integrate=True, position='left'):

        df = df0
        if feature_info_id is None:
            #print('NONE')
            df=df.copy()
            df['FEATURE_ID'] = df.index.to_list()
            feature_info_id = 'FEATURE_ID'


        if curtain is None:
            curtain=np.arange(0,3000).astype(int)
            
        binsize=self.binsize


        if self.iscompressed():
            chr_blocks_refframe=np.hstack([0,np.cumsum(np.ceil(self.chr_len_vec/self.binsize).astype(np.int64))])
            chr_offset_dict={self.chr_LUT[i]:v for i, v in enumerate(chr_blocks_refframe[0:-1])}
            chr_max_pos_dict={self.chr_LUT[i]:v for i, v in enumerate(chr_blocks_refframe[1:])}
        else:
            chr_offset_dict={self.chr_LUT[i]:v for i, v in enumerate(self.chr_blocks[0:-1])}
            chr_max_pos_dict={self.chr_LUT[i]:v for i, v in enumerate(self.chr_blocks[1:])}

        feature_pos=[]
        feature_strand=[]
        feature_name=[]
        feature_chr=[]
        blocks=[]
        ix=0

        this_chr=""
        new_chr_blocks=np.zeros(self.nchr+1)

        data=[]

        pos_offset=1 if integrate else 0

        for ensg, r in df.iterrows():
            blocks.append(ix)
            if r['chr']!=this_chr: #we are starting a new chr
                new_chr_blocks[self.chr_dict[r['chr']]]=ix
                this_chr=r[0]

            this_pos=int(int(r['start'])/binsize)+chr_offset_dict[r['chr']]
            
            if (stranded and r['strand']=="-"):
                this_pos+=pos_offset
                cmin=-chr_max_pos_dict[r[0]]+this_pos  
                cmax=int(int(r['start'])/binsize)+1 #=-chr_offset_dict[r[0]]+this_pos+1
                this_curtain=-curtain[(curtain<cmax) & (curtain>cmin)]+this_pos
                this_curtain_list=this_curtain.tolist()
                this_curtain_list.reverse()
                data+=this_curtain_list
                feature_strand.append(False)
                feature_pos.append(this_pos-pos_offset)
                
            else:
                cmin=-int(int(r['start'])/binsize)-1 #=chr_offset_dict[r[0]]-this_pos-1
                cmax=chr_max_pos_dict[r['chr']]-this_pos
                this_curtain=curtain[(curtain<cmax) & (curtain>cmin)]+this_pos
                data+=this_curtain.tolist()
                feature_strand.append(True)
                feature_pos.append(this_pos)

            feature_name+=[r[feature_info_id]]
            feature_chr+=[self.chr_dict[r['chr']]]
            ix+=max(0,(len(this_curtain)-pos_offset))
                # ix+=len(this_curtain)
                    
        
        blocks.append(ix)

        

        blocks=np.array(blocks)
        
        this_chr_ix = self.chr_dict[this_chr]
        new_chr_blocks[this_chr_ix+1]=ix
        new_chr_blocks[0]=0
        new_chr_blocks[-1]=ix

        for i in range(self.nchr-1,0,-1):
            if (new_chr_blocks[i]==0) and (new_chr_blocks[0:i].sum()>0):
                new_chr_blocks[i]=new_chr_blocks[i+1]

        # new_chr_blocks[0]=0
        # new_chr_blocks[-1]=ix

        # for i in range(self.nchr,0,-1):
        #     if new_chr_blocks[i]==0:
        #         new_chr_blocks[i]=new_chr_blocks[i+1]

        if integrate:
            w=np.diff(curtain)*self.binsize
            x=curtain[0:-1]*self.binsize
        else:
            w=self.binsize
            x=curtain*self.binsize

        return {'data':data, 'blocks':blocks, 'feature_name':feature_name, 'feature_chr':feature_chr, 'feature_pos':feature_pos, 'feature_strand':feature_strand, 'chr_blocks':new_chr_blocks, 'x':x, 'w':w, 'integrate':integrate}


    def bed_to_shelf(self, bedfile, feature_info_id=4, feature_id=8):
        binsize=self.binsize
        print('ok')
        if self.iscompressed():
            chr_blocks_refframe=np.hstack([0,np.cumsum(np.ceil(self.chr_len_vec/self.binsize).astype(np.int64))])
            chr_offset_dict={self.chr_LUT[i]:v for i, v in enumerate(chr_blocks_refframe[0:-1])}
            # chr_max_pos_dict={self.chr_LUT[i]:v for i, v in enumerate(chr_blocks_refframe[1:])}
        else:
            chr_offset_dict={self.chr_LUT[i]:v for i, v in enumerate(self.chr_blocks[0:-1])}
            # chr_max_pos_dict={self.chr_LUT[i]:v for i, v in enumerate(self.chr_blocks[1:])}


        feature_name=[]
        feature_info=[]
        feature_chr=[]

        ix=0

        this_chr=""
        new_chr_blocks=np.zeros(self.nchr+1, int)

        start_pos=[]
        stop_pos=[]
        feature_chr=[]
        this_offset=0
        this_chr_ix=-1

        with open(bedfile) as file:
            reader = csv.reader(file, delimiter="\t")
            iterator = iter(reader)
            
            for r in iterator:
                if r:
                    if r[0]!=this_chr: #we are starting a new chr
                        this_chr=r[0]
                        this_chr_ix=self.chr_dict[r[0]]
                        new_chr_blocks[this_chr_ix]=ix
                        this_offset=chr_offset_dict[r[0]]
                        

                    start_pos.append(int(int(r[1])/binsize)+this_offset)
                    stop_pos.append(int(int(r[2])/binsize)+this_offset)

                    feature_name+=[r[feature_id]]
                    feature_info+=[r[feature_info_id]]
                    feature_chr+=[this_chr_ix]

                    ix+=1

        start_pos_compressed=np.searchsorted(self.x, np.array(start_pos))
        stop_pos_compressed=np.searchsorted(self.x, np.array(stop_pos))

        delta=stop_pos_compressed-start_pos_compressed
#         print(this_offset)
        indptr=np.hstack([0,np.cumsum(delta)])
        n=len(delta)
        indices=np.zeros(indptr[-1])
        for i in range(n):
            if indptr[i+1]>indptr[i]:
                indices[indptr[i]:indptr[i+1]]=np.arange(start_pos_compressed[i],stop_pos_compressed[i])
        data=np.ones(len(indices), int)

        M=sparse.csc_matrix((data,indices,indptr), shape=(len(self.x),n))

        new_chr_blocks[this_chr_ix+1]=ix
        new_chr_blocks[0]=0
        new_chr_blocks[-1]=ix

        for i in range(self.nchr-1,0,-1):
            if (new_chr_blocks[i]==0) and (new_chr_blocks[0:i].sum()>0):
                new_chr_blocks[i]=new_chr_blocks[i+1]

        return {'counting_matrix':M, 'feature_name':feature_name, 'feature_info':feature_info, 'chr_blocks':new_chr_blocks, 'w':delta, 'feature_start':np.array(start_pos), 'feature_stop':np.array(stop_pos), 'feature_chr':np.array(feature_chr), 'binsize':binsize}

        #return {'feature_name':feature_name, 'feature_chr':np.array(feature_chr), 'feature_start':np.array(start_pos), 'feature_stop':np.array(stop_pos), 'chr_blocks':chr_blocks, 'binsize':binsize}


    def df_to_shelf(self, df0, feature_info_id=None): #'name'

        df = df0
        if feature_info_id is None:
            #print('NONE')
            df=df.copy()
            df['FEATURE_ID'] = df.index.to_list()
            feature_info_id = 'FEATURE_ID'

        binsize=self.binsize
        #print('ok')
        if self.iscompressed():
            chr_blocks_refframe=np.hstack([0,np.cumsum(np.ceil(self.chr_len_vec/self.binsize).astype(np.int64))])
            chr_offset_dict={self.chr_LUT[i]:v for i, v in enumerate(chr_blocks_refframe[0:-1])}
            # chr_max_pos_dict={self.chr_LUT[i]:v for i, v in enumerate(chr_blocks_refframe[1:])}
        else:
            chr_offset_dict={self.chr_LUT[i]:v for i, v in enumerate(self.chr_blocks[0:-1])}
            # chr_max_pos_dict={self.chr_LUT[i]:v for i, v in enumerate(self.chr_blocks[1:])}


        feature_name=[]
        feature_info=[]
        feature_chr=[]

        ix=0

        this_chr=""
        new_chr_blocks=np.zeros(self.nchr+1, int)

        start_pos=[]
        stop_pos=[]
        feature_chr=[]
        this_offset=0
        this_chr_ix=-1

        
        for ensg, r in df.iterrows():
            
            if r['chr']!=this_chr: #we are starting a new chr
                if ix>0:
                    new_chr_blocks[this_chr_ix+1]=ix
                    
                this_chr=r['chr']
                this_chr_ix=self.chr_dict[r['chr']]
                
                new_chr_blocks[this_chr_ix]=ix
                this_offset=chr_offset_dict[r['chr']]
                

            start_pos.append(int(int(r['start'])/binsize)+this_offset)
            stop_pos.append(int(int(r['stop'])/binsize)+this_offset)

            feature_name+=[ensg]
            feature_info+=[r[feature_info_id]]
            feature_chr+=[this_chr_ix]

            ix+=1

        start_pos_compressed=np.searchsorted(self.x, np.array(start_pos))
        stop_pos_compressed=np.searchsorted(self.x, np.array(stop_pos))

        delta=stop_pos_compressed-start_pos_compressed
#         print(this_offset)
        indptr=np.hstack([0,np.cumsum(delta)])
        n=len(delta)
        indices=np.zeros(indptr[-1])
        for i in range(n):
            if indptr[i+1]>indptr[i]:
                indices[indptr[i]:indptr[i+1]]=np.arange(start_pos_compressed[i],stop_pos_compressed[i])
        data=np.ones(len(indices), int)

        M=sparse.csc_matrix((data,indices,indptr), shape=(len(self.x),n))

        new_chr_blocks[this_chr_ix+1]=ix
        new_chr_blocks[0]=0
        new_chr_blocks[-1]=ix

        for i in range(self.nchr-1,0,-1):
            if (new_chr_blocks[i]==0) and (new_chr_blocks[0:i].sum()>0):
                new_chr_blocks[i]=new_chr_blocks[i+1]

        return {'counting_matrix':M, 'feature_name':feature_name, 'feature_info':feature_info, 'chr_blocks':new_chr_blocks, 'w':delta, 'feature_start':np.array(start_pos), 'feature_stop':np.array(stop_pos), 'feature_chr':np.array(feature_chr), 'binsize':binsize, 'df':df}
   


    def bed_to_gi(self, bedfile, strand_id=4, feature_id=8):
        binsize=self.binsize
        if self.iscompressed():
            chr_blocks_refframe=np.hstack([0,np.cumsum(np.ceil(self.chr_len_vec/self.binsize).astype(np.int64))])
            chr_offset_dict={self.chr_LUT[i]:v for i, v in enumerate(chr_blocks_refframe[0:-1])}
            # chr_max_pos_dict={self.chr_LUT[i]:v for i, v in enumerate(chr_blocks_refframe[1:])}
        else:
            chr_offset_dict={self.chr_LUT[i]:v for i, v in enumerate(self.chr_blocks[0:-1])}
            # chr_max_pos_dict={self.chr_LUT[i]:v for i, v in enumerate(self.chr_blocks[1:])}

        chr_blocks=np.zeros(len(self.chr_LUT)+1, int)

        feature_name=[]
        feature_chr=[]
        start_pos=[]
        stop_pos=[]

        ix=0

        this_chr=""

        this_chr_ix=-1

        with open(bedfile) as file:
            reader = csv.reader(file, delimiter="\t")
            iterator = iter(reader)
            
            for r in iterator:
                if r:
                    if r[0]!=this_chr: #we are starting a new chr
                        this_chr=r[0]
                        this_chr_ix=self.chr_dict[r[0]]
                        this_offset=chr_offset_dict[r[0]]
                        chr_blocks[this_chr_ix]=ix
                    
                        
                    start_pos.append(int(int(r[1])/binsize)+this_offset)
                    stop_pos.append(int(int(r[2])/binsize)+this_offset)


                    feature_name+=[r[feature_id]]
                    feature_chr+=[this_chr_ix]
                    
                    ix+=1
        chr_blocks[this_chr_ix+1]=ix
        chr_blocks[0]=0
        chr_blocks[-1]=ix

        for i in range(self.nchr-1,0,-1):
            if (chr_blocks[i]==0) and (chr_blocks[0:i].sum()>0):
                chr_blocks[i]=chr_blocks[i+1]

        return {'feature_name':feature_name, 'feature_chr':np.array(feature_chr), 'feature_start':np.array(start_pos), 'feature_stop':np.array(stop_pos), 'chr_blocks':chr_blocks, 'binsize':binsize}

    def gi_to_shelfmatrix(self, gi):
        # binsize=self.binsize

        # if self.iscompressed:
        #     chr_blocks_refframe=np.hstack([0,np.cumsum(np.ceil(self.chr_len_vec/self.binsize).astype(np.int64))])
        #     chr_offset_dict={self.chr_LUT[i]:v for i, v in enumerate(chr_blocks_refframe[0:-1])}
        #     # chr_max_pos_dict={self.chr_LUT[i]:v for i, v in enumerate(chr_blocks_refframe[1:])}
        # else:
        #     chr_offset_dict={self.chr_LUT[i]:v for i, v in enumerate(self.chr_blocks[0:-1])}
        #     # chr_max_pos_dict={self.chr_LUT[i]:v for i, v in enumerate(self.chr_blocks[1:])}


       
        assert gi['binsize']==self.binsize, "Mismatched bin size between data and GI object"
        start_pos_compressed=np.searchsorted(self.x, np.array(gi['feature_start']))
        stop_pos_compressed=np.searchsorted(self.x, np.array(gi['feature_stop']))

        delta=stop_pos_compressed-start_pos_compressed

        indptr=np.hstack([0,np.cumsum(delta)])
        n=len(delta)
        indices=np.zeros(indptr[-1])
        for i in range(n):
            if indptr[i+1]>indptr[i]:
                indices[indptr[i]:indptr[i+1]]=np.arange(start_pos_compressed[i],stop_pos_compressed[i])
        data=np.ones(len(indices), int)

        M=sparse.csc_matrix((data,indices,indptr), shape=(len(self.x),n))

        return M

    #def project_onto_gi(self, gi):
        

    @staticmethod
    def select_features(curtain, features):
        ix_in=np.array([i for i, e in enumerate(curtain['feature_name']) if e in features])
        
        chr_ids=np.searchsorted(curtain['chr_blocks'],curtain['blocks'][ix_in],side='right')-1

        blocks=curtain['blocks']+np.arange(len(curtain['blocks']))
        data_arr=np.array(curtain['data'])
        
        new_data_arr=[]
        new_blocks=np.zeros(len(ix_in)+1,int)

        this_ix=0
        new_chr_blocks=np.zeros(len(curtain['chr_blocks']),int)
        this_chr_id=-1
        for j, i in enumerate(ix_in):
            if chr_ids[j]!=this_chr_id:
                this_chr_id=chr_ids[j]
                new_chr_blocks[this_chr_id]=this_ix

            new_blocks[j]=this_ix
            this_arr=data_arr[blocks[i]:blocks[i+1]]
            new_data_arr+=this_arr.tolist()
            this_ix+=len(this_arr)-1
        
        new_blocks[-1]=this_ix
        
        new_chr_blocks[-1]=this_ix

        for i in range(len(new_chr_blocks)-2,0,-1):
            if new_chr_blocks[i]==0:
                new_chr_blocks[i]=new_chr_blocks[i+1]

        new_curtain={'data':list(new_data_arr), 'blocks':new_blocks, 'feature_name':[f for i,f in enumerate(curtain['feature_name']) if i in ix_in], 'feature_chr':[f for i,f in enumerate(curtain['feature_chr']) if i in ix_in], 'feature_pos':[f for i,f in enumerate(curtain['feature_pos']) if i in ix_in], 'feature_strand':[f for i,f in enumerate(curtain['feature_strand']) if i in ix_in], 'chr_blocks':new_chr_blocks, 'x':curtain['x'].copy(), 'w':curtain['w'].copy(), 'integrate':curtain['integrate']}

        return new_curtain

    def data_curtaining_matrix_(self, curtain):
        #creates a ncompressed x ncurtain matrix
        c=curtain['data']
        blocks=curtain['blocks']
        ncurtain=len(c)
        nblocks=len(blocks)-1
        x=self.x
        index=np.searchsorted(x, c) #size of c

        integrate=curtain['integrate']

        if not integrate:
            data=np.ones(len(index),int)
            data[~np.isin(c, x)]=0 #remove curtain points that have no data to them
            indexptr=np.arange(ncurtain+1)

            m=sparse.csc_matrix((data, index, indexptr), shape=(len(x),ncurtain))
            
        else:
            # idx_blocks=np.searchsorted(self.x,blocks)
            
            #convert back the blocks to old indexing
            blocks=curtain['blocks']+np.arange(nblocks+1)

            d=np.hstack([0,np.diff(index)])
            d[blocks[0:-1]]=0
            pre_idxptr=np.cumsum(d)

            start=blocks[0:-1]
            stop=blocks[1:]-1
            ncols=np.sum(stop-start)
            # print(ncols)
            # print(ncurtain-nblocks)
            ndata=pre_idxptr[stop[-1]]
            # print(ndata)

            data=np.ones(ndata,int)
            
            index_int=np.ones(ndata,int)
            indexptr=np.zeros(ncols+1,int)
            ntot=0
            for i in range(nblocks):
                this_start=start[i]
                this_stop=stop[i]
                n=this_stop-this_start
#                 print("%g_%g"%(pre_idxptr[this_stop]-pre_idxptr[this_start],index[this_stop]-index[this_start]))
                index_int[pre_idxptr[this_start]:pre_idxptr[this_stop]]=np.arange(index[this_start],index[this_stop])
                
                indexptr[ntot:(ntot+n)]=pre_idxptr[this_start:this_stop]
            
                ntot+=n
            indexptr[-1]=pre_idxptr[-1]
            # print(ntot)
            
            m=sparse.csc_matrix((data, index_int, indexptr), shape=(len(x),ncols))
            
            
            blocks=curtain['blocks'] #-np.arange(nblocks+1)
            # print("ncol=%g,blocks_1=%g"%(ncols, new_blocks[-1]))
            
            strandify_index=np.arange(blocks[-1])
            strandify_data=np.ones(blocks[-1], int)
            strandify_indptr=np.arange(blocks[-1]+1)
            
            for i in range(nblocks):
                if not curtain['feature_strand'][i]: #gene on negative strand
                    strandify_index[blocks[i]:blocks[i+1]]=np.arange(blocks[i+1]-1,blocks[i]-1,-1)
            
            strandify_m=sparse.csc_matrix((strandify_data, strandify_index, strandify_indptr), shape=(ncols,ncols))
            
            m=m*strandify_m
        
        return m
            
    
    def meta_matrix(self, curtain, aggregate_features=False):
        curtaining_matrix=self.data_curtaining_matrix_(curtain)
        if aggregate_features:
            aggregating_matrix=Chartable_binned.block_aggregating_matrix(curtain['blocks'])
            return curtaining_matrix*aggregating_matrix
        else:
            return curtaining_matrix


    def make_data_curtain(self, curtain, aggregate_features=False, meta_matrix=None):

        out=Chartable()
        
        out.nannot=self.nannot
        out.annot_chr=self.annot_chr.copy()
        out.annot_dict = {k:v for k, v in self.annot_dict.items()}
        out.annot_LUT=[k for k in self.annot_LUT]

        if self.counts.getformat()=='csr':
            self.counts=self.counts.tocsc()



        
        if aggregate_features:
            
            # aggregating_matrix=Chartable_binned.block_aggregating_matrix(curtain['blocks']-int(curtain['integrate'])*np.arange(len(curtain['blocks'])))
            if meta_matrix is None:
                curtaining_matrix=self.data_curtaining_matrix_(curtain)
                aggregating_matrix=Chartable_binned.block_aggregating_matrix(curtain['blocks'])
                curtaining_matrix=curtaining_matrix*aggregating_matrix
            else:
                curtaining_matrix=meta_matrix
            
            out.counts=(self.counts*curtaining_matrix).tocsr()
            out.nbins=out.counts.shape[1]
            out.bins_LUT=curtain['x'].tolist()

            out.chr_blocks=np.array([0, out.nbins])
            out.chr_dict={'none':0}
            out.chr_LUT=['none']
            out.nchr=1

        else:
            if meta_matrix is None:
                curtaining_matrix=self.data_curtaining_matrix_(curtain)
            else:
                curtaining_matrix=meta_matrix

            out.chr_blocks=curtain['chr_blocks']
            out.chr_dict={k:v for k, v in self.chr_dict.items()}
            out.chr_LUT=[k for k in self.chr_LUT]
            out.nchr=self.nchr

            out.counts=(self.counts*curtaining_matrix).tocsr()
            out.nbins=out.counts.shape[1]
            out.bins_LUT=[] #curtain['feature_name']
        
        return out



    def project(self, shelf, average=False):
        # note that projection count all reads in [start, stop[ with start and stop defined in the bed file
        out=Chartable()

        out.nannot=self.nannot
        out.annot_chr=self.annot_chr.copy()
        out.annot_dict = {k:v for k, v in self.annot_dict.items()}
        out.annot_LUT=[k for k in self.annot_LUT]

        
        new_chr_blocks=shelf['chr_blocks']

        out.chr_blocks=new_chr_blocks
        out.chr_dict={k:v for k, v in self.chr_dict.items()}
        out.chr_LUT=[k for k in self.chr_LUT]
        out.nchr=self.nchr
        out.chr_len_vec=self.chr_len_vec.copy()
        if average:
            counts=self.counts.tocsc()
            counts_norm=counts.copy()
            counts_norm.eliminate_zeros()
            counts_norm.data=1.0*np.ones(len(counts_norm.data))
            #x=(1.0*np.ones(len(counts_norm.data)))/len(counts_norm.data)

            M_norm=counts_norm*shelf['counting_matrix']
            
            M=counts*shelf['counting_matrix']
            M.data=M.data/M_norm.data
        
            out.counts=M.tocsr()
        else:
            out.counts=(self.counts.tocsc()*shelf['counting_matrix']).tocsr()

        out.nbins=out.counts.shape[1]
        out.bins_LUT=shelf['feature_name']

        return out

    def predict_pr(self, bkg, pr):

        out=Chartable()
        out.nannot=self.nannot
        out.annot_chr=self.annot_chr.copy()
        out.annot_dict = {k:v for k, v in self.annot_dict.items()}
        out.annot_LUT=[k for k in self.annot_LUT]


        new_chr_blocks=pr['chr_blocks']

        out.chr_blocks=new_chr_blocks
        out.chr_dict={k:v for k, v in self.chr_dict.items()}
        out.chr_LUT=[k for k in self.chr_LUT]
        out.nchr=self.nchr
        out.chr_len_vec=self.chr_len_vec.copy()
        
        N=self.get_chr_matrix(withcounts=True)
        B=bkg.normalize_row_bychr().counts
        
        out.counts=(N*(B*pr['counting_matrix'])).tocsr()

        out.nbins=out.counts.shape[1]
        out.bins_LUT=pr['feature_name']

        return out
    

    def make_giset_projector(self, gi_set):
        counting_matrix=(self.gi_to_shelfmatrix(gi_set))*gi_set['summing_matrix']

        new_chr_blocks=counting_matrix.shape[1]*np.ones(len(self.chr_blocks),int)
        new_chr_blocks[0]=0
        gi_set['chr_blocks']=new_chr_blocks

        return {'chr_blocks': new_chr_blocks, 'feature_name': gi_set['set_names'], 'counting_matrix':counting_matrix} 

    # def load_bam(self, bamfile):



    #     reader = pysam.AlignmentFile(bamfile)
    #     bamiter = reader.fetch(until_eof=True)

    #     try:
    #         r = bamiter.__next__()
    #     except StopIteration: 
    #         r = None
    #     while not(r==None):

    #     data=np.zeros(len(self.x), float)
    #     with open(bedfile) as file:
    #         reader = csv.reader(file, delimiter="\t")
    #         iterator = iter(reader)
    #         ix=-1
    #         for r in iterator:
    #             if r:
    #                 ix+=1
    #                 data[ix]=float(r[4])

    #     return sparse.csr_matrix(data, shape=(1,len(data)))


    def load_phased_bed(self, bedfile):
        # out=self.copy(copydata=False)
        # out.x=self.x.copy()

        data=np.zeros(len(self.x), float)
        with open(bedfile) as file:
            reader = csv.reader(file, delimiter="\t")
            iterator = iter(reader)
            ix=-1
            for r in iterator:
                if r:
                    ix+=1
                    data[ix]=float(r[4])

        return sparse.csr_matrix(data, shape=(1,len(data)))

    def load_multi_phased_bed(self, bedfiles, labels):
        data=sparse.vstack([self.load_phased_bed(f) for f in bedfiles])
        out=self.copy(copydata=False)
        out.annot_LUT=[l for l in labels]
        out.annot_dict={l:i for i, l in enumerate(labels)}
        out.annot_chr=np.zeros(len(labels))

        out.counts=data
        out.x=self.x.copy()
        #out.iscompressed=self.iscompressed

        return out

    def make_bed(self,f):
        #TODO: hanlde case where not compressed
        out=open(f,"w")
        d=self
        blcks=d.chr_blocks
        binsize=d.binsize
        chr_offset_vec=np.cumsum(np.hstack([0,d.chr_len_vec]))
        ii=-1
        for ch_ix, ch in enumerate(d.chr_LUT):
        
            x=d.x[blcks[ch_ix]:blcks[(ch_ix+1)]].astype(np.int64)
            for i, p in enumerate(x):
                ii+=1
                this_p=p*binsize-chr_offset_vec[ch_ix]
                out.write("%s\t%d\t%d\t%d\n"%(ch,this_p,this_p+binsize,ii))

        out.flush()
        out.close()

    def to_bigwig(self,filename, normalize=False):
        import pyBigWig as pbw
        if self.nannot>1:
            print('Noting done, this dataset has multiple features')
        else:
            

            #chr_off=np.cumsum(np.hstack([0,self.chr_len_vec]))
            chr_blocks=self.chr_blocks

            chr_ids=np.searchsorted(chr_blocks,np.arange(len(self.x)),side='right')-1

            chr_offset_coarseframe=np.hstack([0,np.cumsum(np.ceil(self.chr_len_vec/self.binsize).astype(int))])
            #x_real=np.array((np.array(self.x, dtype=np.int64)*self.binsize)-chr_off[chr_ids], dtype=int)
            x_real=np.array((np.array(self.x, dtype=np.int64)-chr_offset_coarseframe[chr_ids])*self.binsize, dtype=int)
            v=self.counts.toarray().flatten()

            bwheader=[(k,self.chr_len_vec[i]) for i,k in enumerate(self.chr_LUT)]
            bw = pbw.open(filename, "w")
            bw.addHeader(bwheader)

            normfactor=v.sum()*1.0 if normalize else 1.0
            
            for i in range(self.nchr):
                if chr_blocks[i+1]>chr_blocks[i]:
                    bw.addEntries(self.chr_LUT[i],x_real[chr_blocks[i]:chr_blocks[i+1]].tolist(),values=v[chr_blocks[i]:chr_blocks[i+1]]/normfactor,span=1)

            bw.close()
        
    def compress_at(self, others, mincount=1):
        idx = np.unique(sparse.find(self.counts)[1], return_counts=True)
        ix=idx[0]
        y=idx[1]
        if mincount>1:
            ix=ix[y>=mincount]
            y=y[y>=mincount]
        self.counts=self.counts[:,ix]
        self.x=self.x[ix]

        chr_ids=np.searchsorted(self.chr_blocks,ix,side='right')-1
        self.chr_blocks=np.searchsorted(chr_ids,np.arange(self.nchr+1))
        self.nbins=len(self.x)

        for o in others:
            o.x=o.x[ix]
            o.counts=o.counts[:,ix]
            o.chr_blocks=self.chr_blocks.copy()


    def to_bigwig_wide(self,xL, xR, filename, normalize=False):
        import pyBigWig as pbw
        if self.nannot>1:
            print('Noting done, this dataset has multiple features')
        else:
            

            chr_off=np.cumsum(np.hstack([0,self.chr_len_vec]))
            chr_blocks=self.chr_blocks

            chr_ids=np.searchsorted(chr_blocks,np.arange(len(self.x)),side='right')-1
            x_real=np.array((np.array(self.x, dtype=np.int64)*self.binsize)-chr_off[chr_ids], dtype=int)

            v=self.counts.toarray().flatten()
            xL_real=np.array(np.array(xL.counts.toarray().flatten(), dtype=np.int64)*self.binsize-chr_off[chr_ids], dtype=int)
            xR_real=np.array(np.array(xR.counts.toarray().flatten(), dtype=np.int64)*self.binsize-chr_off[chr_ids], dtype=int)

            bwheader=[(k,self.chr_len_vec[i]) for i,k in enumerate(self.chr_LUT)]
            bw = pbw.open(filename, "w")
            bw.addHeader(bwheader)

            normfactor=v.sum()*1.0 if normalize else 1.0
            
            for i in range(self.nchr):
                if chr_blocks[i+1]>chr_blocks[i]:
                    this_chr=self.chr_LUT[i]
                    these_chr=[this_chr]*(chr_blocks[i+1]-chr_blocks[i])
                    #for j in range(chr_blocks[i+1]-chr_blocks[i]):
                    bw.addEntries(these_chr,xL_real[chr_blocks[i]:chr_blocks[i+1]],ends=xR_real[chr_blocks[i]:chr_blocks[i+1]], values=v[chr_blocks[i]:chr_blocks[i+1]]/normfactor)

        bw.close()

    def sum(self, axis=0): #sums over annotations
        out=self.copy()
        if axis==0:
            out.annot_dict={'sum':'0'} # name to ix, chrix  #dictionnary
        
            out.annot_LUT=['sum']
            out.annot_chr=np.array([-1], dtype=int)
            out.nannot=1
       
            out.counts=sparse.csc_matrix(self.counts.sum(axis=0))
        else:
            out=Chartable()
            out.annot_dict={k:v for k,v in self.annot_dict.items()} # name to ix, chrix  #dictionnary
            
            out.annot_LUT=self.annot_LUT.copy()
            out.annot_chr=self.annot_chr.copy()
            out.nannot=self.nannot

            # columns are bins, blocked by chr
            out.chr_dict={'all':0}
            out.chr_blocks=np.array([0,1], int)
            out.nchr=1
            out.nbins=1
            out.bins_LUT=['all']
            
            out.chr_LUT=[0]

            out.counts=sparse.csr_matrix(self.counts.sum(axis=1))
        return out


        # out.counts=self.counts[:,curtain['data']]
        # out.nbins=len(curtain['data'])

        # if aggregate:
        #     aggmatrix, feature_name, feature_chr=Chartable_binned.make_agg_matrix(curtain)
        #     Chartable_binned.aggerate_data_curtain(out, aggmatrix, feature_name, feature_chr)
        #     return out
        # else:
        #     return out
    
    # def make_curtain_uncompress_(self, curtain, aggregate=False):
    #     #make the data curtain 
    #     out=Chartable()
    #     out.chr_blocks=curtain['chr_blocks']
    #     out.chr_dict={k:v for k, v in self.chr_dict.items()}
    #     out.chr_LUT=[k for k in self.chr_LUT]
    #     out.nchr=self.nchr
    #     out.nannot=self.nannot
    #     out.annot_chr=self.annot_chr.copy()
    #     out.annot_dict = {k:v for k, v in self.annot_dict.items()}
    #     out.annot_LUT=[k for k in self.annot_LUT]

    #     out.counts=self.counts[:,curtain['data']]
    #     out.nbins=len(curtain['data'])

    #     if aggregate:
    #         aggmatrix, feature_name, feature_chr=Chartable_binned.make_agg_matrix(curtain)
    #         Chartable_binned.aggerate_data_curtain(out, aggmatrix, feature_name, feature_chr)
    #         return out
    #     else:
    #         return out


    @staticmethod
    def block_aggregating_matrix(blocks): #return a matrix that transforms the curtain matrix to aggregate over features
        bsizes=np.diff(blocks)
        blocks2=blocks[0:-1]
        nbins=blocks[-1]
        aggmatrix_data=np.ones(nbins,int)

        bsize_max=np.max(bsizes)
        ix=0

        aggmatrix_indices=[]
        aggmatrix_indptr=np.zeros(bsize_max+1,int)
        for i in range(bsize_max):
            aggmatrix_indptr[i]=ix
            ith_ind=blocks2[bsizes>i]+i
            aggmatrix_indices+=(ith_ind).tolist()
            ix+=len(ith_ind)

        aggmatrix_indptr[bsize_max]=ix #np.hstack((bsizes,len(aggmatrix_data)))

        aggmatrix=sparse.csc_matrix((aggmatrix_data,aggmatrix_indices, aggmatrix_indptr), shape=(nbins,bsize_max))

        return aggmatrix

    @staticmethod
    def aggerate_data_curtain(out, aggmatrix, feature_name, feature_chr):
        new_chr_blocks=np.zeros(out.nchr+1, int)
        new_chr_blocks[0]=0
        this_chr=-1
        
        for ix, k in enumerate(feature_chr):
            if k!=this_chr:
                new_chr_blocks[this_chr]=ix
                this_chr=k

        new_chr_blocks[-1]=ix

        for i in range(out.nchr,0,-1):
            if new_chr_blocks[i]==0:
                new_chr_blocks[i]=new_chr_blocks[i+1]

        out.chr_blocks=new_chr_blocks
        out.bins_LUT=feature_name
        out.nbins=len(feature_name)

        out.counts=out.counts*aggmatrix

        # return out


    @staticmethod
    def features_aggregating_matrix_(curtain): #return a matrix that transforms the curtain matrix to aggregate over features
        
        bsizes=np.diff(curtain['blocks'])
        blocks=curtain['blocks'][0:-1]
        nbins=curtain['blocks'][-1]
        aggmatrix_data=np.ones(nbins,int)

        bsize_max=np.max(bsizes)
        ix=0

        aggmatrix_indices=[]
        aggmatrix_indptr=np.zeros(bsize_max+1,int)
        for i in range(bsize_max):
            aggmatrix_indptr[i]=ix
            ith_ind=blocks[bsizes>i]+i
            aggmatrix_indices+=(ith_ind).tolist()
            ix+=len(ith_ind)

        aggmatrix_indptr[bsize_max]=ix #np.hstack((bsizes,len(aggmatrix_data)))

        aggmatrix=sparse.csc_matrix((aggmatrix_data,aggmatrix_indices, aggmatrix_indptr), shape=(nbins,bsize_max))

        # featture=[[name, chro, pos, strand] for name, chro, pos, strand in zip(curtain['feature_name'],curtain['feature_chr'],curtain['feature_pos'],curtain['feature_strand'])]

        return aggmatrix

    @staticmethod
    def make_agg_matrix(curtain):
        bsizes=np.diff(curtain['blocks'])
        blocks=curtain['blocks'][0:-1]
        nbins=curtain['blocks'][-1]
        aggmatrix_data=np.ones(nbins,int)

        bsize_max=np.max(bsizes)
        ix=0

        aggmatrix_indices=[]
        aggmatrix_indptr=np.zeros(bsize_max+1,int)
        for i in range(bsize_max):
            aggmatrix_indptr[i]=ix
            ith_ind=blocks[bsizes>i]+i
            aggmatrix_indices+=(ith_ind).tolist()
            ix+=len(ith_ind)

        aggmatrix_indptr[bsize_max]=ix #np.hstack((bsizes,len(aggmatrix_data)))

        aggmatrix=sparse.csc_matrix((aggmatrix_data,aggmatrix_indices, aggmatrix_indptr), shape=(nbins,bsize_max))

        # featture=[[name, chro, pos, strand] for name, chro, pos, strand in zip(curtain['feature_name'],curtain['feature_chr'],curtain['feature_pos'],curtain['feature_strand'])]

        return aggmatrix, curtain['feature_name'],curtain['feature_chr']


    def get_source_target_chr_data(self): #return the chromosome id of the source

        sourcechr_data=np.zeros(len(self.counts.data), int)
        
        indptr=self.counts.indptr
        annot_chr=self.annot_chr
        for i in range(self.nannot):
            sourcechr_data[indptr[i]:indptr[i+1]]=annot_chr[i]
        
        # chr_LUT=self.bin_to_chr_LUT()
        
        targetchr_data=np.searchsorted(self.chr_blocks,self.counts.indices,side='right')-1

        return sourcechr_data, targetchr_data

    def cis_trans(self, format='pandas', annotate=None):
        trans=self.copy(copydata=True)
        sourcechr_data, targetchr_data = self.get_source_target_chr_data()
        
        counts=trans.counts.data
        counts[sourcechr_data==targetchr_data]=0

        cis=self.copy(copydata=False)
        cis.x=self.x.copy()
        
        cis.counts=sparse.csr_matrix((self.counts.data-counts, self.counts.indices.copy(), self.counts.indptr.copy()), shape=(self.nannot, self.nbins))

        cis.counts.eliminate_zeros()
        trans.counts.eliminate_zeros()

        N_cis=np.array(cis.counts.sum(axis=1)).ravel()
        N_trans=np.array(trans.counts.sum(axis=1)).ravel()
        N=np.column_stack([N_cis,N_trans])

        if ((format=='pandas') or (not annotate is None)) :
            N=pd.DataFrame(N, columns=['Ncis','Ntrans'], index=self.annot_LUT)
            if not annotate is None:
                N=pd.concat([annotate.reindex(N.index), N], axis=1)

        return cis, trans, N


    def bed_to_matrix(self, bedfile, curtain=np.arange(0,3000).astype(int), reverse=False, stranded=True, feature_id=8):

        binsize=self.binsize
        chr_offset_dict={self.chr_LUT[i]:v for i, v in enumerate(self.chr_blocks[0:-1])}
        chr_max_pos_dict={self.chr_LUT[i]:v for i, v in enumerate(self.chr_blocks[1:])}

        feature_pos=[]
        feature_strand=[]
        feature_name=[]
        feature_chr=[]
        blocks=[]
        ix=0

        this_chr=""
        new_chr_blocks=np.zeros(self.nchr+1)

        data=[]


        with open(bedfile) as file:
            reader = csv.reader(file, delimiter="\t")
            iterator = iter(reader)
            
            for r in iterator:
                if r:
                    blocks.append(ix)
                    if r[0]!=this_chr: #we are starting a new chr
                        new_chr_blocks[self.chr_dict[r[0]]]=ix


                    this_pos=int(int(r[1])/binsize)+chr_offset_dict[r[0]]
                    
                    if (r[4]=="-" and stranded):
                        cmin=-chr_max_pos_dict[r[0]]+this_pos  
                        cmax=int(int(r[1])/binsize)+1 #=-chr_offset_dict[r[0]]+this_pos+1
                        this_curtain=-curtain[(curtain<cmax) & (curtain>cmin)]+this_pos
                        data+=this_curtain.tolist()
                        feature_strand.append(False)
                        
                    else:
                        cmin=-int(int(r[1])/binsize)-1 #=chr_offset_dict[r[0]]-this_pos-1
                        cmax=chr_max_pos_dict[r[0]]-this_pos
                        this_curtain=curtain[(curtain<cmax) & (curtain>cmin)]+this_pos
                        data+=this_curtain.tolist()
                        feature_strand.append(True)
                    
                    

                    feature_name+=[r[feature_id]]
                    feature_chr+=[self.chr_dict[r[0]]]
                    ix+=len(this_curtain)
                    feature_pos.append(this_pos)
        
        blocks.append(ix)
        
        blocks=np.array(blocks)
        
        
        new_chr_blocks[0]=0
        new_chr_blocks[-1]=ix

        for i in range(self.nchr,0,-1):
            if new_chr_blocks[i]==0:
                new_chr_blocks[i]=new_chr_blocks[i+1]

        # out=Chartable()
        # out.chr_blocks=new_chr_blocks
        # out.chr_dict={k:v for k: v in self.chr_dict.items()}
        # out.chr_LUT=self.chr_LUT.copy()
        # out.nchr=self.nchr

        return {'data':data, 'blocks':blocks, 'feature_name':feature_name, 'feature_chr':feature_chr, 'feature_pos':feature_pos, 'feature_strand':feature_strand, 'chr_blocks':new_chr_blocks}

    def get_block_summator_bychr(self, blocks=None):
        if blocks is None:
            data=np.ones(self.nbins, int)
            indptr=self.chr_blocks.copy() 
            indices=np.arange(self.nbins)
            M=sparse.csc_matrix((data,indices,indptr),shape=(self.nbins,self.nchr))
        #TODO: verify this else statement
        else:
            ndata=blocks[-1]-blocks[0]
            data=np.ones(ndata, int)
            indptr=blocks.copy()
            indices=np.arange(ndata)
            M=sparse.csc_matrix((data,indices,indptr),shape=(self.nbins,len(blocks)-1))
           
        return M

    def N_bychr(self):
        M=self.get_block_summator_bychr()
        
        if self.counts.format=='csr':
            N=self.counts*(M.tocsr())
        else:
            N=self.counts*M

        return N

    

    def tocsc(self):
        if self.counts.getformat()!='csc':
            self.counts=self.counts.tocsc()
        else:
            pass

    def tocsr(self):
        if self.counts.getformat()!='csr':
            self.counts=self.counts.tocsr()
        else:
            pass

    def normalize_row_bychr(self, blocks=None, inplace=False, expand=True):
        # 
        assert self.nannot==1, "More than one row in count matrix!"
        if inplace:
            out=self
        else:
            out=self.copy(copydata=True)
        
        if blocks is None:
            blocks=self.chr_blocks
        
        M=self.get_block_summator_bychr(blocks=blocks)
        nblocks=M.shape[1]
        
        out.tocsc()
        scale=(out.counts*M).toarray().ravel()
        scale[scale==0]=1.0
  
        # print(out.counts.sum())
        # data=out.counts.data
        indptr=out.counts.indptr
       
        data=out.counts.data.copy()*1.0

        for i in range(nblocks):
            data[indptr[blocks[i]]:indptr[blocks[i+1]]]=data[indptr[blocks[i]]:indptr[blocks[i+1]]]/scale[i]

        out.counts=sparse.csc_matrix((data, out.counts.indices,out.counts.indptr), shape=(1,self.nbins))

        if expand:
            out.tocsr()

            new_indptr=np.searchsorted(out.counts.indices, blocks)
            out.counts=sparse.csr_matrix((out.counts.data, out.counts.indices, new_indptr), shape=(nblocks, out.nbins))

        return out

    def normalize_row(self, inplace=False):
        # 
        assert self.nannot==1, "More than one row in count matrix!"
        if inplace:
            out=self
        else:
            out=self.copy(copydata=True)
        
        out.tocsc()
       
        scale=out.counts.sum()
        #scale[scale==0]=1.0
  
        # print(out.counts.sum())
        # data=out.counts.data
        indptr=out.counts.indptr
       
        data=(out.counts.data.copy()*1.0)/scale
        out.counts=sparse.csc_matrix((data, out.counts.indices,out.counts.indptr), shape=(1,self.nbins))

        return out

    def normalize_row_bytrans(self, blocks=None, inplace=False):
        assert self.nannot==1, "More than one row in count matrix!"
        
        if blocks is None:
            blocks=self.chr_blocks

        nblocks=len(blocks)-1

        if inplace:
            out=self
        else:
            out=self.copy(copydata=True)
        
        out.tocsr()
    
        data=1.0*np.tile(out.counts.data,nblocks)
        indices=np.tile(out.counts.indices,nblocks)
        new_indptr=np.zeros(nblocks+1)
        ndata=len(out.counts.data)

        for i in range(nblocks):
            offst=ndata*i
            start=offst+np.searchsorted(out.counts.indices,blocks[i])
            stop=offst+np.searchsorted(out.counts.indices,blocks[i+1])
            data[start:stop]=0
            x=data[offst:(offst+ndata)]
            data[offst:(offst+ndata)]=x/x.sum()

            new_indptr[i]=offst
        
        new_indptr[-1]=ndata*nblocks
            # this_M=out.counts.copy()
            # this_M.data[indptr[blocks[i]]:indptr[blocks[i+1]]]=0
            # this_M.data=(1.0*this_M.data)/this_M.sum()
            # M_to_stack+=[this_M.tocsr()]
        
        out.counts=sparse.csr_matrix((data,indices,new_indptr), shape=(nblocks, self.nbins))
        out.counts.eliminate_zeros()
        # out.counts=sparse.hstack(, format='csr')

        return out

    def get_chr_matrix(self, withcounts=False):
        ngenes=self.counts.shape[0]
        
        if withcounts:
            data=np.array(self.counts.sum(axis=1)).flatten()
        else:
            data=np.ones(ngenes, int)

        M=sparse.csr_matrix((data, self.annot_chr, np.arange(ngenes+1)), shape=((ngenes,self.nchr)))
        return M


    def predict_trans_bytrans(self, bkg, leakage_vector=None):
        #make sure to pass expected bkg_model as a nchr X p matrix block diagonal
        N=self.get_chr_matrix(withcounts=True)
        B=bkg.normalize_row_bytrans().counts
        if not leakage_vector is None:
            L=sparse.diags(leakage_vector, format='csr')
        else: #we use existing trans counts
            L=1
        out=(N*L)*B

        return out

    def predict_cis_withinGI(self, cismodel, gi, gene_strand_dict=None, return_chartable=True):
        start=gi['feature_start']
        stop=gi['feature_stop']
        chrom=gi['feature_chr']
        n_gi=len(start)

        self.tocsr()
        counts=self.counts
        chr_len_vec=self.chr_len_vec
        chr_offset_vec=np.hstack([0,chr_len_vec.cumsum()])[0:-1]

        # this_chr=""
        chrid=-1

        bins=cismodel['bins'] #-1 ?
        #d=np.diff(bins)
        b=bins-1
        #nb=len(bins)
        chi=cismodel['chi'][0]
        rho=cismodel['rho'][0]
        edges=np.hstack([cismodel['bins'],cismodel['bins'][-1]+1])

        out_counts=[]
        
        directional = False if gene_strand_dict is None else True
        dirflag=1

        # if not bkg is None:
        #     bkg_offset_vec={bkg.chr_LUT[i]:v for i, v in enumerate(bkg.chr_blocks[0:-1])}

        #     B=bkg.normalize_row_bychr(expand=False)
        #     if gi_bkg is None:
        #         gi_bkg=gi
        #     M=B.gi_to_shelfmatrix(gi_bkg)
        #     bgk_weights=B.counts*M #csc

        # else:
        #     bkg_weigths=1

            

        
        for i_annot in tqdm(range(self.nannot)):
            y=np.zeros(n_gi,float)
            # if self.annot_chr[i_annot]!=chrid:
            
            # rho=cismodel['rho'][chrid].copy()

            if directional and (gene_strand_dict.get(self.annot_LUT[i_annot], False)): #TRUE means negative strand
                dirflag=-1
            else:
                dirflag=1
            
             #we need to renormalize the densities here to account for finite chromosome size and RNA source position

            
            chrid=self.annot_chr[i_annot]
            chr_len=chr_len_vec[chrid]

            rho=cismodel['rho'][chrid] #.copy()
            chi=cismodel['chi'][chrid]  
            # dchi=cismodel['dchi'][chrid].copy()
            # this_annot_pos=annot_pos_dict[self.annot_LUT[i_annot]]
            # max_flight_pos=(chr_len-this_annot_pos-1) if dirflag>0 else this_annot_pos
            # min_flight_pos=(-this_annot_pos) if dirflag>0 else (-chr_len-this_annot_pos+1) 
            
            # max_flight=np.searchsorted(b,max_flight_pos)
            # min_flight=np.searchsorted(b,min_flight_pos)

            # # dchi[(max_flight+1):(nb+1)]=0
            # # dchi[0:(min_flight)]=0

            # dchi=dchi/dchi.sum()
            # chi=dchi.cumsum()

            # rho=np.hstack([0, dchi[1:-1]/d, 0])


            data=counts.data[counts.indptr[i_annot]:counts.indptr[i_annot+1]]
            pos=self.x[counts.indices[counts.indptr[i_annot]:counts.indptr[i_annot+1]]]

            Pbkg=1
            # if not bkg is None:
            #     pos_bgk_frame=(pos-chr_offset_vec[chrid])/bkg.binsize+bkg_offset_dict[chrid]
            #     Pbkg=np.searchsorted(bkg.x, )
            

            # x1=this_annot_pos+min_flight_pos+1
            # x2=this_annot_pos+max_flight_pos-1
            x1=0
            x2=chr_len-1
            delta1=(x1-pos+chr_offset_vec[chrid]) if dirflag==1 else (pos-x2-chr_offset_vec[chrid])
            delta2=(x2-pos+chr_offset_vec[chrid]) if dirflag==1 else (pos-x1-chr_offset_vec[chrid])
            binid2=np.searchsorted(b,delta2)
            binid1=np.searchsorted(b,delta1)

            Ps=chi[binid2]-rho[binid2]*(edges[binid2]-(delta2)) - chi[binid1]+rho[binid1]*(edges[binid1]-delta1) 
            ntot=np.sum(Ps*data*Pbkg)
            # if ntot==0:
            #     print("NTOT 0")
            #     print(binid1)
            #     print(binid2)
            #     print(chrid)
            corr_factor=np.sum(data)/ntot
            # print(corr_factor)

            for i_gi in range(n_gi):
                if chrom[i_gi]==chrid:
                    x1=start[i_gi]
                    x2=stop[i_gi]
                    delta1=(x1-pos) if dirflag==1 else (pos-x2)
                    delta2=(x2-pos) if dirflag==1 else (pos-x1)
                    binid2=np.searchsorted(b,delta2) #bins?
                    binid1=np.searchsorted(b,delta1)

                    Ps=chi[binid2]-rho[binid2]*(edges[binid2]-(delta2)) - chi[binid1]+rho[binid1]*(edges[binid1]-delta1) 
                    
                    y[i_gi]=np.sum(Ps*data*Pbkg)*corr_factor
            
            out_counts+=[sparse.csr_matrix(y)]

        if return_chartable:
            out=Chartable()
            out.nannot=self.nannot
            out.annot_chr=self.annot_chr.copy()
            out.annot_dict = {k:v for k, v in self.annot_dict.items()}
            out.annot_LUT=[k for k in self.annot_LUT]

            out.chr_dict=self.chr_dict.copy()
            out.chr_LUT=[k for k in self.chr_LUT]
            out.chr_len_vec=self.chr_len_vec.copy()
            out.nchr=self.nchr

            out.counts=sparse.vstack(out_counts)
            out.nbins=out.counts.shape[1]
            out.bins_LUT=gi['feature_name']

            if not 'chr_blocks' in gi:
                new_chr_blocks=out.counts.shape[1]*np.ones(len(self.chr_blocks),int)
                new_chr_blocks[0]=0
                out.chr_blocks=new_chr_blocks
            else:
                out.chr_blocks=gi['chr_blocks']

            return out
        else:
            return sparse.vstack(out_counts)

    # def predict_cis_sparse(self, d, cismodel, gene_strand_dict=None):
    #     self.tocsr()
    #     counts=self.counts
    #     chr_len_vec=self.chr_len_vec
    #     chr_offset_vec=np.hstack([0,chr_len_vec.cumsum()])[0:-1]

    #     chrid=-1
    #     bins=cismodel['bins'] #-1 ?
        
    #     #b=bins-1
        
    #     #chi=cismodel['chi'][0]
    #     #rho=cismodel['rho'][0]
    #     #edges=np.hstack([cismodel['bins'],cismodel['bins'][-1]+1])

    #     #out_counts=[]
        
    #     directional = False if gene_strand_dict is None else True
    #     dirflag=1

    #     chr_blocks_dnaframe=d.chr_blocks
    #     #chr_len_vec=d.chr_len_vec
    #     #chr_offset_vec=np.hstack([0,chr_len_vec.cumsum()])
    #     chr_offset_dnaframe=np.hstack([0,np.cumsum(np.ceil(d.chr_len_vec/d.binsize).astype(np.int64))])


    #     outdata = np.ones(counts.sum())
    #     outindptr = np.hstack([0,np.array(counts.sum(axis=1)).ravel().cumsum()])
    #     outindices = np.zeros(outdata.shape[0], np.int64)

    #     for i_annot in tqdm(range(self.nannot)):

    #         # if self.annot_chr[i_annot]!=chrid:

    #         if directional and (gene_strand_dict.get(self.annot_LUT[i_annot], False)): #TRUE means negative strand
    #             dirflag=-1
    #         else:
    #             dirflag=1
    #          #we need to renormalize the densities here to account for finite chromosome size and RNA source position
    #         chrid=self.annot_chr[i_annot]
    #         #chr_len=chr_len_vec[chrid]
    #         chr_ix_low_dnaframe = chr_blocks_dnaframe[chrid]
    #         chr_ix_high_dnaframe = chr_blocks_dnaframe[chrid+1]
    #         #rho=cismodel['rho'][chrid] #.copy()
    #         chi=cismodel['chi'][chrid][0:-1]
    #         this_chr_rna_offset = chr_offset_vec[chrid]
           
    #         #this_n = counts.indptr[i_annot+1]-counts.indptr[i_annot]
    #         # if this_n>0:

    #         data=counts.data[counts.indptr[i_annot]:counts.indptr[i_annot+1]].copy()
    #         pos=self.x[counts.indices[counts.indptr[i_annot]:counts.indptr[i_annot+1]]]-this_chr_rna_offset
    #         ntotal=data.sum()
    #         ix_sim_dnaframe = np.zeros(ntotal, np.int64)
    #         data_nonzero = data>0
    #         ntosim=data_nonzero.sum()
    #         ncum=0
    #         #outindptr[]
    #         while ntosim>0:
                
    #             deltas_sim = np.interp(random_sample(ntosim), chi, bins)
    #             this_ix_sim_dnaframe = np.searchsorted(d.x,np.ceil((pos[data_nonzero]+dirflag*deltas_sim)/d.binsize) + chr_offset_dnaframe)
                
    #             ixok = (this_ix_sim_dnaframe>=chr_ix_low_dnaframe) & (this_ix_sim_dnaframe<chr_ix_high_dnaframe)
    #             data[data_nonzero*ixok]-=1
                
    #             thisn = ixok.sum()
    #             ix_sim_dnaframe[ncum:(ncum+thisn)] = this_ix_sim_dnaframe[ixok]
                
    #             ncum+=thisn
    #             ntosim = data.sum()
    #             data_nonzero = data>0
                    
                
    #         ix_sim_dnaframe = np.sort(ix_sim_dnaframe)
    #         outindices[outindptr[i_annot]:outindptr[i_annot+1]] = ix_sim_dnaframe

    #     out=d.copy(copydata=False)
    #     out.annot_dict={k:v for k,v in self.annot_dict.items()} # name to ix, chrix  #dictionnary
    #     out.annot_LUT=self.annot_LUT.copy()
    #     out.annot_chr=self.annot_chr.copy()
    #     out.nannot=self.nannot
        
    #     out.x= d.x.copy()
    #     out.counts=sparse.csr_matrix((outdata, outindices, outindptr), shape=d.counts.shape)

    #     return out

    def predict_cis_sparse(self, d, cismodel, bkg=None, wiggle_re=0, gene_strand_dict=None, b_prenormalized=False, isproba=False):
        self.tocsr()
        counts=self.counts
        chr_len_vec=self.chr_len_vec
        chr_offset_vec=np.hstack([0,chr_len_vec.cumsum()])[0:-1]
        
        if not bkg is None:
            if not b_prenormalized:
                B=bkg.normalize_row_bychr(expand=False).counts.toarray().ravel()
            else:
                B=bkg.counts.toarray().ravel().copy()
            
            if not isproba:
                
                for i in range(len(d.chr_blocks)-1):
                    h = B[d.chr_blocks[i]:d.chr_blocks[i+1]]
                    B[d.chr_blocks[i]:d.chr_blocks[i+1]] = h/np.max(h)
            else:
                print("isproba")

        bins=cismodel['bins'] #-1 ?
        
        directional = False if gene_strand_dict is None else True
        dirflag=1

        chr_blocks_dnaframe=d.chr_blocks
        chr_offset_dnaframe=np.hstack([0,np.cumsum(np.ceil(d.chr_len_vec/d.binsize).astype(np.int64))])


        outdata = np.ones(counts.sum())
        outindptr = np.hstack([0,np.array(counts.sum(axis=1)).ravel().cumsum()])
        outindices = np.zeros(outdata.shape[0], np.int64)

        for i_annot in tqdm(range(self.nannot)):

            if directional and (gene_strand_dict.get(self.annot_LUT[i_annot], False)): #TRUE means negative strand
                dirflag=-1
            else:
                dirflag=1
            
            chrid=self.annot_chr[i_annot]
            chr_ix_low_dnaframe = chr_blocks_dnaframe[chrid]
            chr_ix_high_dnaframe = chr_blocks_dnaframe[chrid+1]
            chi=cismodel['chi'][chrid][0:-1]
            this_chr_rna_offset = chr_offset_vec[chrid]
            this_chr_dna_offset = chr_offset_dnaframe[chrid]
            #this_n = counts.indptr[i_annot+1]-counts.indptr[i_annot]
            # if this_n>0:

            data=counts.data[counts.indptr[i_annot]:counts.indptr[i_annot+1]].copy()
            indx=np.arange(len(data))
            pos=self.x[counts.indices[counts.indptr[i_annot]:counts.indptr[i_annot+1]]]-this_chr_rna_offset
            
            ntotal=data.sum()
            ix_sim_dnaframe = np.zeros(ntotal, np.int64)
            data_nonzero = data>0
            ntosim=data_nonzero.sum()
            ncum=0
            #outindptr[]
            while ntosim>0:
                #print(ntotal)
                #print(ntosim)
                deltas_sim = np.interp(random_sample(ntosim), chi, bins)
                #print(pos[data_nonzero].shape)
                this_ix_sim_dnaframe = np.searchsorted(d.x,np.round((pos[data_nonzero]+dirflag*deltas_sim)/d.binsize) + this_chr_dna_offset)
                
                if wiggle_re>0:
                    this_ix_sim_dnaframe += randint(-wiggle_re, wiggle_re+1, ntosim)
                
                
                ixok = (this_ix_sim_dnaframe>=chr_ix_low_dnaframe) & (this_ix_sim_dnaframe<chr_ix_high_dnaframe)
                thisn = ixok.sum()
      
                
                if not bkg is None:
                    accept_bkg = random_sample(thisn)<B[this_ix_sim_dnaframe[ixok]]
                    thisn=accept_bkg.sum()
                    ix_sim_dnaframe[ncum:(ncum+thisn)] = this_ix_sim_dnaframe[ixok][accept_bkg]
                    
                    data[indx[data_nonzero][ixok][accept_bkg]]-=1
                else:
                #data_nonzero[.]
                    data[indx[data_nonzero][ixok]]-=1
                    ix_sim_dnaframe[ncum:(ncum+thisn)] = this_ix_sim_dnaframe[ixok]
                
                
                ncum+=thisn
                ntotal = data.sum()
                data_nonzero = data>0
                ntosim = data_nonzero.sum()
                    
                
            ix_sim_dnaframe = np.sort(ix_sim_dnaframe)
            outindices[outindptr[i_annot]:outindptr[i_annot+1]] = ix_sim_dnaframe

        out=d.copy(copydata=False)

        out.annot_dict={k:v for k,v in self.annot_dict.items()} # name to ix, chrix  #dictionnary
        out.annot_LUT=self.annot_LUT.copy()
        out.annot_chr=self.annot_chr.copy()
        out.nannot=self.nannot

        out.counts=sparse.csr_matrix((outdata, outindices, outindptr), shape=(self.counts.shape[0],d.counts.shape[1]))
        out.x=d.x.copy()

        return out

    def predict_cis_sparse_debunch(self, d, cismodel, wiggle_re=0, gene_strand_dict=None):
        chr_blocks_dnaframe=d.chr_blocks

        p_all = np.zeros(len(d.x))

        for i in range(len(chr_blocks_dnaframe)-1):
            this_x = d.x[chr_blocks_dnaframe[i]:chr_blocks_dnaframe[i+1]]
            deltas = np.diff(this_x)
            p = np.hstack([0,deltas])+np.hstack([deltas,0])

            # p_windowed = np.zeros(p.shape)
            
            # for j in np.arange(wiggle_re, len(p)-wiggle_re):
            #     p_windowed[j] = np.min(p[(j-wiggle_re):(j+wiggle_re)])
           
            # p_windowed[0:wiggle_re] = np.min(p[0:wiggle_re])
            # p_windowed[len(p)-wiggle_re:len(p)] = np.min( p[len(p)-wiggle_re:len(p)])

            p_windowed=pd.Series(p).rolling(wiggle_re*2, center=True, min_periods=1).min().values
            p = p_windowed/p
            #p = p/p.sum()

            p_all[chr_blocks_dnaframe[i]:chr_blocks_dnaframe[i+1]] = p

        self.tocsr()
        counts=self.counts
        chr_len_vec=self.chr_len_vec
        chr_offset_vec=np.hstack([0,chr_len_vec.cumsum()])[0:-1]
        

        bins=cismodel['bins'] #-1 ?
        
        directional = False if gene_strand_dict is None else True
        dirflag=1

        chr_blocks_dnaframe=d.chr_blocks
        chr_offset_dnaframe=np.hstack([0,np.cumsum(np.ceil(d.chr_len_vec/d.binsize).astype(np.int64))])


        outdata = np.ones(counts.sum())
        outindptr = np.hstack([0,np.array(counts.sum(axis=1)).ravel().cumsum()])
        outindices = np.zeros(outdata.shape[0], np.int64)

        for i_annot in tqdm(range(self.nannot)):

            if directional and (gene_strand_dict.get(self.annot_LUT[i_annot], False)): #TRUE means negative strand
                dirflag=-1
            else:
                dirflag=1
            
            chrid=self.annot_chr[i_annot]
            chr_ix_low_dnaframe = chr_blocks_dnaframe[chrid]
            chr_ix_high_dnaframe = chr_blocks_dnaframe[chrid+1]
            chi=cismodel['chi'][chrid][0:-1]
            this_chr_rna_offset = chr_offset_vec[chrid]
            this_chr_dna_offset = chr_offset_dnaframe[chrid]
            #this_n = counts.indptr[i_annot+1]-counts.indptr[i_annot]
            # if this_n>0:

            data=counts.data[counts.indptr[i_annot]:counts.indptr[i_annot+1]].copy()
            indx=np.arange(len(data))
            pos=self.x[counts.indices[counts.indptr[i_annot]:counts.indptr[i_annot+1]]]-this_chr_rna_offset
            
            ntotal=data.sum()
            ix_sim_dnaframe = np.zeros(ntotal, np.int64)
            data_nonzero = data>0
            ntosim=data_nonzero.sum()
            ncum=0
            #outindptr[]
            while ntosim>0:
                #print(ntotal)
                #print(ntosim)
                deltas_sim = np.interp(random_sample(ntosim), chi, bins)
                #print(pos[data_nonzero].shape)
                this_ix_sim_dnaframe = np.searchsorted(d.x,np.round((pos[data_nonzero]+dirflag*deltas_sim)/d.binsize) + this_chr_dna_offset)
                
                
                ixok = (this_ix_sim_dnaframe>=chr_ix_low_dnaframe) & (this_ix_sim_dnaframe<chr_ix_high_dnaframe)
                thisn = ixok.sum()
      
                
                
                accept_bkg = random_sample(thisn)<p_all[this_ix_sim_dnaframe[ixok]]
                thisn=accept_bkg.sum()
                ix_sim_dnaframe[ncum:(ncum+thisn)] = this_ix_sim_dnaframe[ixok][accept_bkg]
                
                data[indx[data_nonzero][ixok][accept_bkg]]-=1
                
                
                
                ncum+=thisn
                ntotal = data.sum()
                data_nonzero = data>0
                ntosim = data_nonzero.sum()
                    
                
            ix_sim_dnaframe = np.sort(ix_sim_dnaframe)
            outindices[outindptr[i_annot]:outindptr[i_annot+1]] = ix_sim_dnaframe

        out=d.copy(copydata=False)

        out.annot_dict={k:v for k,v in self.annot_dict.items()} # name to ix, chrix  #dictionnary
        out.annot_LUT=self.annot_LUT.copy()
        out.annot_chr=self.annot_chr.copy()
        out.nannot=self.nannot

        out.counts=sparse.csr_matrix((outdata, outindices, outindptr), shape=(self.counts.shape[0],d.counts.shape[1]))

        out.counts.sum_duplicates()
        out.x=d.x.copy()

        return out

    def predict_cis_sparse_debunch2(self, d, cismodel, bkg=None, wiggle_re=0, gene_strand_dict=None, b_prenormalized=False, isproba=False):
        
        use_bkg = False
        if not bkg is None:
            use_bkg = True
            if not b_prenormalized:
                bkg_mat=bkg.normalize_row_bychr(expand=False).counts.toarray().ravel()
            else:
                bkg_mat=bkg.counts.toarray().ravel()
        
        self.tocsr()
        counts=self.counts
        chr_len_vec=self.chr_len_vec
        chr_offset_vec=np.hstack([0,chr_len_vec.cumsum()])[0:-1]

        bins=cismodel['bins'][np.newaxis, :] #-1 ?
        
        directional = False if gene_strand_dict is None else True
        dirflag=1

        chr_blocks_dnaframe=d.chr_blocks
        chr_offset_dnaframe=np.hstack([0,np.cumsum(np.ceil(d.chr_len_vec/d.binsize).astype(np.int64))])


        outdata = np.ones(counts.sum())
        n_rnas=np.array(counts.sum(axis=1)).ravel()
        outindptr = np.hstack([0,n_rnas.cumsum()])
        outindices = np.zeros(outdata.shape[0], np.int64)


        this_chr = -1

        for i_annot in tqdm(range(self.nannot)):
            if (n_rnas[i_annot]>0):
                if directional and (gene_strand_dict.get(self.annot_LUT[i_annot], False)): #TRUE means negative strand
                    dirflag=-1
                else:
                    dirflag=1
                
                chrid=self.annot_chr[i_annot]
                

                if this_chr!=chrid:
                    chi=cismodel['chi'][chrid][0:-1]
                    this_chr_rna_offset = chr_offset_vec[chrid]
                    this_chr_dna_offset = chr_offset_dnaframe[chrid]
                    x_dna = d.x[np.newaxis,chr_blocks_dnaframe[chrid]:chr_blocks_dnaframe[chrid+1]]-this_chr_dna_offset
                    n_x = x_dna.shape[1]
                    this_chr = chrid
                    if use_bkg:
                        this_bkg = bkg_mat[chr_blocks_dnaframe[chrid]:chr_blocks_dnaframe[chrid+1]]
                    
                sim_x_pos = np.zeros(n_x, int)
            
                rna_pos=self.x[counts.indices[counts.indptr[i_annot]:counts.indptr[i_annot+1]]]-this_chr_rna_offset

                rna_pos_approx = (np.floor(rna_pos/100)*100).astype(np.int64)
                
                rna_frac_at_pos = counts.data[counts.indptr[i_annot]:counts.indptr[i_annot+1]]
                rna_frac_at_pos = rna_frac_at_pos/rna_frac_at_pos.sum()

                mean_rna_pos = np.sum(rna_pos * rna_frac_at_pos)
                
                #rna_pos=rna_pos[:,np.newaxis]
                n_rna_pos = rna_pos.shape[0]
                

                probas = np.interp(x_dna, mymat.indices+bmin, mymat.data)
                #probas = np.interp(x_dna, mymat.indices+bmin, mymat.data)

                    #ntosim = counts.data[counts.indptr[i_annot]+j]

                #delta = (x_dna-(this_chr_dna_offset+rna_pos))*dirflag #check x_dna binning
                # all_chi = (np.tile(chi,(n_rna_pos,1)) * rna_frac_at_pos).ravel()
                # #if dirflag>0:
                # all_bins = (bins+rna_pos).ravel()
                
                # bmin=np.min(all_bins)
                # bmax=np.max(all_bins)

                # matsummer = sparse.csr_matrix((np.ones(n_rna_pos), np.arange(n_rna_pos), np.array([0, n_rna_pos], dtype=int)), shape=(1, n_rna_pos))

                # mymat = matsummer*sparse.csr_matrix((all_chi, all_bins-bmin,np.hstack([0, np.cumsum(bins.shape[1]*np.ones(n_rna_pos))])), shape=(n_rna_pos, bmax-bmin+1))
                # mymat.sort_indices()
                # mymat.sum_duplicates()
                probas = np.interp(x_dna, mymat.indices+bmin, mymat.data)
                probas = np.diff(np.hstack([0,probas]), axis=1)

        

                    # for j in range(n_rna_pos):
                    #     probas = np.interp(x_dna, bins+rna_pos[j], chi)
                #else:
                    # for j in range(n_rna_pos):
                    #     probas = np.interp(-x_dna, bins-rna_pos[j], chi)
                    # all_bins = (bins-rna_pos).ravel()
                    # probas = np.interp(-x_dna, all_bins, all_chi)
                
                
        

        #         probas = np.reshape(probas, (n_rna_pos, n_x))
        #         probas = np.diff(np.hstack([np.zeros((n_rna_pos,1)),probas]), axis=1)
                

        #         probas = probas/ probas.sum(axis=1)

        #         probas = (probas * rna_frac_at_pos).sum(axis=0)

                probas = probas / probas.sum()

                sim_x_pos = multinomial(n_rnas[i_annot], probas)
                    # if use_bkg:
                    #     probas = probas*this_bkg
                    # probas = probas / probas.sum()
                    
                    
                    # sim_x_pos += multinomial(ntosim, probas)
                    
                this_data=np.ones(n_rnas[i_annot], int)
                this_rowindices=np.arange(n_rnas[i_annot])
                this_colindptr = np.cumsum(np.hstack([0, sim_x_pos]))

                this_colindices =  sparse.csc_matrix((this_data, this_rowindices, this_colindptr), dtype=int, shape=(n_rnas[i_annot], len(x_dna))).tocsr().indices + chr_blocks_dnaframe[chrid]
                np.sort(this_colindices)
                #outindices[ncum:(ncum+ntotal)]=this_colindices
                outindices[outindptr[i_annot]:outindptr[i_annot+1]]=this_colindices

        out=d.copy(copydata=False)

        # out.annot_dict={k:v for k,v in self.annot_dict.items()} # name to ix, chrix  #dictionnary
        # out.annot_LUT=self.annot_LUT.copy()
        # out.annot_chr=self.annot_chr.copy()
        # out.nannot=self.nannot

        # out.counts=sparse.csr_matrix((outdata, outindices, outindptr), shape=(self.counts.shape[0],d.counts.shape[1]))
        # #out.counts.sort_indices()
        # out.counts.sum_duplicates()
        # out.x=d.x.copy()

        return out

    # def predict_cis_sparse_debunch2(self, d, cismodel, bkg=None, wiggle_re=0, gene_strand_dict=None, b_prenormalized=False, isproba=False):
        
    #     use_bkg = False
    #     if not bkg is None:
    #         use_bkg = True
    #         if not b_prenormalized:
    #             bkg_mat=bkg.normalize_row_bychr(expand=False).counts.toarray().ravel()
    #         else:
    #             bkg_mat=bkg.counts.toarray().ravel()
        
    #     self.tocsr()
    #     counts=self.counts
    #     chr_len_vec=self.chr_len_vec
    #     chr_offset_vec=np.hstack([0,chr_len_vec.cumsum()])[0:-1]

    #     bins=cismodel['bins'] #-1 ?
        
    #     directional = False if gene_strand_dict is None else True
    #     dirflag=1

    #     chr_blocks_dnaframe=d.chr_blocks
    #     chr_offset_dnaframe=np.hstack([0,np.cumsum(np.ceil(d.chr_len_vec/d.binsize).astype(np.int64))])


    #     outdata = np.ones(counts.sum())
    #     n_rnas=np.array(counts.sum(axis=1)).ravel()
    #     outindptr = np.hstack([0,n_rnas.cumsum()])
    #     outindices = np.zeros(outdata.shape[0], np.int64)


    #     this_chr = -1

    #     for i_annot in tqdm(range(self.nannot)):
    #         if (n_rnas[i_annot]>0):
    #             if directional and (gene_strand_dict.get(self.annot_LUT[i_annot], False)): #TRUE means negative strand
    #                 dirflag=-1
    #             else:
    #                 dirflag=1
                
    #             chrid=self.annot_chr[i_annot]
                

    #             if this_chr!=chrid:
    #                 chi=cismodel['chi'][chrid][0:-1]
    #                 this_chr_rna_offset = chr_offset_vec[chrid]
    #                 this_chr_dna_offset = chr_offset_dnaframe[chrid]
    #                 x_dna = d.x[chr_blocks_dnaframe[chrid]:chr_blocks_dnaframe[chrid+1]]
    #                 n_x = len(x_dna)
    #                 this_chr = chrid
    #                 if use_bkg:
    #                     this_bkg = bkg_mat[chr_blocks_dnaframe[chrid]:chr_blocks_dnaframe[chrid+1]]
                    
    #             #sim_x_pos = np.zeros(n_x, int)
            
    #             rna_pos=self.x[counts.indices[counts.indptr[i_annot]:counts.indptr[i_annot+1]]]-this_chr_rna_offset


    #             rna_frac_at_pos = counts.data[counts.indptr[i_annot]:counts.indptr[i_annot+1]]
    #             rna_frac_at_pos = rna_frac_at_pos/rna_frac_at_pos.sum()
                
    #             n_rna_pos = rna_pos.shape[0]

    #             probas_cum = np.zeros(n_x)

    #             for j in range(n_rna_pos):

    #                 delta_j = (x_dna-(this_chr_dna_offset+rna_pos[j]))*dirflag #check x_dna binning
    #                 if dirflag<0:
    #                     delta_j = np.flip(delta_j)
                    
    #                     probas = np.diff(np.hstack([0,np.interp(delta_j, bins, chi)]))
    #                     probas = probas
    #                     probas = np.flip(probas)
    #                 else:
    #                     probas = np.diff(np.hstack([0,np.interp(delta_j, bins, chi)]))
    #                     probas = probas

    #                 probas = probas / probas.sum()
                
    #             probas_cum += probas*rna_frac_at_pos[j]
    #             if use_bkg:
    #                 probas = probas*this_bkg
        
    #             probas = probas / probas.sum()
                    
                    
    #             sim_x_pos = multinomial(n_rnas[i_annot], probas)
                    
    #             this_data=np.ones(n_rnas[i_annot], int)
    #             this_rowindices=np.arange(n_rnas[i_annot])
    #             this_colindptr = np.cumsum(np.hstack([0, sim_x_pos]))

    #             this_colindices =  sparse.csc_matrix((this_data, this_rowindices, this_colindptr), dtype=int, shape=(n_rnas[i_annot], len(x_dna))).tocsr().indices + chr_blocks_dnaframe[chrid]
    #             np.sort(this_colindices)
    #             #outindices[ncum:(ncum+ntotal)]=this_colindices
    #             outindices[outindptr[i_annot]:outindptr[i_annot+1]]=this_colindices

    #     out=d.copy(copydata=False)

    #     out.annot_dict={k:v for k,v in self.annot_dict.items()} # name to ix, chrix  #dictionnary
    #     out.annot_LUT=self.annot_LUT.copy()
    #     out.annot_chr=self.annot_chr.copy()
    #     out.nannot=self.nannot

    #     out.counts=sparse.csr_matrix((outdata, outindices, outindptr), shape=(self.counts.shape[0],d.counts.shape[1]))
    #     #out.counts.sort_indices()
    #     out.counts.sum_duplicates()
    #     out.x=d.x.copy()

    #     return out


    # def predict_trans_withinGI(self, bkg, gi, shelfing_matrix=None, leakage_vector=None, return_chartable=True):
    #     #make sure to pass expected the trans matrix as self!
    #     N=self.get_chr_matrix(withcounts=True)
    #     if not (isinstance(bkg, sparse.csc.csc_matrix) or isinstance(bkg, sparse.csr.csr_matrix)): 
    #         B=bkg.normalize_row_bytrans().counts.tocsc() #nchr X p matrix block diagonal
    #     else:
    #         B=bkg
    #     if shelfing_matrix is None:
    #         shelfing_matrix=gi['counting_matrix']
    #     if not leakage_vector is None:
    #         L=sparse.diags(leakage_vector, format='csr')
    #     else: #we use existing trans counts
    #         L=1

    #     M=(B*shelfing_matrix).tocsr()
    #     out_counts=(N*L)*M

    #     if return_chartable:

    #         out=Chartable()
    #         out.nannot=self.nannot
    #         out.annot_chr=self.annot_chr.copy()
    #         out.annot_dict = {k:v for k, v in self.annot_dict.items()}
    #         out.annot_LUT=[k for k in self.annot_LUT]

    #         out.chr_blocks=gi['chr_blocks']
    #         out.chr_dict=self.chr_dict.copy()
    #         out.chr_LUT=[k for k in self.chr_LUT]
    #         out.chr_len_vec=self.chr_len_vec.copy()
    #         out.nchr=self.nchr

    #         out.counts=out_counts
    #         out.nbins=out.counts.shape[1]
    #         out.bins_LUT=gi['feature_name']

    #         return out
    #     else:
    #         return out_counts

    @staticmethod
    def inv_logit(x):
        ex = np.exp(x)
        return ex/(1.0+ex)

    # def predict_trans_leakage(self, betas, shrink_estimates=True):
        
    #     Nct = np.array(self.counts.sum(axis=1)).ravel()[:,np.newaxis]
    #     _, trans, _ = self.cis_trans()
    #     Nbychr = trans.N_bychr()

    #     L_offset = betas.loc[self.chr_LUT][self.chr_LUT].values
    #     L_intercept = betas['(Intercept)'].loc[self.chr_LUT].values
    #     L_alpha = betas['log(Nct)'].loc[self.chr_LUT].values
    #     sigmas = betas['sigma'].loc[self.chr_LUT].values
        
    #     y_offset = L_offset + L_intercept[:, np.newaxis]
    #     y_alpha = np.log(Nct)*L_alpha[self.annot_chr][:, np.newaxis]
        
    #     r_prior = Chartable_binned.inv_logit(y_alpha + y_offset[self.annot_chr])
        
    #     if shrink_estimates:
    #         s = sigmas[self.annot_chr][:, np.newaxis]
    #         alpha_prior=r_prior / s
    #         beta_prior=(1-r_prior)/s
    #         alpha_post=alpha_prior + Nbychr.toarray()
    #         beta_post=beta_prior + Nct
            
    #         r_post=alpha_post/(alpha_post+beta_post)
        
    #         N_pred = Nct*r_post
    #     else:
    #         N_pred = Nct*r_prior
        
    #     return sparse.csr_matrix(N_pred)

    def predict_trans_leakage_bychr(self, betas, shrink_estimates=True, fix_cis=True):
        
        Nct = np.array(self.counts.sum(axis=1)).ravel()[:,np.newaxis]
        _, _, N = self.cis_trans()
        Ncis = N['Ncis'].values.ravel()[:,np.newaxis]
        Nbychr = self.N_bychr().toarray()

        L_offset = betas.loc[self.chr_LUT][self.chr_LUT].values
        L_intercept = betas['(Intercept)'].loc[self.chr_LUT].values
        L_alpha = betas['log(Nct)'].loc[self.chr_LUT].values
        sigmas = betas['sigma'].loc[self.chr_LUT].values
        
        y_offset = L_offset + L_intercept[:, np.newaxis]
        y_alpha = np.log(Nct)*L_alpha[self.annot_chr][:, np.newaxis]
        
        r_prior = Chartable_binned.inv_logit(y_alpha + y_offset[self.annot_chr])
        
        if shrink_estimates:
            s = sigmas[self.annot_chr][:, np.newaxis]
            alpha_prior=r_prior / s
            beta_prior=(1-r_prior)/s
            alpha_post=alpha_prior + Nbychr
            #beta_post=beta_prior + (Nct-Nbychr)
            beta_post=beta_prior + Ncis
            
            r_post=alpha_post/(alpha_post+beta_post)
        
            if fix_cis:
                N_pred = Ncis*(r_post/(1-r_post)) #N/(Ncis+N)=r, N(1-r)=Ncis r,N=r/(1-r)Ncis
            else:
                N_pred = (Nct)*r_post # this seems wrong
        else:
            if fix_cis:
                N_pred = Ncis*(r_prior/(1-r_prior))
            else:
                N_pred = (Nct)*r_prior
        
        return sparse.csr_matrix(N_pred)


    def predict_trans_leakage_bytrans(self, betas, shrink_estimates=True, fix_cis=True):
        
        #Nct = np.array(self.counts.sum(axis=1)).ravel()[:,np.newaxis]
        _, _, N = self.cis_trans()
        Ncis = N['Ncis'].values.ravel()[:,np.newaxis]
        Ntrans = N['Ntrans'].values.ravel()[:,np.newaxis]
        #Nbychr = self.N_bychr().toarray()

        L_offset = betas.loc[self.chr_LUT][self.chr_LUT].values
        L_intercept = betas['(Intercept)'].loc[self.chr_LUT].values
        L_alpha = betas['log(Nct)'].loc[self.chr_LUT].values
        sigmas = betas['sigma'].loc[self.chr_LUT].values
        
        y_offset = L_offset + L_intercept[:, np.newaxis]
        y_alpha = np.log(Ncis+Ntrans)*L_alpha[self.annot_chr][:, np.newaxis]
        
        r_prior = Chartable_binned.inv_logit(y_alpha + y_offset[self.annot_chr])
        
        if shrink_estimates:
            s = sigmas[self.annot_chr][:, np.newaxis]
            alpha_prior=r_prior / s
            beta_prior=(1-r_prior)/s
            alpha_post=alpha_prior + Ntrans
            beta_post=beta_prior + Ncis
            
            r_post=alpha_post/(alpha_post+beta_post)
            if fix_cis:
                N_pred = Ncis*(r_post/(1-r_post))
            else:
                N_pred = (Ncis+Ntrans)*r_post
        else:
            if fix_cis:
                N_pred = Ncis*(r_prior/(1-r_prior))
            else:
                N_pred = (Ncis+Ntrans)*r_prior
        
        return sparse.csr_matrix(N_pred)

    def predict_trans_withinGI_bychr(self, bkg, pr, shelfing_matrix=None, betas=None, shrink_estimates=True, fix_cis=True, return_chartable=True, b_prenormalized=False, include_cis = False):
        #make sure to pass expected the trans matrix as self! --> no longer true
        
        if betas is None:
            if include_cis:
                N = self.N_bychr()
            else:
                _, trans, _ = self.cis_trans()
                N=trans.N_bychr()

        else:
            N = self.predict_trans_leakage_bychr(betas, shrink_estimates=shrink_estimates, fix_cis=fix_cis)

        if not b_prenormalized:
            B=bkg.normalize_row_bychr().counts
        else:
            B=bkg.counts

        if shelfing_matrix is None:
            # assert gi.binsize==bkg.binsize, "Mismatched bin size in background and GI object!" 
            #shelfing_matrix=bkg.gi_to_shelfmatrix(gi)
            shelfing_matrix=pr['counting_matrix']

        M=B*shelfing_matrix
        out_counts=(N*M).tocsr()

        if return_chartable:

            out=Chartable()
            out.nannot=self.nannot
            out.annot_chr=self.annot_chr.copy()
            out.annot_dict = {k:v for k, v in self.annot_dict.items()}
            out.annot_LUT=[k for k in self.annot_LUT]

            out.chr_blocks=pr['chr_blocks']
            out.chr_dict=self.chr_dict.copy()
            out.chr_LUT=[k for k in self.chr_LUT]
            out.chr_len_vec=self.chr_len_vec.copy()
            out.nchr=self.nchr

            out.counts=out_counts
            out.nbins=out.counts.shape[1]
            out.bins_LUT=pr['feature_name']

            return out
        else:
            return out_counts

    # cis: model, no model / + or - bkg
    # trans : true N, predicted N | +/- bkg

    def predict_nomodel_sparse(self, bkg, genes, genes_to_keep=[], betas=None, shrink_estimates=False, fix_cis=True,replace_cis_leak=False):
        if len(genes_to_keep)>0:
            this_dna = self.select_annots2(genes_to_keep, genes)
        else:
            this_dna = self

        Nbychr = this_dna.N_bychr()
        pred_nomodel = this_dna.predict_bychr_bkgonly(bkg=bkg, b_prenormalized=False, Nbychr=Nbychr)
        
       
        if betas is None:
            pred_nomodel_leak=None
        else:

            N = this_dna.predict_trans_leakage_bychr(betas, shrink_estimates=shrink_estimates, fix_cis=fix_cis)
            N.data=np.round(N.data).astype(int)

            pred_nomodel_leak = this_dna.predict_bychr_bkgonly(bkg=bkg, b_prenormalized=False, Nbychr=N)
            if replace_cis_leak:
                pred_nomodel_leak = pred_nomodel_leak.replace_cis(pred_nomodel)
        
        
        return pred_nomodel, pred_nomodel_leak

    def predict_cis_sparse2(self, bkg, genes=None, genes_to_keep=[], cismodel=None, gene_strand_dict=None):
        if len(genes_to_keep)>0:
            this_rna = self.select_annots2(genes_to_keep, genes)
        else:
            this_rna = self

        # no bck fluctiations
        pred_cis = this_rna.predict_cis_sparse(bkg, cismodel, gene_strand_dict=gene_strand_dict, bkg=None, wiggle_re=0)

        # with bck fluctuatations
        pred_cis_flat = pred_cis.sparsify_bychr(bkg)

        return pred_cis, pred_cis_flat
        
    


    # def predict_sparse(self, bkg, genes, rna=None,  predict_cis = True, predict_trans = True, cismodel=None, gene_strand_dict=None, betas=None, shrink_estimates=True, fix_cis=True, genes_to_keep=[]):
       
    #     assert predict_cis or predict_trans

    #     if len(genes_to_keep)>0:
    #             this_dna = self.select_annots2(genes_to_keep, genes)
    #     else:
    #             this_dna = self

    #     if not rna is None and len(genes_to_keep)==0:
    #         assert self.annot_LUT==rna.annot_LUT
    #         this_rna = rna

    #     elif not rna is None:
    #         this_rna = rna.select_annots2(genes_to_keep, genes)
    #     else:
    #         this_rna = rna

        

    #     # start prediction
    #     if predict_cis:
    #         if not this_rna is None: # we use the model
    #             pred_cis = this_rna.predict_cis_sparse(this_dna, cismodel, gene_strand_dict=gene_strand_dict)
    #             if not bkg is None:
    #                 pred_cis = sparsify_bychr(bkg=bkg, b_prenormalized=False)
    #         else: # we 
    #             pass
        
    #     if predict_trans:
    #         this_dna.sparsify_bychr(bkg=bkg, b_prenormalized=False)

    #         if not pr_rna is None: # use the model for predition
    #             cts_cis = this_rna.predict_cis_withinGI(cismodel, pr_rna, gene_strand_dict=gene_strand_dict, return_chartable=True)
                
    #         else: 
    #             if not this_rna is None: #if rna is there use rna itself
    #                 cts_cis, _, _ = this_rna.cis_trans()
    #             else: # otherwise predict without decay model
    #                 if len(genes_to_keep)>0:
    #                     dna_cis, _, _ = self.select_annots(genes_to_keep).cis_trans()
    #                 else:
    #                     dna_cis, _, _ = self.cis_trans()
    #                 cts_cis = dna_cis.predict_trans_withinGI_bychr(bkg, pr, betas=None, include_cis=True)

    #         if len(genes_to_keep)>0:
    #             cts_cis = cts_cis.select_annots2(genes_to_keep, genes)
            
    #         if flatfield and (not cismodel is None): #if cismodel is none flatfielding is automatic
    #             cts_cis, flatv = cts_cis.flat_field(bkg, pr, b_prenormalized=False)
        
    #     #_, dna_trans, _ = self.cis_trans()
        
    #     if predict_trans:
    #         if len(genes_to_keep)>0:
    #             this_dna = self.select_annots(genes_to_keep)
    #         else:
    #             this_dna = self
        
    #         cts_trans = this_dna.predict_trans_withinGI_bychr(bkg, pr, betas=betas, shrink_estimates=shrink_estimates, fix_cis=fix_cis, return_chartable=True, b_prenormalized=b_prenormalized, include_cis=False)
    #         _, cts_trans, _ = cts_trans.cis_trans()

    #         if len(genes_to_keep)>0:
    #             cts_trans = cts_trans.select_annots2(genes_to_keep, genes)
            
    #     if predict_cis and predict_trans:
    #         cts = cts_trans
    #         cts.counts = cts.counts + cts_cis.counts
    #         return cts, flatv
    #     elif predict_cis:
    #         return cts_cis, flatv
    #     else:
    #         return cts_trans, flatv

    def predict(self, bkg, pr, rna, genes, predict_cis = True, predict_trans = True, pr_rna=None, b_prenormalized=False, cismodel=None, gene_strand_dict=None, betas=None, shrink_estimates=True, fix_cis=True,flatfield=True, genes_to_keep=[]):
        flatv = None
        assert predict_cis or predict_trans

        if predict_cis and predict_trans and len(genes_to_keep)==0:
            if not rna is None:
                assert self.annot_LUT==rna.annot_LUT

        if predict_cis:
            #
            if (not rna is None) and (len(genes_to_keep)>0):
                this_rna = rna.select_annots(genes_to_keep)
            else:
                this_rna = rna

            if not pr_rna is None: # use the model for predition
                cts_cis = this_rna.predict_cis_withinGI(cismodel, pr_rna, gene_strand_dict=gene_strand_dict, return_chartable=True)
                
            else: 
                if not this_rna is None: #if rna is there use rna itself
                    cts_cis, _, _ = this_rna.cis_trans()
                else: # otherwise predict without decay model
                    if len(genes_to_keep)>0:
                        dna_cis, _, _ = self.select_annots(genes_to_keep).cis_trans()
                    else:
                        dna_cis, _, _ = self.cis_trans()
                    cts_cis = dna_cis.predict_trans_withinGI_bychr(bkg, pr, betas=None, include_cis=True)

            if len(genes_to_keep)>0:
                cts_cis = cts_cis.select_annots2(genes_to_keep, genes)
            
            if flatfield and (not cismodel is None): #if cismodel is none flatfielding is automatic
                cts_cis, flatv = cts_cis.flat_field(bkg, pr, b_prenormalized=False)
        
        #_, dna_trans, _ = self.cis_trans()
        
        if predict_trans:
            if len(genes_to_keep)>0:
                this_dna = self.select_annots(genes_to_keep)
            else:
                this_dna = self
        
            cts_trans = this_dna.predict_trans_withinGI_bychr(bkg, pr, betas=betas, shrink_estimates=shrink_estimates, fix_cis=fix_cis, return_chartable=True, b_prenormalized=b_prenormalized, include_cis=False)
            _, cts_trans, _ = cts_trans.cis_trans()

            if len(genes_to_keep)>0:
                cts_trans = cts_trans.select_annots2(genes_to_keep, genes)
            
        if predict_cis and predict_trans:
            cts = cts_trans
            cts.counts = cts.counts + cts_cis.counts
            return cts, flatv
        elif predict_cis:
            return cts_cis, flatv
        else:
            return cts_trans, flatv

    def replace_cis(self, new_cis, annots=None, genes=None, how='outter',  reorder_genes=True):
        _, out, _ = self.cis_trans()
        new_cis_cis, _, _ = new_cis.cis_trans()

        if (annots is None) and (out.annot_LUT==new_cis_cis.annot_LUT):
            annots_new = None
        else:
            if annots is None:
                if self.annot_LUT == new_cis_cis.annot_LUT:
                    annots_new = None
                else:
                    if how=='inner':
                        annots_new = list((set(self.annots).intersect(set(new_cis_cis.annots))))
                    else:
                        annots_new = list((set(self.annots).union(set(new_cis_cis.annots))))
                    annots_new = genes.loc[genes.index.isin(annots_new)].index.to_list()
            else:
                if (out.annot_LUT==annots) and (new_cis_cis.annot_LUT==annots):
                    annots_new = None
                else:
                    annots_new = annots

        if reorder_genes and (not annots_new is None):
            annots_new = genes.loc[genes.index.isin(annots_new)].index.to_list()

        if not annots_new is None:
            out = out.select_annots2(annots_new, genes=genes)
            new_cis_cis = new_cis_cis.select_annots2(annots_new, genes=genes)

        out.counts = out.counts + new_cis_cis.counts

        return out

    def replace_trans(self, new_trans, annots=None, genes=None, how='outter', reorder_genes=True):
        out, _, _ = self.cis_trans()
        _, new_trans_trans, _ = new_trans.cis_trans()

        if (annots is None) and (out.annot_LUT==new_trans_trans.annot_LUT):
            annots_new = None
        else:
            if annots is None:
                if self.annot_LUT == new_trans_trans.annot_LUT:
                    annots_new = None
                else:
                    if how=='inner':
                        annots_new = list((set(self.annots).intersect(set(new_trans_trans.annots))))
                    else:
                        annots_new = list((set(self.annots).union(set(new_trans_trans.annots))))
                    
            else:
                if (out.annot_LUT==annots) and (new_trans_trans.annot_LUT==annots):
                    annots_new = None
                else:
                    annots_new = annots

        if reorder_genes and (not annots_new is None):
            annots_new = genes.loc[genes.index.isin(annots_new)].index.to_list()

        
        if not annots_new is None:
            out = out.select_annots2(annots_new, genes=genes)
            new_trans_trans = new_trans_trans.select_annots2(annots_new, genes=genes)

        out.counts = out.counts + new_trans_trans.counts

        return out


    # @staticmethod
    # def get_smoothening_intervals_(x, y, N_thr, w_thr):
   
    #     n=len(y)

    #     ycum=np.cumsum(y)
    #     # xcum=np.cumsum(x)

    #     iR=np.minimum(n-1,np.searchsorted(ycum, ycum-y+N_thr))

    #     nR=ycum[iR]-ycum+y
    #     xR=x[iR]

    #     iL=np.searchsorted(ycum, ycum-y-N_thr)
    #     nL=ycum-ycum[iL]+y[iL]
    #     xL=x[iL]
        
    #     ysmooth=nL+nR-y
    #     flagout=(nL<N_thr) | (nR<N_thr) | ((xR-x)>w_thr) | ((x-xL)>w_thr)
        
    #     #ysmooth_norm=

    #     ysmooth[flagout]=0
    #     #dsmooth=ysmooth/(xR-xL) 
    #     fout=(n-np.sum(flagout))/n
    #     return xL, xR, ysmooth, fout, ysmooth_norm #dsmooth

    @staticmethod
    def get_smoothening_intervals_(x, y, N_thr, w_thr):
   
        n=len(y)

        ycum=np.cumsum(y)
        # xcum=np.cumsum(x)

        iR=np.minimum(n-1,np.searchsorted(ycum-y, ycum-y+N_thr)) #added -y to not include iR in the bound
        #iR such that sum([i,iR[)>=Nthr

        nR=ycum[iR]-ycum+y -y[iR] # [i,iR[ added -y[iR] so as not to include the right bound
        xR=x[iR]

        iL=np.searchsorted(ycum, ycum-y-N_thr+1) #added -1, 
        #iL such that sum([iL, i[)>=Nthr

        nL=ycum-ycum[iL]+y[iL]-y #[iL, i[ #added -y to not include right boundary
        xL=x[iL]
        
        #ysmooth=nL+nR-y 
        ysmooth=nL+nR #[iL, iR[
        flagout=(nL<N_thr) | (nR<N_thr) | ((xR-x)>w_thr) | ((x-xL)>w_thr)
        
        ysmooth[flagout]=0
        #dsmooth=ysmooth/(xR-xL) 

        ysmooth_norm=ysmooth/(iR-iL)
        ysmooth_norm[flagout]=0

        fout=(n-np.sum(flagout))/n
        return xL, xR, ysmooth, fout, ysmooth_norm


    def get_smoothening_intervals_2(x, y, N_thr, w_thr):
   
        n=len(y)

        ycum=np.cumsum(y)
        # xcum=np.cumsum(x)

        iR=np.minimum(n-1,np.searchsorted(ycum-y, ycum-y+N_thr)) #added -y to not include iR in the bound
        #iR such that sum([i,iR[)>=Nthr

        nR=ycum[iR]-ycum+y -y[iR] # [i,iR[ added -y[iR] so as not to include the right bound
        xR=x[iR]

        iL=np.searchsorted(ycum, ycum-y-N_thr+1) #added -1, 
        #iL such that sum([iL, i[)>=Nthr

        nL=ycum-ycum[iL]+y[iL]-y #[iL, i[ #added -y to not include right boundary
        xL=x[iL]
        
        #ysmooth=nL+nR-y 
        ysmooth=nL+nR #[iL, iR[
        flagout=(nL<N_thr) | (nR<N_thr) | ((xR-x)>w_thr) | ((x-xL)>w_thr)
        
        ysmooth[flagout]=0
        #dsmooth=ysmooth/(xR-xL) 
        fout=(n-np.sum(flagout))/n
        return xL, xR, ysmooth, fout, iL, iR

    def smoothen_counts(self, N_thr=50, w_thr=250000, return_chartable=True):
        
        #self is DNA data, rna is for where the counts are coming from

        self.tocsr()
        counts=self.counts


        smooth_counts=np.zeros(len(counts.data))
        smooth_density=np.zeros(len(counts.data))
        edgesL=np.zeros(len(counts.data))
        edgesR=np.zeros(len(counts.data))

        chr_blocks=self.chr_blocks

        for i in range(self.nannot):

            yall=counts.data[counts.indptr[i]:counts.indptr[i+1]]
            targetchr=np.searchsorted(chr_blocks,counts.indices[counts.indptr[i]:counts.indptr[i+1]],side='right')-1

            cb=np.searchsorted(targetchr, np.arange(len(chr_blocks)))
            xall=self.x[counts.indices[counts.indptr[i]:counts.indptr[i+1]]]

            p=np.arange(counts.indptr[i],counts.indptr[i+1])

            for j in range(self.nchr):
                if cb[j+1]-cb[j]>2*N_thr:
                    x=xall[cb[j]:cb[j+1]]
                    y=yall[cb[j]:cb[j+1]]

                    xL, xR, ysmooth, _, dsmooth=Chartable_binned.get_smoothening_intervals_(x, y, N_thr, w_thr)
                    smooth_counts[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=ysmooth
                    smooth_density[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=dsmooth*self.binsize
                    edgesL[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=xL
                    edgesR[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=xR

        if return_chartable:
            out_smooth=self.copy()
            out_smooth.counts=sparse.csr_matrix((smooth_counts, counts.indices, counts.indptr), shape=counts.shape)

            out_density=self.copy()
            out_density.counts=sparse.csr_matrix((smooth_density, counts.indices, counts.indptr), shape=counts.shape)

            return out_density, out_smooth, edgesL, edgesR
        else:
            return smooth_density, smooth_counts, edgesL, edgesR

    def predict_cis_smooth(self, cismodel, rna, N_thr=50, w_thr=250000, gene_strand_dict=None, return_chartable=True):
        self.tocsr()
        counts=self.counts


        smooth_counts=np.zeros(len(counts.data))
        smooth_counts_pred=np.zeros(len(counts.data))
        smooth_density=np.zeros(len(counts.data))
        smooth_density_pred=np.zeros(len(counts.data))
        # edgesL=np.zeros(len(counts.data))
        # edgesR=np.zeros(len(counts.data))

        chr_blocks=self.chr_blocks
        chr_len_vec=self.chr_len_vec
        chr_offset_vec=np.hstack([0,chr_len_vec.cumsum()])


        chr_offset_refframe=np.hstack([0,np.cumsum(np.ceil(self.chr_len_vec/self.binsize).astype(np.int64))])
    


        rna_annot_ixs=[rna.annot_dict.get(k,-1) for k in self.annot_LUT]

        chrid=-1

        bins=cismodel['bins'] #-1
        b=bins-1
        chi=cismodel['chi'][0]
        rho=cismodel['rho'][0]
        edges=np.hstack([cismodel['bins'],cismodel['bins'][-1]+1])

        directional = False if gene_strand_dict is None else True
        dirflag=1

        for i in range(self.nannot):
            
            yall=counts.data[counts.indptr[i]:counts.indptr[i+1]]
            targetchr=np.searchsorted(chr_blocks,counts.indices[counts.indptr[i]:counts.indptr[i+1]],side='right')-1
            #targetchr=np.searchsorted(chr_blocks,counts.indices[counts.indptr[i]:counts.indptr[i+1]]) #this does not work

            cb=np.searchsorted(targetchr, np.arange(len(chr_blocks)))
            xall=self.x[counts.indices[counts.indptr[i]:counts.indptr[i+1]]]

            #p=np.arange(counts.indptr[i],counts.indptr[i+1])

            rna_ix=rna_annot_ixs[i]
            
            if rna_ix>-1:
                ix1=rna.counts.indptr[rna_ix]
                ix2=rna.counts.indptr[rna_ix+1]
                rna_counts=rna.counts.data[ix1:ix2]
                rna_pos=rna.x[rna.counts.indices[ix1:ix2]]
                
                chrid=self.annot_chr[i]
                chr_len=chr_len_vec[chrid]

                chi=cismodel['chi'][chrid]
                rho=cismodel['rho'][chrid]

                if directional and (gene_strand_dict.get(self.annot_LUT[i], False)): #TRUE means negative strand
                    dirflag=-1
                else:
                    dirflag=1

                x1=0
                x2=chr_len-1
                delta1=(x1-rna_pos+chr_offset_vec[chrid]) if dirflag==1 else (rna_pos-x2-chr_offset_vec[chrid])
                delta2=(x2-rna_pos+chr_offset_vec[chrid]) if dirflag==1 else (rna_pos-x1-chr_offset_vec[chrid])
                binid2=np.searchsorted(b,delta2)
                binid1=np.searchsorted(b,delta1)

                Ps=chi[binid2]-rho[binid2]*(edges[binid2]-(delta2)) - chi[binid1]+rho[binid1]*(edges[binid1]-delta1) 
                ntot=np.sum(Ps*rna_counts)

                corr_factor=np.sum(rna_counts)/ntot
                print(rna_counts.sum())
                print(corr_factor)

                for j in range(self.nchr):
                    if cb[j+1]-cb[j]>0: #2*N_thr:
                        x=xall[cb[j]:cb[j+1]]
                        y=yall[cb[j]:cb[j+1]]

                        xL, xR, ysmooth, fout=Chartable_binned.get_smoothening_intervals_(x, y, N_thr, w_thr)
                        smooth_counts[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=ysmooth
                        
                        # edgesL[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=xL
                        # edgesR[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=xR
                        n_gi=len(xL)
                        
                        xL=(xL-chr_offset_refframe[j])*self.binsize+chr_offset_vec[j]
                        xR=(xR-chr_offset_refframe[j])*self.binsize+chr_offset_vec[j]

                        smooth_density[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=ysmooth/(xR-xL)


                        if (j==self.annot_chr[i]):
                            print(chr_offset_refframe[j])
                            print(chr_offset_vec[j])
                            
                            #print(n_gi)
                            for i_gi in range(n_gi):
                                if ysmooth[i_gi]>0:
                                    #print('ok')
                                    x1=xL[i_gi]
                                    x2=xR[i_gi]
                                    delta1=(x1-rna_pos) if dirflag==1 else (rna_pos-x2)
                                    delta2=(x2-rna_pos) if dirflag==1 else (rna_pos-x1)
                                    binid2=np.searchsorted(b,delta2)
                                    binid1=np.searchsorted(b,delta1)

                                    Ps=chi[binid2]-rho[binid2]*(edges[binid2]-(delta2)) - chi[binid1]+rho[binid1]*(edges[binid1]-delta1) 
                                    tt=np.sum(Ps*rna_counts)*corr_factor
                                    smooth_counts_pred[cb[j]+counts.indptr[i]+i_gi]=tt

                                    #smooth_density_pred[cb[j]+counts.indptr[i]+i_gi]=tt/(xR[i_gi]-xL[i_gi])
                                
                            # corr_factor2=smooth_counts[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])].sum()/smooth_counts_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])].sum()

                            # print(corr_factor2)
                            # smooth_counts_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=smooth_counts_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]*corr_factor2

                            smooth_counts_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=smooth_counts[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]*fout
                            print(fout)
                            smooth_density_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=smooth_counts_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]/(xR-xL)
            
            if return_chartable:
                out_smooth_density=self.copy()
                out_smooth_density_pred=self.copy()

                out_smooth_counts=self.copy()
                out_smooth_counts_pred=self.copy()


                out_smooth_density.counts=sparse.csr_matrix((smooth_density, counts.indices, counts.indptr), shape=counts.shape)
                out_smooth_density_pred.counts=sparse.csr_matrix((smooth_density_pred, counts.indices, counts.indptr), shape=counts.shape)

                out_smooth_counts.counts=sparse.csr_matrix((smooth_counts, counts.indices, counts.indptr), shape=counts.shape)

                out_smooth_counts_pred.counts=sparse.csr_matrix((smooth_counts_pred, counts.indices, counts.indptr), shape=counts.shape)
        

                return out_smooth_density, out_smooth_density_pred, out_smooth_counts, out_smooth_counts_pred
            else:
                return smooth_density, smooth_density_pred, smooth_counts, smooth_counts_pred

# if not (isinstance(bkg, sparse.csc.csc_matrix) or isinstance(bkg, sparse.csr.csr_matrix)): 
#             B=bkg.normalize_row_bytrans().counts.tocsc()

    def predict_smooth(self, cismodel, rna, bkg, N_thr=50, w_thr=250000, gene_strand_dict=None, return_chartable=True):
        _, _, N_cistrans=self.cis_trans(format='array')
        
        self.tocsr()
        counts=self.counts


        smooth_counts=np.zeros(len(counts.data))
        smooth_counts_pred=np.zeros(len(counts.data))
        smooth_density=np.zeros(len(counts.data))
        smooth_density_pred=np.zeros(len(counts.data))
        smooth_denisty_ratio=np.zeros(len(counts.data))
        # edgesL=np.zeros(len(counts.data))
        # edgesR=np.zeros(len(counts.data))

        chr_blocks=self.chr_blocks
        chr_len_vec=self.chr_len_vec
        chr_offset_vec=np.hstack([0,chr_len_vec.cumsum()])


        chr_offset_refframe=np.hstack([0,np.cumsum(np.ceil(self.chr_len_vec/self.binsize).astype(np.int64))])
    
        rna_annot_ixs=[rna.annot_dict.get(k,-1) for k in self.annot_LUT]

        if not (isinstance(bkg, sparse.csc.csc_matrix) or isinstance(bkg, sparse.csr.csr_matrix)): 
            B=bkg.normalize_row_bytrans().counts.tocsr() #nchr X p matrix block diagonal
        else:
            B=bkg.tocsr()

        chrid=-1

        bins=cismodel['bins'] #-1
        b=bins-1
        chi=cismodel['chi'][0]
        rho=cismodel['rho'][0]
        edges=np.hstack([cismodel['bins'],cismodel['bins'][-1]+1])

        directional = False if gene_strand_dict is None else True
        dirflag=1

        # all_xL=np.zeros(len(self.x),int)
        # all_xR=np.zeros(len(self.x),int)

        for i in range(self.nannot):
            
            yall=counts.data[counts.indptr[i]:counts.indptr[i+1]]
            targetchr=np.searchsorted(chr_blocks,counts.indices[counts.indptr[i]:counts.indptr[i+1]],side='right')-1
            #targetchr=np.searchsorted(chr_blocks,counts.indices[counts.indptr[i]:counts.indptr[i+1]]) #this does not work

            cb=np.searchsorted(targetchr, np.arange(len(chr_blocks))) #index in local compressed space of where each chromosome start (length nchr+1)
            kall=counts.indices[counts.indptr[i]:counts.indptr[i+1]] #position (indices) in compressed space
            xall=self.x[kall] #positions in genome space

            #p=np.arange(counts.indptr[i],counts.indptr[i+1])

            rna_ix=rna_annot_ixs[i]
            
            if rna_ix>-1:
                ix1=rna.counts.indptr[rna_ix]
                ix2=rna.counts.indptr[rna_ix+1]
                rna_counts=rna.counts.data[ix1:ix2]
                rna_pos=rna.x[rna.counts.indices[ix1:ix2]]
                
                chrid=self.annot_chr[i]
                chr_len=chr_len_vec[chrid]

                chi=cismodel['chi'][chrid]
                rho=cismodel['rho'][chrid]

                if directional and (gene_strand_dict.get(self.annot_LUT[i], False)): #TRUE means negative strand
                    dirflag=-1
                else:
                    dirflag=1

                x1=0
                x2=chr_len-1
                delta1=(x1-rna_pos+chr_offset_vec[chrid]) if dirflag==1 else (rna_pos-x2-chr_offset_vec[chrid])
                delta2=(x2-rna_pos+chr_offset_vec[chrid]) if dirflag==1 else (rna_pos-x1-chr_offset_vec[chrid])
                binid2=np.searchsorted(b,delta2)
                binid1=np.searchsorted(b,delta1)

                Ps=chi[binid2]-rho[binid2]*(edges[binid2]-(delta2)) - chi[binid1]+rho[binid1]*(edges[binid1]-delta1) 
                ntot=np.sum(Ps*rna_counts)

                corr_factor=np.sum(rna_counts)/ntot
                print(rna_counts.sum())
                print(corr_factor)


                this_B=B[chrid,:].toarray().ravel().cumsum()
                
                for j in range(self.nchr):
                    if cb[j+1]-cb[j]>0: #2*N_thr: if there are any data on chr j
                        x=xall[cb[j]:cb[j+1]]
                        y=yall[cb[j]:cb[j+1]]
                        k=kall[cb[j]:cb[j+1]] # compressed genomic positions of the data points on chr j

                        xL0, xR0, ysmooth, fout, iL, iR=Chartable_binned.get_smoothening_intervals_2(x, y, N_thr, w_thr) #position in binned genome space
                        smooth_counts[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=ysmooth
                        #all_xL[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=xL0
                        #all_xR[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=xR0
                        # edgesL[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=xL
                        # edgesR[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=xR
                        n_gi=len(xL0)
                        


                        xL=(xL0-chr_offset_refframe[j])*self.binsize+chr_offset_vec[j] #position in true genome space
                        xR=(xR0-chr_offset_refframe[j])*self.binsize+chr_offset_vec[j]

                        smooth_density[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=ysmooth/(xR-xL)

                        isz=ysmooth==0

                        if (j==self.annot_chr[i]):
                            for i_gi in range(n_gi):
                                if ysmooth[i_gi]>0:
                                    x1=xL[i_gi]
                                    x2=xR[i_gi]
                                    delta1=(x1-rna_pos) if dirflag==1 else (rna_pos-x2)
                                    delta2=(x2-rna_pos) if dirflag==1 else (rna_pos-x1)
                                    binid2=np.searchsorted(b,delta2)
                                    binid1=np.searchsorted(b,delta1)

                                    Ps=chi[binid2]-rho[binid2]*(edges[binid2]-(delta2)) - chi[binid1]+rho[binid1]*(edges[binid1]-delta1) 
                                    tt=np.sum(Ps*rna_counts)*corr_factor
                                    smooth_counts_pred[cb[j]+counts.indptr[i]+i_gi]=tt*fout


                            # smooth_counts_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=smooth_counts[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]*fout
                            print(fout)
                            smooth_density_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=smooth_counts_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]/(xR-xL)
                        else:
                            trans_pred=(this_B[k[iR]]-this_B[k[iL]])*N_cistrans[i,1]
                            trans_pred[isz]=0
                            smooth_counts_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=trans_pred

                            smooth_density_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=trans_pred/(xR0-xL0)/self.binsize

                            # start_pos_compressed=np.searchsorted(x, xL0)
                            # stop_pos_compressed=np.searchsorted(x, xR0)

                            # delta=stop_pos_compressed-start_pos_compressed

                            # indptr=np.hstack([0,np.cumsum(delta)])
                            # n=len(delta)
                            # indices=np.zeros(indptr[-1])
                            # for i in range(n):
                            #     if indptr[i+1]>indptr[i]:
                            #         indices[indptr[i]:indptr[i+1]]=np.arange(start_pos_compressed[i],stop_pos_compressed[i])
                            # data=np.ones(len(indices), int)


                            # y_trans=B.data[B.indptr[j]:B.indptr[j+1]]*N_cistrans[i, 1]
                            # x_trans=self.x[B.indices[B.indptr[j]:B.indptr[j+1]]]
                            
                            # xL_pred, xR_pred, ysmooth_trans_pred, fout_trans=Chartable_binned.get_smoothening_intervals_(x_trans, y_trans, N_thr, w_thr)
                            # print(fout_trans)
                            # smooth_counts_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=ysmooth_trans_pred*fout_trans

                            # smooth_density_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=smooth_counts_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]/(xR_pred-xL_pred)
            
            if return_chartable:
                out_smooth_density=self.copy()
                out_smooth_density_pred=self.copy()

                out_smooth_counts=self.copy()
                out_smooth_counts_pred=self.copy()

                out_smooth_oe=self.copy()

                out_smooth_density.counts=sparse.csr_matrix((smooth_density, counts.indices, counts.indptr), shape=counts.shape)
                out_smooth_density_pred.counts=sparse.csr_matrix((smooth_density_pred, counts.indices, counts.indptr), shape=counts.shape)

                out_smooth_counts.counts=sparse.csr_matrix((smooth_counts, counts.indices, counts.indptr), shape=counts.shape)

                out_smooth_counts_pred.counts=sparse.csr_matrix((smooth_counts_pred, counts.indices, counts.indptr), shape=counts.shape)
                xx=smooth_density/smooth_density_pred
                xx[smooth_density==0]=0
                out_smooth_oe.counts=sparse.csr_matrix((xx, counts.indices, counts.indptr), shape=counts.shape)

                return out_smooth_density, out_smooth_density_pred, out_smooth_counts, out_smooth_counts_pred, out_smooth_oe
            else:
                return smooth_density, smooth_density_pred, smooth_counts, smooth_counts_pred

    def predict_smooth_2(self, cismodel, rna, bkg, N_thr=50, w_thr=250000, sig_thr=30, w_max=1000, w_max_trans=1000, gene_strand_dict=None, return_chartable=True, logger=None, normalize=True):
        

        _, _, N_cistrans=self.cis_trans(format='array')
        


        self.tocsr()
        counts=self.counts


        smooth_counts=np.zeros(len(counts.data))
        smooth_counts_pred=np.zeros(len(counts.data))
        smooth_density=np.zeros(len(counts.data))
        smooth_density_pred=np.zeros(len(counts.data))
        
        counts_sig_obs=np.zeros(len(counts.data))
        counts_sig_pred=np.zeros(len(counts.data))
        counts_sig_xL=np.zeros(len(counts.data))
        counts_sig_xR=np.zeros(len(counts.data))
        
        chr_blocks=self.chr_blocks
        chr_len_vec=self.chr_len_vec
        chr_offset_vec=np.hstack([0,chr_len_vec.cumsum()])


        chr_offset_refframe=np.hstack([0,np.cumsum(np.ceil(self.chr_len_vec/self.binsize).astype(np.int64))])
    
        rna_annot_ixs=[rna.annot_dict.get(k,-1) for k in self.annot_LUT]

        if not (isinstance(bkg, sparse.csc.csc_matrix) or isinstance(bkg, sparse.csr.csr_matrix)): 
            B=bkg.normalize_row_bytrans().counts.tocsr() #nchr X p matrix block diagonal
        else:
            B=bkg.tocsr()

        chrid=-1

        bins=cismodel['bins'] #-1
        b=bins-1
        chi=cismodel['chi'][0]
        rho=cismodel['rho'][0]
        edges=np.hstack([cismodel['bins'],cismodel['bins'][-1]+1])

        directional = False if gene_strand_dict is None else True
        dirflag=1

        # all_xL=np.zeros(len(self.x),int)
        # all_xR=np.zeros(len(self.x),int)

        for i in range(self.nannot):
            info="Started RNA %g out of %g : %s"%(i+1, self.nannot, self.annot_LUT[i])
            print(info)
            if not logger is None:
                logger.info(info)
            yall=counts.data[counts.indptr[i]:counts.indptr[i+1]]
            targetchr=np.searchsorted(chr_blocks,counts.indices[counts.indptr[i]:counts.indptr[i+1]],side='right')-1
            #targetchr=np.searchsorted(chr_blocks,counts.indices[counts.indptr[i]:counts.indptr[i+1]]) #this does not work

            cb=np.searchsorted(targetchr, np.arange(len(chr_blocks)))
            kall=counts.indices[counts.indptr[i]:counts.indptr[i+1]]
            xall=self.x[kall]

            #p=np.arange(counts.indptr[i],counts.indptr[i+1])

            rna_ix=rna_annot_ixs[i]
            #print(rna_ix)        
            if rna_ix>-1:
                ix1=rna.counts.indptr[rna_ix]
                ix2=rna.counts.indptr[rna_ix+1]
                rna_counts=rna.counts.data[ix1:ix2]
                rna_pos=rna.x[rna.counts.indices[ix1:ix2]]
                
                chrid=self.annot_chr[i]
                chr_len=chr_len_vec[chrid]

                chi=cismodel['chi'][chrid]
                rho=cismodel['rho'][chrid]

                if directional and (gene_strand_dict.get(self.annot_LUT[i], False)): #TRUE means negative strand
                    dirflag=-1
                else:
                    dirflag=1

                x1=0
                x2=chr_len-1
                delta1=(x1-rna_pos+chr_offset_vec[chrid]) if dirflag==1 else (rna_pos-x2-chr_offset_vec[chrid])
                delta2=(x2-rna_pos+chr_offset_vec[chrid]) if dirflag==1 else (rna_pos-x1-chr_offset_vec[chrid])
                binid2=np.searchsorted(b,delta2)
                binid1=np.searchsorted(b,delta1)

                Ps=chi[binid2]-rho[binid2]*(edges[binid2]-(delta2)) - chi[binid1]+rho[binid1]*(edges[binid1]-delta1) 
                ntot=np.sum(Ps*rna_counts)

                corr_factor=np.sum(rna_counts)/ntot
                #print(rna_counts.sum())
                #print(corr_factor)


                this_B=B[chrid,:].toarray().ravel().cumsum()
                
                for j in range(self.nchr):
                    if cb[j+1]-cb[j]>0: #2*N_thr:
                        x=xall[cb[j]:cb[j+1]]
                        y=yall[cb[j]:cb[j+1]]
                        ycum=np.cumsum(y)
                        k=kall[cb[j]:cb[j+1]] # compressed genomic positions of the data points on chr j

                        xL0, xR0, ysmooth, fout, iL, iR=Chartable_binned.get_smoothening_intervals_2(x, y, N_thr, w_thr)
        
                        if normalize:
                            
                            ysmooth_norm=ysmooth/(iR-iL)
                            ysmooth_norm[ysmooth==0]=0

                        else:
                            ysmooth_norm=ysmooth

                        smooth_counts[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=ysmooth_norm
                       
                        n_gi=len(xL0)


                        xL=(xL0-chr_offset_refframe[j])*self.binsize+chr_offset_vec[j]
                        xR=(xR0-chr_offset_refframe[j])*self.binsize+chr_offset_vec[j]

                        smooth_density[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=ysmooth/(xR-xL)

                        isz=ysmooth==0

                        if (j==self.annot_chr[i]):
                            for i_gi in range(n_gi):
                                if ((ysmooth[i_gi]>0) and (y[i_gi]>0)):
                                    x1=xL[i_gi]
                                    x2=xR[i_gi]
                                    delta1=(x1-rna_pos) if dirflag==1 else (rna_pos-x2)
                                    delta2=(x2-rna_pos) if dirflag==1 else (rna_pos-x1)
                                    binid2=np.searchsorted(b,delta2)
                                    binid1=np.searchsorted(b,delta1)

                                    Ps=chi[binid2]-rho[binid2]*(edges[binid2]-(delta2)) - chi[binid1]+rho[binid1]*(edges[binid1]-delta1) 
                                    tt=np.sum(Ps*rna_counts)*corr_factor
                                    
                                    smooth_counts_pred[cb[j]+counts.indptr[i]+i_gi]=tt*fout

                                



                            print(fout)
                            ypred=smooth_counts_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]
                            
                            smooth_density_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=ypred/(xR-xL)

                            if normalize:
                                ypred_norm=ypred/(iR-iL)
                                ypred_norm[ypred==0]=0
                                smooth_counts_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=ypred_norm

                            #2nd pass on significant
                            sig=(ypred<sig_thr)*(ypred>0)
                            if np.sum(sig):
                                
                                
                                if w_max>0:
                                    xL0_sig, xR0_sig=Chartable_binned.merge_intervals(xL0[sig],xR0[sig])
                                    xL0_sig_float, xR0_sig_float = Chartable_binned.split_intervals(xL0_sig, xR0_sig, w_max)
                                    ixL=np.minimum(np.searchsorted(x, xL0_sig_float), len(x)-1)
                                    ixR=np.minimum(np.searchsorted(x, xR0_sig_float), len(x)-1)
                                    xL0_sig=x[ixL]
                                    xR0_sig=x[ixR]
                                else:
                                    xL0_sig, xR0_sig=Chartable_binned.merge_intervals(xL0[sig],xR0[sig])
                                    ixL=np.minimum(np.searchsorted(x, xL0_sig), len(x)-1)
                                    ixR=np.minimum(np.searchsorted(x, xR0_sig), len(x)-1)


                                ixOK=ixR>ixL
                                ixL=ixL[ixOK]
                                ixR=ixR[ixOK]
                                xL0_sig=xL0_sig[ixOK]
                                xR0_sig=xR0_sig[ixOK]

                                ixData=np.searchsorted(x, (xL0_sig+xR0_sig)/2)
                               
                                xL_sig=(xL0_sig-chr_offset_refframe[j])*self.binsize+chr_offset_vec[j]
                                xR_sig=(xR0_sig-chr_offset_refframe[j])*self.binsize+chr_offset_vec[j]
                                n_gi_sig=len(xL0_sig)
                                for i_gi in range(n_gi_sig):
                                    
                                    x1=xL_sig[i_gi]
                                    x2=xR_sig[i_gi]
                                    delta1=(x1-rna_pos) if dirflag==1 else (rna_pos-x2)
                                    delta2=(x2-rna_pos) if dirflag==1 else (rna_pos-x1)
                                    binid2=np.searchsorted(b,delta2)
                                    binid1=np.searchsorted(b,delta1)

                                    Ps=chi[binid2]-rho[binid2]*(edges[binid2]-(delta2)) - chi[binid1]+rho[binid1]*(edges[binid1]-delta1) 
                                    tt=np.sum(Ps*rna_counts)*corr_factor
                                    
                                    counts_sig_pred[cb[j]+counts.indptr[i]+ixData[i_gi]]=tt*fout

                            
                                #counts_sig_obs[cb[j]+counts.indptr[i]+ixData_translated]=ycum[iR[ix_translator[ixR]]]-ycum[iL[ix_translator[ixL]]]+yall[iL[ix_translator[ixL]]]
                                counts_sig_obs[cb[j]+counts.indptr[i]+ixData]=ycum[ixR]-ycum[ixL]+y[ixL]-y[ixR]
                                counts_sig_xL[cb[j]+counts.indptr[i]+ixData]=xL0_sig
                                counts_sig_xR[cb[j]+counts.indptr[i]+ixData]=xR0_sig

                        else:
                            trans_pred=(this_B[k[iR]]-this_B[k[iL]])*N_cistrans[i,1]
                            
                            if normalize:
                                trans_pred_norm=trans_pred/(iR-iL)
                                trans_pred[isz]=0
                                trans_pred_norm[isz]=0
                            else:
                                trans_pred[isz]=0
                                trans_pred_norm=trans_pred

                            

                            smooth_counts_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=trans_pred_norm


                            smooth_density_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=trans_pred/(xR0-xL0)/self.binsize
                            

                            sig=(trans_pred<sig_thr)*(trans_pred>0)
                 
                            if np.sum(sig):
                                if w_max_trans>0:
                                    xL0_sig, xR0_sig=Chartable_binned.merge_intervals(xL0[sig],xR0[sig])
                                    xL0_sig_float, xR0_sig_float = Chartable_binned.split_intervals(xL0_sig, xR0_sig, w_max_trans)
                                    ixL=np.minimum(np.searchsorted(x, xL0_sig_float), len(x)-1)
                                    ixR=np.minimum(np.searchsorted(x, xR0_sig_float), len(x)-1)
                                    xL0_sig=x[ixL]
                                    xR0_sig=x[ixR]
                                else:
                                    xL0_sig, xR0_sig=Chartable_binned.merge_intervals(xL0[sig],xR0[sig])
                                    ixL=np.minimum(np.searchsorted(x, xL0_sig), len(x)-1)
                                    ixR=np.minimum(np.searchsorted(x, xR0_sig), len(x)-1)

                                ixOK=ixR>ixL
                                ixL=ixL[ixOK]
                                ixR=ixR[ixOK]
                                xL0_sig=xL0_sig[ixOK]
                                xR0_sig=xR0_sig[ixOK]

                                ixData=np.searchsorted(x, (xL0_sig+xR0_sig)/2)

                                counts_sig_pred[cb[j]+counts.indptr[i]+ixData]=(this_B[k[ixR]]-this_B[k[ixL]])*N_cistrans[i,1]

                                counts_sig_obs[cb[j]+counts.indptr[i]+ixData]=ycum[ixR]-ycum[ixL]+y[ixL]-y[ixR]
                                counts_sig_xL[cb[j]+counts.indptr[i]+ixData]=xL0_sig
                                counts_sig_xR[cb[j]+counts.indptr[i]+ixData]=xR0_sig
            
            info="Finished RNA %g out of %g : %s"%(i+1, self.nannot, self.annot_LUT[i])
            print(info)
            if not logger is None:
                logger.info(info)
            
        if return_chartable:
            out_smooth_density=self.copy()
            out_smooth_density_pred=self.copy()

            out_smooth_counts=self.copy()
            out_smooth_counts_pred=self.copy()

            out_smooth_oe=self.copy()

            out_sig_obs=self.copy()
            out_sig_pred=self.copy()
            out_sig_oe=self.copy()
            out_sig_xL=self.copy()
            out_sig_xR=self.copy()

            out_smooth_density.counts=sparse.csr_matrix((smooth_density, counts.indices, counts.indptr), shape=counts.shape)
            out_smooth_density_pred.counts=sparse.csr_matrix((smooth_density_pred, counts.indices, counts.indptr), shape=counts.shape)

            out_smooth_counts.counts=sparse.csr_matrix((smooth_counts, counts.indices, counts.indptr), shape=counts.shape)

            out_smooth_counts_pred.counts=sparse.csr_matrix((smooth_counts_pred, counts.indices, counts.indptr), shape=counts.shape)
            xx=smooth_density/smooth_density_pred
            xx[smooth_density==0]=0
            out_smooth_oe.counts=sparse.csr_matrix((xx, counts.indices, counts.indptr), shape=counts.shape)

            out_sig_obs.counts=sparse.csr_matrix((counts_sig_obs, counts.indices, counts.indptr), shape=counts.shape)
            out_sig_pred.counts=sparse.csr_matrix((counts_sig_pred, counts.indices, counts.indptr), shape=counts.shape)
            yy=counts_sig_obs/counts_sig_pred
            yy[counts_sig_obs==0]=0
            out_sig_oe.counts=sparse.csr_matrix((yy, counts.indices, counts.indptr), shape=counts.shape)
            out_sig_xL.counts=sparse.csr_matrix((counts_sig_xL, counts.indices, counts.indptr), shape=counts.shape)
            out_sig_xR.counts=sparse.csr_matrix((counts_sig_xR, counts.indices, counts.indptr), shape=counts.shape)

            out_sig_obs.compress_at([out_sig_pred, out_sig_oe, out_sig_xL, out_sig_xR])

            return out_smooth_density, out_smooth_density_pred, out_smooth_counts, out_smooth_counts_pred, out_smooth_oe, out_sig_obs,out_sig_pred, out_sig_oe, out_sig_xL, out_sig_xR
        
        else:
            return smooth_density, smooth_density_pred, smooth_counts, smooth_counts_pred
    
    def predict_smooth_pval(self, cismodel, rna, bkg, N_thr=50, w_thr=250000, gene_strand_dict=None, return_chartable=True):
        _, _, N_cistrans=self.cis_trans(format='array')
        
        self.tocsr()
        counts=self.counts

        nbychr=self.N_bychr().toarray()

        smooth_counts=np.zeros(len(counts.data))
        smooth_counts_pred=np.zeros(len(counts.data))
        smooth_density=np.zeros(len(counts.data))
        smooth_density_pred=np.zeros(len(counts.data))
        #smooth_denisty_ratio=np.zeros(len(counts.data))


        flagout=np.zeros(len(counts.data), dtype=bool)
        pvals=np.ones(len(counts.data))
        edgesL=np.zeros(len(counts.data), dtype=int)
        edgesR=np.zeros(len(counts.data), dtype=int)

        chr_blocks=self.chr_blocks
        chr_len_vec=self.chr_len_vec
        chr_offset_vec=np.hstack([0,chr_len_vec.cumsum()])


        chr_offset_refframe=np.hstack([0,np.cumsum(np.ceil(self.chr_len_vec/self.binsize).astype(np.int64))])

    
        rna_annot_ixs=[rna.annot_dict.get(k,-1) for k in self.annot_LUT]


        B=bkg.normalize_row_bychr().sum().counts.tocsr()
        this_B=B.toarray().ravel().cumsum()
        # if not (isinstance(bkg, sparse.csc.csc_matrix) or isinstance(bkg, sparse.csr.csr_matrix)): 
        #     B=bkg.normalize_row_bytrans().counts.tocsr() #nchr X p matrix block diagonal

        # else:
        #     B=bkg.tocsr()

        chrid=-1

        bins=cismodel['bins'] #-1
        b=bins-1
        chi=cismodel['chi'][0]
        rho=cismodel['rho'][0]
        edges=np.hstack([cismodel['bins'],cismodel['bins'][-1]+1])

        directional = False if gene_strand_dict is None else True
        dirflag=1

        # all_xL=np.zeros(len(self.x),int)
        # all_xR=np.zeros(len(self.x),int)

        for i in range(self.nannot):
            yall=counts.data[counts.indptr[i]:counts.indptr[i+1]]
            targetchr=np.searchsorted(chr_blocks,counts.indices[counts.indptr[i]:counts.indptr[i+1]],side='right')-1
            #targetchr=np.searchsorted(chr_blocks,counts.indices[counts.indptr[i]:counts.indptr[i+1]]) #this does not work
            cb=np.searchsorted(targetchr, np.arange(len(chr_blocks))) #index in local compressed space of where each chromosome start (length nchr+1)
            kall=counts.indices[counts.indptr[i]:counts.indptr[i+1]] #position (indices) in compressed space
            xall=self.x[kall] #positions in genome space
            #p=np.arange(counts.indptr[i],counts.indptr[i+1])
            rna_ix=rna_annot_ixs[i]
            if rna_ix>-1:
                ix1=rna.counts.indptr[rna_ix]
                ix2=rna.counts.indptr[rna_ix+1]
                rna_counts=rna.counts.data[ix1:ix2]
                rna_pos=rna.x[rna.counts.indices[ix1:ix2]]
                
                chrid=self.annot_chr[i]
                chr_len=chr_len_vec[chrid]

                chi=cismodel['chi'][chrid]
                rho=cismodel['rho'][chrid]

                if directional and (gene_strand_dict.get(self.annot_LUT[i], False)): #TRUE means negative strand
                    dirflag=-1
                else:
                    dirflag=1

                x1=0
                x2=chr_len-1
                delta1=(x1-rna_pos+chr_offset_vec[chrid]) if dirflag==1 else (rna_pos-x2-chr_offset_vec[chrid])
                delta2=(x2-rna_pos+chr_offset_vec[chrid]) if dirflag==1 else (rna_pos-x1-chr_offset_vec[chrid])
                binid2=np.searchsorted(b,delta2)
                binid1=np.searchsorted(b,delta1)

                Ps=chi[binid2]-rho[binid2]*(edges[binid2]-(delta2)) - chi[binid1]+rho[binid1]*(edges[binid1]-delta1) 
                ntot=np.sum(Ps*rna_counts)

                corr_factor=np.sum(rna_counts)/ntot
                print(rna_counts.sum())
                print(corr_factor)


                #this_B=B[chrid,:].toarray().ravel().cumsum()
                
                for j in range(self.nchr):
                    if cb[j+1]-cb[j]>0: #2*N_thr: if there are any data on chr j
                        x=xall[cb[j]:cb[j+1]]
                        y=yall[cb[j]:cb[j+1]]
                        k=kall[cb[j]:cb[j+1]] # compressed genomic positions of the data points on chr j

                        n_in_chr=y.sum()
                        xL0, xR0, ysmooth, fout, iL, iR=Chartable_binned.get_smoothening_intervals_2(x, y, N_thr, w_thr) #position in binned genome space
                        smooth_counts[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=ysmooth
                        #all_xL[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=xL0
                        #all_xR[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=xR0
                        # edgesL[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=xL
                        # edgesR[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=xR
                        n_gi=len(xL0)

                        xL=(xL0-chr_offset_refframe[j])*self.binsize+chr_offset_vec[j] #position in true genome space
                        xR=(xR0-chr_offset_refframe[j])*self.binsize+chr_offset_vec[j]

                        sd=ysmooth/(xR-xL)
                        smooth_density[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=sd

                        isz=ysmooth==0

                        sd[isz]=0


                        if (j==self.annot_chr[i]):
                            for i_gi in range(n_gi):
                                if ysmooth[i_gi]>0:
                                    x1=xL[i_gi]
                                    x2=xR[i_gi]
                                    delta1=(x1-rna_pos) if dirflag==1 else (rna_pos-x2)
                                    delta2=(x2-rna_pos) if dirflag==1 else (rna_pos-x1)
                                    binid2=np.searchsorted(b,delta2)
                                    binid1=np.searchsorted(b,delta1)

                                    Ps=chi[binid2]-rho[binid2]*(edges[binid2]-(delta2)) - chi[binid1]+rho[binid1]*(edges[binid1]-delta1) 
                                    tt=np.sum(Ps*rna_counts)*corr_factor
                                    smooth_counts_pred[cb[j]+counts.indptr[i]+i_gi]=tt #*fout


                            # smooth_counts_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=smooth_counts[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]*fout
                            print(fout)
                            scp=smooth_counts_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]
                            sdp=scp/(xR-xL)
                            smooth_density_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=sdp
                        else:
                            #print(nbychr[i,self.annot_chr[i]])
                            #print(this_B[k[iR]]-this_B[k[iL]])
                            trans_pred=(this_B[k[iR]]-this_B[k[iL]])*nbychr[i,self.annot_chr[i]] #N_cistrans[i,1]
                            trans_pred[isz]=0
                            scp=trans_pred
                            smooth_counts_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=trans_pred

                            sdp=trans_pred/(xR0-xL0)/self.binsize
                            smooth_density_pred[(cb[j]+counts.indptr[i]):(cb[j+1]+counts.indptr[i])]=sdp

                            #select significant ones!
                            # sort by density, select ones where counts measured > counts bkg + 2, do binomial test with p_success=nbk/ntot, n_try=counts rna -1 , nsuccess=counts oberved in bin -1 
                            
                        density_rank=np.argsort(-sd)

                        i_gi=0
                        goon=True
                        while goon and (i_gi<len(density_rank)):
                            ix=density_rank[i_gi]
                            if not(flagout[cb[j]+counts.indptr[i]+ix]):
                                if (ysmooth[ix]>(scp[ix]+2)):
                                    ntot=np.sum(y)
                                    pval=stats.binom_test(ysmooth[ix]-1, ntot-1, scp[ix]/ntot)
                                    if pval<0.05:
                                        pvals[cb[j]+counts.indptr[i]+ix]=pval
                                        flagout[(cb[j]+counts.indptr[i]+iL[ix]):(cb[j]+counts.indptr[i]+iR[ix])]=True
                                        edgesL[cb[j]+counts.indptr[i]+ix]=k[iL[ix]]
                                        edgesR[cb[j]+counts.indptr[i]+ix]=k[iR[ix]]
                                    else:
                                        goon=True
                                else:
                                    goon=True
                            i_gi+=1
            sigix=pvals<0.05
            sigintL=edgesL[sigix]
            sigintR=edgesR[sigix]
            sigintP=pvals[sigix]

            chr_ids=np.searchsorted(chr_blocks,sigintL,side='right')-1
            chr_ids2=np.searchsorted(chr_blocks,sigintR,side='right')-1
            print(sigintL)
            print(sigintR)
            print(chr_ids)
            print(chr_ids2)
            xL_real=np.array((np.array(self.x[sigintL], dtype=np.int64)-chr_offset_refframe[chr_ids])*self.binsize, dtype=int)
            xR_real=np.array((np.array(self.x[sigintR], dtype=np.int64)-chr_offset_refframe[chr_ids])*self.binsize, dtype=int)

            chr_dict_rev={i:k for i,k in enumerate(self.chr_LUT)}

            sigint={'chr':[chr_dict_rev[kk] for kk in chr_ids], 'chr2':[chr_dict_rev[kk] for kk in chr_ids2], 'start':xL_real,'stop':xR_real,'p':sigintP, 'observed':smooth_counts[sigix], 'model':smooth_counts_pred[sigix]}
            if return_chartable:
                out_smooth_density=self.copy()
                out_smooth_density_pred=self.copy()

                out_smooth_counts=self.copy()
                out_smooth_counts_pred=self.copy()

                out_smooth_oe=self.copy()

                out_smooth_density.counts=sparse.csr_matrix((smooth_density, counts.indices, counts.indptr), shape=counts.shape)
                out_smooth_density_pred.counts=sparse.csr_matrix((smooth_density_pred, counts.indices, counts.indptr), shape=counts.shape)

                out_smooth_counts.counts=sparse.csr_matrix((smooth_counts, counts.indices, counts.indptr), shape=counts.shape)

                out_smooth_counts_pred.counts=sparse.csr_matrix((smooth_counts_pred, counts.indices, counts.indptr), shape=counts.shape)
                xx=smooth_density/smooth_density_pred
                xx[smooth_density==0]=0
                out_smooth_oe.counts=sparse.csr_matrix((xx, counts.indices, counts.indptr), shape=counts.shape)

                return out_smooth_density, out_smooth_density_pred, out_smooth_counts, out_smooth_counts_pred, out_smooth_oe, sigint
            else:
                return smooth_density, smooth_density_pred, smooth_counts, smooth_counts_pred, sigint


    @staticmethod
    def merge_intervals(xL, xR):
        ixs=np.argsort(xL)
        merged=[]
        ixL=[]
        ixR=[]
        for j, i in enumerate(ixs):
            if j==0:
                merged.append((xL[i],xR[i]))
                ixL.append(i)
                ixR.append(i)
                
            else:
                lower=merged[-1]
                if lower[1]<xL[i]:
                    merged.append((xL[i],xR[i]))
                    ixL.append(i)
                    ixR.append(i)
                else:
                    merged[-1]=(lower[0], max(xR[i],lower[1]))
                    if xR[i]>lower[1]:
                        ixR[-1]=i

        xL_merged=np.array([k[0] for k in merged])
        xR_merged=np.array([k[1] for k in merged])
        return xL_merged, xR_merged

    @staticmethod
    def merge_intervals_2(xL, xR, chrom):
        ixs=np.argsort(xL)
        merged=[]
        ixL=[]
        ixR=[]
        chrom_merged=[]
        for j, i in enumerate(ixs):
            if j==0:
                merged.append((xL[i],xR[i]))
                ixL.append(i)
                ixR.append(i)
                chrom_merged.append(chrom[i])
            else:
                lower=merged[-1]
                if lower[1]<xL[i]:
                    merged.append((xL[i],xR[i]))
                    ixL.append(i)
                    ixR.append(i)
                    chrom_merged.append(chrom[i])
                else:
                    merged[-1]=(lower[0], max(xR[i],lower[1]))
                    if xR[i]>lower[1]:
                        ixR[-1]=i

        xL_merged=np.array([k[0] for k in merged])
        xR_merged=np.array([k[1] for k in merged])
        chrom_merged=np.array(chrom_merged)
        return xL_merged, xR_merged, chrom_merged
    
    @staticmethod
    def split_intervals(xL, xR, w_max):

        split=[]
        for i in range(len(xL)):
            x1=xL[i]
            x2=xR[i]
            d=x2-x1
            if d<w_max:
                split.append((x1,x2))
            else:
                n=int(np.floor(d/w_max))
                l=x1+(d-n*w_max)/2
                #r=m[1]-(d-n*w_thr)/2
                p=np.arange(l, x2+1, w_max)

                #ps=x[np.minimum(np.searchsorted(x, p), len(x)-1)]
                for i in range(len(p)-1):
                    split.append((p[i],p[i+1]))
                    #split.append((x[np.minimum(np.searchsorted(x, p[i]), len(x)-1)],x[np.minimum(np.searchsorted(x, p[i+1]), len(x)-1)]))


        xL_split=np.array([k[0] for k in split])
        xR_split=np.array([k[1] for k in split])
        
        return xL_split, xR_split



    def make_metagene(self, bedfile, bins, bkg, meta_matrix=None, feature_id=5, strand_id=4, feature_info_id = None, stranded = True):
    
        data=self
        #data.iscompressed=True
        
        if type(bedfile)==str:
            curtain=data.bed_to_curtain(bedfile, curtain=bins, feature_id=feature_id, integrate=True, strand_id=strand_id)
        else:
            curtain=data.df_to_curtain(bedfile, curtain=bins, feature_info_id=feature_info_id, integrate=True, stranded=stranded)
        
        if meta_matrix is None:
            meta_matrix=data.meta_matrix(curtain, aggregate_features=True)

        out=data.make_data_curtain(curtain, aggregate_features=True, meta_matrix=meta_matrix)
        out.counts=out.counts*sparse.diags(1000/curtain['w'], format='csr')
        
        N=np.array(data.counts.sum(axis=1)).flatten()
        
        out_bkg=bkg.make_data_curtain(curtain, aggregate_features=True, meta_matrix=meta_matrix).counts.toarray().flatten()
        
        N_bkg=bkg.counts.sum()
        
        out_bkg=out_bkg*(1000/curtain['w'])*(1e6/N_bkg)
        
        curtainized_data={'data':out, 'N':N, 'bkg':out_bkg, 'curtain':curtain}
        return curtainized_data

    def mask_gene_body(self, genes):
        binsize=self.binsize
        if self.iscompressed():
            chr_blocks_refframe=np.hstack([0,np.cumsum(np.ceil(self.chr_len_vec/self.binsize).astype(np.int64))])
            chr_offset_dict={self.chr_LUT[i]:v for i, v in enumerate(chr_blocks_refframe[0:-1])}
            # chr_max_pos_dict={self.chr_LUT[i]:v for i, v in enumerate(chr_blocks_refframe[1:])}
        else:
            chr_offset_dict={self.chr_LUT[i]:v for i, v in enumerate(self.chr_blocks[0:-1])}
            # chr_max_pos_dict={self.chr_LUT[i]:v for i, v in enumerate(self.chr_blocks[1:])}


        n=len(self.annot_LUT)
            
        start_pos_abs=genes['start'].loc[self.annot_LUT].values
        stop_pos_abs=genes['stop'].loc[self.annot_LUT].values
        start_pos=np.zeros(n, np.int64)
        stop_pos=np.zeros(n, np.int64)
            
        for i in range(n):
            this_offset=chr_offset_dict[self.chr_LUT[self.annot_chr[i]]]
            start_pos[i]=int(start_pos_abs[i]/binsize)+this_offset
            stop_pos[i]=int(stop_pos_abs[i]/binsize)+this_offset
            

        start_pos_compressed=np.searchsorted(self.x, start_pos)
        stop_pos_compressed=np.searchsorted(self.x, stop_pos)

        
        out=self.copy(copydata=True)
        c=out.counts
        
        out_in=self.copy(copydata=True)
        c_in=out_in.counts
        
        #print(len(c.data))
        is_in_body=np.zeros(len(c.data), bool)
        for i in range(n):
            ixs=c.indices[c.indptr[i]:c.indptr[i+1]]
            is_in_body[c.indptr[i]:c.indptr[i+1]]=((ixs>=start_pos_compressed[i]) & (ixs<=stop_pos_compressed[i]))
    #        print(sum(is_in_body))
        c.data[is_in_body]=0
        c_in.data[~is_in_body]=0
        
        #c.eliminate_zeros()
        
        return out, out_in

    def sparsify_bychr(self, bkg=None, b_prenormalized=False):

        if bkg is None:
            bkg_flat = self.copy(copydata=False)
            bkg_flat.x = self.x.copy()
            bkg_flat.counts=sparse.csr_matrix(np.ones((1,len(bkg_flat.x)),  dtype=int))
            bkg_mat=bkg_flat.normalize_row_bychr(expand=False).counts.toarray().ravel()

        else:
            if not b_prenormalized:
                bkg_mat=bkg.normalize_row_bychr(expand=False).counts.toarray().ravel()
            else:
                bkg_mat=bkg.counts.toarray().ravel()

        out = self.copy(copydata = False)
        counts = self.counts
        out_data = np.zeros(counts.data.shape, int)
        for i in range(counts.shape[0]):
            this_data = counts.data[counts.indptr[i]:counts.indptr[i+1]]
            this_idx = counts.indices[counts.indptr[i]:counts.indptr[i+1]]
            this_bkg = bkg_mat[this_idx]

            chr_ids=np.searchsorted(self.chr_blocks,this_idx,side='right')-1
            chr_blocks=np.searchsorted(chr_ids,np.arange(self.nchr+1))


            for j in range(len(chr_blocks)-1):

                if chr_blocks[j+1]>chr_blocks[j]:
                    this_data_chr = this_data[chr_blocks[j]:chr_blocks[j+1]]
                    S = this_data_chr.sum()
                    
                    if S>0:
                        Sbkg = this_bkg[chr_blocks[j]:chr_blocks[j+1]]
                        if Sbkg.sum()>0:
                            probas = this_data_chr * this_bkg[chr_blocks[j]:chr_blocks[j+1]]
                            probas = probas / probas.sum()
                            #print(probas)
                        else:
                            probas = this_data_chr / S
                    
                    
                        nsample = int(np.ceil(S))
                    
                        new_data = multinomial(nsample, probas)

                        out_data[(counts.indptr[i]+chr_blocks[j]):(counts.indptr[i]+chr_blocks[j+1])]=new_data
                
        out.counts = sparse.csr_matrix((out_data, counts.indices.copy(), counts.indptr.copy()), shape = counts.shape)
        out.x = self.x.copy()
        return out

    def relevel_bychr(self, ref, genes=None):

        if not genes is None:
            ref2 = ref.select_annots2(self.annot_LUT.copy(), genes)
        else:
            ref2 = ref

        out = self.copy(copydata = False)
        counts = self.counts
        counts_ref = ref2.counts
        out_data = np.zeros(counts.data.shape, int)
        for i in range(counts.shape[0]):
            this_data = counts.data[counts.indptr[i]:counts.indptr[i+1]]
            this_idx = counts.indices[counts.indptr[i]:counts.indptr[i+1]]
          
            chr_ids=np.searchsorted(self.chr_blocks,this_idx,side='right')-1
            chr_blocks=np.searchsorted(chr_ids,np.arange(self.nchr+1))


            this_data_ref = counts_ref.data[counts_ref.indptr[i]:counts_ref.indptr[i+1]]
            Sref = this_data_ref.sum()
            # this_idx_ref = counts_ref.indices[counts_ref.indptr[i]:counts_ref.indptr[i+1]]
          
            # chr_ids_ref=np.searchsorted(ref2.chr_blocks,this_idx_ref,side='right')-1
            # chr_blocks_ref=np.searchsorted(chr_ids_ref,np.arange(ref2.nchr+1))


            for j in range(len(chr_blocks)-1):

                if chr_blocks[j+1]>chr_blocks[j]:
                    this_data_chr = this_data[chr_blocks[j]:chr_blocks[j+1]]
                    #this_data_chr_ref = this_data_ref[chr_blocks_ref[j]:chr_blocks_ref[j+1]]
                    Sx = 1.0*(this_data_chr.sum())
                    #Sref = this_data_chr_ref.sum()

                    if (Sx>0) and (Sref>0):
                        probas = this_data_chr / Sx
                    
                        nsample = int(np.ceil(Sref))
                    
                        new_data = multinomial(nsample, probas)

                        out_data[(counts.indptr[i]+chr_blocks[j]):(counts.indptr[i]+chr_blocks[j+1])]=new_data
                
        out.counts = sparse.csr_matrix((out_data, counts.indices.copy(), counts.indptr.copy()), shape = counts.shape)
        out.x = self.x.copy()
        return out

    def predict_bychr_bkgonly(self, bkg=None, b_prenormalized=False, Nbychr = None):
        if bkg is None:
            bkg_flat = self.copy(copydata=False)
            bkg_flat.x = self.x.copy()
            bkg_flat.counts=sparse.csr_matrix(np.ones((1,len(bkg_flat.x)),  dtype=int))
            bkg_mat=bkg_flat.normalize_row_bychr(expand=False).counts.toarray().ravel()

        else:
            if not b_prenormalized:
                bkg_mat=bkg.normalize_row_bychr(expand=False).counts.toarray().ravel()
            else:
                bkg_mat=bkg.counts.toarray().ravel()

        if Nbychr is None:
            Nbychr = self.N_bychr()

       
        chr_blocks = self.chr_blocks
        Nj = np.array(Nbychr.sum(axis=0)).ravel()
            
        alldata = []
        for j in range(len(chr_blocks)-1):
            
            
            p = bkg_mat[chr_blocks[j]:chr_blocks[j+1]]
            p = p / p.sum()
            
            if Nj[j]>0:
                y = multinomial(Nj[j], p)

                this_data=np.ones(Nj[j], int)
                this_rowindices=np.arange(Nj[j])
                this_colindptr = np.cumsum(np.hstack([0, y]))

                this_colindices = sparse.csc_matrix((this_data, this_rowindices, this_colindptr), dtype=int, shape=(Nj[j], len(p))).tocsr().indices

                np.random.shuffle(this_colindices)
                this_rowindptr =  np.cumsum(np.hstack([0, Nbychr[:,j].toarray().ravel()]))

                alldata += [sparse.csr_matrix((this_data, this_colindices, this_rowindptr), shape=(Nbychr.shape[0],p.shape[0]), dtype=int)]
            else:
                alldata += [sparse.csr_matrix(shape=(Nbychr.shape[0],p.shape[0]), dtype=int)]

        out = self.copy(copydata = False)
        out.counts= sparse.hstack(alldata).tocsr()
        out.x = self.x.copy()
        return out

def get_metagene_values_raw_(cd, select_annot_ixs=None):
    data=cd['data']
    if not select_annot_ixs is None:
        
        tdc=data.annot_dict
        annots_in=[k for k in select_annot_ixs if (k in tdc)]

        ixs=np.array([tdc[k] for k in annots_in])
#         ixs=np.array([i for i, k in enumerate(data.annot_LUT) if k in select_annot_ixs],int)
        
        N=cd['N'][ixs].sum()
#         print(N)
        y=np.array(data.select_annots(select_annot_ixs).counts.sum(axis=0)).flatten()
        
    else:
        N=cd['N'].sum()
        y=np.array(data.counts.sum(axis=0)).flatten()
        
    
    y=y*(1e6/N)
    
    x=cd['curtain']['x']
    
    y_bkg=cd['bkg']
    
    return x, y, y_bkg


def get_metagene_values_normalized_(cd, cd_coarse, correct_local=True, select_annot_ixs=None):

    x, y, y_bkg = get_metagene_values_raw_(cd, select_annot_ixs=select_annot_ixs)
    
    if not cd_coarse is None:
        _, y_coarse, y_bkg_coarse = get_metagene_values_raw_(cd_coarse, select_annot_ixs=select_annot_ixs)
    
        m=np.median(y_coarse)
        m_bkg=np.median(y_bkg_coarse)

        if correct_local:
            y=y/m
        else:
            y=y/m_bkg
        
        y_bkg=y_bkg/m_bkg

    return x, y, y_bkg
    

def get_metagene_values(mg, mg_coarse=None, select_annot=None, correct_local=False):
    if mg_coarse is None:
        x, y, y_bkg = get_metagene_values_raw_(mg, select_annot_ixs=select_annot)
    else:
        x, y, y_bkg = get_metagene_values_normalized_(mg, mg_coarse, select_annot_ixs=select_annot, correct_local=correct_local)
    return x, y, y_bkg

def group_gi_by_sets(gi, features, no_double_count=True):
    lud={k: i for i, k in enumerate(gi['feature_name'])}
    
    feature_start=np.zeros(0, int)
    feature_stop=np.zeros(0, int)
    feature_chr=np.zeros(0, int)
    feature_name=[]

    transfer_data=np.zeros(0,int)
    # transfer_indices=np.zeros(0,int)
    transfer_indptr=np.zeros(len(features)+1,int)
    ni=0
    i=0

    features_order=[]
    for features_set_name, features_set_members in features.items():
        features_order.append(features_set_name)
        ix_gi=[lud.get(e,-1) for e in features_set_members if e in lud]
        if len(ix_gi)>0:
            l=gi['feature_start'][ix_gi]
            r=gi['feature_stop'][ix_gi]
            chrom=gi['feature_chr'][ix_gi]

            l_merged, r_merged, chrom_merged=Chartable_binned.merge_intervals_2(l, r, chrom)
            feature_start=np.append(feature_start,l_merged)
            feature_stop=np.append(feature_stop,r_merged)
            feature_chr=np.append(feature_chr,chrom_merged)
            feature_name+=[features_set_name]*len(l_merged)
            transfer_data=np.append(transfer_data,np.ones(len(l_merged)))
            # transfer_indices=np.append(transfer_indices,np.arange(ix_gi))
            ni+=len(l_merged)

        i+=1
        transfer_indptr[i]=ni
    transfer_indices=np.arange(len(transfer_data))
    M=sparse.csc_matrix((transfer_data, transfer_indices, transfer_indptr), shape=(len(feature_start),len(features)))


    gi_out={'feature_start':feature_start, 'feature_stop':feature_stop, 'feature_name':feature_name, 'feature_chr':feature_chr, 'summing_matrix':M, 'set_names':features_order, 'binsize':gi['binsize']}
    return gi_out


def projector_subset(proj, features, no_double_count=True):
    lud={k: i for i, k in enumerate(proj['feature_name'])}
    transfer_data=np.ones(0, int)
    transfer_indices=np.ones(0, int)
    transfer_indptr=np.zeros(1,int)
    n=0
    for features_set_name, features_set_members in features.items():
        ix_in=np.array([i for i, e in enumerate(features_set_members) if e in lud])
        if len(ix_in)>0:
            transfer_data=np.append(transfer_data,np.ones(len(ix_in), int))
            transfer_indices=np.append(transfer_indices, ix_in)
            n+=len(ix_in)
        transfer_indptr=np.append(transfer_indptr,n)
        

    transfer_matrix=sparse.csc_matrix((transfer_data,transfer_indices, transfer_indptr), shape=(proj['counting_matrix'].shape[1], len(features)))
    
    new_counting_matrix=proj['counting_matrix']*transfer_matrix
    new_counting_matrix.eliminate_zeros()
    if no_double_count:
        new_counting_matrix=sparse.csc_matrix((np.ones(len(new_counting_matrix.data)), new_counting_matrix.indices, new_counting_matrix.indptr), shape=(proj['counting_matrix'].shape[0], len(features)))
#         out.counts=(self.counts*transfer_matrix).tocsr()
    
    new_chr_blocks=len(features)*np.ones(len(proj['chr_blocks']),int)
    new_chr_blocks[0]=0
    new_proj={'counting_matrix':new_counting_matrix, 'feature_name':[k for k in features.keys()], 'feature_info':[k for k in features.keys()], 'chr_blocks':new_chr_blocks, 'w':np.zeros(len(features),int), 'feature_start':np.zeros(len(features),int), 'feature_stop':np.zeros(len(features),int), 'feature_chr':np.zeros(len(features),int)}
    return new_proj


    
    

    # def predict_trans_bychr(self, bkg, leakage_matrix=None):
    #     #make sure to pass expected bkg_model as a nchr X p matrix block diagonal
    #     N=self.get_chr_matrix()
    #     B=bkg.normalize_row_bychr(expand=True)
    #     if not leakage_vector is None:   
    #         L=sparse.diags(leakage_vector, format='csr')
    #     else: #we use existing trans counts
    #         L=1
    #     out=(N*L)*B 

    #     return out

        #takes the data and predict the count


        
        

        
            
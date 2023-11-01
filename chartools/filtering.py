from tagtools import utils as tagutils 
from .bamutils import fetch_ID
import sys
import pysam

class AmbivReader_from_pairs():
    def __init__(self, pairsfile, bygene=True, tx_dict=None):
        self.bygene=bygene
        if bygene:
            self.agtag = 24
            self.tx_dict=tx_dict
           
        else:
            self.agtag = 23

        self.reader=open(pairsfile,'r')

    def __iter__(self):
        return self

    def __next__(self):
        line=self.reader.readline()
        
        if len(line)<2:
            self.reader.close()
            raise StopIteration
        else:
            read_data=line.strip().split("\t")
            # read_data=read_data[0,17,self.agtag]
            ag=int(read_data[self.agtag])
            chro=read_data[3]
            pos=int(read_data[4])
            chroR=read_data[1][2:]
            posR=int(read_data[2])
            if self.bygene:
                name=self.tx_dict.get(read_data[17],"*")
            else:
                name=read_data[17]
         
            

            return read_data[0], name, ag, line, chro, pos, chroR, posR

    def close(self):
        self.reader.close()


def _make_tx_LUD(tx_list_list):
    tx_LUD={}
    # n=0
    for ix, tx_list in enumerate(tx_list_list):
        for tx in tx_list:
            if not tx in tx_LUD:
                tx_LUD[tx]=ix
                # n+=1
    return tx_LUD #tx_name to ix for fanout dictionnary

class ReadsFanout():
    def __init__(self, tx_list_list, ambivalence_LUT, bygene=False, tx_dict=None, group_names=None, min_flight=0, locus_LUD=None, noambig=False):
        # if isinstance(gps, list):
        #     pass
        self.bygene=bygene
        self.tx_LUD=_make_tx_LUD(tx_list_list)
        self.tx_dict=tx_dict #ENST: ENSG
        self.min_flight=min_flight
        self.locus_LUD=locus_LUD
        self.noambig=noambig
        nmax=max((max(v) for _, v in ambivalence_LUT.items()))+1
        self.ag_truth_table=[tagutils.get_ambiv_truthtable(ambivalence_LUT, nmax, tx_list) for tx_list in tx_list_list]

        
        if group_names is None:
            group_names=[str(i) for i in range(len(tx_list_list))]
        self.group_names=group_names
    
    def _classify(self, name, ag): #note that name is either tx name or gene name depending on bygene flag
        if self.noambig and (ag>1):
            return -1
        else:
            c=self.tx_LUD.get(name,-1)
            if c>-1:
                return c
    
            for i in range(len(self.ag_truth_table)):
                if self.ag_truth_table[i][ag]: #the first group to which the ag has a non zero intersection is returned
                    return i
            return -1

    def _classify_with_flight(self, name, ag, chro, pos, chroR, posR): #this one does not allow ambivalence
        if ag>1:
            return -1
        else:
            c=self.tx_LUD.get(name,-1)
            if c>-1:
                if ((not(chro==chroR)) or (pos<(posR-self.min_flight)) or (pos>(int(posR+self.min_flight)))):
                    return c
                else:
                    return -1 #return (c+len(self.ag_truth_table))

                # locus=self.locus_LUD.get(name,None)

                # if not locus is None:
                #     # print(chro)
                    #  if ((not(chro==locus[0])) or (pos<(int(locus[1])-self.min_flight)) or (pos>(int(locus[2])+self.min_flight))):
                        
                #         return c
                #     else:
                #         return (c+len(self.ag_truth_table))
            
            # return c
        return -1
        
    def _fanout_reads_stream(self, ambiviter, readid_out_streams, ft_out_streams, nmax=0):
        ngps=len(self.ag_truth_table)
        # if self.min_flight>0:
        #     n=[0]*(2*ngps+1)
        # else:
        #     n=[0]*(ngps+1)
        n=[0]*(ngps+1)
        ntot=0
        while ((ntot<nmax) or nmax==0):
            try:
                readid, name, ag, ft, chro, pos, chroR, posR = ambiviter.__next__()
                
            except StopIteration:
                break
            # if self.noambig:
            #     ag=int(0)
            ntot+=1
            if self.min_flight>0:
                c=self._classify_with_flight(name, abs(ag), chro, pos, chroR, posR)
            else:
                c=self._classify(name, abs(ag))
            if c>-1:
                if not readid_out_streams is None:
                    readid_out_streams[c].write(readid)
                    readid_out_streams[c].write("\n")
                if not ft_out_streams is None:
                    ft_out_streams[c].write(ft)
            n[c+1]+=1
            
        return n

    def _fanout_reads_bams_stream(self, ambiviter, rnabamiter, dnabamiter, readid_out_streams, ft_out_streams, rna_out_streams, dna_out_streams, nmax=0):
        ngps=len(self.ag_truth_table)
        # if self.min_flight>0:
        #     n=[0]*(2*ngps+1)
        # else:
        #     n=[0]*(ngps+1)
        n=[0]*(ngps+1)
        ntot=0
        try:
            rna=rnabamiter.__next__()
        except StopIteration:
            rna=None
        try:
            dna=dnabamiter.__next__()
        except StopIteration:
            dna=None

        while ((ntot<nmax) or nmax==0):
            try:
                readid, name, ag, ft, chro, pos, chroR, posR = ambiviter.__next__()
                
            except StopIteration:
                break
            # if self.noambig:
            #     ag=int(0)
            ntot+=1
            if self.min_flight>0:
                c=self._classify_with_flight(name, abs(ag), chro, pos, chroR, posR)
            else:
                c=self._classify(name, abs(ag))
            
            if c>-1:
                if not readid_out_streams is None:
                    readid_out_streams[c].write(readid)
                    readid_out_streams[c].write("\n")
                if not ft_out_streams is None:
                    ft_out_streams[c].write(ft)
            
                #seek the rna
                rec_rna, rna= fetch_ID(rnabamiter, rna, readid)
                rec_dna, dna= fetch_ID(dnabamiter, dna, readid)
                if len(rec_rna)>0:
                    for r in rec_rna:
                        rna_out_streams[c].write(r)
                if len(rec_dna)>0:
                    for r in rec_dna:
                        dna_out_streams[c].write(r)

            n[c+1]+=1
            
        return n

    def fanout(self, pairs_in, rna_bam_in=None, dna_bam_in=None, readid_out_prfx=None, pairs_out_prfx=None, rnadna_out_prfx=None, to_stdout=None, indexbams=False, nmax=0):
        ambiviter=AmbivReader_from_pairs(pairs_in, bygene=self.bygene, tx_dict=self.tx_dict)
        
        if to_stdout=="readids":
            readid_out_streams=[sys.stdout]*len(self.group_names)
        else:
            if readid_out_prfx is None:           
                readid_out_streams=None
            else:
                readid_out_streams=[open(readid_out_prfx+name+".txt","w") for name in self.group_names]

        
        if to_stdout=="pairs":
            ft_out_streams=[sys.stdout]*len(self.group_names)
        else:
            if pairs_out_prfx is None:
                ft_out_streams=None
            else:
                ft_out_streams=[open(pairs_out_prfx+name+".pairs","w") for name in self.group_names]
            
        if rna_bam_in is None:
            n=self._fanout_reads_stream(ambiviter, readid_out_streams, ft_out_streams, nmax=nmax)
        else:

            inr = pysam.AlignmentFile(rna_bam_in)
            ind = pysam.AlignmentFile(dna_bam_in)
            rnabamiter = inr.fetch(until_eof=True)
            dnabamiter = ind.fetch(until_eof=True)

            if to_stdout=="rna":
                rna_out_streams=[pysam.AlignmentFile("/dev/stdout","w", template=inr)]*len(self.group_names)
            else:
                rna_out_streams=[pysam.AlignmentFile(rnadna_out_prfx+name+".rna.bam","wb", template=inr) for name in self.group_names]

            if to_stdout=="dna":
                dna_out_streams=[pysam.AlignmentFile("/dev/stdout","w", template=ind)]*len(self.group_names)
            else:
                dna_out_streams=[pysam.AlignmentFile(rnadna_out_prfx+name+".dna.bam","wb", template=ind) for name in self.group_names]
        
            n=self._fanout_reads_bams_stream(ambiviter, rnabamiter, dnabamiter, readid_out_streams, ft_out_streams, rna_out_streams, dna_out_streams, nmax=0)

        # close everything
        if not readid_out_prfx is None:
            for stream in readid_out_streams:
                stream.close()
        
        if not pairs_out_prfx is None:
            for stream in ft_out_streams:
                stream.close()

        ambiviter.close()

        if not rna_bam_in is None:
            inr.close()
            if not (to_stdout=="rna"):
                for stream in rna_out_streams:
                    stream.close()
        
        if not dna_bam_in is None:
            ind.close()
            if not (to_stdout=="dna"):
                for stream in dna_out_streams:
                    stream.close()

        if indexbams:
            
            if (not(to_stdout=="rna")) and (not rna_bam_in is None):
                print('Sorting RNA BAM files')
                for name in self.group_names:
                    rna_file_in=rnadna_out_prfx+name+".rna.bam"
                    rna_file_out=rnadna_out_prfx+name+".rna.sorted.bam"
                    pysam.sort("-o", rna_file_out, rna_file_in, catch_stdout=False)
                    pysam.index(rna_file_out, catch_stdout=False)

            if (not(to_stdout=="dna")) and (not dna_bam_in is None):
                print('Sorting DNA BAM files')
                for name in self.group_names:
                    dna_file_in=rnadna_out_prfx+name+".dna.bam"
                    dna_file_out=rnadna_out_prfx+name+".dna.sorted.bam"
                    pysam.sort("-o", dna_file_out, dna_file_in, catch_stdout=False)
                    pysam.index(dna_file_out, catch_stdout=False)



        return n
# def make_tx_list_list(s, byname=False, name_to_ENSG=None):
#     pass

# def 
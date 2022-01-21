import os
import sys, errno
import pysam
import itertools
from .bamutils import fetch_upto_next_ID, _lt_natural, nextAligned, fetch_upto_next_ID_aligned


def bam2pairs_record(rR,rD,nr,nd,invertFlag,triangular_mode):

    is_reverseD=rD.is_reverse^invertFlag
    is_reverseR=rR.is_reverse^invertFlag
    # s+=" "+"-"*x+ "+"*(not(x))
    # mapq, gap2bridge, AS, XS, secondary_flag
    if is_reverseR:
        posR=rR.reference_start #5'end of RNA
        sR="-"

    else:
        posR=rR.reference_end-1 #5'end of DNA
        sR="+"
        # lenR=str(-rR.reference_length)
    if is_reverseD:
        posD=rD.reference_end-1 #5'end of RNA
        sD="-"
    else:
        posD=rD.reference_start #5'end of DNA
        sD="+"

    gap2bridgeD= str((rD.query_length-rD.query_alignment_end) if invertFlag else rD.query_alignment_start)
    gap2bridgeR= str(rR.query_alignment_start if invertFlag else (rR.query_length-rR.query_alignment_end))

    
    annot_ref=(rR.get_tag('ar') if rR.has_tag('ar') else "*")
    annot_pos=str(rR.get_tag('ap') if rR.has_tag('ap') else 0)
    annot_name=(rR.get_tag('an') if rR.has_tag('an') else "*")
    annot_type=(rR.get_tag('at') if rR.has_tag('at') else "*")
    # annot_strand=(rR.get_tag('as') if rR.has_tag('as') else "*")
    # annot_length=(str(rR.get_tag('al')) if rR.has_tag('al') else "0")

    
    nR=str(rR.get_tag('NH') if rR.has_tag('NH') else nr)
    nD=str(rD.get_tag('NH') if rD.has_tag('NH') else nd)

    annot_gS=str(rR.get_tag('gS') if rR.has_tag('aS') else 0)
    annot_aS=str(rR.get_tag('aS') if rR.has_tag('aS') else 0)
    annot_ai=str(rR.get_tag('ai') if rR.has_tag('ai') else 0)
    annot_aI=str(rR.get_tag('aI') if rR.has_tag('aI') else 0)
    annot_aG=str(rR.get_tag('ag') if rR.has_tag('ag') else "*")

    # is_lowtrig=_lt_chrpos(rD.reference_name,rR.reference_name,int(posD),int(posR))
    # TODO : CATCH bad reference ids when unknown target
    is_lowtrig=_lt_chrpos(rD.reference_id,rR.reference_id,int(posD),int(posR))
    # refD=rD.reference_name
    refR=(triangular_mode==1)*"R_"+rR.reference_name #add a prefix to the RNA which makes it always be smaller than the DNA
    s=""
    if ((triangular_mode==2) & is_lowtrig): #reverse the R and D, do not add prefix
        s="\t".join([
            rD.query_name,rD.reference_name,str(posD),rR.reference_name,str(posR),sD,sR, 
            str(rR.mapping_quality),str(rD.mapping_quality), str(rR.flag), str(rD.flag), gap2bridgeR,gap2bridgeD, str(rR.query_alignment_length), 
            str(rD.query_alignment_length), nR, nD, annot_ref, annot_pos, annot_name, annot_type, annot_gS, annot_aS, annot_ai, annot_aI, annot_aG])
    else:
        #NEW OUTPUT 0=readID, 1=chrR, 2=posR, 3=chrD, 4=posD, 5=strandR, 6=strandD, 
        # 7=QR, 8=QD, 9=flagR, 10=flagD, 11=gapR, 12=gapD, 13=query alignment length RNA, 14=query aligment length DNA, 
        # 15=nR, 16=nD, 17=annot_ref(ENST), 18=annot_pos, 19=annot_name(ex: gapdh), 20=annot_type, 21=gS (#of genomic aligmnet with compatible annotation RNA side), 22=aS (#of annotations), 23=ai, 24=aI, 25=ENSG 
        s="\t".join([
            rD.query_name,refR,str(posR),rD.reference_name,str(posD),sR,sD,
            str(rR.mapping_quality),str(rD.mapping_quality), str(rR.flag), str(rD.flag), gap2bridgeR,gap2bridgeD, str(rR.query_alignment_length), 
            str(rD.query_alignment_length), nR, nD, annot_ref, annot_pos, annot_name, annot_type, annot_gS, annot_aS, annot_ai, annot_aI, annot_aG])

       
    return s, is_lowtrig

# def annotate_read_RNA(rR,invertFlag):
#     is_reverseR=rR.is_reverse^invertFlag
#     posR=rR.reference_start if is_reverseR else (rR.reference_end-1) #5'end of RNA
#     gap2bridgeR= rR.query_alignment_start if invertFlag else (rR.query_length-rR.query_alignment_end)
#     alg_len=rR.reference_length
#     rR.set_tag('ct','R',value_type='A')
#     rR.set_tag('cp',posR,value_type='i')
#     rR.set_tag('cg',gap2bridgeR,value_type='i')
#     rR.set_tag('ca',alg_len,value_type='i')
#     rR.set_tag('co',int(invertFlag),value_type='i')

# def annotate_read_DNA(rD,invertFlag):
#     is_reverseD=rD.is_reverse^invertFlag
#     posD=(rD.reference_end-1) if is_reverseD else rD.reference_start #5'end of RNA
#     gap2bridgeD= (rD.query_length-rD.query_alignment_end) if invertFlag else rD.query_alignment_start
#     alg_len=rD.reference_length
#     rD.set_tag('ct','D',value_type='A')
#     rD.set_tag('cp',posD,value_type='i')
#     rD.set_tag('cg',gap2bridgeD,value_type='i')
#     rD.set_tag('ca',alg_len,value_type='i')
#     rD.set_tag('co',int(invertFlag),value_type='i')

# def annotate_read_pair(rR,rD,iR,iD,nR,nD,invertFlag):
#     is_reverseD=rD.is_reverse^invertFlag
#     posD=(rD.reference_end-1) if is_reverseD else rD.reference_start #5'end of RNA
#     gap2bridgeD= (rD.query_length-rD.query_alignment_end) if invertFlag else rD.query_alignment_start
#     alg_len=rD.reference_length
#     rD.set_tag('dp',posD,value_type='i')
#     rD.set_tag('dg',gap2bridgeD,value_type='i')
#     rD.set_tag('da',alg_len,value_type='i')
#     rD.set_tag('do',int(invertFlag),value_type='i')
#     rD.set_tag('dn',int(nD),value_type='i')  
#     rR.set_tag('rp',posD,value_type='i')
#     rR.set_tag('rg',gap2bridgeD,value_type='i')
#     rR.set_tag('ra',alg_len,value_type='i')
#     rR.set_tag('ro',int(invertFlag),value_type='i')
#     rR.set_tag('dn',int(nD),value_type='i')

def bampipe_RNA(bamiterR,out,invertFlag):
    nR=0
    if not out is None:
        try:
            while True:
                rR = nextAligned(bamiterR)
                nR+=1
                # annotate_read_RNA(rR,invertFlag)
                out.write(rR)
        except StopIteration: # no RNA, move everything to DNA
            return nR
    else:
        try:
            while True:
                rR = nextAligned(bamiterR)
                nR+=1
                # annotate_read_RNA(rR,invertFlag)
        except StopIteration: # no RNA, move everything to DNA
            return nR

def bampipe_DNA(bamiterD,out,invertFlag):
    nD=0
    if not out is None:
        try:
            while True:
                rD = nextAligned(bamiterD)
                nD+=1
                # annotate_read_DNA(rD,invertFlag)
                out.write(rD)
        except StopIteration: # no RNA, move everything to DNA
            return nD
    else:
        try:
            while True:
                rD = nextAligned(bamiterD)
                nD+=1
                # annotate_read_DNA(rD,invertFlag)
        except StopIteration: # no RNA, move everything to DNA
            return nD

# def _lt_chrpos(x,y,xpos,ypos):
#     if _lt_natural(x,y):
#         return True
#     elif x==y:
#         (return True) if x<y else (return False)
#     else:
#         return False

def _lt_chrpos(x,y,xpos,ypos):
    return ( (x<y) | ((x==y) & (xpos<ypos)) )

def bams2pairs_stream(bamiterR,bamiterD,OUT_pair,OUT_pairedRNA, OUT_pairedDNA, OUT_r, OUT_d, invert=False, triangular_mode=0, OUT_pair_lowtrig=None, naturalSort=False, reduceMultimapper=1, filter_with_annots=False): #0, keep everything intact, 1=enforce upper triangular by prefixng the rna with R_, 2=enforce upper triangular by splitting into to two files (file 2 has RNA and DNA exchanged and therefore stays upper triang )
    # IMPORTANT: the streams need to generate records sorted by readIDs
    npairs=0
    npairs_lowtrig=0
    nD=0
    nR=0
    try:
        rR = nextAligned(bamiterR)
    except StopIteration: # no RNA, move everything to DNA and done
        nD=bampipe_DNA(bamiterD,OUT_d,invert)
        return (npairs,nR,nD,npairs_lowtrig)

    try:
        rD = nextAligned(bamiterD)
    except StopIteration: # no DNA, move everything to RNA and done
        nR=bampipe_RNA(bamiterR,OUT_r,invert)
        return (npairs,nR,nD,npairs_lowtrig)

    while not((rD==None) & (rR==None)):
        if rD==None: # no more DNA data, finish writing everything to RNA file and done
            # annotate_read_RNA(rR, invert) #we don't annotate these reads anymore, while waste time
            if not OUT_r is None:
                OUT_r.write(rR)
            nRpipe=bampipe_RNA(bamiterR,OUT_r,invert)
            nR+=(1+nRpipe)
            return (npairs,nR,nD,npairs_lowtrig)
        elif rR==None: # no more RNA data, finish writing everything to DNA file and done
            # annotate_read_DNA(rD, invert)
            if not OUT_d is None:
                OUT_d.write(rD)
            nDpipe=bampipe_DNA(bamiterD,OUT_d,invert)
            nD+=(1+nDpipe)
            return (npairs,nR,nD,npairs_lowtrig)
        else:
            while not(rR.query_name==rD.query_name): # no match
                if  _lt_natural(rD.query_name,rR.query_name, naturalSort): #no match and the DNA stream is behind
                    # annotate_read_DNA(rD, invert)
                    if not OUT_d is None:
                        OUT_d.write(rD)
                    nD+=1
                    try:
                        rD=nextAligned(bamiterD)
                    except StopIteration: # no more DNA data, finish writing everything to RNA file and done
                        # annotate_read_RNA(rR,invert)
                        if not OUT_r is None:
                            OUT_r.write(rR)
                        nRpipe=bampipe_RNA(bamiterR,OUT_r,invert)
                        nR+=(1+nRpipe)
                        return (npairs,nR,nD,npairs_lowtrig)

                else: # no match, and the RNA stream is behind
                    #annotate_read_RNA(rR, invert)
                    if not OUT_r is None:
                        OUT_r.write(rR)
                    nR+=1
                    try:
                        rR=nextAligned(bamiterR)
                    except StopIteration: # no more RNA data, finish writing everything to RNA file and done
                        #annotate_read_DNA(rD,invert)
                        if not OUT_d is None:
                            OUT_d.write(rD)
                        nDpipe=bampipe_RNA(bamiterD,OUT_d,invert)
                        nD+=(1+nDpipe)
                        return (npairs,nR,nD,npairs_lowtrig)

            # at this point, the R and D streams are pointing to unanalzyed records with matching IDs, so we gather them along with possible duplicates

            recordsD, rD = fetch_upto_next_ID_aligned(bamiterD,rD)
            recordsR, rR = fetch_upto_next_ID_aligned(bamiterR,rR)
#             n=bam2pairs_record_list(recordsR,recordsD,OUT_pair,invert)
            
            if not(filter_with_annots and recordsR[0].has_tag('aS') and (recordsR[0].get_tag('aS')>0)): 
                npairs+=1


                # no more extra bam files
                # if not(OUT_fastq_D==None):
                #     OUT_fastq_D.write("@%s\n" % recordsD[0].query_name)
                #     OUT_fastq_D.write("%s\n+\n" % recordsD[0].get_forward_sequence())
                #     xx=recordsD[0].get_forward_qualities()
                #     for i in range(len(xx)):
                #         xx[i]+=33
                #     OUT_fastq_D.write("%s\n" % xx.tostring().decode('utf-8'))

                #     OUT_fastq_R.write("@%s\n" % recordsR[0].query_name)
                #     OUT_fastq_R.write("%s\n+\n" % recordsR[0].get_forward_sequence())
                #     xx=recordsR[0].get_forward_qualities()
                #     for i in range(len(xx)):
                #         xx[i]+=33
                #     OUT_fastq_R.write("%s\n" % xx.tostring().decode('utf-8'))

                n_recR=len(recordsR)
                n_recD=len(recordsD)
                if reduceMultimapper==2:
                    recordsR=[recordsR[0]]
                    recordsD=[recordsD[0]]
                if not OUT_pairedRNA is None:
                    for recR in recordsR:
                        # annotate_read_RNA(recR, invert) # TO DO: keep track of duplication here for a duplication stats
                        OUT_pairedRNA.write(recR)

                if not OUT_pairedDNA is None:
                    for recD in recordsD:
                        # annotate_read_DNA(recD, invert)
                        OUT_pairedDNA.write(recD)
                
                if reduceMultimapper==1: # the RNA and DNA bam file are not squished, but the pairs file is, meaning that only one of the mapping option is given
                    recordsR=[recordsR[0]]
                    recordsD=[recordsD[0]]
                #TO DO: offer other modes of operation where we don't take combinatorial matches
                if triangular_mode==2:
                    for (_, elem) in enumerate(itertools.product(recordsR,recordsD)):
                        s, is_lowtrig = bam2pairs_record(elem[0],elem[1],n_recR,n_recD,invert,triangular_mode)
                        if is_lowtrig:
                            npairs_lowtrig+=1
                            OUT_pair_lowtrig.write(s)
                            OUT_pair_lowtrig.write("\n")
                        else:
                            OUT_pair.write(s)
                            OUT_pair.write("\n")
                else:
                    for (_, elem) in enumerate(itertools.product(recordsR,recordsD)):
                        s, is_lowtrig = bam2pairs_record(elem[0],elem[1],n_recR,n_recD,invert,triangular_mode)
                        OUT_pair.write(s)
                        OUT_pair.write("\n")
                        if is_lowtrig:
                            npairs_lowtrig+=1
    return npairs,nR,nD,npairs_lowtrig


def bams2pairs(infile_rna,infile_dna,outfile_pairs=None, outfile_pairedRNA=None, outfile_pairedDNA=None, outfile_unpairedRNA=None,outfile_unpairedDNA=None,pairingmode="prefix", outfile_pairs_lowtrig=None, invert=False, naturalSort=False, outfile_stats=None, reduceMultimapper=1, filter_with_annots=False):
    triangular_mode=0
    if pairingmode=="prefix":
        triangular_mode=1
    elif pairingmode=="triangular":
        triangular_mode=2
    elif pairingmode=="rd":
        triangular_mode=0


    if invert:
        print("Running in reverse bridge mode")

    inr = pysam.AlignmentFile(infile_rna)
    ind = pysam.AlignmentFile(infile_dna)
    bamiterR = inr.fetch(until_eof=True)
    bamiterD = ind.fetch(until_eof=True)

    if outfile_pairs is None:
        OUT_pair = sys.stdout
    else:
        OUT_pair = open(outfile_pairs, "w", encoding="utf-8")

    if triangular_mode==2:
        if outfile_pairs_lowtrig is None:
            OUT_pair_lowtrig = sys.stdout
        else:
            OUT_pair_lowtrig= open(outfile_pairs_lowtrig, "w", encoding="utf-8")
    else:
        OUT_pair_lowtrig = None

    if not outfile_pairedRNA is None:
        outmode_pairedRNA="wb" if outfile_pairedRNA.endswith(".bam") else "w"
        OUT_pairedRNA= pysam.AlignmentFile(outfile_pairedRNA, outmode_pairedRNA, template=inr)
    else:
        OUT_pairedRNA=None

    if not outfile_pairedDNA is None:
        outmode_pairedDNA="wb" if outfile_pairedDNA.endswith(".bam") else "w"
        OUT_pairedDNA= pysam.AlignmentFile(outfile_pairedDNA, outmode_pairedDNA, template=ind)
    else:
        OUT_pairedDNA=None


    if outfile_unpairedRNA is None:
        OUT_r=None
    else:
        outmode_r="wb" if outfile_unpairedRNA.endswith(".bam") else "w"
        OUT_r = pysam.AlignmentFile(outfile_unpairedRNA,outmode_r, template=inr) #add header

    if outfile_unpairedDNA is None:
        OUT_d=None
    else:
        outmode_d="wb" if outfile_unpairedDNA.endswith(".bam") else "w"
        OUT_d = pysam.AlignmentFile(outfile_unpairedDNA,outmode_d, template=ind) #add header
        
    bpipe=False
    # try:
    n=bams2pairs_stream(bamiterR,bamiterD,OUT_pair,OUT_pairedRNA, OUT_pairedDNA, OUT_r, OUT_d, invert=invert, triangular_mode=triangular_mode, OUT_pair_lowtrig=OUT_pair_lowtrig, naturalSort=naturalSort, reduceMultimapper=reduceMultimapper, filter_with_annots=filter_with_annots)
    # except IOError as e:
    #     print("here")
    #     print(e.errno)
    #     if e.errno == errno.EPIPE:
    #         n=[0,0,0,0]
    #         bpipe=True
    #     else:
    #         # pass
    #         raise(e)

    print("Closing streams")
    # write stats
    if (not outfile_stats is None) and (not bpipe):
        OUT_stats=open(outfile_stats, "w")
    
        OUT_stats.write("n,nRNAsOnly,nDNAsOnly,nlowtrig\n")
        OUT_stats.write(",".join([str(n[0]),str(n[1]),str(n[2]),str(n[3])]))
        OUT_stats.write("\n")

        OUT_stats.flush()
        OUT_stats.close()
  
    #close all streams
    if not bpipe:
        OUT_pair.flush()
    
    if not OUT_pair is None and (not bpipe):
        OUT_pair.close()
    if triangular_mode==2:
        if not bpipe:
            OUT_pair_lowtrig.flush()
        if (not outfile_pairs_lowtrig is None) and (not bpipe):
            OUT_pair_lowtrig.close()
    
    if not outfile_pairedRNA is None:
        # OUT_pairedRNA.flush()
        OUT_pairedRNA.close()
    if not outfile_pairedDNA is None:
        # OUT_pairedDNA.flush()
        OUT_pairedDNA.close()

    if not outfile_unpairedRNA is None:
        # OUT_r.flush()
        OUT_r.close()
    if not outfile_unpairedDNA is None:
        # OUT_d.flush()
        OUT_d.close()
  
    
    return n

# def triangularize(pairsfile, ltfile, utfile):
#     reader=open(pairsfile,'r')
#     try:
#         line=reader.readline()
#     except StopIteration:
#         break
    
#     if len(line)<2:
#         raise StopIteration
#     else:
#         read_data=line.strip().split("\t")
#         rR=read_data[1][2:]
#         rD=read_data[3]
#         if 

#             # read_data=read_data[0,17,self.agtag]
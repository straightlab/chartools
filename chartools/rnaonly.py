import os
import sys, errno
import pysam
import itertools
from .bamutils import fetch_upto_next_ID, _lt_natural, nextAligned, fetch_upto_next_ID_aligned


def bam2single_record(rR, nr, invertFlag):

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


    
    annot_ref=(rR.get_tag('ar') if rR.has_tag('ar') else "*")
    annot_pos=str(rR.get_tag('ap') if rR.has_tag('ap') else 0)
    annot_name=(rR.get_tag('an') if rR.has_tag('an') else "*")
    annot_type=(rR.get_tag('at') if rR.has_tag('at') else "*")
    # annot_strand=(rR.get_tag('as') if rR.has_tag('as') else "*")
    # annot_length=(str(rR.get_tag('al')) if rR.has_tag('al') else "0")

    
    nR=str(rR.get_tag('NH') if rR.has_tag('NH') else nr)

    annot_gS=str(rR.get_tag('gS') if rR.has_tag('aS') else 0)
    annot_aS=str(rR.get_tag('aS') if rR.has_tag('aS') else 0)
    annot_ai=str(rR.get_tag('ai') if rR.has_tag('ai') else 0)
    annot_aI=str(rR.get_tag('aI') if rR.has_tag('aI') else 0)
    annot_aG=str(rR.get_tag('ag') if rR.has_tag('ag') else "*")
    #1=readID, 2=chrR, 3=posR, 4=strandR, 5=QR, 6=flagR, 7=query alignment length RNA, 8=nR, 9=annot_ref(ENST), 10=annot_pos, 11=annot_name(ex: gapdh), 12=annot_type, 13=gS (#of genomic aligmnet with compatible annotation RNA side), 14=aS (#of annotations), 15=ai, 16=aI, 17=ENSG 

    s="\t".join([rR.query_name, rR.reference_name,str(posR),sR, str(rR.mapping_quality), str(rR.flag), str(rR.query_alignment_length), nR, annot_ref, annot_pos, annot_name, annot_type, annot_gS, annot_aS, annot_ai, annot_aI, annot_aG])
    
    return s


def _lt_chrpos(x,y,xpos,ypos):
    return ( (x<y) | ((x==y) & (xpos<ypos)) )

def bam2single_stream(bamiterR, OUT_single, invert=False, naturalSort=False, reduceMultimapper=1, filter_with_annots=False): #0, keep everything intact, 1=enforce upper triangular by prefixng the rna with R_, 2=enforce upper triangular by splitting into to two files (file 2 has RNA and DNA exchanged and therefore stays upper triang )
    # IMPORTANT: the streams need to generate records sorted by readIDs
    npairs=0
    try:
        rR = nextAligned(bamiterR)
    except StopIteration: # no RNA, move everything to DNA and done
        return npairs


    while not(rR==None):
        

        recordsR, rR = fetch_upto_next_ID_aligned(bamiterR,rR)
#             n=bam2pairs_record_list(recordsR,recordsD,OUT_pair,invert)
            
        if not(filter_with_annots and recordsR[0].has_tag('aS') and (recordsR[0].get_tag('aS')>0)): 
            npairs+=1

            n_recR=len(recordsR)
                
            if reduceMultimapper==2:
                recordsR=[recordsR[0]]

            
            if reduceMultimapper==1: # the RNA and DNA bam file are not squished, but the pairs file is, meaning that only one of the mapping option is given
                recordsR=[recordsR[0]]
        
            for (_, elem) in enumerate(recordsR):
                s = bam2single_record(elem,n_recR,invert)
                OUT_single.write(s)
                OUT_single.write("\n")
    
    return npairs


def bam2single(infile_rna,outfile_single=None, invert=False, naturalSort=False, outfile_stats=None, reduceMultimapper=1, filter_with_annots=False):

    if invert:
        print("Running in reverse bridge mode")

    inr = pysam.AlignmentFile(infile_rna)
    bamiterR = inr.fetch(until_eof=True)

    if outfile_single is None:
        OUT_single = sys.stdout
    else:
        OUT_single = open(outfile_single, "w", encoding="utf-8")

    bpipe=False
    # try:
    n=bam2single_stream(bamiterR,OUT_single, invert=invert, naturalSort=naturalSort, reduceMultimapper=reduceMultimapper, filter_with_annots=filter_with_annots)
    
    print("Closing streams")
    # write stats
  
    #close all streams
    if not bpipe:
        OUT_single.flush()
    
    if not OUT_single is None and (not bpipe):
        OUT_single.close()
    
    return n

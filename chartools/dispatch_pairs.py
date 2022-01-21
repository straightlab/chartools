import os
import sys, errno
import pysam


def cistrans_stream(bamiterR,bamiterD,OUT_pairedRNA_CIS, OUT_pairedDNA_CIS,OUT_pairedRNA_TRANS, OUT_pairedDNA_TRANS, OUT_pairedRNA_ambig, OUT_pairedDNA_ambig, qminRNA=0, qminDNA=0, nmax=0):
    # IMPORTANT: the streams need to generate records sorted by readIDs
    isRNA=True
    isDNA=True
    n=0
    try:
        rR = bamiterR.__next__()
    except StopIteration: # no RNA, we're done
        isRNA=False
    try:
        rD = bamiterD.__next__()
    except StopIteration: # no DNA, move everything to RNA and done
        isDNA=False

    assert (isRNA - isDNA)==0, 'Discordant number of RNA and DNA reads'
    
    while (isRNA and isDNA and ((nmax==0) or (n<nmax))):
        assert rR.query_name==rD.query_name, 'Discordant RNA and DNA read IDs'
        n+=1
        if ((rR.mapping_quality>=qminRNA) and (rD.mapping_quality>=qminDNA)):
            if rR.reference_id==rD.reference_id:
                OUT_pairedRNA_CIS.write(rR)
                OUT_pairedDNA_CIS.write(rD)
            else:
                OUT_pairedRNA_TRANS.write(rR)
                OUT_pairedDNA_TRANS.write(rD)
        else:
            OUT_pairedRNA_ambig.write(rR)
            OUT_pairedDNA_ambig.write(rD)

        try:
            rR = bamiterR.__next__()
        except StopIteration: # no RNA, we're done
            isRNA=False
        try:
            rD = bamiterD.__next__()
        except StopIteration: # no DNA, move everything to RNA and done
            isDNA=False
        assert (isRNA - isDNA)==0, 'Discordant number of RNA and DNA reads'

def cistrans(infile_rna,infile_dna, rcis=None, dcis=None, rtrans=None, dtrans=None, rambig=None, dambig=None, qminRNA=0, qminDNA=0, nmax=0):

    inr = pysam.AlignmentFile(infile_rna)
    ind = pysam.AlignmentFile(infile_dna)
    bamiterR = inr.fetch(until_eof=True)
    bamiterD = ind.fetch(until_eof=True)

    OUT_pairedRNA_CIS= pysam.AlignmentFile(rcis, "wb", template=inr)
    OUT_pairedDNA_CIS= pysam.AlignmentFile(dcis, "wb", template=ind)

    OUT_pairedRNA_TRANS= pysam.AlignmentFile(rtrans, "wb", template=inr)
    OUT_pairedDNA_TRANS= pysam.AlignmentFile(dtrans, "wb", template=ind)

    OUT_pairedRNA_ambig= pysam.AlignmentFile(rambig, "wb", template=inr)
    OUT_pairedDNA_ambig= pysam.AlignmentFile(dambig, "wb", template=ind)

    cistrans_stream(bamiterR,bamiterD,OUT_pairedRNA_CIS, OUT_pairedDNA_CIS,OUT_pairedRNA_TRANS, OUT_pairedDNA_TRANS, OUT_pairedRNA_ambig, OUT_pairedDNA_ambig, qminRNA=qminRNA, qminDNA=qminDNA, nmax=nmax)
   
    
    OUT_pairedRNA_CIS.close()
    OUT_pairedDNA_CIS.close()
    OUT_pairedRNA_TRANS.close()
    OUT_pairedDNA_TRANS.close()
    OUT_pairedRNA_ambig.close()
    OUT_pairedDNA_ambig.close()



def within_bp_stream(bamiterR,bamiterD,OUT_pairedRNA_within, OUT_pairedDNA_within, qminRNA=0, qminDNA=0, mind=-1000, maxd=1000, nmax=0):
    # IMPORTANT: the streams need to generate records sorted by readIDs
    isRNA=True
    isDNA=True
    n=0
    try:
        rR = bamiterR.__next__()
    except StopIteration: # no RNA, we're done
        isRNA=False
    try:
        rD = bamiterD.__next__()
    except StopIteration: # no DNA, move everything to RNA and done
        isDNA=False

    assert (isRNA - isDNA)==0, 'Discordant number of RNA and DNA reads'
    
    while (isRNA and isDNA and ((nmax==0) or (n<nmax))):
        assert rR.query_name==rD.query_name, 'Discordant RNA and DNA read IDs'
        n+=1
        if ((rR.mapping_quality>=qminRNA) and (rD.mapping_quality>=qminDNA)):
            if rR.is_reverse:
                posR=rR.reference_start #5'end of RNA

            else:
                posR=rR.reference_end-1 #5'end of DNA
            
            if rD.is_reverse:
                posD=rD.reference_end-1 #5'end of RNA
            else:
                posD=rD.reference_start #5'end of DNA

            if ((rR.reference_id==rD.reference_id) and ((posD-posR)>mind) and ((posD-posR)<maxd)):
                OUT_pairedRNA_within.write(rR)
                OUT_pairedDNA_within.write(rD)

        try:
            rR = bamiterR.__next__()
        except StopIteration: # no RNA, we're done
            isRNA=False
        try:
            rD = bamiterD.__next__()
        except StopIteration: # no DNA, move everything to RNA and done
            isDNA=False
        assert (isRNA - isDNA)==0, 'Discordant number of RNA and DNA reads'

def within_bp(infile_rna,infile_dna, rwithin, dwithin, qminRNA=0, qminDNA=0, mind=-1000, maxd=1000, nmax=0):

    inr = pysam.AlignmentFile(infile_rna)
    ind = pysam.AlignmentFile(infile_dna)
    bamiterR = inr.fetch(until_eof=True)
    bamiterD = ind.fetch(until_eof=True)

    OUT_pairedRNA_within= pysam.AlignmentFile(rwithin, "wb", template=inr)
    OUT_pairedDNA_within= pysam.AlignmentFile(dwithin, "wb", template=ind)


    within_bp_stream(bamiterR,bamiterD,OUT_pairedRNA_within, OUT_pairedDNA_within, qminRNA=qminRNA, qminDNA=qminDNA,mind=mind, maxd=maxd, nmax=nmax)
   
    
    OUT_pairedRNA_within.close()
    OUT_pairedDNA_within.close()






def deambig_stream(bamiterR,bamiterD,OUT_pairedRNA_DNAmultimap, OUT_pairedDNA_DNAmultimap, qminRNA=0, nmax=0):
    # IMPORTANT: the streams need to generate records sorted by readIDs
    isRNA=True
    isDNA=True
    n=0
    try:
        rR = bamiterR.__next__()
    except StopIteration: # no RNA, we're done
        isRNA=False
    try:
        rD = bamiterD.__next__()
    except StopIteration: # no DNA, move everything to RNA and done
        isDNA=False

    assert (isRNA - isDNA)==0, 'Discordant number of RNA and DNA reads'
    
    while (isRNA and isDNA and ((nmax==0) or (n<nmax))):
        assert rR.query_name==rD.query_name, 'Discordant RNA and DNA read IDs'
        n+=1
        if ((rR.mapping_quality>=qminRNA)):
            OUT_pairedRNA_DNAmultimap.write(rR)
            OUT_pairedDNA_DNAmultimap.write(rD)


        try:
            rR = bamiterR.__next__()
        except StopIteration: # no RNA, we're done
            isRNA=False
        try:
            rD = bamiterD.__next__()
        except StopIteration: # no DNA, move everything to RNA and done
            isDNA=False
        assert (isRNA - isDNA)==0, 'Discordant number of RNA and DNA reads'


def deambig(infile_rna,infile_dna, rout=None, dout=None, qminRNA=255,nmax=0):

    inr = pysam.AlignmentFile(infile_rna)
    ind = pysam.AlignmentFile(infile_dna)
    bamiterR = inr.fetch(until_eof=True)
    bamiterD = ind.fetch(until_eof=True)

    OUT_pairedRNA_DNAmultimap= pysam.AlignmentFile(rout, "wb", template=inr)
    OUT_pairedDNA_DNAmultimap= pysam.AlignmentFile(dout, "wb", template=ind)

    

    deambig_stream(bamiterR,bamiterD,OUT_pairedRNA_DNAmultimap, OUT_pairedDNA_DNAmultimap, qminRNA=qminRNA,nmax=nmax)
   
    
    OUT_pairedRNA_DNAmultimap.close()
    OUT_pairedDNA_DNAmultimap.close()
    
def bamsToPairsSimple_stream(bamiterR,bamiterD, out_pairs, qminRNA=0, qminDNA=0, nmax=0, chr_to_keep=None):
    # IMPORTANT: the streams need to generate records sorted by readIDs
    isRNA=True
    isDNA=True
    n=0
    non_spliced=[-1]
    tab="\t"
    newline="\n"
    ensu="ENSU"
    try:
        rR = bamiterR.__next__()
    except StopIteration: # no RNA, we're done
        isRNA=False
    try:
        rD = bamiterD.__next__()
    except StopIteration: # no DNA, move everything to RNA and done
        isDNA=False

    assert (isRNA - isDNA)==0, 'Discordant number of RNA and DNA reads'
    
    print(chr_to_keep)
    while (isRNA and isDNA and ((nmax==0) or (n<nmax))):
        assert rR.query_name==rD.query_name, 'Discordant RNA and DNA read IDs'
        n+=1
        if ((rR.mapping_quality>=qminRNA) and (rD.mapping_quality>=qminDNA) and (rR.reference_name in chr_to_keep) and (rD.reference_name in chr_to_keep)):
        
            posR = rR.reference_start if rR.is_reverse else (rR.reference_end-1)
            posD = (rD.reference_end-1) if rD.is_reverse else (rD.reference_start)
            
            cisVsTrans=str(int(rR.reference_id==rD.reference_id))
            flight=str(posD-posR)

            splice_tag = rR.get_tag('jM') if rR.has_tag('jM') else non_spliced
            isspliced = 0 if ((len(splice_tag)==1) and (splice_tag[0]==-1)) else 1

            exon_tag=rR.get_tag('ar') if rR.has_tag('ar') else "*"
            if exon_tag=="*":
                exonType='unannotated'
            else:
                exonType='exon' if exon_tag.startswith('ENST') else 'intron'

            gene_name=rR.get_tag('ag') if rR.has_tag('ag') else "*"

            
            annot_aI=str(rR.get_tag('aI') if rR.has_tag('aI') else 1)

            out_pairs.write(rR.reference_name)
            out_pairs.write(tab)
            out_pairs.write(str(posR))
            out_pairs.write(tab)
            out_pairs.write(rD.reference_name)
            out_pairs.write(tab)
            out_pairs.write(str(posD))
            out_pairs.write(tab)
            out_pairs.write(gene_name)
            out_pairs.write(tab)
            out_pairs.write(annot_aI)
            out_pairs.write(tab)
            out_pairs.write(exonType)
            out_pairs.write(tab)
            out_pairs.write(str(isspliced))
            out_pairs.write(tab)
            out_pairs.write(cisVsTrans)
            out_pairs.write(tab)
            out_pairs.write(flight)
            out_pairs.write(newline)

        try:
            rR = bamiterR.__next__()
        except StopIteration: # no RNA, we're done
            isRNA=False
        try:
            rD = bamiterD.__next__()
        except StopIteration: # no DNA, move everything to RNA and done
            isDNA=False
        assert (isRNA - isDNA)==0, 'Discordant number of RNA and DNA reads'

def bamsToPairsSimple(infile_rna,infile_dna, outfile_pairs, chr_file, qminRNA=0, qminDNA=0, nmax=0):

    inr = pysam.AlignmentFile(infile_rna)
    ind = pysam.AlignmentFile(infile_dna)
    bamiterR = inr.fetch(until_eof=True)
    bamiterD = ind.fetch(until_eof=True)

    OUT_pair = open(outfile_pairs, "w", encoding="utf-8")

    with open(chr_file) as f:
        chr_to_keep={k.split("\t")[0] : True for k in f.read().splitlines()}

    
    bamsToPairsSimple_stream(bamiterR,bamiterD, OUT_pair, qminRNA=qminRNA, qminDNA=qminDNA, nmax=nmax, chr_to_keep=chr_to_keep)
   
    OUT_pair.flush()
    OUT_pair.close()


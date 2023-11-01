import sys
# def _tobed_stream_both(pairs_stream,rna_out_stream, dna_out_stream, nmax=0):
    
#     ntot=0
#     dna_fields=[0,8,6,10]+[7,5,9]+[17,18,19,20,24,25]
#     rna_fields=[0,7,5,9]+[8,6,10]+[17,18,19,20,24,25]
#     while ((ntot<nmax) or nmax==0):
#         try:
#             line=pairs_stream.readline()
#         except StopIteration:
#             break
#         ntot+=1
    
        
#         if len(line)<2:
#             break
#         else:
#             read_data=line.strip().split("\t")
#             if len(read_data[1])>2:
#                 read_data[1]=read_data[1][2:]
            
#             #0=readID, 1=chrR, 2=posR, 3=chrD, 4=posD, 5=strandR, 6=strandD, 
#         # 7=QR, 8=QD, 9=flagR, 10=flagD, 11=gapR, 12=gapD, 13=query alignment length RNA, 14=query aligment length DNA, 
#         # 15=nR, 16=nD, 17=annot_ref(ENST), 18=annot_pos, 19=annot_name(ex: gapdh), 20=annot_type, 21=gS (#of genomic aligmnet with compatible annotation RNA side), 22=aS (#of annotations), 23=ai, 24=aI, 25=ENSG
            
#             if read_data[5]=="+":
#                 rna_start=int(read_data[2])
#                 rna_end=rna_start+int(read_data[13])
#             else:
#                 rna_end=int(read_data[2])
#                 rna_start=rna_end-int(read_data[13])

#             if read_data[6]=="+":
#                 dna_start=int(read_data[4])
#                 dna_end=dna_start+int(read_data[14])
#             else:
#                 dna_end=int(read_data[4])
#                 dna_start=dna_end-int(read_data[14])

#             rna_out_stream.write(read_data[1]+"\t")
#             rna_out_stream.write(str(rna_start)+"\t")
#             rna_out_stream.write(str(rna_end)+"\t")

#             rna_out_stream.write("\t".join([read_data[i] for i in rna_fields]))
#             rna_out_stream.write("\n")

#             dna_out_stream.write(read_data[3]+"\t")
#             dna_out_stream.write(str(dna_start)+"\t")
#             dna_out_stream.write(str(dna_end)+"\t")
#             dna_out_stream.write("\t".join([read_data[i] for i in dna_fields]))
#             dna_out_stream.write("\n")
            
def _tobed_stream_rna(pairs_stream,rna_out_stream, nmax=0):
    
    ntot=0
    rna_fields_bed=[0,7,5]
    rna_fields=[17,18,19,20,24,25]
    dna_fields_bed=[8,6]
    while ((ntot<nmax) or nmax==0):
        try:
            line=pairs_stream.readline()
        except StopIteration:
            break
        ntot+=1
    
        
        if len(line)<2:

            break
        else:
            read_data=line.strip().split("\t")
            if len(read_data[1])>2:
                read_data[1]=read_data[1][2:]
            
            #0=readID, 1=chrR, 2=posR, 3=chrD, 4=posD, 5=strandR, 6=strandD, 
        # 7=QR, 8=QD, 9=flagR, 10=flagD, 11=gapR, 12=gapD, 13=query alignment length RNA, 14=query aligment length DNA, 
        # 15=nR, 16=nD, 17=annot_ref(ENST), 18=annot_pos, 19=annot_name(ex: gapdh), 20=annot_type, 21=gS (#of genomic aligmnet with compatible annotation RNA side), 22=aS (#of annotations), 23=ai, 24=aI, 25=ENSG
            
            if read_data[5]=="+":
                rna_start=int(read_data[2])
                rna_end=rna_start+int(read_data[13])
            else:
                rna_end=int(read_data[2])
                rna_start=rna_end-int(read_data[13])

            if read_data[6]=="+":
                dna_start=int(read_data[4])
                dna_end=dna_start+int(read_data[14])
            else:
                dna_end=int(read_data[4])
                dna_start=dna_end-int(read_data[14])


            rna_out_stream.write(read_data[1]+"\t")
            rna_out_stream.write(str(rna_start)+"\t")
            rna_out_stream.write(str(rna_end)+"\t")
    
            rna_out_stream.write("\t".join([read_data[i] for i in rna_fields_bed]))
            rna_out_stream.write("\t")
            rna_out_stream.write(read_data[3]+"\t")
            rna_out_stream.write(str(dna_start)+"\t")
            rna_out_stream.write(str(dna_end)+"\t")
            rna_out_stream.write("\t".join([read_data[i] for i in dna_fields_bed]))
            rna_out_stream.write("\t")
            rna_out_stream.write("\t".join([read_data[i] for i in rna_fields]))
            rna_out_stream.write("\n")

def _tobed_stream_dna(pairs_stream,dna_out_stream, nmax=0):
    
    ntot=0
    # dna_fields=[0,8,6,10]
    rna_fields_bed=[7,5]
    rna_fields=[17,18,19,20,24,25]
    dna_fields_bed=[0,8,6]

    while ((ntot<nmax) or nmax==0):
        try:
            line=pairs_stream.readline()
        except StopIteration:
            break
        ntot+=1
    
        
        if len(line)<2:
            break
        else:
            read_data=line.strip().split("\t")
            
            
            #0=readID, 1=chrR, 2=posR, 3=chrD, 4=posD, 5=strandR, 6=strandD, 
        # 7=QR, 8=QD, 9=flagR, 10=flagD, 11=gapR, 12=gapD, 13=query alignment length RNA, 14=query aligment length DNA, 
        # 15=nR, 16=nD, 17=annot_ref(ENST), 18=annot_pos, 19=annot_name(ex: gapdh), 20=annot_type, 21=gS (#of genomic aligmnet with compatible annotation RNA side), 22=aS (#of annotations), 23=ai, 24=aI, 25=ENSG
            

            if read_data[6]=="+":
                dna_start=int(read_data[4])
                dna_end=dna_start+int(read_data[14])
            else:
                dna_end=int(read_data[4])
                dna_start=dna_end-int(read_data[14])

            if read_data[5]=="+":
                rna_start=int(read_data[2])
                rna_end=rna_start+int(read_data[13])
            else:
                rna_end=int(read_data[2])
                rna_start=rna_end-int(read_data[13])


            dna_out_stream.write(read_data[3]+"\t")
            dna_out_stream.write(str(dna_start)+"\t")
            dna_out_stream.write(str(dna_end)+"\t")

            dna_out_stream.write("\t".join([read_data[i] for i in dna_fields_bed]))
            dna_out_stream.write("\t")
            dna_out_stream.write(read_data[1][2:]+"\t")
            dna_out_stream.write(str(rna_start)+"\t")
            dna_out_stream.write(str(rna_end)+"\t")
            dna_out_stream.write("\t".join([read_data[i] for i in rna_fields_bed]))
            dna_out_stream.write("\t")
            dna_out_stream.write("\t".join([read_data[i] for i in rna_fields]))
            dna_out_stream.write("\n")

        

def tobed(pairsfile, to_stdout=None, nmax=0):
    pairs_stream=open(pairsfile,'r')

    if to_stdout=="rna":
        rna_out_stream=sys.stdout
        _tobed_stream_rna(pairs_stream, rna_out_stream, nmax=nmax)
    
    elif to_stdout=="dna":
        dna_out_stream=sys.stdout
        _tobed_stream_dna(pairs_stream, dna_out_stream, nmax=nmax)

    # else:
    #     if to_stdout=="both":
    #         rna_out_stream=sys.stdout
    #         dna_out_stream=sys.stdout
    #     else:
    #         dna_out_stream=open(outprfx+"dna.bed","w")
    #         rna_out_stream=open(outprfx+"rna.bed","w")

    #     _tobed_stream_both(pairs_stream, rna_out_stream, dna_out_stream, nmax=nmax)
    

    pairs_stream.close()

    # if to_stdout is None:
    #    rna_out_stream.close()
    #    dna_out_stream.close()
    
    
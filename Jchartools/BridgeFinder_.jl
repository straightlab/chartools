using CodecZlib
using BioSequences
using TranscodingStreams
using FuzzySearch


function write_record(stream::IO,record::BioSequences.FASTQ.Record,rng::UnitRange{Int64},rc::Bool)
    # pline=[0x0a,0x2b,0x0a]
    if rc
        revcomplement_record!(record,rng)
    end
    if length(rng)>0 && length(rng)<length(record.sequence)
        seq_range=rng+first(record.sequence)-1
        qual_range=rng+first(record.quality)-1
        id_range=1:(first(record.sequence)-1)
        unsafe_write(stream,pointer(record.data,first(id_range)),length(id_range))
        unsafe_write(stream,pointer(record.data,first(seq_range)),length(seq_range))
        # the code below is slower
        # for i=id_range
        #     write(stream,record.data[i])
        # end
        # for i=seq_range
        #     write(stream,record.data[i])
        # end
        write(stream,0x0a,0x2b,0x0a)
        unsafe_write(stream,pointer(record.data,first(qual_range)),length(qual_range))
        # for i=qual_range
        #     write(stream,record.data[i])
        # end
        write(stream,0x0a)
    elseif length(rng)>0
        write(stream,record.data)
        write(stream,0x0a)
    else
        id_range=1:(first(record.sequence)-1)
        unsafe_write(stream,pointer(record.data,first(id_range)),length(id_range))
        write(stream,0x4e,0x0a,0x2b,0x0a,0x21)
        write(stream,0x0a);
    end
end

function revcomplement_record!(record::BioSequences.FASTQ.Record,rng::UnitRange{Int64})
    seq_range=rng+first(record.sequence)-1
    qual_range=rng+first(record.quality)-1
    for i=seq_range
        if !(record.data[i]==0x4e)
            (record.data[i] & 0x02)>0 ? (record.data[i] ⊻= 0x04) : (record.data[i] ⊻= 0x15)
        end
    end
    reverse!(record.data,seq_range.start,seq_range.stop)
    reverse!(record.data,qual_range.start,qual_range.stop)
end




function suffix_mismatch_(data_short::Vector{UInt8},data_long::Vector{UInt8},rng_short::UnitRange{Int64},rng_long::UnitRange{Int64},maxerr::Int64,maxlu::Int64)
    nerr=0
    nmatch=0
    lu=0
    i=rng_short.start
    imax=rng_short.stop+1
    j=rng_long.stop
    jmin=rng_short.start-1
    while i<imax && j>jmin && (maxerr<0 || nerr<=maxerr) && (maxlu<0 || lu<maxlu)
        lu+=1
        x=data_short[i] | 0x20 #lower case it
        y=data_long[j] | 0x20
        if (x!==0x6e)
            (x & 0x02)>0 ? (x ⊻= 0x04) : (x ⊻= 0x15) # complement
            if x!==y
                nerr+=1
            else
                nmatch+=1
            end
        # else
            # nmatch+=1 # do not match N
        end
        j-=1
        i+=1
    end
    return nerr, nmatch
end

function is_suffix_(read1::BioSequences.FASTQ.Record,read2::BioSequences.FASTQ.Record,maxerr::Int64,minmatch::Int64,maxlu::Int64)
    #RNA1 is long, DNA1 is short
    #RNA2 is short, DNA2 is long
    (nerr,nmatch)=suffix_mismatch_(read1.data,read2.data,read1.sequence,read2.sequence,maxerr,maxlu)
    flag=(nerr<=maxerr && nmatch>=minmatch)
    return flag
end

function revcomp_stream(streamIN::IO,streamOUT::IO, verbose::Bool )
    reader = FASTQ.Reader(streamIN)
    record = FASTQ.Record()
    nread=0::Int64
    verbose_period=100000

    cumtime=0.0
    tic()
    while !eof(reader)
        nread+=1
        if verbose && mod(nread,verbose_period)==0
            cumtime+=toq()
            print(STDERR,"Reverse complemented ", nread/(1000000.0), " million reads [", nread/cumtime/1000000*60, " million reads/min]\n")
            tic()
        end
        read!(reader, record)
        revcomplement_record!(record,1:length(record.sequence))
        write(streamOUT,record.data)
        write(streamOUT,0x0a)
    end
    cumtime+=toq()
    print(STDERR,"DONE. Reverse complemented ", nread/(1000000.0), " million reads [", nread/cumtime/1000000*60, " million reads/min]\n")
    return nread
end

function revcompFQ(fqIN::String,fqOUT::String, verbose::Bool=true )
    if length(fqOUT)==0
        streamOUT=STDOUT
    else
        if endswith(fqOUT, ".gz")
            codec = GzipCompressor()
            streamOUT=TranscodingStream(codec, open(fqOUT,"w"))
        else
            codec = Noop()
            TranscodingStream(codec, open(fqOUT,"w"))
        end
    end
    streamIN=makestream_in(fqIN)
    revcomp_stream(streamIN,streamOUT,verbose)
    !(streamOUT==STDOUT) ? close(streamOUT) : true
    close(streamIN)
end

function readlength_summary_stream_(stream::IO,blen::Int64,maxreadl::Int64, verbose::Bool)
    M_F0=zeros(Int64,maxreadl,maxreadl)
    M_0R=zeros(Int64,maxreadl,maxreadl)
    M_F02=zeros(Int64,maxreadl,maxreadl)
    M_0R2=zeros(Int64,maxreadl,maxreadl)
    lR=0::Int64
    lD=0::Int64
    lR2=0::Int64
    lD2=0::Int64

    nread=0::Int64
    verbose_period=100000

    cumtime=0.0
    tic()
    while !eof(stream)
        nread+=1
        # print(nread,"\n")
        if verbose && mod(nread,verbose_period)==0
            cumtime+=toq()
            print(STDERR,"Processed ", nread/(1000000.0), " million reads [", nread/cumtime/1000000*60, " million reads/min]\n")
            tic()
        end

        s=readline(stream)
        x=parse.(split(s,'\t'))
        # print(x)
        if x[1]==1 && x[2]==2 #OF
            lR=x[6]-1
            lD=x[8]-(x[6]+blen)+1
            lR2=lR
            lD2=lD+x[7]
            M_F0[lR+1,lD+1]+=1
            M_F02[lR2+1,lD2+1]+=1
        elseif x[2]==1 && x[1]==2 #F0
            lR=x[5]-1
            lD=x[7]-(x[5]+blen)+1
            lR2=lR
            lD2=lD+x[8]
            M_F0[lR+1,lD+1]+=1
            M_F02[lR2+1,lD2+1]+=1
        elseif x[1]==3 && x[2]==1 #OR
            lD=x[5]-1
            lR=x[7]-(x[5]+blen)+1
            lD2=lD
            lR2=lR+x[8]
            M_0R[lR+1,lD+1]+=1
            M_0R2[lR2+1,lD2+1]+=1
        elseif x[2]==3 && x[1]==1 #OR
            lD=x[6]-1
            lR=x[8]-(x[6]+blen)+1
            lD2=lD
            lR2=lR+x[7]
            M_0R[lR+1,lD+1]+=1
            M_0R2[lR2+1,lD2+1]+=1
        end
    end
    return M_F02, M_0R2, M_F0, M_0R
end

function readlength_SE_summary_stream_(stream::IO,maxreadl::Int64, verbose::Bool)
    M=zeros(Int64,maxreadl,maxreadl)

    nread=0::Int64
    verbose_period=100000

    cumtime=0.0
    tic()
    while !eof(stream)
        nread+=1
        # print(nread,"\n")
        if verbose && mod(nread,verbose_period)==0
            cumtime+=toq()
            print(STDERR,"Processed ", nread/(1000000.0), " million reads [", nread/cumtime/1000000*60, " million reads/min]\n")
            tic()
        end

        s=readline(stream)
        x=parse.(split(s,' '))
        M[x[1]+1,x[2]+1]+=1
        # print(x)

    end
    return M
end

# function readlength_summary_stream_FAST_(stream::IO,blen::Int64,maxreadl::Int64, verbose::Bool)
#     M_F0=zeros(Int64,maxreadl,maxreadl)
#     M_0R=zeros(Int64,maxreadl,maxreadl)
#     M_F02=zeros(Int64,maxreadl,maxreadl)
#     M_0R2=zeros(Int64,maxreadl,maxreadl)
#     lR=0::Int64
#     lD=0::Int64
#     lR2=0::Int64
#     lD2=0::Int64
#
#     nread=0::Int64
#     verbose_period=100000
#
#     cumtime=0.0
#     tic()
#     while !eof(stream)
#         nread+=1
#         # print(nread,"\n")
#         if verbose && mod(nread,verbose_period)==0
#             cumtime+=toq()
#             print(STDERR,"Processed ", nread/(1000000.0), " million reads [", nread/cumtime/1000000*60, " million reads/min]\n")
#             tic()
#         end
#
#         x=read(stream,UInt8)
#         s=readline(stream)
#         x=parse.(split(s,'\t'))
#         # print(x)
#         if x[1]==1 && x[2]==2 #OF
#             lR=x[6]-1
#             lD=x[8]-(x[6]+blen)+1
#             lR2=lR
#             lD2=lD+x[7]
#             M_F0[lR+1,lD+1]+=1
#             M_F02[lR2+1,lD2+1]+=1
#         elseif x[2]==1 && x[1]==2 #F0
#             lR=x[5]-1
#             lD=x[7]-(x[5]+blen)+1
#             lR2=lR
#             lD2=lD+x[8]
#             M_F0[lR+1,lD+1]+=1
#             M_F02[lR2+1,lD2+1]+=1
#         elseif x[1]==3 && x[2]==1 #OR
#             lD=x[5]-1
#             lR=x[7]-(x[5]+blen)+1
#             lD2=lD
#             lR2=lR+x[8]
#             M_0R[lR+1,lD+1]+=1
#             M_0R2[lR2+1,lD2+1]+=1
#         elseif x[2]==3 && x[1]==1 #OR
#             lD=x[6]-1
#             lR=x[8]-(x[6]+blen)+1
#             lD2=lD
#             lR2=lR+x[7]
#             M_0R[lR+1,lD+1]+=1
#             M_0R2[lR2+1,lD2+1]+=1
#         end
#     end
#     return M_F02, M_0R2, M_F0, M_0R
# end

function readlength_summary(fname::String,fout::String,blen::Int64,maxreadl::Int64,verbose::Bool=true)
    in=makestream_in(fname)
    stream=open(fout,"w")
    M_F02, M_0R2, M_F0, M_0R = readlength_summary_stream_(in,blen,maxreadl,verbose)
    println(stream, "\n>>RNA/DNA length distribution")
    print(stream,"MF02"*"\t")
    Base.showarray(stream,M_F02,true)
    println(stream,"")
    print(stream,"M0R2"*"\t")
    Base.showarray(stream,M_0R2,true)
    println(stream,"")
    print(stream,"MF0"*"\t")
    Base.showarray(stream,M_F0,true)
    println(stream,"")
    print(stream,"M0R"*"\t")
    Base.showarray(stream,M_0R,true)
    println(stream,"")
    close(stream)
    close(in)
    return M_F02, M_0R2, M_F0, M_0R
end

function readlength_SE_summary(fname::String,fout::String,maxreadl::Int64,verbose::Bool=true)
    in=makestream_in(fname)
    stream=open(fout,"w")
    M = readlength_SE_summary_stream_(in,maxreadl,verbose)
    println(stream, "\n>>RNA/DNA length distribution")
    print(stream,"M"*"\t")
    Base.showarray(stream,M,true)
    println(stream,"")
    close(stream)
    close(in)
    return M
end

function splice_mates_stream_(streamR1::IO,streamR2::IO,streamD1::IO,streamD2::IO,stream_splicedR::IO,stream_splicedD::IO,stream_unsplicedR1::IO,stream_unsplicedR2::IO,stream_unsplicedD1::IO,stream_unsplicedD2::IO, err_match_lu_R::Tuple{Int64,Int64,Int64},err_match_lu_D::Tuple{Int64,Int64,Int64}, maxreadl::Int64, checkpairs::Bool, verbose::Bool )
    readerR1 = FASTQ.Reader(streamR1)
    recordR1 = FASTQ.Record()
    readerR2 = FASTQ.Reader(streamR2)
    recordR2 = FASTQ.Record()
    readerD1 = FASTQ.Reader(streamD1)
    recordD1 = FASTQ.Record()
    readerD2 = FASTQ.Reader(streamD2)
    recordD2 = FASTQ.Record()


    rd_len=zeros(Int64,maxreadl,maxreadl)
    rd_len_unspliced1=zeros(Int64,maxreadl,maxreadl)
    rd_len_unspliced2=zeros(Int64,maxreadl,maxreadl)

    nread=0::Int64
    nspliced=0::Int64
    verbose_period=100000

    cumtime=0.0
    tic()
    while !eof(readerR1)
        nread+=1
        # print(nread,"\n")
        if verbose && mod(nread,verbose_period)==0
            cumtime+=toq()
            print(STDERR,"Spliced ", nread/(1000000.0), " million reads [", nread/cumtime/1000000*60, " million reads/min]\n")
            tic()
        end
        read!(readerR1, recordR1)
        read!(readerR2, recordR2)
        read!(readerD1, recordD1)
        read!(readerD2, recordD2)

        if checkpairs
        end

        if (is_suffix_(recordR2,recordR1,err_match_lu_R[1],err_match_lu_R[2],err_match_lu_R[3]) && is_suffix_(recordD1,recordD2,err_match_lu_D[1],err_match_lu_D[2],err_match_lu_D[3]))
            nspliced+=1
            write(stream_splicedR,recordR1.data)
            write(stream_splicedR,0x0a)

            revcomplement_record!(recordD2,1:length(recordD2.sequence))
            write(stream_splicedD,recordD2.data)
            write(stream_splicedD,0x0a)

            rd_len[min(maxreadl+1,length(recordR1.sequence)+1),min(maxreadl+1,length(recordD2.sequence)+1)]+=1
        else
            write(stream_unsplicedR1,recordR1.data)
            write(stream_unsplicedR1,0x0a)
            write(stream_unsplicedR2,recordR2.data)
            write(stream_unsplicedR2,0x0a)
            write(stream_unsplicedD1,recordD1.data)
            write(stream_unsplicedD1,0x0a)
            write(stream_unsplicedD2,recordD2.data)
            write(stream_unsplicedD2,0x0a)
            rd_len_unspliced1[min(maxreadl+1,length(recordR1.sequence)+1),min(maxreadl+1,length(recordD1.sequence)+1)]+=1
            rd_len_unspliced2[min(maxreadl+1,length(recordR1.sequence)+1),min(maxreadl+1,length(recordD1.sequence)+1)]+=1

        end
    end
    return (nspliced,nread,rd_len,rd_len_unspliced1,rd_len_unspliced2)
end

function summarize_splicing_(read_counts, parsable::Bool, stream::IO)
    println(stream, ">>Total counts")
    println(stream,read_counts[2])
    println(stream, ">>Spliced counts")
    println(stream,read_counts[1])
    println(stream, "\n>>RNA/DNA length distribution")
    print(stream,"spliced"*"\t")
    Base.showarray(stream,read_counts[3],parsable)
    println(stream,"")
    print(stream,"unspliced_1"*"\t")
    Base.showarray(stream,read_counts[4],parsable)
    println(stream,"")
    print(stream,"unspliced_2"*"\t")
    Base.showarray(stream,read_counts[5],parsable)
    println(stream,"")

    # println("")
end

function splice_mates(fq_prefix::String,fq_spliced_prefix::String,fq_unspliced_prefix::String,err_match_lu_R_string::String,err_match_lu_D_string::String, maxreadl::Int64, checkpairs::Bool=false, verbose::Bool=true)

    streamSUMMARY=makestream_out(fq_spliced_prefix,".summary",false,".txt")

    streamR1=makestream_in(fq_prefix*".rna.1.fastq.gz")
    streamR2=makestream_in(fq_prefix*".rna.2.fastq.gz")
    streamD1=makestream_in(fq_prefix*".dna.1.fastq.gz")
    streamD2=makestream_in(fq_prefix*".dna.2.fastq.gz")
    stream_splicedR=makestream_out(fq_spliced_prefix,".rna",true,".fastq")
    stream_splicedD=makestream_out(fq_spliced_prefix,".dna",true,".fastq")
    stream_unsplicedR1=makestream_out(fq_unspliced_prefix,".rna.1",true,".fastq")
    stream_unsplicedR2=makestream_out(fq_unspliced_prefix,".rna.2",true,".fastq")
    stream_unsplicedD1=makestream_out(fq_unspliced_prefix,".dna.1",true,".fastq")
    stream_unsplicedD2=makestream_out(fq_unspliced_prefix,".dna.2",true,".fastq")

    err_match_lu_R=tuple(map(x->parse(Int64,x),split(err_match_lu_R_string,":"))...)
    err_match_lu_D=tuple(map(x->parse(Int64,x),split(err_match_lu_D_string,":"))...)

    read_counts = splice_mates_stream_(streamR1,streamR2,streamD1,streamD2,stream_splicedR,stream_splicedD,stream_unsplicedR1,stream_unsplicedR2,stream_unsplicedD1,stream_unsplicedD2,err_match_lu_R,err_match_lu_D,maxreadl,checkpairs, verbose)

    # println(streamSUMMARY, "total,spliced")
    # print(streamSUMMARY,nread)
    # print(streamSUMMARY,",")
    # println(streamSUMMARY, nspliced)

    summarize_splicing_(read_counts,true,streamSUMMARY)

    close(streamSUMMARY)
    close(stream_splicedR)
    close(stream_splicedD)
    close(stream_unsplicedR1)
    close(stream_unsplicedR2)
    close(stream_unsplicedD1)
    close(stream_unsplicedD2)
    close(streamR1)
    close(streamR2)
    close(streamD1)
    close(streamD2)

    return read_counts

end



function write_record_skipempty(stream::IO,record::BioSequences.FASTQ.Record,rng::UnitRange{Int64},rc::Bool)
    # pline=[0x0a,0x2b,0x0a]
    if length(rng)>0 && length(rng)<length(record.sequence)
        seq_range=rng+first(record.sequence)-1
        qual_range=rng+first(record.quality)-1
        id_range=1:(first(record.sequence)-1)
        if rc
            revcomplement_record!(record,rng)
        end
        unsafe_write(stream,pointer(record.data,first(id_range)),length(id_range))
        unsafe_write(stream,pointer(record.data,first(seq_range)),length(seq_range))
        # the code below is slower
        # for i=id_range
        #     write(stream,record.data[i])
        # end
        # for i=seq_range
        #     write(stream,record.data[i])
        # end
        write(stream,0x0a,0x2b,0x0a)
        unsafe_write(stream,pointer(record.data,first(qual_range)),length(qual_range))
        # for i=qual_range
        #     write(stream,record.data[i])
        # end
        write(stream,0x0a)
    elseif length(rng)>0
        if rc
            revcomplement_record!(record,rng)
        end
        write(stream,record.data)
        write(stream,0x0a)
    end
end

function debridgeSE_stream_(bridgeF::String, bridgeR::String, streamX::IO, streamSTATS::Union{IO,Bool},streamsOUT::Union{Array{Array{IO,1},1},Bool}, verbose::Bool,maxreadlen::Int64,revcomp_bRreads::Bool)
    wF=Vector{UInt8}(bridgeF)
    wR=Vector{UInt8}(bridgeR)
    wF_length=length(wF)
    wR_length=length(wR)

    read_counts=zeros(Int64,4,1)
    bridge_pos_mean_X=zeros(Int64,4,1)
    maxl=maxreadlen-min(wF_length,wR_length)+1 #+1 because 1st element is length 0
    rd_len=vcat([zeros(Int64,maxreadlen+1,1)],[zeros(Int64,maxl,maxl) for i in range(1,2)],[zeros(Int64,maxreadlen+1,1)]) #rna length / dna length when both, otherwise just read length

    readerX = FASTQ.Reader(streamX)
    recordX = FASTQ.Record()
    #flag1,flag2,N1,N2,position1,position2,l1,l2
    read_stats=zeros(Int64,1,4)
    writeDebridged=!(streamsOUT==false)
    writePositions=!(streamSTATS==false)
    if writeDebridged
        stream_array_DNA=streamsOUT[1]
        stream_array_RNA=streamsOUT[2]

        flip_read=convert(Vector{Bool},[0;0;1;0]) .& revcomp_bRreads
    end

    nread=0::Int64
    cumtime=0.0
    tic()
    while !eof(readerX)
        nread+=1
        # print(nread,"\n")
        if verbose && mod(nread,100000)==0
            cumtime+=toq()
            print(STDERR,"Processed ", nread/1000000.0, " million reads [", round(nread/cumtime), " reads/s]\n")
            tic()
        end
        read!(readerX, recordX)
        n1X=first(recordX.sequence)
        bridge_pos_X_F=brute_force_search_wildcard(wF,recordX.data,n1X,last(recordX.sequence))
        bridge_pos_X_R=brute_force_search_wildcard(wR,recordX.data,n1X,last(recordX.sequence))

        readlength_X=length(recordX.sequence)
        flagX=0x01+0x01*(bridge_pos_X_F[2]>0)+0x02*(bridge_pos_X_R[2]>0)
        bridge_pos_X_F[2]+bridge_pos_X_R[2]>1 ? flagX=0x04 : true

        read_counts[flagX]+=1
        bridge_pos_mean_X[flagX]+=bridge_pos_X_F[1]+bridge_pos_X_R[1]

        if writePositions
            read_stats[1]=flagX
            read_stats[2]=bridge_pos_X_F[2]+bridge_pos_X_R[2]
            read_stats[3]=max(bridge_pos_X_F[1],bridge_pos_X_R[1])
            read_stats[4]=readlength_X
            writedlm(streamSTATS,read_stats)
        end


        dnaX=1:readlength_X
        rnaX=0:-1
        if flagX==0x02 #unique bridgeF
            rnaX=1:(bridge_pos_X_F[1]-1)
            dnaX=(bridge_pos_X_F[1]+wF_length):readlength_X
            rd_len[flagX][min(maxl,length(rnaX)+1),min(maxl,length(dnaX)+1)]+=1
        elseif flagX==0x03 #unique bridgeR
            dnaX=1:(bridge_pos_X_R[1]-1)
            rnaX=(bridge_pos_X_R[1]+wR_length):readlength_X
            rd_len[flagX][min(maxl,length(rnaX)+1),min(maxl,length(dnaX)+1)]+=1
        else
            rd_len[flagX][min(maxreadlen,readlength_X)+1]+=1
        end


        if writeDebridged

            write_record(stream_array_DNA[flagX],recordX,dnaX,flip_read[flagX])
            write_record(stream_array_RNA[flagX],recordX,rnaX,flip_read[flagX])
        end


    end
    bridge_pos_mean_X_ok=bridge_pos_mean_X./read_counts
    return read_counts, bridge_pos_mean_X_ok, rd_len
end

function debridgeSE_stream_fuzzy_(bridgeF::String, bridgeR::String, streamX::IO, maxerr::Int64, streamSTATS::Union{IO,Bool}, streamsOUT::Union{Array{Array{IO,1},1},Bool}, verbose::Bool, maxreadlen::Int64,revcomp_bRreads::Bool) #note: using that method with maxerr=0 gives same result as debridgePE_stream_ but is 1.5x slower
    wF=Vector{UInt8}(bridgeF)
    wR=Vector{UInt8}(bridgeR)
    wF_length=length(wF)
    wR_length=length(wR)

    read_counts=zeros(Int64,4,1)
    bridge_pos_mean_X=zeros(Int64,4,1)
    maxl=maxreadlen-min(wF_length,wR_length)+1 #+1 because 1st element is length 0
    rd_len=vcat([zeros(Int64,maxreadlen+1,1)],[zeros(Int64,maxl,maxl) for i in range(1,2)],[zeros(Int64,maxreadlen+1,1)]) #rna length / dna length when both, otherwise just read length

    readerX = FASTQ.Reader(streamX)
    recordX = FASTQ.Record()

    #flag1,flag2,N1,N2,position1,position2,l1,l2
    read_stats=zeros(Int64,1,6)
    writeDebridged=!(streamsOUT==false)
    writePositions=!(streamSTATS==false)
    if writeDebridged
        stream_array_DNA=streamsOUT[1]
        stream_array_RNA=streamsOUT[2]
        flip_read=convert(Array{Bool,1},[0 0 1 0]) .& revcomp_bRreads
    end
    nread=0
    cumtime=0.0
    tic()
    while !eof(readerX)
        nread+=1
        # print(nread,"\n")
        if verbose && mod(nread,100000)==0
            cumtime+=toq()
            print(STDERR,"Processed ", nread/1000000.0, " million reads [", round(nread/cumtime), " reads/s]\n")
            tic()
        end
        read!(readerX, recordX)
        n1X=first(recordX.sequence)
        bridge_pos_X_F=brute_force_search_bestfuzzy(wF,recordX.data,n1X,last(recordX.sequence),maxerr)
        bridge_pos_X_R=brute_force_search_bestfuzzy(wR,recordX.data,n1X,last(recordX.sequence),maxerr)

        bridge_pos_X=bestof2searches(bridge_pos_X_F,bridge_pos_X_R)

        readlength_X=length(recordX.sequence)

        flagX=bridge_pos_X[1]

        read_counts[flagX]+=1
        bridge_pos_mean_X[flagX]+=bridge_pos_X_F[1]+bridge_pos_X_R[1]

        if writePositions
            read_stats[1]=flagX
            read_stats[2]=bridge_pos_X[3] #total including confusing ones
            read_stats[3]=bridge_pos_X[1]
            read_stats[4]=readlength_X
            read_stats[5]=bridge_pos_X[2]
            read_stats[6]=bridge_pos_X[4]
            writedlm(streamSTATS,read_stats)
        end


        dnaX=1:readlength_X
        rnaX=0:-1

        if flagX==0x02 #unique bridgeF
            rnaX=1:(bridge_pos_X_F[1]-1)
            dnaX=(bridge_pos_X_F[1]+wF_length):readlength_X
            rd_len[flagX][min(maxl,length(rnaX)+1),min(maxl,length(dnaX)+1)]+=1
        elseif flagX==0x03 #unique bridgeR
            dnaX=1:(bridge_pos_X_R[1]-1)
            rnaX=(bridge_pos_X_R[1]+wR_length):readlength_X
            rd_len[flagX][min(maxl,length(rnaX)+1),min(maxl,length(dnaX)+1)]+=1
        else
            rd_len[flagX][min(maxreadlen,readlength_X)+1]+=1
        end

        if writeDebridged
            write_record_skipempty(stream_array_DNA[flagX],recordX,dnaX,flip_read[flagX])
            write_record_skipempty(stream_array_RNA[flagX],recordX,rnaX,flip_read[flagX])
        end
    end
    bridge_pos_mean_X_ok=bridge_pos_mean_X./read_counts
    return read_counts, bridge_pos_mean_X_ok, rd_len
end

function debridgePE_stream_(bridgeF::String, bridgeR::String, streamX::IO, streamY::IO, streamSTATS::Union{IO,Bool},streamsOUT::Union{Array{Array{IO,2},1},Bool}, verbose::Bool, maxreadlen::Int64, revcomp_bRreads::Bool)
    wF=Vector{UInt8}(bridgeF)
    wR=Vector{UInt8}(bridgeR)
    wF_length=length(wF)
    wR_length=length(wR)

    read_counts=zeros(Int64,4,4)
    bridge_pos_mean_X=zeros(Int64,4,4)
    bridge_pos_mean_Y=zeros(Int64,4,4)
    maxl=maxreadlen-min(wF_length,wR_length)+1 #+1 because 1st element is length 0
    rd_lenX=vcat([zeros(Int64,maxreadlen+1,1)],[zeros(Int64,maxl,maxl) for i in range(1,2)],[zeros(Int64,maxreadlen+1,1)])
    rd_lenY=vcat([zeros(Int64,maxreadlen+1,1)],[zeros(Int64,maxl,maxl) for i in range(1,2)],[zeros(Int64,maxreadlen+1,1)])

    readerX = FASTQ.Reader(streamX)
    readerY = FASTQ.Reader(streamY)
    recordX = FASTQ.Record()
    recordY = FASTQ.Record()


    #flag1,flag2,N1,N2,position1,position2,l1,l2
    read_stats=zeros(Int64,1,8)
    writeDebridged=!(streamsOUT==false)
    writePositions=!(streamSTATS==false)
    if writeDebridged
        stream_array_DNA_PE1=streamsOUT[1]
        stream_array_DNA_PE2=streamsOUT[2]
        stream_array_RNA_PE1=streamsOUT[3]
        stream_array_RNA_PE2=streamsOUT[4]

        # RC in flip are :
        # dnaY in 0R
        # dnaX in R0
        # all in RM or RR or MR
        flip_DNA_PE1=convert(Array{Bool,2},[0 0 0 0 ; 0 0 0 0 ; 1 0 1 1 ; 0 0 1 0]) .& revcomp_bRreads
        flip_DNA_PE2=convert(Array{Bool,2},[0 0 1 0 ; 0 0 0 0 ; 0 0 1 1 ; 0 0 1 0]) .& revcomp_bRreads
        flip_RNA_PE1=convert(Array{Bool,2},[0 0 0 0 ; 0 0 0 0 ; 0 0 1 1 ; 0 0 1 0]) .& revcomp_bRreads
        flip_RNA_PE2=convert(Array{Bool,2},[0 0 0 0 ; 0 0 0 0 ; 0 0 1 1 ; 0 0 1 0]) .& revcomp_bRreads

    end

    nread=0::Int64
    cumtime=0.0
    tic()
    while !eof(readerX)
        nread+=1
        # print(nread,"\n")
        if verbose && mod(nread,100000)==0
            cumtime+=toq()
            print(STDERR,"Processed ", nread/1000000.0, " million reads [", round(nread/cumtime), " reads/s]\n")
            tic()
        end
        read!(readerX, recordX)
        read!(readerY, recordY)
        n1X=first(recordX.sequence)
        n1Y=first(recordY.sequence)
        bridge_pos_X_F=brute_force_search_wildcard(wF,recordX.data,n1X,last(recordX.sequence))
        bridge_pos_X_R=brute_force_search_wildcard(wR,recordX.data,n1X,last(recordX.sequence))
        bridge_pos_Y_F=brute_force_search_wildcard(wF,recordY.data,n1Y,last(recordY.sequence))
        bridge_pos_Y_R=brute_force_search_wildcard(wR,recordY.data,n1Y,last(recordY.sequence))

        readlength_X=length(recordX.sequence)
        readlength_Y=length(recordY.sequence)
        flagX=0x01+0x01*(bridge_pos_X_F[2]>0)+0x02*(bridge_pos_X_R[2]>0)
        flagY=0x01+0x01*(bridge_pos_Y_F[2]>0)+0x02*(bridge_pos_Y_R[2]>0)
        bridge_pos_X_F[2]+bridge_pos_X_R[2]>1 ? flagX=0x04 : true
        bridge_pos_Y_F[2]+bridge_pos_Y_R[2]>1 ? flagY=0x04 : true

        read_counts[flagX,flagY]+=1
        bridge_pos_mean_X[flagX,flagY]+=bridge_pos_X_F[1]+bridge_pos_X_R[1]
        bridge_pos_mean_Y[flagX,flagY]+=bridge_pos_Y_F[1]+bridge_pos_Y_R[1]

        if writePositions
            read_stats[1]=flagX
            read_stats[2]=flagY
            read_stats[3]=bridge_pos_X_F[2]+bridge_pos_X_R[2]
            read_stats[4]=bridge_pos_Y_F[2]+bridge_pos_Y_R[2]
            read_stats[5]=max(bridge_pos_X_F[1],bridge_pos_X_R[1])
            read_stats[6]=max(bridge_pos_Y_F[1],bridge_pos_Y_R[1])
            read_stats[7]=readlength_X
            read_stats[8]=readlength_Y
            writedlm(streamSTATS,read_stats)
        end


        dnaX=1:readlength_X
        dnaY=1:readlength_Y
        rnaX=0:-1
        rnaY=0:-1

        if flagX==0x02 #unique bridgeF
            rnaX=1:(bridge_pos_X_F[1]-1)
            dnaX=(bridge_pos_X_F[1]+wF_length):readlength_X
            rd_lenX[flagX][min(maxl,length(rnaX)+1),min(maxl,length(dnaX)+1)]+=1
        elseif flagX==0x03 #unique bridgeR
            dnaX=1:(bridge_pos_X_R[1]-1)
            rnaX=(bridge_pos_X_R[1]+wR_length):readlength_X
            rd_lenX[flagX][min(maxl,length(rnaX)+1),min(maxl,length(dnaX)+1)]+=1
        else
            rd_lenX[flagX][min(maxreadlen,readlength_X)+1]+=1
        end

        if flagY==0x02 #unique bridgeF
            rnaY=1:(bridge_pos_Y_F[1]-1)
            dnaY=(bridge_pos_Y_F[1]+wF_length):readlength_Y
            rd_lenY[flagY][min(maxl,length(rnaY)+1),min(maxl,length(dnaY)+1)]+=1
        elseif flagY==0x03 #unique bridgeR
            dnaY=1:(bridge_pos_Y_R[1]-1)
            rnaY=(bridge_pos_Y_R[1]+wR_length):readlength_Y
            rd_lenY[flagY][min(maxl,length(rnaY)+1),min(maxl,length(dnaY)+1)]+=1
        else
            rd_lenY[flagY][min(maxreadlen,readlength_Y)+1]+=1
        end



        if writeDebridged
            write_record(stream_array_DNA_PE1[flagX,flagY],recordX,dnaX,flip_DNA_PE1[flagX,flagY])
            write_record(stream_array_DNA_PE2[flagX,flagY],recordY,dnaY,flip_DNA_PE2[flagX,flagY])
            write_record(stream_array_RNA_PE1[flagX,flagY],recordX,rnaX,flip_RNA_PE1[flagX,flagY])
            write_record(stream_array_RNA_PE2[flagX,flagY],recordY,rnaY,flip_RNA_PE2[flagX,flagY])
        end
    end
    bridge_pos_mean_X_ok=bridge_pos_mean_X./read_counts
    bridge_pos_mean_Y_ok=bridge_pos_mean_Y./read_counts
    return read_counts, bridge_pos_mean_X_ok, bridge_pos_mean_Y_ok, rd_lenX, rd_lenY
end

function debridgePE_stream_fuzzy_(bridgeF::String, bridgeR::String, streamX::IO, streamY::IO, maxerr::Int64, streamSTATS::Union{IO,Bool}, streamsOUT::Union{Array{Array{IO,2},1},Bool}, verbose::Bool, maxreadlen::Int64, revcomp_bRreads::Bool) #note: using that method with maxerr=0 gives same result as debridgePE_stream_ but is 1.5x slower
    wF=Vector{UInt8}(bridgeF)
    wR=Vector{UInt8}(bridgeR)
    wF_length=length(wF)
    wR_length=length(wR)

    read_counts=zeros(Int64,4,4)
    bridge_pos_mean_X=zeros(Int64,4,4)
    bridge_pos_mean_Y=zeros(Int64,4,4)
    maxl=maxreadlen-min(wF_length,wR_length)+1 #+1 because 1st element is length 0
    rd_lenX=vcat([zeros(Int64,maxreadlen+1,1)],[zeros(Int64,maxl,maxl) for i in range(1,2)],[zeros(Int64,maxreadlen+1,1)])
    rd_lenY=vcat([zeros(Int64,maxreadlen+1,1)],[zeros(Int64,maxl,maxl) for i in range(1,2)],[zeros(Int64,maxreadlen+1,1)])


    readerX = FASTQ.Reader(streamX)
    readerY = FASTQ.Reader(streamY)
    recordX = FASTQ.Record()
    recordY = FASTQ.Record()


    #flag1,flag2,N1,N2,position1,position2,l1,l2
    read_stats=zeros(Int64,1,12)
    writeDebridged=!(streamsOUT==false)
    writePositions=!(streamSTATS==false)
    if writeDebridged
        stream_array_DNA_PE1=streamsOUT[1]
        stream_array_DNA_PE2=streamsOUT[2]
        stream_array_RNA_PE1=streamsOUT[3]
        stream_array_RNA_PE2=streamsOUT[4]

        # RC in flip are :
        # dnaY in 0R
        # dnaX in R0
        # all in RM or RR or MR
        flip_DNA_PE1=convert(Array{Bool,2},[0 0 0 0 ; 0 0 0 0 ; 1 0 1 1 ; 0 0 1 0]) .& revcomp_bRreads
        flip_DNA_PE2=convert(Array{Bool,2},[0 0 1 0 ; 0 0 0 0 ; 0 0 1 1 ; 0 0 1 0]) .& revcomp_bRreads
        flip_RNA_PE1=convert(Array{Bool,2},[0 0 0 0 ; 0 0 0 0 ; 0 0 1 1 ; 0 0 1 0]) .& revcomp_bRreads
        flip_RNA_PE2=convert(Array{Bool,2},[0 0 0 0 ; 0 0 0 0 ; 0 0 1 1 ; 0 0 1 0]) .& revcomp_bRreads
    end
    nread=0
    cumtime=0.0
    tic()
    while !eof(readerX)
        nread+=1
        # print(nread,"\n")
        if verbose && mod(nread,100000)==0
            cumtime+=toq()
            print(STDERR,"Processed ", nread/1000000.0, " million reads [", round(nread/cumtime), " reads/s]\n")
            tic()
        end
        read!(readerX, recordX)
        read!(readerY, recordY)
        n1X=first(recordX.sequence)
        n1Y=first(recordY.sequence)
        bridge_pos_X_F=brute_force_search_bestfuzzy(wF,recordX.data,n1X,last(recordX.sequence),maxerr)
        bridge_pos_X_R=brute_force_search_bestfuzzy(wR,recordX.data,n1X,last(recordX.sequence),maxerr)
        bridge_pos_Y_F=brute_force_search_bestfuzzy(wF,recordY.data,n1Y,last(recordY.sequence),maxerr)
        bridge_pos_Y_R=brute_force_search_bestfuzzy(wR,recordY.data,n1Y,last(recordY.sequence),maxerr)

        bridge_pos_X=bestof2searches(bridge_pos_X_F,bridge_pos_X_R)
        bridge_pos_Y=bestof2searches(bridge_pos_Y_F,bridge_pos_Y_R)

        readlength_X=length(recordX.sequence)
        readlength_Y=length(recordY.sequence)

        flagX=bridge_pos_X[1]
        flagY=bridge_pos_Y[1]

        read_counts[flagX,flagY]+=1
        bridge_pos_mean_X[flagX,flagY]+=bridge_pos_X_F[1]+bridge_pos_X_R[1]
        bridge_pos_mean_Y[flagX,flagY]+=bridge_pos_Y_F[1]+bridge_pos_Y_R[1]

        if writePositions
            read_stats[1]=flagX
            read_stats[2]=flagY
            read_stats[3]=bridge_pos_X[3] #total including confusing ones
            read_stats[4]=bridge_pos_Y[3]
            read_stats[5]=bridge_pos_X[1]
            read_stats[6]=bridge_pos_Y[1]
            read_stats[7]=readlength_X
            read_stats[8]=readlength_Y
            read_stats[9]=bridge_pos_X[2]
            read_stats[10]=bridge_pos_Y[2]
            read_stats[11]=bridge_pos_X[4]
            read_stats[12]=bridge_pos_Y[4]
            writedlm(streamSTATS,read_stats)
        end


        dnaX=1:readlength_X
        dnaY=1:readlength_Y
        rnaX=0:-1
        rnaY=0:-1

        if flagX==0x02 #unique bridgeF
            rnaX=1:(bridge_pos_X_F[1]-1)
            dnaX=(bridge_pos_X_F[1]+wF_length):readlength_X
            rd_lenX[flagX][min(maxl,length(rnaX)+1),min(maxl,length(dnaX)+1)]+=1
        elseif flagX==0x03 #unique bridgeR
            dnaX=1:(bridge_pos_X_R[1]-1)
            rnaX=(bridge_pos_X_R[1]+wR_length):readlength_X
            rd_lenX[flagX][min(maxl,length(rnaX)+1),min(maxl,length(dnaX)+1)]+=1
        else
            rd_lenX[flagX][min(maxreadlen,readlength_X)+1]+=1
        end

        if flagY==0x02 #unique bridgeF
            rnaY=1:(bridge_pos_Y_F[1]-1)
            dnaY=(bridge_pos_Y_F[1]+wF_length):readlength_Y
            rd_lenY[flagY][min(maxl,length(rnaY)+1),min(maxl,length(dnaY)+1)]+=1
        elseif flagY==0x03 #unique bridgeR
            dnaY=1:(bridge_pos_Y_R[1]-1)
            rnaY=(bridge_pos_Y_R[1]+wR_length):readlength_Y
            rd_lenY[flagY][min(maxl,length(rnaY)+1),min(maxl,length(dnaY)+1)]+=1
        else
            rd_lenY[flagY][min(maxreadl,readlength_Y)+1]+=1
        end



        if writeDebridged
            # write_record(stream_array_DNA_PE1[flagX,flagY],recordX,dnaX)
            # write_record(stream_array_DNA_PE2[flagX,flagY],recordY,dnaY)
            # write_record(stream_array_RNA_PE1[flagX,flagY],recordX,rnaX)
            # write_record(stream_array_RNA_PE2[flagX,flagY],recordY,rnaY)
            write_record(stream_array_DNA_PE1[flagX,flagY],recordX,dnaX,flip_DNA_PE1[flagX,flagY])
            write_record(stream_array_DNA_PE2[flagX,flagY],recordY,dnaY,flip_DNA_PE2[flagX,flagY])
            write_record(stream_array_RNA_PE1[flagX,flagY],recordX,rnaX,flip_RNA_PE1[flagX,flagY])
            write_record(stream_array_RNA_PE2[flagX,flagY],recordY,rnaY,flip_RNA_PE2[flagX,flagY])
        end
    end
    bridge_pos_mean_X_ok=bridge_pos_mean_X./read_counts
    bridge_pos_mean_Y_ok=bridge_pos_mean_Y./read_counts
    return read_counts, bridge_pos_mean_X_ok, bridge_pos_mean_Y_ok, rd_lenX, rd_lenY
end


function makestream_out(fileprefix::String,filename::String, compress::Bool=true, ext::String=".fastq")
    if isempty(fileprefix)
        return STDOUT
    else
        if compress
            codec=GzipCompressor()
            return TranscodingStream(codec, open(fileprefix * filename * ext * ".gz","w"))
        else
            codec=Noop()
            return TranscodingStream(codec, open(fileprefix * filename * ext,"w"))
        end
    end
end

function makestream_in(filepath::String)
    if endswith(filepath, ".gz")
        codec = GzipDecompressor()
        return TranscodingStream(codec, open(filepath))
    else
        codec = Noop()
        return TranscodingStream(codec, open(filepath))
    end
end

function summarize_debridgingPE_(read_counts, parsable::Bool, stream::IO)
    println(stream, ">>Total counts")
    println(stream, sum(read_counts[1]))
    println(stream,">>Count matrix")
    Base.showarray(stream,read_counts[1],parsable)
    println(stream, "\n>>Count %")
    Base.showarray(stream,read_counts[1]*100.0/sum(read_counts[1]),parsable)
    println(stream, "\n>>Bridge mean position read 1")
    Base.showarray(stream,read_counts[2],parsable)
    println(stream, "\n>>Bridge mean position read 2")
    Base.showarray(stream,read_counts[3],parsable)
    println(stream, "\n>>RNA/DNA length distribution read 1/2")
    bridgecodes=["0","F","R","M"]
    for i=1:4
        print(stream,"bridge_" * bridgecodes[i] * "_1\t")
        Base.showarray(stream,read_counts[4][i],parsable)
        print(stream,"\nbridge_" * bridgecodes[i] * "_2\t")
        Base.showarray(stream,read_counts[5][i],parsable)
        println(stream,"")
    end
    # println("")
end

function summarize_debridgingSE_(read_counts, parsable::Bool, stream::IO)
    println(stream, ">>Total counts")
    println(stream,sum(read_counts[1]))
    println(stream,">>Count matrix")
    Base.showarray(stream,read_counts[1],parsable)
    println(stream, "\n>>Count %")
    Base.showarray(stream,read_counts[1]*100.0/sum(read_counts[1]),parsable)
    println(stream, "\n>>Bridge mean position")
    Base.showarray(stream,read_counts[2],parsable)
    println(stream, "\n>>RNA/DNA length distribution")
    bridgecodes=["0","F","R","M"]
    for i=1:4
        print(stream,"bridge_" * bridgecodes[i] * "\t")
        Base.showarray(stream,read_counts[3][i],parsable)
        println(stream,"")
    end
    # println("")
end

function debridgePE(bridgeF::String, bridgeR::String, fqX::String, fqY::String, fqout_prefix::String, maxerr::Int64=0, writePositions::Int64=2, writeDebridged::Int64=0,  verbose::Bool=true, maxreadl::Int=151, rc::Bool=true)

    streamX=makestream_in(fqX)
    streamY=makestream_in(fqY)

    streamSUMMARY=makestream_out(fqout_prefix,"summary.PE",false,".txt")
    if rc
        println(streamSUMMARY, "Running debridge PE with reverse complement mode ON")
    end
    if writePositions>0
        streamSTATS=makestream_out(fqout_prefix,"positions.PE",writePositions>1,".txt") #2 stands for PE
    else
        streamSTATS=false
    end
    if writeDebridged>0
        compressOut=writeDebridged>1
        #00 reads

        # RC in flip are :
        # dnaY in 0R
        # dnaX in R0
        # all in RM or RR
        NULL1=makestream_out(fqout_prefix,"00.xxx.1",compressOut)
        NULL2=makestream_out(fqout_prefix,"00.xxx.2",compressOut)
        #F0 reads
        F0rna1=makestream_out(fqout_prefix,"F0.rna.1",compressOut)#sense
        F0dna1=makestream_out(fqout_prefix,"F0.dna.1",compressOut)
        F0dna2=makestream_out(fqout_prefix,"F0.dna.2",compressOut)
        #R0 reads (read2 more likely RNA?)

        if rc
            R0rna1=makestream_out(fqout_prefix, "0R.rna.2",compressOut) #comes from the R read, which is mate 2 in flip mode. This RNA is antisense
            R0dna1=makestream_out(fqout_prefix, "0R.dna",compressOut) #THIS WILL GET RC
            R0rna2=makestream_out(fqout_prefix, "0R.rna.1",compressOut) # comes from the zero read which is .1 in flip version. This RNA is sense
        else
            R0rna1=makestream_out(fqout_prefix, "R0.rna.1",compressOut) # comes from the R read, which is mate 1 in original mode. This RNA is antisense.
            R0dna1=makestream_out(fqout_prefix, "R0.dna.1",compressOut)
            R0rna2=makestream_out(fqout_prefix, "R0.rna.2",compressOut) # comes from the zero read which is .2 in the initial no flip version. This RNA is sense
        end

        #FM reads
        FMrna1=makestream_out(fqout_prefix, "FM.rna.1",compressOut) #sense
        FMdna1=makestream_out(fqout_prefix, "FM.dna.1",compressOut)
        FMdna2=makestream_out(fqout_prefix, "FM.xxx.2",compressOut)
        #RM reads
        RMrna1=makestream_out(fqout_prefix, "RM.rna.1",compressOut) #antisense, or gets RC in flip mode so sense
        RMdna1=makestream_out(fqout_prefix, "RM.dna.1",compressOut) #antisense, or gets RC in flip mode so sense
        RMdna2=makestream_out(fqout_prefix, "RM.xxx.2",compressOut) # gets RC in flip mode
        #M0 reads
        M01=makestream_out(fqout_prefix, "M0.xxx.1",compressOut)
        M02=makestream_out(fqout_prefix, "M0.xxx.2",compressOut)
        #MM reads
        MM1=makestream_out(fqout_prefix, "MM.xxx.1",compressOut)
        MM2=makestream_out(fqout_prefix, "MM.xxx.2",compressOut)
        #FF reads
        FFrna1=makestream_out(fqout_prefix, "FF.rna.1",compressOut)#sense
        FFdna1=makestream_out(fqout_prefix, "FF.dna.1",compressOut)
        FFrna2=makestream_out(fqout_prefix, "FF.rna.2",compressOut)#sense
        FFdna2=makestream_out(fqout_prefix, "FF.dna.2",compressOut)
        #RR reads
        RRrna1=makestream_out(fqout_prefix, "RR.rna.1",compressOut) #antisense, get Rc in flip
        RRdna1=makestream_out(fqout_prefix, "RR.dna.1",compressOut) # get RC in flip
        RRrna2=makestream_out(fqout_prefix, "RR.rna.2",compressOut)#antisense #gets RC in flip
        RRdna2=makestream_out(fqout_prefix, "RR.dna.2",compressOut)#gets RC in flip
        #FR reads
        FRrna1=makestream_out(fqout_prefix, "FR.rna.1",compressOut) #sense
        FRdna1=makestream_out(fqout_prefix, "FR.dna.1",compressOut)
        FRrna2=makestream_out(fqout_prefix, "FR.rna.2",compressOut) #antisense
        FRdna2=makestream_out(fqout_prefix, "FR.dna.2",compressOut)

        NULLstream=STDOUT #

        # pritories for reversing reads: F, R, M, 0
        #[00, 0F*, 0R*&, 0M*; F0, FF, FR, FM; R0&, RF*, RR, RM; M0, MF*, MR*, MM] * we reverse the reads, & 2nd read (after reverse) contains no DNA

        #if rc
        #    R0rna1=makestream_out(fqout_prefix, "0R.rna.2",compressOut) #comes from the R read, which is mate 2 in flip mode. This RNA is antisense
    #        R0dna1=makestream_out(fqout_prefix, "0R.dna",compressOut) #THIS WILL GET RC
#            R0rna2=makestream_out(fqout_prefix, "0R.rna.1",compressOut) # comes from the zero read which is .1 in flip version. This RNA is sense

        stream_array_DNA_PE1=[NULL1 F0dna2 R0rna2 M02; F0dna1 FFdna1 FRdna1 FMdna1; R0dna1 FRdna2 RRdna1 RMdna1; M01 FMdna2 RMdna2 MM1]
        stream_array_DNA_PE2=[NULL2 F0dna1 R0dna1 M01; F0dna2 FFdna2 FRdna2 FMdna2; R0rna2 FRdna1 RRdna2 RMdna2; M02 FMdna1 RMdna1 MM2]
        #[00&&, 0F*&, 0R*, 0M*&&; F0&, FF, FR, FM&; R0, RF*, RR, RM&; M0&&, MF*&, MR*&, MM&&] * we reverse the reads, & 2nd read contains no RNA
        stream_array_RNA_PE1=[NULLstream NULLstream NULLstream NULLstream; F0rna1 FFrna1 FRrna1 FMrna1; R0rna1 FRrna2 RRrna1 RMrna1; NULLstream NULLstream NULLstream NULLstream]
        stream_array_RNA_PE2=[NULLstream F0rna1 R0rna1 NULLstream; NULLstream FFrna2 FRrna2 NULLstream; NULLstream FRrna1 RRrna2 NULLstream; NULLstream FMrna1 RMrna1 NULLstream]

        streamsOUT=[stream_array_DNA_PE1,stream_array_DNA_PE2,stream_array_RNA_PE1,stream_array_RNA_PE2]
    else
        streamsOUT=false
    end
    if maxerr>0
        read_counts=debridgePE_stream_fuzzy_(bridgeF, bridgeR, streamX, streamY, maxerr,streamSTATS,streamsOUT, verbose, maxreadl, rc)
    else
        read_counts=debridgePE_stream_(bridgeF, bridgeR, streamX, streamY,streamSTATS,streamsOUT, verbose, maxreadl, rc)
    end
    if writeDebridged>0
        for i=1:4
            for j=1:4
                !(stream_array_DNA_PE1[i,j]==STDOUT) ? close(stream_array_DNA_PE1[i,j]) : true
                !(stream_array_DNA_PE2[i,j]==STDOUT) ? close(stream_array_DNA_PE2[i,j]) : true
                !(stream_array_RNA_PE1[i,j]==STDOUT) ? close(stream_array_RNA_PE1[i,j]) : true
                !(stream_array_RNA_PE1[i,j]==STDOUT) ? close(stream_array_RNA_PE1[i,j]) : true
            end
        end
    end
    if writePositions>0
        !(streamSTATS==STDOUT) ? close(streamSTATS) : true
    end

    close(streamX)
    close(streamY)

    summarize_debridgingPE_(read_counts,true,streamSUMMARY)
    !(streamSUMMARY==STDOUT) ? close(streamSUMMARY) : true
    return read_counts

end


function debridgeSE(bridgeF::String, bridgeR::String, fqX::String, fqout_prefix::String, maxerr::Int64=0, writePositions::Int64=1, writeDebridged::Int64=0,  verbose::Bool=true, maxreadl::Int=151, rc::Bool=true)
    streamX=makestream_in(fqX)
    streamSUMMARY=makestream_out(fqout_prefix,"summary.SE",false,".txt")

    if rc
        println(streamSUMMARY, "Running debridge SE with reverse complement mode ON")
    end

    if writePositions>0
        streamSTATS=makestream_out(fqout_prefix,"positions.SE",writePositions>1,".txt")
    else
        streamSTATS=false
    end

    if writeDebridged>0
        compressOut=writeDebridged>1
        #0 reads
        NULL=makestream_out(fqout_prefix,"0.xxx",compressOut)
        #F reads
        Frna=makestream_out(fqout_prefix,"F.rna",compressOut)#sense
        Fdna=makestream_out(fqout_prefix,"F.dna",compressOut)

        #R reads (read2 more likely RNA?)
        Rrna=makestream_out(fqout_prefix, "R.rna",compressOut) # RNA is antisense in nomrmal or sense in FLIP mode
        Rdna=makestream_out(fqout_prefix, "R.dna",compressOut)

        #M reads
        M=makestream_out(fqout_prefix, "M.xxx",compressOut)

        NULLstream=STDOUT #

        # pritories for reversing reads: F, R, M, 0
        #[00, 0F*, 0R*&, 0M*; F0, FF, FR, FM; R0&, RF*, RR, RM; M0, MF*, MR*, MM] * we reverse the reads, & 2nd read (after reverse) contains no DNA
        stream_array_DNA=[NULL, Fdna, Rdna, M]
        stream_array_RNA=[NULLstream, Frna, Rrna, NULLstream]
        streamsOUT=[stream_array_DNA,stream_array_RNA]
    else
        streamsOUT=false
    end

    if maxerr>0
        read_counts=debridgeSE_stream_fuzzy_(bridgeF, bridgeR, streamX, maxerr,streamSTATS,streamsOUT, verbose, maxreadl, rc)
    else
        read_counts=debridgeSE_stream_(bridgeF, bridgeR, streamX ,streamSTATS,streamsOUT, verbose, maxreadl, rc)
    end
    if writeDebridged>0
        for i=1:4
            !(stream_array_DNA[i]==STDOUT) ? close(stream_array_DNA[i]) : true
            !(stream_array_RNA[i]==STDOUT) ? close(stream_array_RNA[i]) : true
        end
    end
    if writePositions>0
        !(streamSTATS==STDOUT) ? close(streamSTATS) : true
    end

    close(streamX)
    summarize_debridgingSE_(read_counts,true,streamSUMMARY)
    !(streamSUMMARY==STDOUT) ? close(streamSUMMARY) : true
    return read_counts

end


#


# bF="AAAAAAACCGGCGTCCAAG"
# bR="CTTGGACGCCGGTTTTTTT"
# fqX="s1head_R1.fastq.gz"
# fqY="s1head_R2.fastq.gz"
# fqOutPrefix="out_"

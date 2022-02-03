#!/usr/bin/env julia
using ArgParse
using BridgeFinder
#

function parse_commandline()
    s = ArgParseSettings(description= "Finds bridge and split reads in paired end CHarSeq data")
    @add_arg_table s begin
        "--writedebridged", "-d"
            help = "output mode: 0=don't write out, 1=write fastq, 2=compress"
            # nargs = 1
            arg_type = Int
            default = 0
        "--writepositions", "-p"
            help = "output mode: 0=don't write out, 1=write fastq, 2=compress"
            # nargs = 1
            arg_type = Int
            default = 1
        "--maxerr", "-e"
            help = "Maximum number of mismatches tolerated during bridge search"
            arg_type = Int
            default = 0
        "--maxreadl", "-l"
            help = "Maximum read length"
            arg_type = Int
            default = 151
        "--verbose", "-v"
            help = "verbose, prints progress"
            action = :store_true
        "--singleend", "-s"
            help = "single end reads"
            action = :store_true
        "--revcomp", "-r"
            help = "reverse complement R type reads"
            action = :store_true
        "bridgeF"
            help = "bridge forward sequence"
            required = true
        "bridgeR"
            help = "bridge reverse sequence"
            required = true
        "fq_out_prefix"
            help = "prefix of output files"
            required = true
        "fqX"
            help = "fastq file with reads1"
            required = true
        "fqY"
            help = "fastq file with reads2 (used only in paired end mode)"
            required = false


    end
    parsed_args = parse_args(s)
end

function main()
    dict = parse_commandline()
    if dict["singleend"]
        debridgeSE(dict["bridgeF"],dict["bridgeR"],dict["fqX"],dict["fq_out_prefix"],dict["maxerr"],dict["writepositions"],dict["writedebridged"],dict["verbose"],dict["maxreadl"], dict["revcomp"])
    else
        debridgePE(dict["bridgeF"],dict["bridgeR"],dict["fqX"],dict["fqY"],dict["fq_out_prefix"],dict["maxerr"],dict["writepositions"],dict["writedebridged"],dict["verbose"],dict["maxreadl"], dict["revcomp"])
    end
end

main()

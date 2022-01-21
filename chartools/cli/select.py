from ..filtering import ReadsFanout
import json
from tagtools import utils

def add_subcommand_select(subparsers):
    parser = subparsers.add_parser('select')
    parser.set_defaults(func=select)
    parser.add_argument('--version', action='version', version='1.0.0')
    parser.add_argument('gps', help='groups')
    parser.add_argument('pairs', help='pairs file (MUST be sorted by readID)')
    parser.add_argument('tagfile', help='yaml tagfile')
    parser.add_argument('outprfx', help='output pairs file prfx')
    parser.add_argument('--txdict','-t', help='txtogene conversion')
    parser.add_argument('--gpnames','-n', help='gp1,gp2,etc...')
    parser.add_argument('--genes_by_name', '-e', action='store_true', help='gene names in english, rather that ENSG') #english
    parser.add_argument('--bytx', '-g', action='store_true', help='switches to tx level mode') 
    #switches to gen
    parser.add_argument('--noambig', '-s', action='store_true', help='switches to STRICT mode')
    #  rna_bam_in=rna, dna_bam_in=dna, readid_out_prfx='k_', pairs_out_prfx='k_', rnadna_out_prfx='k_', to_stdout="rna", nmax=0
    parser.add_argument('--rna','-R', help='rna bam in')
    parser.add_argument('--dna','-D', help='dna bam in')
    parser.add_argument('--indexbams', '-i', action='store_true', help='sort and indexes out bams')
    parser.add_argument('--stdout','-o', help='choose what goes to stdout', choices=['pairs','readids','rna','dna'])
    parser.add_argument('--minflight', help='min flight distance', type=int, default=0)
    parser.add_argument('--nmax', help='limits parsing to that many reads')
   
    # parser.set_defaults(func=tag)

def select(args):
    print('pairup', args)
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE,SIG_DFL) #allows to do head operations without throwing errors

    print("loading tag file")
    with open(args.tagfile, "r") as ambiv_file:
        ambiv_data = json.load(ambiv_file)
    if args.bytx:
        # ambiv_LUT=ambiv_data[1][0]
        ambiv_LUT=utils.make_ambiv_LUT(ambiv_data[0])[0]
    else:
        ambiv_LUT=utils.make_ambiv_LUT(ambiv_data[1])[0]
        # ambiv_LUT=ambiv_data[0][0]

    print("loading annotations file")
    with open(args.txdict, "r") as read_file:
        annot_data=json.load(read_file)

    tx_dict=annot_data[0]
    name_to_ENSG=annot_data[1]
    bygene=not(args.bytx)

    tx_pre_list=args.gps.split(",")
    tx_list_names=[pre.split("+") for pre in tx_pre_list]

    if args.genes_by_name:
        tx_list_list=[[name_to_ENSG.get(s,"*") for s in g] for g in tx_list_names]
    else:
        tx_list_list=tx_list_names

    if args.gpnames is None:
        gpnames=[str(i+1) for i in range(len(tx_list_list))]
    else:
        gpnames=args.gpnames.split(",")
    
    # if args.minflight is None:
    #     myminflight=0
    # else:
    #     myminflight=args.minflight

    # if (args.minflight>0) and (len(gpnames)<2*len(tx_list_list)):
    #     gpnames2=gpnames+[x+"_selflocus" for x in gpnames]
    #     gpnames=gpnames2 

    print("started!")
    fanout=ReadsFanout(tx_list_list, ambiv_LUT, bygene=bygene, tx_dict=tx_dict, group_names=gpnames, noambig=args.noambig, min_flight=args.minflight)
    fanout.fanout(args.pairs, rna_bam_in=args.rna, dna_bam_in=args.dna, readid_out_prfx=args.outprfx, pairs_out_prfx=args.outprfx, rnadna_out_prfx=args.outprfx, to_stdout=args.stdout, indexbams=args.indexbams, nmax=args.nmax)
    

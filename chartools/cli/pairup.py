from ..pairup import bams2pairs
from ..pairup_R2 import bams2pairs as bams2pairs_R2

def add_subcommand_pairup(subparsers):
    parser = subparsers.add_parser('pairup')
    parser.set_defaults(func=pairup)
    parser.add_argument('--version', action='version', version='1.0.0')
    parser.add_argument('rna', help='RNA bam (MUST be sorted by readID)')
    parser.add_argument('dna', help='DNA bam (MUST be sorted by readID)')

    parser.add_argument('--outpair','-o', help='write pairs to OUTPAIRS instead of stdout')
    parser.add_argument('--outpair_lowtrig','-l', help='write lowtrig pairs to OUTPAIRS2 instead of stdout')
    parser.add_argument('--outpairedRNA', help='write paired RNA to bam/sam file')
    parser.add_argument('--outpairedDNA', help='write paired DNA to bam/sam file')
    parser.add_argument('--outunpairedRNA', help='write unpaired RNA to bam/sam file')
    parser.add_argument('--outunpairedDNA', help='write unpaired DNA to bam/sam file')
    parser.add_argument('--pairingmode','-m', help='pairing mode', type=str, choices=["triangular","prefix","rd"])
    parser.add_argument('--outstats','-s', help='write statistics to file rather than STDERR')
    parser.add_argument('--reducemmap',help='keep only one for each multimapper', type=int, choices=[0,1,2])
    parser.add_argument('--removeannots',help='reads with annotations do not go into pairs file', action='store_true')
    parser.add_argument('--r2',help='RNA-RNA mode', action='store_true')
    # parser.set_defaults(func=tag)

def pairup(args):
    print('pairup', args)
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE,SIG_DFL)
    if args.r2:
        bams2pairs_R2(args.rna,args.dna,outfile_pairs=args.outpair, outfile_pairedRNA=args.outpairedRNA, outfile_pairedDNA=args.outpairedDNA, outfile_unpairedRNA=args.outunpairedRNA,outfile_unpairedDNA=args.outunpairedDNA, pairingmode=args.pairingmode, outfile_pairs_lowtrig=args.outpair_lowtrig, invert=False, naturalSort=False, outfile_stats=args.outstats, reduceMultimapper=args.reducemmap, filter_with_annots=args.removeannots)
    else:
        bams2pairs(args.rna,args.dna,outfile_pairs=args.outpair, outfile_pairedRNA=args.outpairedRNA, outfile_pairedDNA=args.outpairedDNA, outfile_unpairedRNA=args.outunpairedRNA,outfile_unpairedDNA=args.outunpairedDNA, pairingmode=args.pairingmode, outfile_pairs_lowtrig=args.outpair_lowtrig, invert=False, naturalSort=False, outfile_stats=args.outstats, reduceMultimapper=args.reducemmap, filter_with_annots=args.removeannots)

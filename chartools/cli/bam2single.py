from ..rnaonly import bam2single

def add_subcommand_bam2single(subparsers):
    parser = subparsers.add_parser('bam2single')
    parser.set_defaults(func=bam2single_cmd)
    parser.add_argument('--version', action='version', version='1.0.0')
    parser.add_argument('rna', help='RNA bam (MUST be sorted by readID)')

    parser.add_argument('--outsingle','-o', help='write pairs to OUTPAIRS instead of stdout')
    parser.add_argument('--reducemmap',help='keep only one for each multimapper', type=int, choices=[0,1,2])
    parser.add_argument('--removeannots',help='reads with annotations do not go into pairs file', action='store_true')
    # parser.set_defaults(func=tag)

def bam2single_cmd(args):
    print('bam2single', args)
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE,SIG_DFL)
    bam2single(args.rna,outfile_single=args.outsingle, invert=False, naturalSort=False, reduceMultimapper=args.reducemmap, filter_with_annots=args.removeannots)

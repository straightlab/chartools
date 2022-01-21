from ..dispatch_pairs import cistrans, deambig, within_bp


def add_subcommand_cistrans(subparsers):
    parser = subparsers.add_parser('cistrans')
    parser.set_defaults(func=cistrans_rd)
    parser.add_argument('--version', action='version', version='1.0.0')
    parser.add_argument('rna', help='RNA bam (MUST be sorted by readID)')
    parser.add_argument('dna', help='DNA bam (MUST be sorted by readID)')

    parser.add_argument('--rcis','-r', help='CIS rna')
    parser.add_argument('--dcis','-d', help='CIS dna')
    parser.add_argument('--rtrans','-R', help='TRANS rna')
    parser.add_argument('--dtrans','-D', help='TRANS dna')
    parser.add_argument('--rambig','-s', help='AMBIG rna')
    parser.add_argument('--dambig','-e', help='AMBIG dna')
    parser.add_argument('--nmax','-n', type=int, help='up to n reads', default=0)
    parser.add_argument('--qrna', type=int, help='Min q for RNA', default=0)
    parser.add_argument('--qdna', type=int, help='Min q for DNA', default=15)

    # parser.set_defaults(func=tag)

def add_subcommand_deambig(subparsers):
    parser = subparsers.add_parser('deambig')
    parser.set_defaults(func=deambig_rd)
    parser.add_argument('--version', action='version', version='1.0.0')
    parser.add_argument('rna', help='RNA bam (MUST be sorted by readID)')
    parser.add_argument('dna', help='DNA bam (MUST be sorted by readID)')

    parser.add_argument('--rout','-r', help='CIS rna')
    parser.add_argument('--dout','-d', help='CIS dna')
    
    parser.add_argument('--nmax','-n', type=int, help='up to n reads', default=0)
    parser.add_argument('--qrna', type=int, help='Min q for RNA', default=255)

def add_subcommand_withinbp(subparsers):
    parser = subparsers.add_parser('withinbp')
    parser.set_defaults(func=withinbp_rd)
    parser.add_argument('--version', action='version', version='1.0.0')
    parser.add_argument('rna', help='RNA bam (MUST be sorted by readID)')
    parser.add_argument('dna', help='DNA bam (MUST be sorted by readID)')

    parser.add_argument('--rnaout','-r', help='CIS rna')
    parser.add_argument('--dnaout','-d', help='CIS dna')
    
    parser.add_argument('--nmax','-n', type=int, help='up to n reads', default=0)
    parser.add_argument('--qrna', type=int, help='Min q for RNA', default=0)
    parser.add_argument('--qdna', type=int, help='Min q for DNA', default=15)
    parser.add_argument('--dmin', type=int, help='Min q for RNA', default=-1000)
    parser.add_argument('--dmax', type=int, help='Min q for DNA', default=1000)

def cistrans_rd(args):
    print('cistrans', args)
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE,SIG_DFL)
    cistrans(args.rna,args.dna,rcis=args.rcis, dcis=args.dcis, rtrans=args.rtrans, dtrans=args.dtrans, rambig=args.rambig, dambig=args.dambig, qminRNA=args.qrna, qminDNA=args.qdna, nmax=args.nmax)

def deambig_rd(args):
    print('deambig', args)
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE,SIG_DFL)
    deambig(args.rna,args.dna,rout=args.rout, dout=args.dout, qminRNA=args.qrna,nmax=args.nmax)

def withinbp_rd(args):
    print('withinbp', args)
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE,SIG_DFL)
    within_bp(args.rna,args.dna, rwithin=args.rnaout, dwithin=args.dnaout, qminRNA=args.qrna, qminDNA=args.qdna, mind=args.dmin, maxd=args.dmax, nmax=args.nmax)
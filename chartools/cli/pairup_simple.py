from ..dispatch_pairs import bamsToPairsSimple


def add_subcommand_pairup_simple(subparsers):
    parser = subparsers.add_parser('pairup_simple')
    parser.set_defaults(func=pairup_simple)
    parser.add_argument('--version', action='version', version='1.0.0')
    parser.add_argument('rna', help='RNA bam (MUST be sorted by readID)')
    parser.add_argument('dna', help='DNA bam (MUST be sorted by readID)')
    parser.add_argument('out', help='name of out pairs file')
    parser.add_argument('chrtokeep', help='name of file with chromosome to keep')
    parser.add_argument('--nmax','-n', type=int, help='up to n reads', default=0)
    parser.add_argument('--qrna', type=int, help='Min q for RNA', default=0)
    parser.add_argument('--qdna', type=int, help='Min q for DNA', default=15)

    # parser.set_defaults(func=tag)

def pairup_simple(args):
    print('pairup_simple', args)
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE,SIG_DFL)
    bamsToPairsSimple(args.rna,args.dna,args.out, args.chrtokeep, qminRNA=args.qrna, qminDNA=args.qdna, nmax=args.nmax)
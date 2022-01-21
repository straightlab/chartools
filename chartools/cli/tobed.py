from ..convert import tobed

def add_subcommand_tobed(subparsers):
    parser = subparsers.add_parser('tobed')
    parser.set_defaults(func=tobeds)
    parser.add_argument('--version', action='version', version='1.0.0')
    parser.add_argument('pairs', help='pairs file (MUST be sorted by readID)')
    # parser.add_argument('--out', '-o', help='output file prfx')
    parser.add_argument('--rna','-R', help='rna to stdout', action="store_true")
    parser.add_argument('--dna','-D', help='dna to stdout', action="store_true")
    parser.add_argument('--nmax', help='limits parsing to that many reads')
   
    # parser.set_defaults(func=tag)

def tobeds(args):
    # print('tobed', args)
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE,SIG_DFL) #allows to do head operations without throwing errors

    to_stdout="rna"
    # outprfx=args.out
    if (args.rna and args.dna):
        # print("Asked for both RNA and DNA heavy output, switched to default RNA")
        to_stdout="rna"
        # outprfx=None
    elif args.rna:
        to_stdout="rna"
    elif args.dna:
        to_stdout="dna"

    # if (args.out is None) and (to_stdout is None):
    #     to_stdout="both"

    # print(to_stdout)
    tobed(args.pairs, to_stdout=to_stdout, nmax=0)
from ..genes import get_annotations_db


def add_subcommand_mkannotdb(subparsers):
    parser = subparsers.add_parser('make_annot_db')
    parser.set_defaults(func=make_annot_db)
    parser.add_argument('--version', action='version', version='1.0.0')
    parser.add_argument('chrNameLengthFile', help='chr file')
    parser.add_argument('txdbFile', help='tx db file')
    parser.add_argument('genedbFile', help='gene db file')

    # parser.add_argument('--outpair','-o', help='write pairs to OUTPAIRS instead of stdout')
    # parser.add_argument('--outpair_lowtrig','-l', help='write lowtrig pairs to OUTPAIRS2 instead of stdout')
    # parser.add_argument('--outpairedRNA', help='write paired RNA to bam/sam file')
    # parser.add_argument('--outpairedDNA', help='write paired DNA to bam/sam file')
    # parser.add_argument('--outunpairedRNA', help='write unpaired RNA to bam/sam file')
    # parser.add_argument('--outunpairedDNA', help='write unpaired DNA to bam/sam file')
    # parser.add_argument('--pairingmode','-m', help='pairing mode', type=str, choices=["triangular","prefix","rd"])
    # parser.add_argument('--outstats','-s', help='write statistics to file rather than STDERR')
    # parser.add_argument('--reducemmap',help='keep only one for each multimapper', type=int, choices=[0,1,2])
    # parser.add_argument('--removeannots',help='reads with annotations do not go into pairs file', action='store_true')
    # parser.add_argument('--r2',help='RNA-RNA mode', action='store_true')
    # # parser.set_defaults(func=tag)

def make_annot_db(args):
    print('make_annot_db', args)
    genes_df, tx_df, chr_dict, chr_norm = get_annotations_db(chrNameLengthFile, txdbFile, genedbFile)
    
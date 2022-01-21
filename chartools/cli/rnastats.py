from ..genes import get_annotations_db, rna_stats, add_gene_info, load_salmon, rna_stats2, rna_stats2noDNA, summarize_counts, get_annotations_db2
import pyarrow.parquet as pq
import pyarrow as pa

def add_subcommand_rnastats(subparsers):
    parser = subparsers.add_parser('genestats')
    parser.set_defaults(func=make_rna_stats)
    parser.add_argument('--version', action='version', version='1.0.0')
    #exons_stats_file, , intergene_stats_file
    parser.add_argument('exons_stats_file', help='rna stats file')
    parser.add_argument('introns_stats_file', help='rna stats file')
    parser.add_argument('intergene_stats_file', help='rna stats file')
    parser.add_argument('chrNameLengthFile', help='chr file')
    parser.add_argument('txdbFile', help='tx db file')
    parser.add_argument('genedbFile', help='gene db file')
    parser.add_argument('out_prefix', help='output file prfx')
    
   
    # parser.set_defaults(func=tag)

def make_rna_stats(args):
    # print('tobed', args)
    outprefix=args.out_prefix
    genes_df, tx_df, chr_dict, chr_norm = get_annotations_db(args.chrNameLengthFile, args.txdbFile, args.genedbFile)
    data_unq, data_all, data_intergenic = rna_stats(args.exons_stats_file, args.introns_stats_file, args.intergene_stats_file, tx_df, genes_df, chr_norm)
    data_unq2=add_gene_info(data_unq, genes_df, join='inner').sort_values(('all','NumContacts'), ascending=False)
    data_unq2.to_excel(outprefix+"nomultimap.xlsx")
    data_unq2.to_csv(outprefix+"nomultimap.csv")
    table_unq = pa.Table.from_pandas(data_unq2)
    pq.write_table(table_unq, outprefix+'nomultimap.parquet')

    data_all2=add_gene_info(data_all, genes_df, join='inner').sort_values(('all','NumContacts'), ascending=False)
    data_all2.to_excel(outprefix+"all.xlsx")
    data_all2.to_csv(outprefix+"all.csv")
    table_all = pa.Table.from_pandas(data_all2)
    pq.write_table(table_all, outprefix+'all.parquet')


def add_subcommand_rnastats2(subparsers):
    parser = subparsers.add_parser('genestats2')
    parser.set_defaults(func=make_rna_stats2)
    parser.add_argument('--version', action='version', version='1.0.0')
    #exons_stats_file, , intergene_stats_file
    parser.add_argument('stats_file', help='rna stats file')
    parser.add_argument('chrNameLengthFile', help='chr file')
    parser.add_argument('txdbFile', help='tx db file')
    parser.add_argument('genedbFile', help='gene db file')
    parser.add_argument('out_prefix', help='output file prfx')
    parser.add_argument('--nodna', '-n', action='store_true', help='RNA only')
    parser.add_argument('--tx_start', default='ENST', help='ENST?')
    parser.add_argument('--gene_start', default='ENSG', help='ENSG?')

def add_subcommand_rnastats3(subparsers):
    parser = subparsers.add_parser('genestats3')
    parser.set_defaults(func=make_rna_stats3)
    parser.add_argument('--version', action='version', version='1.0.0')
    #exons_stats_file, , intergene_stats_file
    parser.add_argument('stats_file', help='rna stats file')
    parser.add_argument('chrNameLengthFile', help='chr file')
    parser.add_argument('txdbFile', help='tx db file')
    parser.add_argument('genedbFile', help='gene db file')
    parser.add_argument('out_prefix', help='output file prfx')
    # parser.add_argument('--nodna', '-n', action='store_true', help='RNA only')
    parser.add_argument('--tx_start', default='ENST', help='ENST?')
    parser.add_argument('--gene_start', default='ENSG', help='ENSG?')

    # parser.set_defaults(func=tag)

def make_rna_stats3(args):
    # print('tobed', args)
    outprefix=args.out_prefix


    genes_df, tx_df, chr_df = get_annotations_db2(args.chrNameLengthFile, args.txdbFile, args.genedbFile)


    sufx=['multimapRNA','NOmultimapRNA']
    unqRNA=[False, True]

    for s,unq in zip(sufx, unqRNA):
        ct, ct_name, ct_intergenic = summarize_counts(args.stats_file, tx_df, genes_df, chr_df, unq=unq, tx_start=args.tx_start, gene_start=args.gene_start)
        ct=ct.sort_values(('all','N'), ascending=False)    
        ct_name=ct_name.sort_values(('all','N'), ascending=False)
        ct.to_excel(outprefix+s+".xlsx")
        ct_name.to_excel(outprefix+s+"_byname.xlsx")
        ct.to_csv(outprefix+s+".csv", sep="\t")
        ct_name.to_csv(outprefix+s+"_byname.csv", sep="\t")
        ct_tb=pa.Table.from_pandas(ct)
        ct_name_tb=pa.Table.from_pandas(ct_name)
        pq.write_table(ct_tb, outprefix+s+'.parquet')
        pq.write_table(ct_name_tb, outprefix+s+'_byname.parquet')

        if not ct_intergenic is None:
            ct_intergenic.to_csv(outprefix+s+"_intergenic.csv", sep="\t")
            ct_intergenic.to_excel(outprefix+s+"_intergenic.xlsx")
            ct_intergenic_tb=pa.Table.from_pandas(ct_intergenic)
            pq.write_table(ct_intergenic_tb, outprefix+s+'_intergenic.parquet')

def make_rna_stats2(args):
    # print('tobed', args)
    outprefix=args.out_prefix

    if args.nodna:
        genes_df, tx_df, chr_dict, chr_norm = get_annotations_db2(args.chrNameLengthFile, args.txdbFile, args.genedbFile)



        data_unq, data_all, data_intergenic = rna_stats2noDNA(args.stats_file, tx_df, genes_df, chr_norm)
        data_unq2=add_gene_info(data_unq, genes_df, join='inner').sort_values(('all','NumContacts'), ascending=False)
        data_unq2.to_excel(outprefix+"nomultimap.xlsx")
        data_unq2.to_csv(outprefix+"nomultimap.csv")
        table_unq = pa.Table.from_pandas(data_unq2)
        pq.write_table(table_unq, outprefix+'nomultimap.parquet')

        data_all2=add_gene_info(data_all, genes_df, join='inner').sort_values(('all','NumContacts'), ascending=False)
        data_all2.to_excel(outprefix+"all.xlsx")
        data_all2.to_csv(outprefix+"all.csv")
        table_all = pa.Table.from_pandas(data_all2)
        pq.write_table(table_all, outprefix+'all.parquet')
    else:
        genes_df, tx_df, chr_dict, chr_norm = get_annotations_db(args.chrNameLengthFile, args.txdbFile, args.genedbFile)
        data_unq, data_all, data_intergenic = rna_stats2(args.stats_file, tx_df, genes_df, chr_norm, tx_start=args.tx_start, gene_start=args.gene_start)
        data_unq2=add_gene_info(data_unq, genes_df, join='inner').sort_values(('all','NumContacts'), ascending=False)
        data_unq2.to_excel(outprefix+"nomultimap.xlsx")
        data_unq2.to_csv(outprefix+"nomultimap.csv")
        table_unq = pa.Table.from_pandas(data_unq2)
        pq.write_table(table_unq, outprefix+'nomultimap.parquet')

        data_all2=add_gene_info(data_all, genes_df, join='inner').sort_values(('all','NumContacts'), ascending=False)
        data_all2.to_excel(outprefix+"all.xlsx")
        data_all2.to_csv(outprefix+"all.csv")
        table_all = pa.Table.from_pandas(data_all2)
        pq.write_table(table_all, outprefix+'all.parquet')

def add_subcommand_salmonstats(subparsers):
    parser = subparsers.add_parser('salmonstats')
    parser.set_defaults(func=make_salmon_stats)
    parser.add_argument('--version', action='version', version='1.0.0')
    #exons_stats_file, , intergene_stats_file
    parser.add_argument('salmon_file', help='salmon file')
    parser.add_argument('chrNameLengthFile', help='chr file')
    parser.add_argument('txdbFile', help='tx db file')
    parser.add_argument('genedbFile', help='gene db file')
    parser.add_argument('out_prefix', help='output file prfx')
   
    # parser.set_defaults(func=tag)

def make_salmon_stats(args):
    # print('tobed', args)
    outprefix=args.out_prefix
    # genes_df, tx_df, chr_dict, chr_norm = get_annotations_db(args.chrNameLengthFile, args.txdbFile, args.genedbFile)
    genes_df, tx_df, chr_df = get_annotations_db2(args.chrNameLengthFile, args.txdbFile, args.genedbFile)
    data=load_salmon(args.salmon_file, tx_df)
    data2=add_gene_info(data, genes_df, join='inner').sort_values('TPM', ascending=False)
    data2.to_excel(outprefix+".xlsx")
    data2.to_csv(outprefix+".csv", sep="\t")
    table_unq = pa.Table.from_pandas(data2)
    pq.write_table(table_unq, outprefix+'.parquet')

    
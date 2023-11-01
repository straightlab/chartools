import os
import pandas as pd
import subprocess
import io
import numpy as np
import re

import pyarrow.parquet as pq
import pyarrow as pa

#TODO: move to separate module
#TODO: make annotation class for portability

def get_annotations_db(chrNameLengthFile, genedbFile):
    genes_df=pd.read_csv(genedbFile, header=None, sep='\t', names=['ENSG','type','name','chr','strand','position','L'], index_col=0, usecols=['ENSG','type','name','chr','strand','L'])
    chr_df=pd.read_csv(chrNameLengthFile, sep='\t', header=None, names=['chr','L'], index_col=0)


    return genes_df, chr_df

def parse_clumpify_report (file):
    with io.open(file) as ifile:
        raw_data=ifile.read()
        regexes = {
            'Ninput':                  r"Reads In:\s+(\d+)",
            'Nout':        r"Reads Out:\s+(\d+)",
#             'time':        r"Total time:\s+(\d+\)\s+.+",
        }
        parsed_data = {}
        for k, r in regexes.items():
            r_search = re.search(r, raw_data, re.MULTILINE)
            if r_search:
                parsed_data[k] = int(r_search.group(1))
    return parsed_data

def get_stats_duplication(sample_root_folder):        
    fname=os.path.join(sample_root_folder,'chimeras/deduped/deduping.log')
    parsed_data=parse_clumpify_report(fname)
    Nin=parsed_data['Ninput']/2
    Nout=parsed_data['Nout']/2
    stats={'step input': Nin, 'step output': Nout, 'discard at step %': 100-Nout/Nin*100}
    df=pd.DataFrame.from_dict(stats,orient='index').reindex(["step input","discard at step %","step output"])
    df.index.names=['stats']
    return df 

def parse_trimmomatic_report(file):
    with io.open(file) as ifile:
        raw_data=ifile.read()
        regexes = {
            'Ninput':    r"Input Read Pairs: (\d+)\s+",
            'N2':        r"Both Surviving: (\d+)\s+",
            'NF':        r"Forward Only Surviving: (\d+)\s+",
            'NR':        r"Reverse Only Surviving: (\d+)\s+",
            'Ndropped':  r"Dropped: (\d+)\s+"
#             'time':        r"Total time:\s+(\d+\)\s+.+",
        }
        parsed_data = {}
        for k, r in regexes.items():
            r_search = re.search(r, raw_data, re.MULTILINE)
            if r_search:
                parsed_data[k] = int(r_search.group(1))
    return parsed_data

def get_stats_trimming(sample_root_folder):        
    fname=os.path.join(sample_root_folder,'chimeras/deduped.trimmed/trimming.log')
    parsed_data=parse_trimmomatic_report(fname)
    Nin=parsed_data['Ninput']
    Nout=Nin-parsed_data['Ndropped']

    stats={'step input': Nin, 'step output': Nout, 'discard at step %': 100-Nout/Nin*100, '%F':parsed_data['NF']/Nout*100, '%R':parsed_data['NR']/Nout*100, '%FR':parsed_data['N2']/Nout*100}
    df=pd.DataFrame.from_dict(stats,orient='index').reindex(["step input",'%F','%R','%FR', "discard at step %","step output"])
    df.index.names=['stats']
    return df 

def parse_bridge_count (file, sectionHeader=">>Count matrix"):
    with io.open(file) as ifile:
        raw_data=ifile.read()
        regexes = {
            'counts':                  r"^"+sectionHeader+".*\n(.+)\n"
        }
        parsed_data = {}
        for k, r in regexes.items():
            r_search = re.search(r, raw_data, re.MULTILINE)
            if r_search:
                parsed_data[k] = r_search.group(1)
    return parsed_data

def get_bridge_counts(fname):
    parsed_data=parse_bridge_count(fname,">>Count matrix")
    # M=np.fromstring(parsed_data['counts'][1:-1], sep=';', dtype=int)
    M=np.array(np.matrix(parsed_data['counts']))
    return M

def get_COMBINEDstats_trim_pear(sample_root_folder, fname={'se2':'split_chimeras/_debridging/mergedPairs_pear/bridge.stats.txt', 'se1':'split_chimeras/_debridging/readthrough/bridge.stats.txt', 'pe':'split_chimeras/_debridging/unmergedPairs_pear/bridge.stats.txt'}):
    df_trim=get_stats_trimming(sample_root_folder)
    Nin=df_trim[0].at['step input']
    # use bridge stats to get pear out
    M_se2=get_bridge_counts(os.path.join(sample_root_folder,fname['se2']))
    M_se1=get_bridge_counts(os.path.join(sample_root_folder,fname['se1']))
    M_pe=get_bridge_counts(os.path.join(sample_root_folder,fname['pe']))

    N_merged=np.sum(M_se2)
    N_readthrough=np.sum(M_se1)
    N_nomerge=np.sum(M_pe)
    Nout=N_merged+N_readthrough+N_nomerge

    stats={'step input': Nin, 'step output': Nout, 'discard at step %': 100-Nout/Nin*100}
    df=pd.DataFrame.from_dict(stats,orient='index').reindex(["step input","discard at step %","step output"])
    df.index.names=['stats']
    return df 

def get_stats_bridge(sample_root_folder, fname={'se2':'split_chimeras/_debridging/mergedPairs_pear/bridge.stats.txt', 'se1':'split_chimeras/_debridging/readthrough/bridge.stats.txt', 'pe':'split_chimeras/_debridging/unmergedPairs_pear/bridge.stats.txt'}):
    M_se2=get_bridge_counts(os.path.join(sample_root_folder,fname['se2'])).flatten()
    M_se1=get_bridge_counts(os.path.join(sample_root_folder,fname['se1'])).flatten()
    M_pe=get_bridge_counts(os.path.join(sample_root_folder,fname['pe']))
    N=np.sum(M_pe)+np.sum(M_se1)+np.sum(M_se2)

    Nambiguous=(M_pe[1,2]+M_pe[2,1])
    Ndiscordant=M_pe[0,3]+M_pe[1,3]+M_pe[2,3]+M_pe[3,3]+M_pe[3,0]+M_pe[3,1]+M_pe[3,2]+M_pe[1,1]+M_pe[2,2]+M_se1[3]+M_se2[3]
#     N_discordant_nomerge=M_pe[0,3]+M_pe[1,3]+M_pe[2,3]+M_pe[3,3]+M_pe[3,0]+M_pe[3,1]+M_pe[3,2]+M_pe[2,1]+M_pe[1,1]+M_pe[2,2]
    
    Nconcordant=(M_pe[1,0]+M_pe[0,1]+M_pe[2,0]+M_pe[0,2]+M_se1[1]+M_se1[2]+M_se2[1]+M_se2[2]) #+M_pe[1,2] note that this is now discordant because 
    Nnobridge=M_pe[0,0]+M_se1[0]+M_se2[0]
    
    stats={'step input':N, '%NObridge':Nnobridge/N*100, '%discordant':Ndiscordant/N*100, '%ambiguous':Nambiguous/N*100,'%concordant':Nconcordant/N*100,  'step output':Nconcordant, 'discard at step %': 100-(Nconcordant/N*100)}
    
    df=pd.DataFrame.from_dict(stats,orient='index').reindex(["step input","%NObridge","%discordant","%ambiguous","%concordant","discard at step %","step output"])

    df.index.names=['stats']
    return df



def get_stats_filtering(sample_root_folder, mates_pairing_mode="SE_merge_pear"):
        
    
    df_T=pd.read_csv(os.path.join(sample_root_folder,"split_chimeras", mates_pairing_mode, 'filtering.stats.txt'),delimiter=',', dtype=int).transpose()

    N=df_T.at["total",0]
    RNAonlyShort=(df_T.at["shortR",0])/N*100 
    DNAonlyShort=(df_T.at["shortD",0])/N*100
    RNAandDNAShort=(df_T.at["shortRD",0])/N*100
    final=df_T.at["long",0]
    
    stats={'step input': N, '%RNA too short':RNAonlyShort , '%DNA too short':DNAonlyShort, '%RD too short': RNAandDNAShort, 'discard at step %': 100-final/N*100, 'step output':final}
    
    df=pd.DataFrame.from_dict(stats,orient='index').reindex(["step input","%RNA too short","%DNA too short","%RD too short","discard at step %","step output"])

    df.index.names=['stats']
    return df


def parse_bowtie(fname):
    f=open(fname,'r')
    regexes = {
                'unpaired': {
                    'unpaired_aligned_none': r"(\d+) \([\d\.]+%\) aligned 0 times",
                    'unpaired_aligned_one': r"(\d+) \([\d\.]+%\) aligned exactly 1 time",
                    'unpaired_aligned_multi': r"(\d+) \([\d\.]+%\) aligned >1 times"
                },
                'paired': {
                    'paired_aligned_none': r"(\d+) \([\d\.]+%\) aligned concordantly 0 times",
                    'paired_aligned_one': r"(\d+) \([\d\.]+%\) aligned concordantly exactly 1 time",
                    'paired_aligned_multi': r"(\d+) \([\d\.]+%\) aligned concordantly >1 times",
                    'paired_aligned_discord_one': r"(\d+) \([\d\.]+%\) aligned discordantly 1 time",
                    'paired_aligned_discord_multi': r"(\d+) \([\d\.]+%\) aligned discordantly >1 times",
                    'paired_aligned_mate_one': r"(\d+) \([\d\.]+%\) aligned exactly 1 time",
                    'paired_aligned_mate_multi': r"(\d+) \([\d\.]+%\) aligned >1 times",
                    'paired_aligned_mate_none': r"(\d+) \([\d\.]+%\) aligned 0 times"
                }
            }

    # Go through log file line by line
    parsed_data = {}

    for l in f:
        # Attempt in vain to find original bowtie2 command, logged by another program

        # Total reads
        total = re.search(r"(\d+) reads; of these:", l)
        if total:
            parsed_data['total_reads'] = int(total.group(1))

        # Single end reads
        unpaired = re.search(r"(\d+) \([\d\.]+%\) were unpaired; of these:", l)
        if unpaired:
            parsed_data['unpaired_total'] = int(unpaired.group(1))

            # Do nested loop whilst we have this level of indentation
            l = f.readline()
            while l.startswith('    '):
                for k, r in regexes['unpaired'].items():
                    match = re.search(r, l)
                    if match:
                        parsed_data[k] = int(match.group(1))
                l = f.readline()

        # Paired end reads
        paired = re.search(r"(\d+) \([\d\.]+%\) were paired; of these:", l)
        if paired:
            parsed_data['paired_total'] = int(paired.group(1))


            # Do nested loop whilst we have this level of indentation
            l = f.readline()
            while l.startswith('    '):
                for k, r in regexes['paired'].items():
                    match = re.search(r, l)
                    if match:
                        parsed_data[k] = int(match.group(1))
                l = f['f'].readline()

        # Overall alignment rate
        overall = re.search(r"([\d\.]+)% overall alignment rate", l)
        if overall:
            parsed_data['overall_alignment_rate'] = float(overall.group(1))

            # End of log section
            # Save half 'pairs' of mate counts
            m_keys = ['paired_aligned_mate_multi', 'paired_aligned_mate_none', 'paired_aligned_mate_one']
            for k in m_keys:
                if k in parsed_data:
                    parsed_data['{}_halved'.format(k)] = float(parsed_data[k]) / 2.0
            # Save parsed data

    f.close()
    return parsed_data

def get_stats_decon(sample_root_folder,  mates_paring_mode="SE_merge_pear"):
    parsed_data=parse_bowtie(os.path.join(sample_root_folder, 'split_chimeras',mates_paring_mode,'long.contaminants/rna.stats.txt'))
    N=parsed_data['unpaired_aligned_none']+parsed_data['unpaired_aligned_one']+parsed_data['unpaired_aligned_multi']
    NrRNA=N-parsed_data['unpaired_aligned_none']
    final=N-NrRNA
    stats={'step input': N, '%rRNA':NrRNA/max(1,N)*100, '%rRNA %single align': parsed_data['unpaired_aligned_one']/max(1,NrRNA), 'discard at step %': NrRNA/max(1,N)*100, 'step output':final}
    
    df=pd.DataFrame.from_dict(stats,orient='index').reindex(["step input","%rRNA","%rRNA %single align","discard at step %","step output"])
    df.index.names=['stats']

    return df

def parse_star_report (file):
    with io.open(file) as ifile:
        raw_data=ifile.read()
        regexes = {
            'total_reads':                  r"Number of input reads \|\s+(\d+)",
            'avg_input_read_length':        r"Average input read length \|\s+([\d\.]+)",
            'uniquely_mapped':              r"Uniquely mapped reads number \|\s+(\d+)",
            'uniquely_mapped_percent':      r"Uniquely mapped reads % \|\s+([\d\.]+)",
            'avg_mapped_read_length':       r"Average mapped length \|\s+([\d\.]+)",
            'num_splices':                  r"Number of splices: Total \|\s+(\d+)",
            'num_annotated_splices':        r"Number of splices: Annotated \(sjdb\) \|\s+(\d+)",
            'num_GTAG_splices':             r"Number of splices: GT/AG \|\s+(\d+)",
            'num_GCAG_splices':             r"Number of splices: GC/AG \|\s+(\d+)",
            'num_ATAC_splices':             r"Number of splices: AT/AC \|\s+(\d+)",
            'num_noncanonical_splices':     r"Number of splices: Non-canonical \|\s+(\d+)",
            'mismatch_rate':                r"Mismatch rate per base, % \|\s+([\d\.]+)",
            'deletion_rate':                r"Deletion rate per base \|\s+([\d\.]+)",
            'deletion_length':              r"Deletion average length \|\s+([\d\.]+)",
            'insertion_rate':               r"Insertion rate per base \|\s+([\d\.]+)",
            'insertion_length':             r"Insertion average length \|\s+([\d\.]+)",
            'multimapped':                  r"Number of reads mapped to multiple loci \|\s+(\d+)",
            'multimapped_percent':          r"% of reads mapped to multiple loci \|\s+([\d\.]+)",
            'multimapped_toomany':          r"Number of reads mapped to too many loci \|\s+(\d+)",
            'multimapped_toomany_percent':  r"% of reads mapped to too many loci \|\s+([\d\.]+)",
            'unmapped_mismatches_percent':  r"% of reads unmapped: too many mismatches \|\s+([\d\.]+)",
            'unmapped_tooshort_percent':    r"% of reads unmapped: too short \|\s+([\d\.]+)",
            'unmapped_other_percent':       r"% of reads unmapped: other \|\s+([\d\.]+)",
        }
        parsed_data = {}
        for k, r in regexes.items():
            r_search = re.search(r, raw_data, re.MULTILINE)
            if r_search:
                parsed_data[k] = float(r_search.group(1))
        # Figure out the numbers for unmapped as for some reason only the percentages are given
        try:
            total_mapped = parsed_data['uniquely_mapped'] + parsed_data['multimapped'] + parsed_data['multimapped_toomany']
            unmapped_count = parsed_data['total_reads'] - total_mapped
            total_unmapped_percent = parsed_data['unmapped_mismatches_percent'] + parsed_data['unmapped_tooshort_percent'] + parsed_data['unmapped_other_percent']
            try:
                parsed_data['unmapped_mismatches'] = int(round(unmapped_count * (parsed_data['unmapped_mismatches_percent'] / total_unmapped_percent), 0))
                parsed_data['unmapped_tooshort'] = int(round(unmapped_count * (parsed_data['unmapped_tooshort_percent'] / total_unmapped_percent), 0))
                parsed_data['unmapped_other'] = int(round(unmapped_count * (parsed_data['unmapped_other_percent'] / total_unmapped_percent), 0))
            except ZeroDivisionError:
                parsed_data['unmapped_mismatches'] = 0
                parsed_data['unmapped_tooshort'] = 0
                parsed_data['unmapped_other'] = 0
        except KeyError:
            pass

        if len(parsed_data) == 0: return None
    return parsed_data


def get_stats_alignment_star(sample_root_folder, rna_alignment_mode='star_gencodeV29', allow_multimap=True):
    parsed_data=parse_star_report(os.path.join(sample_root_folder, 'alignments/rna/',rna_alignment_mode,'rna.Log.final.out'))
    total=pd.DataFrame.from_dict(parsed_data,orient='index')
    
    N_multimapped=total.loc['multimapped']
    N_uniquely_mapped=total.loc['uniquely_mapped']
    N_unmapped_mismatches=total.loc['unmapped_mismatches']
    N_unmapped_other=total.loc['unmapped_other']
    N_unmapped_tooshort=total.loc['unmapped_tooshort']
    N_multimapped_toomany=total.loc['multimapped_toomany']

    
    N=total.loc['total_reads']
    
    if allow_multimap:
        N_final=N_multimapped+N_uniquely_mapped
    else:
        N_final=N_uniquely_mapped

    stats={'step input': N, '%Uniquely mapped':N_uniquely_mapped/N*100, '%Multimapped':N_multimapped/N*100, '%Unmapped_tooshort':N_unmapped_tooshort/N*100,  '%Multimapped_toomany':N_multimapped_toomany/N*100, '%Unmapped_mismatches': N_unmapped_mismatches/N*100,'% ':N_unmapped_other/N*100,'avg_length': total.at['avg_mapped_read_length'], 'N annottated splices':total.at['num_annotated_splices'], 'N unannottated splices':total.at['num_noncanonical_splices'], 'discard at step %': (N-N_final)/N*100, 'step output':N_final}
    
    df=pd.DataFrame.from_dict(stats,orient='index').reindex(["step input","%Uniquely mapped",'%Multimapped','%Multimapped_toomany','%Unmapped_tooshort','%Unmapped_mismatches','%Unmapped_other','N annottated splices','N unannottated splices', "discard at step %","step output"])
 
    df.index.names=['stats']
    return df 

def get_stats_alignment_DNA(sample_root_folder, dna_alignment_mode='bowtie_hg38', allow_multimap=True):
    parsed_data=parse_bowtie(os.path.join(sample_root_folder, 'alignments/dna', dna_alignment_mode, 'dna.stats.txt'))
    total=pd.DataFrame.from_dict(parsed_data,orient='index')

    N=total.loc['total_reads']
    N_aligned_one=total.loc['unpaired_aligned_one']
    N_aligned_multi=total.loc['unpaired_aligned_multi']
    N_aligned_none=total.loc['unpaired_aligned_none']

    if allow_multimap:
        final=N-N_aligned_none
    else:
        final=N-N_aligned_none-N_aligned_multi
    stats={'step input': N, '%Uniquely mapped': N_aligned_one/N*100, '%Multimapped': N_aligned_multi/N*100, '%Unmapped': N_aligned_none/N*100, 'discard at step %': (N-final)/N*100, 'step output':final}
    df=pd.DataFrame.from_dict(stats,orient='index').reindex(["step input","%Uniquely mapped",'%Multimapped','%Unmapped',"discard at step %","step output"])

    df.index.names=['stats']
    return df


def parse_pairing_reports(sample_root_folder, pairing_mode='gencondeV29_hg38', rna_alignment_mode='star_gencodeV29'):
    ptypes=['exons','introns','intergenic','novel']
    pairing_mode='gencondeV29_hg38'
    df=pd.concat([pd.read_csv(os.path.join(sample_root_folder,'pairs',pairing_mode,pairing_type,'rd.stats.txt'), sep=',', usecols=[0,1]) for pairing_type in ptypes], axis=0)
    df.index=ptypes

    df_noalign=pd.read_csv(os.path.join(sample_root_folder,'alignments/rna',rna_alignment_mode,'noalign.stats.txt'), sep=',').transpose()
    df_noalign.columns=['nRNAsOnly']
    df_noalign['n']=np.zeros((2,1), dtype=int)

    df=pd.concat([df, df_noalign], axis=0, sort=False)

    df['nRNAsOnly']=df.sum(axis=1)
    df.columns=['paired','total']
    return df

def get_stats_alignment_RNA(sample_root_folder, pairing_mode='gencondeV29_hg38', rna_alignment_mode='star_gencodeV29'):
    data=parse_pairing_reports(sample_root_folder, pairing_mode=pairing_mode)['total']
    N=data.sum()
    Nmapped=data.loc[['exons','introns','intergenic','novel']].sum()
    Nexons=data.at['exons']

    stats={'step input': N, '%Mapped': Nmapped/N*100, '%TooMany':  data.at['toomany']/N*100, '%Unmapped':  data.at['unmapped']/N*100, '%Exons': data.at['exons']/Nmapped*100, '%Introns': data.at['introns']/Nmapped*100, '%Intergenic': data.at['intergenic']/Nmapped*100, '%Ambiguous': data.at['novel']/Nmapped*100, 'discard at step %': (N-Nmapped)/N*100, 'step output':Nmapped}


    stardata=parse_star_report(os.path.join(sample_root_folder,'alignments/rna', rna_alignment_mode, 'rna.Log.final.out'))

    stats['Annotated Splice Rate']=stardata['num_annotated_splices']/Nexons
    stats['Non Canonical Splice Rate']=stardata['num_noncanonical_splices']/Nexons

    df=pd.DataFrame.from_dict(stats,orient='index').reindex(["step input","%Mapped",'%TooMany','%Unmapped','%Exons', '%Introns', '%Intergenic', '%Ambiguous', 'Annotated Splice Rate', 'Non Canonical Splice Rate',"discard at step %","step output"])
 
    df.index.names=['stats']

    return df

def get_stats_pairing(sample_root_folder, pairing_mode='gencondeV29_hg38'):
    data=parse_pairing_reports(sample_root_folder, pairing_mode=pairing_mode)
    N=data.loc[['exons','introns','intergenic','novel'],'total'].sum()
    Npaired=data['paired'].sum()

    stats={'step input': N, '%Paired': Npaired/N*100, '%Exons': data.at['exons','paired']/Npaired*100, '%Introns': data.at['introns','paired']/Npaired*100, '%Intergenic': data.at['intergenic','paired']/Npaired*100, '%Ambiguous': data.at['novel','paired']/Npaired*100, 'discard at step %': (N-Npaired)/N*100, 'step output':Npaired}


    df=pd.DataFrame.from_dict(stats,orient='index').reindex(["step input","%Paired",'%Exons', '%Introns', '%Intergenic', '%Ambiguous', "discard at step %","step output"])
 
    df.index.names=['stats']

    return df

def get_stats_annotating(sample_root_folder, pairing_mode='gencondeV29_hg38'): #just remove the intergenic and ambiguous stuff
    df_all=pd.read_parquet(os.path.join(sample_root_folder, 'pairs', pairing_mode, 'all/stats/gene_expression_multimapRNA.parquet')).xs('N', axis=1, level=1).sum()
    #N_intergenic=pd.read_parquet(os.path.join(sample_root_folder, 'pairs', pairing_mode, 'all/stats/gene_expression_multimapRNA_intergenic.parquet')).iat[0,0]

    data2=parse_pairing_reports(sample_root_folder, pairing_mode=pairing_mode)
    N_intergenic = data2.at['intergenic','paired']
    #N=data.loc[['exons','introns','intergenic','novel'],'total'].sum()
    #Npaired=data['paired'].sum()

    Nin=df_all['all']+N_intergenic

    stats={'step input': Nin, 'step output': Nin-N_intergenic, 'discard at step %': 100-(Nin-N_intergenic)/Nin*100}
    df=pd.DataFrame.from_dict(stats,orient='index').reindex(["step input","discard at step %","step output"])
    df.index.names=['stats']
    return df 

def get_stats_unqmappers_filtering(sample_root_folder, pairing_mode='gencondeV29_hg38'):
    df_all=pd.read_parquet(os.path.join(sample_root_folder, 'pairs', pairing_mode, 'all/stats/gene_expression_multimapRNA.parquet')).xs('N', axis=1, level=1).sum()
    df_nommap=pd.read_parquet(os.path.join(sample_root_folder, 'pairs', pairing_mode, 'all/stats/gene_expression_NOmultimapRNA.parquet')).xs('N', axis=1, level=1).sum()

    Nin=df_all['all']
    Nout=df_nommap['all']
    stats={'step input': Nin, 'step output': Nout, '%Exons': df_nommap['exons']/Nout*100, '%Introns': df_nommap['introns']/Nout*100, 'discard at step %': 100-Nout/Nin*100}
    df=pd.DataFrame.from_dict(stats,orient='index').reindex(["step input",'%Exons', '%Introns', "discard at step %","step output"])
    df.index.names=['stats']
    return df 
    

def get_all_stats(sample_root_folder, pairing_mode='gencondeV29_hg38', rna_alignment_mode='star_gencodeV29', dna_alignment_mode='hg38_unmasked', mates_pairing_mode='SE_merge'):
    stats2get={
        'duplication':get_stats_duplication,
        'trim and splice':get_COMBINEDstats_trim_pear,
        'bridge':get_stats_bridge,
        'length filter':get_stats_filtering,
        'rRNA decon': get_stats_decon, 
        'DNA alignment' :get_stats_alignment_DNA,
        'RNA alignment' :get_stats_alignment_RNA,
        'pairing': get_stats_pairing,
        'annotating': get_stats_annotating,
        'unique annotating': get_stats_unqmappers_filtering
    }
    
    stats2get_order=['duplication','trim and splice','bridge','length filter','rRNA decon','DNA alignment','RNA alignment','pairing', 'annotating', 'unique annotating']
    # allstats_array=[]
    

    allstats=[]
    N=stats2get['duplication'](sample_root_folder).at['step input',0]
    for this_stat in stats2get_order:
        this_df=stats2get[this_stat](sample_root_folder)

        remain_in=pd.DataFrame([this_df.at['step output',0]/N*100], index=['%library remains'])
        this_df=this_df.append(remain_in, sort=False)

        allstats+=[this_df]

    df=pd.concat(allstats, axis=0, keys=stats2get_order) #.fillna(-1)
    
    df.index.names=['step','stat']
    
    return df

def get_all_stats_simple(sample_root_folder, pairing_mode='gencondeV29_hg38', rna_alignment_mode='star_gencodeV29', dna_alignment_mode='hg38_unmasked', mates_pairing_mode='SE_merge'):
    stats2get={
        'duplication':get_stats_duplication,
        'trim and splice':get_COMBINEDstats_trim_pear,
        'bridge':get_stats_bridge,
        'length filter':get_stats_filtering,
        'rRNA decon': get_stats_decon, 
        'DNA alignment' :get_stats_alignment_DNA,
        'RNA alignment' :get_stats_alignment_RNA,
        'pairing': get_stats_pairing
    }
    
    stats2get_order=['duplication','trim and splice','bridge','length filter','rRNA decon','DNA alignment','RNA alignment','pairing']
    # allstats_array=[]
    

    allstats=[]
    N=stats2get['duplication'](sample_root_folder).at['step input',0]
    for this_stat in stats2get_order:
        this_df=stats2get[this_stat](sample_root_folder)

        remain_in=pd.DataFrame([this_df.at['step output',0]/N*100], index=['%library remains'])
        this_df=this_df.append(remain_in, sort=False)

        allstats+=[this_df]

    df=pd.concat(allstats, axis=0, keys=stats2get_order) #.fillna(-1)
    
    df.index.names=['step','stat']
    
    return df

def get_biotype_census(sample_root_folder, byN='CPKM', pairing_mode='gencondeV29_hg38', types_coarsening_def='/home/groups/astraigh/CHARSEQ-analysis/analysis/ESCs_nova/CL/mygenetypes4.csv'):
    df_nommap=pd.read_parquet(os.path.join(sample_root_folder, 'pairs', pairing_mode, 'all/stats/gene_expression_NOmultimapRNA.parquet'))

    if not types_coarsening_def is None:
        types_coarsening=pd.read_csv(types_coarsening_def,index_col=0)
        alltypes=list(set(df_nommap.index.get_level_values(2)))
        alltypes_df=pd.DataFrame(np.array(['other' for _ in range(len(alltypes))]), index=alltypes)

        # types_coarsening=types_coarsening.reindex(list(alltypes_df.index), fill_value='other')
        types_coarsening2=types_coarsening.loc[types_coarsening.index.isin(alltypes_df.index),:]
        alltypes_df.loc[types_coarsening2.index,0]=types_coarsening2['type2'].values
        alltypes_df.index.name='type'
        alltypes_df.columns=['type2']
        types_coarsening=alltypes_df

        census=pd.concat([types_coarsening,df_nommap.xs(byN, axis=1, level=1).groupby('type').sum()], axis=1, sort=False).fillna(0).groupby('type2').sum().sort_values('all',ascending=False)
        census.index.name='type'
    else:
        census=df_nommap.xs(byN, axis=1, level=1).groupby('type').sum().sort_values('all',ascending=False)

    census=census/census.sum(axis=0)*100
    
    return census

def get_stability_ratios(sample_root_folder, pairing_mode='gencondeV29_hg38', types_coarsening_def='/home/groups/astraigh/CHARSEQ-analysis/analysis/ESCs_nova/CL/mygenetypes4.csv'):
    df_nommap=pd.read_parquet(os.path.join(sample_root_folder, 'pairs', pairing_mode, 'all/stats/gene_expression_NOmultimapRNA.parquet'))

    if not types_coarsening_def is None:
        types_coarsening=pd.read_csv(types_coarsening_def,index_col=0)
        N_cis=pd.concat([types_coarsening, df_nommap.xs('N_cis', axis=1, level=1).groupby('type').sum()], axis=1, sort=False).groupby('type2').sum()
        N_trans=pd.concat([types_coarsening,df_nommap.xs('N_trans', axis=1, level=1).groupby('type').sum()], axis=1, sort=False).groupby('type2').sum()
        
        N_total=N_trans+N_cis.reindex(N_trans.index)
        
        ratios=(N_trans/N_total*100).sort_values('all', ascending=False)
        ratios.index.name='type'
    else:
        N_cis=df_nommap.xs('N_cis', axis=1, level=1).groupby('type').sum()
        N_trans=df_nommap.xs('N_cis', axis=1, level=1).groupby('type').sum()
        
        N_total=N_trans+N_cis.reindex(N_trans.index)
        
        ratios=(N_trans/N_total*100).sort_values('all', ascending=False)

    return ratios



def get_tscores_bytype(sample_root_folder, pairing_mode='gencondeV29_hg38', types_coarsening_def='/home/groups/astraigh/CHARSEQ-analysis/analysis/ESCs_nova/CL/mygenetypes4.csv', chrNameLengthFile='/oak/stanford/groups/astraigh/charseq2.0/genomes/hsapiens/grch38_foralign/chrNameLength_ALL.txt', genedbFile='/oak/stanford/groups/astraigh/charseq2.0/genomes/hsapiens/grch38_foralign/gene_table_complete_NOPARY_withLength.txt', drop_chr=['chrM','chrY']):

    genes_df, chr_df = get_annotations_db(chrNameLengthFile, genedbFile)
    alldata=pd.read_parquet(os.path.join(sample_root_folder, 'pairs', pairing_mode, 'all/stats/gene_expression_NOmultimapRNA.parquet'))

    chr_len_df=chr_df.copy()   
    chr_len_df['Ltrans']=chr_len_df['L'].sum()-chr_len_df['L']
    if not drop_chr is None:
        chr_len_df['Ltrans']=chr_len_df['Ltrans']-chr_len_df['L'].loc[drop_chr].sum()

    chr_len_cis_dict=chr_len_df['L'].to_dict()
    chr_len_trans_dict=chr_len_df['Ltrans'].to_dict()

    
    types_coarsening=pd.read_csv(types_coarsening_def,index_col=0)
    data_to_combine=[]
    for k in ['all','exons','introns']:
        data=pd.concat([alldata.xs(k, level=0, axis=1).loc[:,['N_cis','N_trans']].reset_index().set_index('ENSN'),genes_df['chr']], axis=1, join='inner').drop('name',axis=1)

        data['cis_density']=data['N_cis']/(data['chr'].map(chr_len_cis_dict))
        data['trans_density']=data['N_trans']/(data['chr'].map(chr_len_trans_dict))

        if not drop_chr is None:
            data=data.loc[~(data['chr'].isin(drop_chr))]

        data=data.groupby(['type']).sum()

        if not types_coarsening_def is None:
            data=pd.concat([types_coarsening, data], axis=1, join='inner').groupby('type2').sum()
            data.index.name='type'

        data['tscore']=data['trans_density']/data['cis_density']

        data_to_combine+=[data['tscore']]

    data=pd.concat(data_to_combine, axis=1, keys=['all','exons','introns']).sort_values('all',ascending=False)
    return data

def get_stability_summary(sample_root_folder, pairing_mode='gencondeV29_hg38', types_coarsening_def='/home/groups/astraigh/CHARSEQ-analysis/analysis/ESCs_nova/CL/mygenetypes4.csv', chrNameLengthFile='/oak/stanford/groups/astraigh/charseq2.0/genomes/hsapiens/grch38_foralign/chrNameLength_ALL.txt', genedbFile='/oak/stanford/groups/astraigh/charseq2.0/genomes/hsapiens/grch38_foralign/gene_table_complete_NOPARY_withLength.txt'):

    ratios=get_stability_ratios(sample_root_folder, pairing_mode=pairing_mode, types_coarsening_def=types_coarsening_def)
    tscores=get_tscores_bytype(sample_root_folder, pairing_mode=pairing_mode, types_coarsening_def=types_coarsening_def, chrNameLengthFile=chrNameLengthFile, genedbFile=genedbFile)

    df=pd.concat([ratios,tscores], axis=1, sort=False, keys=['TRANS%','tscore'])
    df.index.name='type'
    return df


def get_top_hits(sample_root_folder, pairing_mode='gencondeV29_hg38', types_coarsening_def='/home/groups/astraigh/CHARSEQ-analysis/analysis/ESCs_nova/CL/mygenetypes4.csv', n=5, nALL=10, byN='CPKM'):
    counts=pd.read_parquet(os.path.join(sample_root_folder, 'pairs', pairing_mode, 'all/stats/gene_expression_NOmultimapRNA_byname.parquet'))[('all',byN)]

    if not types_coarsening_def is None:
        types_coarsening=pd.read_csv(types_coarsening_def,index_col=0)
        counts.index=pd.MultiIndex.from_arrays([counts.index.get_level_values(0), types_coarsening['type2'].reindex(counts.index.get_level_values(1)).values])
        counts.index.names=['name','type']
        counts.name=byN
        # counts=counts.reset_index().set_index('type')
        
        
    # else:
    #     counts_gp=counts.groupby('type').apply(lambda x: x.sort_values(ascending=False).head(n))
    #     counts_gp.index=counts_gp.index.droplevel(2)
    counts_gp=pd.DataFrame(counts).groupby('type').apply(lambda x: pd.concat([pd.DataFrame(np.arange(1,n+1), columns=['rank'], dtype=int), x.sort_values(byN,ascending=False).head(n).reset_index()], axis=1)).drop('type', axis=1)
    counts_gp.index=counts_gp.index.droplevel(1)
    counts_gp=counts_gp.reset_index().set_index(['type','rank'])


    # counts_gp.name=byN
    counts_global=counts.sort_values(ascending=False).head(nALL)
    counts_global.name=byN
    counts_global=pd.concat([pd.DataFrame(np.arange(1,nALL+1), columns=['rank']), pd.DataFrame(counts_global).reset_index()], axis=1).set_index('rank')
    
    
    return counts_gp, counts_global


def run_qc(sample_root_folder, savefolder=""):
    if (not savefolder is None) and (len(savefolder)==0):
        savefolder='pairs/gencondeV29_hg38/all/QC'
    top_hits_bytype, top_hits_global=get_top_hits(sample_root_folder,nALL=1000, n=10, byN='CPKM')
    top_hits_bytype_CPM, top_hits_global_CPM=get_top_hits(sample_root_folder,nALL=1000, n=10, byN='CPM')
    eff=get_all_stats(sample_root_folder)
    eff_loss=eff.xs('discard at step %', axis=0, level=1)
    eff_remains=pd.DataFrame([eff.at[('duplication','step input'),0]], index=['Nreads input'])
    eff_remains=eff_remains.append(eff.xs('%library remains', axis=0, level=1))
    eff_remains=eff_remains.append(pd.DataFrame([eff.at[('unique annotating','step output'),0]], index=['Nreads output']))
    qc_data={
        'efficiency_complete_stats':eff,
        'efficiency_LOSS':eff_loss,
        'efficiency_REMAINS':eff_remains,
        'biotype_census_CPM':get_biotype_census(sample_root_folder, byN='CPM'),
        'biotype_census_CPKM':get_biotype_census(sample_root_folder, byN='CPKM'),
        'trans_statistics':get_stability_summary(sample_root_folder),
        'top_hits_bytype_CPKM':top_hits_bytype,
        'top_hits_CPKM':top_hits_global,
        'top_hits_bytype_CPM':top_hits_bytype_CPM,
        'top_hits_CPM':top_hits_global_CPM
    }
    if not savefolder is None:
        with pd.ExcelWriter(os.path.join(sample_root_folder,savefolder,'qc.xlsx')) as writer:
            for k, v in qc_data.items():
                v.to_excel(writer, sheet_name=k)
        for k, v in qc_data.items():
            v.to_csv(os.path.join(sample_root_folder,savefolder,k+'.csv'), sep='\t', header=True)
            print("Writing %s QC for sample %s"%(k,sample_root_folder))
            pq.write_table(pa.Table.from_pandas(v), os.path.join(sample_root_folder,savefolder,'_parquet',k+'.parquet'))
    return qc_data

def run_qc_simple(sample_root_folder, savefolder=""):
    if (not savefolder is None) and (len(savefolder)==0):
        savefolder='pairs/gencondeV29_hg38/all/QC'
   
    eff=get_all_stats_simple(sample_root_folder)
    eff_loss=eff.xs('discard at step %', axis=0, level=1)
    eff_remains=pd.DataFrame([eff.at[('duplication','step input'),0]], index=['Nreads input'])
    eff_remains=eff_remains.append(eff.xs('%library remains', axis=0, level=1))
    eff_remains=eff_remains.append(pd.DataFrame([eff.at[('pairing','step output'),0]], index=['Nreads output']))
    qc_data={
        'efficiency_complete_stats':eff,
        'efficiency_LOSS':eff_loss,
        'efficiency_REMAINS':eff_remains
    }
    if not savefolder is None:
        with pd.ExcelWriter(os.path.join(sample_root_folder,savefolder,'qc.xlsx')) as writer:
            for k, v in qc_data.items():
                v.to_excel(writer, sheet_name=k)
        for k, v in qc_data.items():
            v.to_csv(os.path.join(sample_root_folder,savefolder,k+'.simple.csv'), sep='\t', header=True)
            print("Writing %s QC for sample %s"%(k,sample_root_folder))
            pq.write_table(pa.Table.from_pandas(v), os.path.join(sample_root_folder,savefolder,'_parquet',k+'.simple.parquet'))
    return qc_data

def load_qc(smpls, smpls_names, qc_root_folder='pairs/gencondeV29_hg38/all/QC/_parquet'):
    qc_to_load=['efficiency_REMAINS','efficiency_LOSS','efficiency_complete_stats','biotype_census_CPM','biotype_census_CPKM','trans_statistics','top_hits_bytype_CPM','top_hits_bytype_CPKM','top_hits_CPM','top_hits_CPKM']
    qc_data={}
    
    for this_qc in qc_to_load:
        qc_data[this_qc]=pd.concat([pd.read_parquet(os.path.join(s, qc_root_folder, this_qc+'.parquet')) for s in smpls], axis=1, keys=smpls_names)
        if this_qc in ['efficiency_complete_stats','efficiency_LOSS','efficiency_REMAINS']:
            qc_data[this_qc].columns=qc_data[this_qc].columns.droplevel(1)
            

    return qc_data

def load_qc_simple(smpls, smpls_names, qc_root_folder='pairs/gencondeV29_hg38/all/QC/_parquet'):
    qc_to_load=['efficiency_REMAINS','efficiency_LOSS','efficiency_complete_stats']
    qc_data={}
    
    for this_qc in qc_to_load:
        qc_data[this_qc]=pd.concat([pd.read_parquet(os.path.join(s, qc_root_folder, this_qc+'.simple.parquet')) for s in smpls], axis=1, keys=smpls_names)
        if this_qc in ['efficiency_complete_stats','efficiency_LOSS','efficiency_REMAINS']:
            qc_data[this_qc].columns=qc_data[this_qc].columns.droplevel(1)
            

    return qc_data

def export_qc(qc_data, fname):
     with pd.ExcelWriter(fname) as writer:
        for k, v in qc_data.items():
            v.to_excel(writer, sheet_name=k)

# def run_qc_multisamples(smpls):
# test
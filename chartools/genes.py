from . import analysis
import numpy as np
import pandas as pd

def geneify(df, tx_df):
    if df.columns.nlevels>1:
        tx_df_multiindex=tx_df.copy()
        tx_df_multiindex.columns=pd.MultiIndex.from_product([tx_df.columns]+[['']]*(len(df.columns.levels)-1))
    else:
        tx_df_multiindex=tx_df
    
    df_withinfo=pd.concat([tx_df_multiindex,df], axis=1, join='inner').groupby('ENSG').sum()
    
    return df_withinfo

def add_gene_info(df, genes_df, join='outer'):
    if df.columns.nlevels>1:
        genes_df_multiindex=genes_df.copy()
        genes_df_multiindex.columns=pd.MultiIndex.from_product([genes_df.columns]+[['']]*(len(df.columns.levels)-1))
    else:
        genes_df_multiindex=genes_df
    
    df_withinfo=pd.concat([genes_df_multiindex,df], axis=1, join=join)
    return df_withinfo


def rna_stats(exons_stats_file, introns_stats_file, intergene_stats_file, tx_df, genes_df, chr_norm):
    np.seterr(divide='ignore', invalid='ignore')
    data_exons=pd.read_csv(exons_stats_file, delimiter=" ", header=None, names=['Count','ENST','uniqmapping','istrans'])
    
    data_tx=data_exons.join(tx_df, on='ENST').drop('ENST',axis=1).groupby(['ENSG','uniqmapping','istrans']).sum().unstack(level=[1,2]).fillna(0)
    
    data_intergenic=pd.read_csv(intergene_stats_file, delimiter=" ", header=None, names=['Count','ENSG','uniqmapping','istrans']).groupby(['ENSG','uniqmapping','istrans']).sum().unstack(level=[1,2]).fillna(0)
    
    data_gb=pd.read_csv(introns_stats_file, delimiter=" ", header=None, names=['Count','ENSG','uniqmapping','istrans']).groupby(['ENSG','uniqmapping','istrans']).sum().unstack(level=[1,2]).fillna(0)
        
    mydf={}
    n=['exons','introns']
    for nn, data_bygene in zip(n,[data_tx, data_gb]):
        unq=data_bygene.loc[:, (slice(None),[1],slice(None))].sum(axis=1,level=2)
        alls=data_bygene.sum(axis=1,level=2)
        unq['NumContacts']=unq.sum(axis=1)
        alls['NumContacts']=alls.sum(axis=1)
        
        mydf[nn]=pd.concat([unq,alls], axis=1, keys=['unq','all']).fillna(0)
#     x=add_gene_info(mydf, genes_df).sort_values(('unq','tot'), ascending=False)
    
    x0=pd.concat([mydf[k] for k in n], keys=n, axis=1).fillna(0).reorder_levels([1,0,2],axis=1) #.sort_index(level=[0,1], axis=1)

    x=add_gene_info(x0,genes_df[['chr']], join='inner')
    mynorm=x['chr'].map(chr_norm).values
    for k in ['unq','all']:
        
        
        for kk in [-1,0,1,'NumContacts']:
            x[(k,"all",kk)]=x[(k,'exons',kk)]+x[(k,'introns',kk)]
#         for kk in [-1,0,1]:    
        for hh in ['exons','introns','all']:
            x[(k,hh,-1)]=x[(k,hh,-1)].values/x[(k,hh,'NumContacts')].values*100 #ambiguous
            Ntrans=x[(k,hh,1)].values
            Ncis=x[(k,hh,0)].values
            x[(k,hh,0)]=Ntrans/(Ntrans+Ncis)*100 #trans
            x[(k,hh,1)]=(Ntrans/Ncis)*mynorm #trans
            Nreads=x[(k,hh,'NumContacts')].sum()/1000000
            x[(k,hh,'CPM')]=x[(k,hh,'NumContacts')]/Nreads
        x[(k,'all','%exons')]=x[(k,'exons','NumContacts')].values/(x[(k,'exons','NumContacts')].values+x[(k,'introns','NumContacts')].values)*100
            
        
        
    x.rename(columns={-1:'%ambiguousDNA',0:'%trans',1:'trans_chr_score'}, inplace=True, level=2)      
    x.sort_index(level=[1,2], axis=1, inplace=True)
    
    data_all=x.xs('all', level=0, axis=1).sort_values(('all','NumContacts'), ascending=False).reindex(['NumContacts','CPM','%exons','%ambiguousDNA','%trans','trans_chr_score'], axis=1, level=1)
    data_unq=x.xs('all', level=0, axis=1).sort_values(('all','NumContacts'), ascending=False).reindex(['NumContacts','CPM','%exons','%ambiguousDNA','%trans','trans_chr_score'], axis=1, level=1)
    data_intergenic.columns=data_intergenic.columns.droplevel([0,1])
    ig_tot=data_intergenic.sum(axis=1)
    

    # # for c in [-1,0,1]:
    # # data_intergenic[c]=data_intergenic[c].values/data_intergenic['tot'].values*100
    ig_Ntrans=data_intergenic[1].values
    ig_Ncis=data_intergenic[0].values
    ig_Nambig=data_intergenic[-1].values

    
    data_intergenic[0]=ig_Ntrans/(ig_Ntrans+ig_Ncis)*100 #trans
    data_intergenic[-1]=ig_Nambig/ig_tot #trans
    data_intergenic[1]=ig_tot
    # # data_intergenic.drop([1], axis=1)
    # #         Nreads=x[(k,hh,'tot')].sum()/1000000
    data_intergenic.rename(columns={-1:'%ambiguousDNA',0:'%trans',1:'N'}, inplace=True)
    data_intergenic.columns.names=[None]
    
    return data_unq, data_all, data_intergenic


def rna_stats2(stats_file, tx_df, genes_df, chr_norm, tx_start='ENST', gene_start='ENSG'):
    np.seterr(divide='ignore', invalid='ignore')
    data=pd.read_csv(stats_file, delimiter=" ", header=None, names=['Count','ENSN','uniqmapping','istrans'])
    data_exons=data.loc[data['ENSN'].str.startswith(tx_start)].copy()
    data_exons.rename(columns={"ENSN": "ENST"}, inplace=True)
    # data_exons=pd.read_csv(exons_stats_file, delimiter=" ", header=None, names=['Count','ENST','uniqmapping','istrans'])
    
    data_tx=data_exons.join(tx_df, on='ENST').drop('ENST',axis=1).groupby(['ENSG','uniqmapping','istrans']).sum().unstack(level=[1,2]).fillna(0)
    
    data_intergenic=data.loc[data['ENSN'].str.startswith('*')].copy()
    data_intergenic.rename(columns={"ENSN": "ENSG"}, inplace=True)
    data_intergenic=data_intergenic.groupby(['ENSG','uniqmapping','istrans']).sum().unstack(level=[1,2]).fillna(0)
    
    # pd.read_csv(intergene_stats_file, delimiter=" ", header=None, names=['Count','ENSG','uniqmapping','istrans']).groupby(['ENSG','uniqmapping','istrans']).sum().unstack(level=[1,2]).fillna(0)
    
    data_gb=data.loc[data['ENSN'].str.startswith(gene_start)].copy()
    data_gb.rename(columns={"ENSN": "ENSG"}, inplace=True)
    data_gb=data_gb.groupby(['ENSG','uniqmapping','istrans']).sum().unstack(level=[1,2]).fillna(0)
    # pd.read_csv(introns_stats_file, delimiter=" ", header=None, names=['Count','ENSG','uniqmapping','istrans'])
    mydf={}
    n=['exons','introns']
    for nn, data_bygene in zip(n,[data_tx, data_gb]):
        unq=data_bygene.loc[:, (slice(None),[1],slice(None))].sum(axis=1,level=2)
        alls=data_bygene.sum(axis=1,level=2)
        unq['NumContacts']=unq.sum(axis=1)
        alls['NumContacts']=alls.sum(axis=1)
        
        mydf[nn]=pd.concat([unq,alls], axis=1, keys=['unq','all']).fillna(0)
#     x=add_gene_info(mydf, genes_df).sort_values(('unq','tot'), ascending=False)
    
    x0=pd.concat([mydf[k] for k in n], keys=n, axis=1).fillna(0).reorder_levels([1,0,2],axis=1) #.sort_index(level=[0,1], axis=1)

    x=add_gene_info(x0,genes_df[['chr']], join='inner')
    mynorm=x['chr'].map(chr_norm).values
    for k in ['unq','all']:
        
        
        for kk in [-1,0,1,'NumContacts']:
            x[(k,"all",kk)]=x[(k,'exons',kk)]+x[(k,'introns',kk)]
#         for kk in [-1,0,1]:    
        for hh in ['exons','introns','all']:
            x[(k,hh,-1)]=x[(k,hh,-1)].values/x[(k,hh,'NumContacts')].values*100 #ambiguous
            Ntrans=x[(k,hh,1)].values
            Ncis=x[(k,hh,0)].values
            x[(k,hh,0)]=Ntrans/(Ntrans+Ncis)*100 #trans
            x[(k,hh,1)]=(Ntrans/Ncis)*mynorm #trans
            Nreads=x[(k,hh,'NumContacts')].sum()/1000000
            x[(k,hh,'CPM')]=x[(k,hh,'NumContacts')]/Nreads
        x[(k,'all','%exons')]=x[(k,'exons','NumContacts')].values/(x[(k,'exons','NumContacts')].values+x[(k,'introns','NumContacts')].values)*100
            
        
        
    x.rename(columns={-1:'%ambiguousDNA',0:'%trans',1:'trans_chr_score'}, inplace=True, level=2)      
    x.sort_index(level=[1,2], axis=1, inplace=True)
    
    data_all=x.xs('all', level=0, axis=1).sort_values(('all','NumContacts'), ascending=False).reindex(['NumContacts','CPM','%exons','%ambiguousDNA','%trans','trans_chr_score'], axis=1, level=1)
    data_unq=x.xs('unq', level=0, axis=1).sort_values(('all','NumContacts'), ascending=False).reindex(['NumContacts','CPM','%exons','%ambiguousDNA','%trans','trans_chr_score'], axis=1, level=1)
    data_intergenic.columns=data_intergenic.columns.droplevel([0,1])
    ig_tot=data_intergenic.sum(axis=1)
    

    # # for c in [-1,0,1]:
    # # data_intergenic[c]=data_intergenic[c].values/data_intergenic['tot'].values*100
    ig_Ntrans=data_intergenic[1].values
    ig_Ncis=data_intergenic[0].values
    ig_Nambig=data_intergenic[-1].values

    
    data_intergenic[0]=ig_Ntrans/(ig_Ntrans+ig_Ncis)*100 #trans
    data_intergenic[-1]=ig_Nambig/ig_tot #trans
    data_intergenic[1]=ig_tot
    # # data_intergenic.drop([1], axis=1)
    # #         Nreads=x[(k,hh,'tot')].sum()/1000000
    data_intergenic.rename(columns={-1:'%ambiguousDNA',0:'%trans',1:'N'}, inplace=True)
    data_intergenic.columns.names=[None]
    
    return data_unq, data_all, data_intergenic

def count_table(data, tx_df, genes_df, chr_df=None):
    tx_dict=tx_df['ENSG'].to_dict()
    for g in set(tx_df['ENSG'].values):
        tx_dict[g]=g
        tx_dict['*']='*'

    
    len_df=pd.concat([tx_df['L'],genes_df['L']])
    d=data.copy()
    d.index=pd.MultiIndex.from_arrays([d.index.get_level_values(0).map(tx_dict), d.index.get_level_values(1)])
    d=d.groupby(['ENSN','istrans']).sum()
    
    d_total=d.groupby('ENSN').sum()
    try:
        d_cis=d.xs(0, level=1)
    except KeyError:
        d_cis=pd.DataFrame(np.zeros(len(d_total), int), index=d_total.index)
    try:
        d_trans=d.xs(1, level=1)
    except KeyError:
        d_trans=pd.DataFrame(np.zeros(len(d_total), int), index=d_total.index)


    N=pd.concat([d_total,d_cis,d_trans], axis=1, sort=False).fillna(0).astype(int)
    
    N.columns=['N','N_cis','N_trans']
    N['CPM']=N['N']/N['N'].sum()*1000000
    
    N['CPKM']=N['CPM']*1000/len_df.reindex(N.index)

    if not chr_df is None:
        chr_len_df=chr_df.copy()
        chr_len_df['Ltrans']=chr_len_df['L'].sum()-chr_len_df['L']

        # chr_len_cis_dict=chr_len_df['L'].to_dict()
        # chr_len_trans_dict=chr_len_df['Ltrans'].to_dict()
        chr_norm=(chr_len_df['L']/chr_len_df['Ltrans']).to_dict()

        Lcis_over_Ltrans=genes_df['chr'].map(chr_norm)
        
     
        N['tscore']=(N['N_trans']/N['N_cis']*Lcis_over_Ltrans.reindex(N.index)).values
    
    N.index.name='ENSN'
    N['name']=genes_df['name'].reindex(N.index)
    N['type']=genes_df['type'].reindex(N.index)
    N=N.reset_index().set_index(['ENSN','name','type'])
    return N

def namify_count_table(ct, genes_df, chr_df):
    # expects:
    # chr_df=pd.read_csv(chrNameLengthFile, sep='\t', header=None, names=['chr','L'], index_col=0)
    chr_len_df=chr_df.copy()
   
    chr_len_df['Ltrans']=chr_len_df['L'].sum()-chr_len_df['L']

    chr_len_cis_dict=chr_len_df['L'].to_dict()
    chr_len_trans_dict=chr_len_df['Ltrans'].to_dict()
    
    ct_byname=pd.concat([genes_df[['chr']],ct[['N','N_cis','N_trans','CPM','CPKM']].reset_index().set_index(['ENSN'])], axis=1, join='inner').set_index(['name','type'])
    
    ct_byname['cis_density']=ct_byname['N_cis']/ct_byname['chr'].map(chr_len_cis_dict)
    ct_byname['trans_density']=ct_byname['N_trans']/ct_byname['chr'].map(chr_len_trans_dict)
    
#     ct_byname['NLoci']=np.ones(len(ct_byname), dtype=int)
    ct_byname=ct_byname.groupby(['name','type']).sum()
    
    
    ct_byname['tscore']=ct_byname['trans_density']/ct_byname['cis_density']
    ct_byname=ct_byname.drop(['cis_density','trans_density'], axis=1)
    
    return ct_byname

def summarize_counts(stats_file, tx_df, genes_df, chr_df, unq=True, tx_start='ENST', gene_start='ENSG'):
    
    data=pd.read_csv(stats_file, delimiter=" ", header=None, names=['Count','ENSN','uniqmapping','istrans'], index_col=['ENSN','uniqmapping','istrans'])
    
    if '*' in data.index.get_level_values(0): 
        data_intergenic=data.loc['*'].xs(-1, level=0).reindex([-1,0,1]).fillna(0)
        d_intergenic=pd.DataFrame([data_intergenic.sum().iat[0],data_intergenic.iat[1,0],data_intergenic.iat[2,0]]).transpose()
        d_intergenic.columns=['N','N_cis','N_trans']
        d_intergenic.index=pd.Index(['*'], name='ENSN')
    else:
        d_intergenic=None

    if unq:
        data_unq=data.xs(1, level=1) #.drop('*', axis=0, level=1)
        data_unq_exons=data_unq.loc[data_unq.index.get_level_values(0).str.startswith('ENST')]
        data_unq_introns=data_unq.loc[data_unq.index.get_level_values(0).str.startswith('ENSG')]
    else:
        data_unq=data.groupby(['ENSN','istrans']).sum()
        if '*' in data_unq.index.get_level_values(0):
            data_unq=data_unq.drop('*', axis=0, level=0)
        data_unq_exons=data_unq.loc[data_unq.index.get_level_values(0).str.startswith('ENST')]
        data_unq_introns=data_unq.loc[data_unq.index.get_level_values(0).str.startswith('ENSG')]
    
    if (len(data_unq_exons)>0) and (len(data_unq_introns)>0):
        data_unq_all=[data_unq,data_unq_exons, data_unq_introns]
        names=['all','exons','introns']
        
    
    elif (len(data_unq_exons)>0):
        data_unq_all=[data_unq_exons]
        names=['all']
    else:
        data_unq_all=[data_unq_introns]
        names=['all']

    cts=[count_table(d, tx_df, genes_df, chr_df) for d in data_unq_all]
    cts_name=[namify_count_table(ct,genes_df, chr_df) for ct in cts]
    out={'ct':pd.concat([ct for ct in cts], axis=1, keys=names, sort=False), 'ct_name': pd.concat([ct_name for ct_name in cts_name], axis=1, keys=names, sort=False)}

    out['ct'].loc[:,(slice(None),['N','N_cis','N_trans'])]=out['ct'].loc[:,(slice(None),['N','N_cis','N_trans','CPM'])].fillna(0).astype(int)
    out['ct'].loc[:,(slice(None),['CPM','CPKM'])]=out['ct'].loc[:,(slice(None),['CPM','CPKM'])].fillna(0)

    out['ct_name'].loc[:,(slice(None),['N','N_cis','N_trans'])]=out['ct_name'].loc[:,(slice(None),['N','N_cis','N_trans'])].fillna(0).astype(int)
    out['ct_name'].loc[:,(slice(None),['CPM','CPKM'])]=out['ct_name'].loc[:,(slice(None),['CPM','CPKM'])].fillna(0)
    
#     out['ct_name']=out['ct_name'].drop('')
    return out['ct'], out['ct_name'], d_intergenic




def rna_stats2noDNA(stats_file, tx_df, genes_df, chr_norm):
    np.seterr(divide='ignore', invalid='ignore')
    data=pd.read_csv(stats_file, delimiter=" ", header=None, names=['Count','ENSN','uniqmapping','istrans'])
    data_exons=data.loc[data['ENSN'].str.startswith('ENST')].copy()
    data_exons.rename(columns={"ENSN": "ENST"}, inplace=True)
    # data_exons=pd.read_csv(exons_stats_file, delimiter=" ", header=None, names=['Count','ENST','uniqmapping','istrans'])
    
    data_tx=data_exons.join(tx_df, on='ENST').drop('ENST',axis=1).groupby(['ENSG','uniqmapping','istrans']).sum().unstack(level=[1,2]).fillna(0)
    
    data_intergenic=data.loc[data['ENSN'].str.startswith('*')].copy()
    data_intergenic.rename(columns={"ENSN": "ENSG"}, inplace=True)
    data_intergenic=data_intergenic.groupby(['ENSG','uniqmapping','istrans']).sum().unstack(level=[1,2]).fillna(0)
    
    # pd.read_csv(intergene_stats_file, delimiter=" ", header=None, names=['Count','ENSG','uniqmapping','istrans']).groupby(['ENSG','uniqmapping','istrans']).sum().unstack(level=[1,2]).fillna(0)
    
    data_gb=data.loc[data['ENSN'].str.startswith('ENSG')].copy()
    data_gb.rename(columns={"ENSN": "ENSG"}, inplace=True)
    data_gb=data_gb.groupby(['ENSG','uniqmapping','istrans']).sum().unstack(level=[1,2]).fillna(0)
    # pd.read_csv(introns_stats_file, delimiter=" ", header=None, names=['Count','ENSG','uniqmapping','istrans'])
    mydf={}
    n=['exons','introns']
    for nn, data_bygene in zip(n,[data_tx, data_gb]):
        unq=data_bygene.loc[:, (slice(None),[1],slice(None))].sum(axis=1,level=2)
        alls=data_bygene.sum(axis=1,level=2)
        unq['NumContacts']=unq.sum(axis=1)
        alls['NumContacts']=alls.sum(axis=1)
        
        mydf[nn]=pd.concat([unq,alls], axis=1, keys=['unq','all']).fillna(0)
#     x=add_gene_info(mydf, genes_df).sort_values(('unq','tot'), ascending=False)
    
    x0=pd.concat([mydf[k] for k in n], keys=n, axis=1).fillna(0).reorder_levels([1,0,2],axis=1) #.sort_index(level=[0,1], axis=1)

    x=add_gene_info(x0,genes_df[['chr']], join='inner')
    mynorm=x['chr'].map(chr_norm).values
    for k in ['unq','all']:
        
        
        for kk in ['NumContacts']:
            x[(k,"all",kk)]=x[(k,'exons',kk)]+x[(k,'introns',kk)]
#         for kk in [-1,0,1]:    
        for hh in ['exons','introns','all']:
            Nreads=x[(k,hh,'NumContacts')].sum()/1000000
            x[(k,hh,'CPM')]=x[(k,hh,'NumContacts')]/Nreads
        x[(k,'all','%exons')]=x[(k,'exons','NumContacts')].values/(x[(k,'exons','NumContacts')].values+x[(k,'introns','NumContacts')].values)*100
            
        
        
    # x.rename(columns={-1:'%ambiguousDNA',0:'%trans',1:'trans_chr_score'}, inplace=True, level=2)      
    # x.sort_index(level=[1,2], axis=1, inplace=True)
    
    data_all=x.xs('all', level=0, axis=1).sort_values(('all','NumContacts'), ascending=False).reindex(['NumContacts','CPM','%exons'], axis=1, level=1)
    data_unq=x.xs('unq', level=0, axis=1).sort_values(('all','NumContacts'), ascending=False).reindex(['NumContacts','CPM','%exons'], axis=1, level=1)
    data_intergenic.columns=data_intergenic.columns.droplevel([0,1])
    ig_tot=data_intergenic.sum(axis=1)
    

    # # for c in [-1,0,1]:
    # # data_intergenic[c]=data_intergenic[c].values/data_intergenic['tot'].values*100
    

    
    
    data_intergenic[1]=ig_tot
    # # data_intergenic.drop([1], axis=1)
    # #         Nreads=x[(k,hh,'tot')].sum()/1000000
    data_intergenic.rename(columns={1:'N'}, inplace=True)
    data_intergenic.columns.names=[None]
    
    return data_unq, data_all, data_intergenic


def rna_stats_noDNA(exons_stats_file, introns_stats_file, tx_df, genes_df, chr_norm):
    np.seterr(divide='ignore', invalid='ignore')
    data_exons=pd.read_csv(exons_stats_file, delimiter=" ", header=None, names=['Count','ENST','uniqmapping'])
    data_exons['istrans']=-1
    
    data_tx=data_exons.join(tx_df, on='ENST').drop('ENST',axis=1)
    # data_tx['istrans']=-1
    data_tx=data_tx.groupby(['ENSG','uniqmapping','istrans']).sum().unstack(level=[1,2]).fillna(0)
    
    # data_intergenic=pd.read_csv(intergene_stats_file, delimiter=" ", header=None, names=['Count','ENSG','uniqmapping'])
    # data_intergenic['istrans']=-1
    # data_intergenic=data_intergenic.groupby(['ENSG','uniqmapping','istrans']).sum().unstack(level=[1,2]).fillna(0)
    
    data_gb=pd.read_csv(introns_stats_file, delimiter=" ", header=None, names=['Count','ENSG','uniqmapping'])
    data_gb['istrans']=-1
    data_gb=data_gb.groupby(['ENSG','uniqmapping','istrans']).sum().unstack(level=[1,2]).fillna(0)
        
    mydf={}
    n=['exons','introns']
    for nn, data_bygene in zip(n,[data_tx, data_gb]):
        unq=data_bygene.loc[:, (slice(None),[1],slice(None))].sum(axis=1,level=2)
        alls=data_bygene.sum(axis=1,level=2)
        unq['NumContacts']=unq.sum(axis=1)
        alls['NumContacts']=alls.sum(axis=1)
        
        mydf[nn]=pd.concat([unq,alls], axis=1, keys=['unq','all']).fillna(0)
#     x=add_gene_info(mydf, genes_df).sort_values(('unq','tot'), ascending=False)
    
    x0=pd.concat([mydf[k] for k in n], keys=n, axis=1).fillna(0).reorder_levels([1,0,2],axis=1) #.sort_index(level=[0,1], axis=1)

    x=add_gene_info(x0,genes_df[['chr']], join='inner')
    mynorm=x['chr'].map(chr_norm).values
    for k in ['unq','all']:
        
        
        for kk in [-1,'NumContacts']:
            x[(k,"all",kk)]=x[(k,'exons',kk)]+x[(k,'introns',kk)]
#         for kk in [-1,0,1]:    
        for hh in ['exons','introns','all']:
            x[(k,hh,-1)]=x[(k,hh,-1)].values/x[(k,hh,'NumContacts')].values*100 #ambiguous
            # Ntrans=x[(k,hh,1)].values
            # Ncis=x[(k,hh,0)].values
            # x[(k,hh,0)]=Ntrans/(Ntrans+Ncis)*100 #trans
            # x[(k,hh,1)]=(Ntrans/Ncis)*mynorm #trans
            Nreads=x[(k,hh,'NumContacts')].sum()/1000000
            x[(k,hh,'CPM')]=x[(k,hh,'NumContacts')]/Nreads
        x[(k,'all','%exons')]=x[(k,'exons','NumContacts')].values/(x[(k,'exons','NumContacts')].values+x[(k,'introns','NumContacts')].values)*100
            
        
        
    x.rename(columns={-1:'%ambiguousDNA'}, inplace=True, level=2)      
    x.sort_index(level=[1,2], axis=1, inplace=True)
    
    data_all=x.xs('all', level=0, axis=1).sort_values(('all','NumContacts'), ascending=False).reindex(['NumContacts','CPM','%exons','%ambiguousDNA'], axis=1, level=1).drop([('all','%ambiguousDNA'),('exons','%ambiguousDNA'),('introns','%ambiguousDNA')], axis=1)
    data_unq=x.xs('unq', level=0, axis=1).sort_values(('all','NumContacts'), ascending=False).reindex(['NumContacts','CPM','%exons','%ambiguousDNA'], axis=1, level=1).drop([('all','%ambiguousDNA'),('exons','%ambiguousDNA'),('introns','%ambiguousDNA')], axis=1)
    # data_intergenic.columns=data_intergenic.columns.droplevel([0,1])
    # ig_tot=data_intergenic.sum(axis=1)
    

    # # for c in [-1,0,1]:
    # # data_intergenic[c]=data_intergenic[c].values/data_intergenic['tot'].values*100
    
    # ig_Nambig=data_intergenic[-1].values

    
    # data_intergenic[0]=ig_Ntrans/(ig_Ntrans+ig_Ncis)*100 #trans
    # data_intergenic[-1]=ig_Nambig/ig_tot #trans
    # data_intergenic[1]=ig_tot
    # # # data_intergenic.drop([1], axis=1)
    # # #         Nreads=x[(k,hh,'tot')].sum()/1000000
    # data_intergenic.rename(columns={-1:'%ambiguousDNA',0:'%trans',1:'N'}, inplace=True)
    # data_intergenic.columns.names=[None]
    
    return data_unq, data_all #, data_intergenic


def rna_stats_sample(s, tx_df, genes_df, chr_norm):
    interegene_stats_file(s+'/pairs/SE_merge_pear_decon/starbowtie_resolved/rnaTable.stats.txt',s+'/pairs/SE_merge_pear_decon/starbowtie_genebodies/rnaTable.stats.txt',s+'/pairs/SE_merge_pear_decon/starbowtie_NOgenebodies/rnaTable.stats.txt')

    data_unq, data_all, data_intergenic = rna_stats(exons_stats_file, introns_stats_file, intergene_stats_file, tx_df, genes_df, chr_norm)
    return data_unq, data_all, data_intergenic

def get_annotations_db(chrNameLengthFile, txdbFile, genedbFile):

    chr_dict=analysis.load_index_fromCSV(chrNameLengthFile)
    chr_len=analysis.make_chr_length_dict_fromCSV(chrNameLengthFile) 
    chr_len_vec=np.zeros(len(chr_dict))
    for k, v in chr_dict.items():
        chr_len_vec[v]=chr_len[k]
    
    tx_dict=analysis.make_tx_dict_fromCSV(txdbFile)

    genes_df=pd.read_csv(genedbFile, header=None, sep='\t', names=['ENSG','type','name','chr','strand','position'], index_col=0, usecols=['ENSG','type','name','chr','strand'])

    tx_dict_simple={k: v[0] for k, v in tx_dict.items()}
    tx_df=pd.DataFrame.from_dict(tx_dict_simple, orient='index')
    tx_df.index.names=['ENST']
    tx_df.columns=['ENSG']

    Ltot=np.sum([L for _,L in chr_len.items()])
    chr_norm={k:L/(Ltot-L) for k, L in chr_len.items()}

    return genes_df, tx_df, chr_dict, chr_norm

def get_annotations_db2(chrNameLengthFile, txdbFile, genedbFile):
    tx_df=pd.read_csv(txdbFile, header=None, sep='\t', names=['ENST','ENSG','strand','position','L'], index_col=0, usecols=['ENST','ENSG','L'])
    genes_df=pd.read_csv(genedbFile, header=None, sep='\t', names=['ENSG','type','name','chr','strand','position','L'], index_col=0, usecols=['ENSG','type','name','chr','strand','L'])

    chr_df=pd.read_csv(chrNameLengthFile, sep='\t', header=None, names=['chr','L'], index_col=0)


    return genes_df, tx_df, chr_df

def load_salmon(salmonFile, tx_df):
    sdf=pd.read_csv(salmonFile, header=0, sep='\t', index_col=0)
    N=sdf['NumReads'].sum()
    sdf['RPM']=sdf['NumReads']/N*1000000
    sdf_genes=geneify(sdf, tx_df[['ENSG']])

    tx_len=pd.concat([tx_df[['ENSG']],sdf['Length']], axis=1, join='outer').groupby('ENSG').max()
    tx_len.columns=['MaxLength']

    sdf_genes_withlen=pd.concat([sdf_genes,tx_len], axis=1, join='inner')
    sdf_genes_withlen.drop(['Length','EffectiveLength'], axis=1, inplace=True)
    sdf_genes_withlen['RPKM_maxlen']=sdf_genes_withlen['RPM']/sdf_genes_withlen['MaxLength']*1000
    # sdf['RPKM']=sdf['RPM']/sdf['Length']*1000

    return sdf_genes_withlen

import pandas as pd
import numpy as np

# in this module, df is a dataframe with samples at level 0 of columns, outcome variables at level -1, and possible measurement type at level 1 and all intermediate levels

# row index is gene/ID at level 0, and type at level -1

# def load_star

# def load_salmon

# def geneify

# def recompute_tscores

def renormalize(df, v=['CPM'], bytype=False):
    out=df.copy()
    slicer=(slice(None),)+tuple([out.columns.get_level_values(i)[0] for i in range(1,out.columns.nlevels-1)])+(v,)
    x=out.loc[:,slicer]
    if bytype:
        out.loc[:,slicer]=x.groupby('type').transform(lambda x: x/np.sum(x)*1000000)
    else:
        out.loc[:,slicer]=x/np.sum(x, axis=1)*1000000
    # for s in set(df.columns.get_level_values(0)):
    #     for t in set(df.columns.get_level_values(1)):
    #         for v in var:
    #             if (s,t,v) in df.columns:
    #                 if bytype:
    #                     out.loc[:,(s,t,v)]=df.loc[:,(s,t,v)].groupby('type').transform(lambda x: x/np.sum(x)*1000000)
    #                 else:
    #                     out.loc[:,(s,t,v)]=df.loc[:,(s,t,v)]/np.sum(df.loc[:,(s,t,v)])*1000000
                        
    return out

def pool(df, sample_groups, sample_order=None, normalization='average', normalization_var='N'):
    if sample_order is None:
        sample_order=list(sample_groups.items())
    
    
    pool=[]
    
    for g in sample_order:
        these_samples=sample_groups[g]
        this_df=df.loc[:,(these_samples, slice(None), slice(None))] #.copy()
        
        grouped=this_df.groupby(level=list(range(1,this_df.columns.nlevels)), axis=1)
        
        if normalization=='average':
                pool+=[grouped.mean()]

        elif normalization=='weighted_global':
            def wm(x):
                
                slicer=(slice(None),)+tuple([x.columns.get_level_values(i)[0] for i in range(1,this_df.columns.nlevels-1)])+(normalization_var,)
                
                weights=this_df.loc[:,slicer].sum(axis=0).values
                weights=weights/weights.sum()
                
                return (x*weights[np.newaxis,:]).sum(axis=1)
        
            pool+=[grouped.apply(wm)]
            
        elif normalization=='weighted_gene':

            
            def wm(x):
                
                slicer=(slice(None),)+tuple([x.columns.get_level_values(i)[0] for i in range(1,this_df.columns.nlevels-1)])+(normalization_var,)
                
                weights_norm=this_df.loc[:,slicer].sum(axis=1).values[:,np.newaxis]
                weights=this_df.loc[:,slicer].values/weights_norm
                return (x*weights).sum(axis=1)
        
            pool+=[grouped.apply(wm)]
            
    return pd.concat(pool, axis=1, keys=sample_order)

def simplify(df,level=None,vars=['CPM','tscore'], types_subset=None, threshold_var='CPM', threshold=0):
    
    if types_subset is None:
        slicer=slice(None)
    else:
        slicer=types_subset
        
    if level is None:
        out=df
    else:
        out=df.xs(level, axis=1, level=1)

    sliced=out.loc[(slice(None),slice(None),slicer),:]
    
    m=sliced.loc[:,(slice(None),threshold_var)].mean(axis=1)>threshold
    
    out=sliced.loc[m,(slice(None),vars)]
    out_null=out.sum(axis=0)==0
    out=out.drop(out_null.loc[out_null].index, axis=1)
    return out
        
def pairwise_comparison_table(mt, sample_groups, sample_order, comps, comp_vars, type_subset=None, level='all', normalization='average', normalization_var='N', threshold_var='CPM', threshold=0.1):
    #mt is a master table

    pt=pool(mt, sample_groups, sample_order=sample_order, normalization=normalization, normalization_var=normalization_var)
#     pt is pulled table
    
    ct=simplify(pt,level=level, vars=comp_vars, types_subset=type_subset, threshold_var=threshold_var, threshold=threshold)
    # ct is caomparision table
    
    index_order=sample_order.copy()

    for this_comp in comps:
        c0=this_comp[0]
        c1=this_comp[1]
        index_order+=[c0+"_"+c1,c0+"+"+c1]
        
        for v in comp_vars:
            if (((c0,v) in ct.columns) and ((c1,v) in ct.columns)):
                ct[c0+"_"+c1, v]=ct[c0,v]/ct[c1,v]
                ct[c0+"+"+c1, v]=(ct[c0,v]+ct[c1,v])/2

    ct=ct.reindex(index_order, axis=1, level=0)
    return ct
        
def coarsen_type(df, d): # d is a coarsening dictionnary
    out=df.copy()
    coarse_type=[d.get(k,'other') for k in df.index.get_level_values(2)]
    out.index=pd.MultiIndex.from_arrays([out.index.get_level_values(0),out.index.get_level_values(1),coarse_type])
    out.index.names=['ENSG','name','type']
    return out

def get_scale_factors(df, threshold=0):
    included_data=df.loc[df.sum(axis=1)>threshold]
    x=np.maximum(1.0,included_data.fillna(0))
    scale=1/np.median(x/np.exp(np.mean(np.log(x), axis=1))[:,np.newaxis], axis=0)
    return scale

def scale(df, scale_columns=['N'], types_subset=None, threshold=0):
    # geometric scaling of sample counts as standard in DEseq
    grouped=df.groupby(level=list(range(1,df.columns.nlevels)), axis=1)
    
    if types_subset is None:
        row_slicer=(slice(None),slice(None), slice(None))
    else:
        row_slicer=(slice(None),slice(None), types_subset)
        
    col_slicer=(slice(None),)+tuple([df.columns.get_level_values(i)[0] for i in range(1,df.columns.nlevels-1)])+('N',)
    
    
    def scalefun(x):
        if x.columns.get_level_values(-1)[0] in scale_columns:
           
            included_data=df.loc[row_slicer,col_slicer].fillna(0)
            included_data=included_data.loc[included_data.sum(axis=1)>threshold]
            Nvalues=np.maximum(included_data,1)
            
            
            scale=1/np.median(Nvalues/np.exp(np.mean(np.log(Nvalues), axis=1))[:,np.newaxis], axis=0)
            print([x.columns.get_level_values(i)[0] for i in range(1,df.columns.nlevels)])
            print(scale)
            return x*scale[np.newaxis,:]
        else:
            return x
                                      
                                      
    return grouped.apply(scalefun)

# remove variables and drop intermediate level to desired if any

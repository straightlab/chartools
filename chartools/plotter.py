import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
from bioframe import fetch_chromsizes, binnify
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import pyBigWig as pbw
import bioframe as bf


def get_intervals(chro, start=0, stop=-1, binsize=10000):
    
    if (chro is None) or (type(chro)==list):
        intervals = get_full_chromsomes(chros=chro, binsize=binsize)
    else: 
    
    
        chr_l = fetch_chromsizes('hg38')
        if stop==-1:
            stop=chr_l[chro]
        else:
            stop=min(chr_l[chro], stop)

        breaks = np.arange(start, stop+1, binsize)
        intervals = pd.DataFrame({'chr':[chro]*(len(breaks)-1), 
        'start':breaks[0:-1], 
        'stop':breaks[1:],
        'name':np.arange(len(breaks)-1)})

        intervals = intervals.astype({'start':np.int64,'stop':np.int64})
        

    return intervals

def get_full_chromsomes(chros=None, binsize=int(1e6)):
    chr_l = fetch_chromsizes('hg38')
    if chros is None:
        chros=['chr%g'%k for k in range(1,23)]+['chrX']
    intervals = binnify(chr_l, binsize)
    intervals = intervals.loc[intervals['chrom'].isin(chros)].rename({'chrom':'chr','end':'stop'}, axis=1, inplace=False)
    intervals['name']=np.arange(intervals.shape[0])
    intervals['chr_tmp'] = intervals['chr'].map({k:i for i, k in enumerate(chros)})
    intervals=intervals.sort_values(['chr_tmp','start']).drop('chr_tmp', axis=1)
    
    return intervals


def get_tracks(d, genes, intervals): 

    d_slice = d.select_annots(genes.index.to_list())

    chr_indexes = np.unique(intervals['chr'], return_index=True)[1]
    mychr = [intervals['chr'].values[index] for index in sorted(chr_indexes)]

    if len(mychr)>1:
        intervals_sorted=intervals.sort_values(['chr','start'])
        pr = d.df_to_shelf(intervals_sorted, feature_info_id='name')
        d_slice_proj = d_slice.project(pr)
        d_slice_proj, _ = d_slice_proj.reindex_chr(mychr, drop_annots=False)
    else:
        pr = d.df_to_shelf(intervals, feature_info_id='name')
        d_slice_proj = d_slice.project(pr)

    tracks = d_slice_proj.counts.toarray()

    nreads = np.array(d_slice.counts.sum(axis=1))

    ntot = d.counts.sum()
    return tracks, nreads, ntot


def readBigWig(path,intervals, oper="mean"):
    bw = pbw.open(path)
    vals = np.zeros(intervals.shape[0])
    for i in range(intervals.shape[0]):
        vals[i] = bw.stats(intervals.iat[i,0], intervals.iat[i,1], intervals.iat[i,2], type=oper)[0]
    bw.close()
    vals = np.nan_to_num(vals, copy=False)
    return vals

def get_group_(x0, x1, chros):
    n = len(x0)
    ngroups = 1
    groups_max = -1*np.ones(1)
    groups = np.zeros(n, int)
    this_chro = -1
    for i in range(n):
        if chros[i]!=this_chro:
            groups_max = -1*np.ones(ngroups)
            this_chro = chros[i]

        groups_ok = np.flatnonzero(x0[i]>groups_max)
        if len(groups_ok) == 0:
            ngroups +=1 
            groups[i] = (ngroups - 1)
            groups_max = np.append(groups_max, x1[i])
        else:
            groups[i] = groups_ok[0]
            groups_max[groups_ok[0]] = x1[i]

    return groups


def read_genes(genes_df, intervals, nmax = 20):

    genes_to_plot = bf.closest(genes_df.reset_index(), intervals, cols1=("chr","start","stop"), cols2=("chr","start","stop"))
    genes_to_plot_ensg = genes_to_plot['ENSG_1'].loc[genes_to_plot['distance']==0].values

    genes_to_plot = genes_df.loc[genes_to_plot_ensg].sort_values(['chr','start'])

    #genes_to_plot['group']=0
    #genes_to_plot['start']=np.maximum(genes_to_plot['start'], intervals.iloc[0,1])
    #genes_to_plot['stop']=np.minimum(genes_to_plot['stop'], intervals.iloc[-1,2])


    genes_to_plot['group'] = get_group_(genes_to_plot.iloc[:,1], genes_to_plot.iloc[:,2], genes_to_plot.iloc[:,0])
    
    return genes_to_plot

def normalize_tracks(tracks, scale_factor=1.0, nreads=None, intervals=None):
    normalization_factor = scale_factor
    if not nreads is None:
        normalization_factor = normalization_factor * (nreads / 1000.0)
    if not intervals is None:
        normalization_factor = normalization_factor * (intervals['stop']-intervals['start']).values.reshape(1, tracks.shape[1])/ 1000.0
    
    tracks_normalized = tracks / normalization_factor
    
    return tracks_normalized


    
def drawTrack(ax, sig, intervals, color, nticks = 4, unit_size=1e6, zero_offset=True):

    chr_indexes = np.unique(intervals['chr'], return_index=True)[1]
    mychr = [intervals['chr'].values[index] for index in sorted(chr_indexes)]

    if len(mychr)>1:
        start = 0
        stop = intervals.shape[0]+1
        Xs = np.arange(stop-1)

        chr_width=intervals[['chr','start']].groupby('chr').count().reset_index().sort_values('chr', key= lambda h: h.map({k: i for i, k in enumerate(mychr)}))['start'].values

        ticks = np.cumsum(np.hstack([0, np.array(chr_width)]))

    else:
        start = intervals.iat[0,1]
        stop = intervals.iat[-1,2]
        
        delta = 0
        if zero_offset:
            delta = start
            start=0
            stop=stop-start
        
        start=start/unit_size
        stop = stop / unit_size
        Xs = (intervals['start'].values*1.0 - delta) / unit_size
    
        ticks=np.linspace(start, stop, nticks)
    
    ax.set_xlim([start,stop])
    ax.xaxis.set_ticks(ticks)
    ax.tick_params(which='both', direction='out',
       right=False,left=True,top=False,
       labelbottom=False,labelleft=True,
       length=4)
    
    ax.fill_between(Xs,0, sig,color=color) #range(n),0,seqDepth,color=color


def draw_chr_patches(ax, chr_width, ym):
    chr_x = np.cumsum(np.hstack([0, np.array(chr_width)]))[0:-1]
    alternate_chr_x = chr_x[1::2]
    alternate_chr_width = chr_width[1::2]
    #ym=ax.get_ylim()
    alternate_chr_boxes = [Rectangle((x, 0), xw, ym)
                  for x, xw in zip(alternate_chr_x, alternate_chr_width)]
    pc = PatchCollection(alternate_chr_boxes, facecolor='lightgrey', alpha=0.1)
    ax.add_collection(pc)
    ax.vlines(chr_x, 0, ym, color='black',linestyle=":", linewidth=0.5)


def draw_genes(ax, genes_to_plot, intervals, zero_offset = False, unit_size = 1e6):

    ax.tick_params(which='both',
                right=False,left=False,top=False,bottom=False,
                labelbottom=False,labelleft=False)
    cmap=plt.get_cmap("tab20b")

    cmap = {"+":"firebrick", "-":"navy"}

    chr_indexes = np.unique(intervals['chr'], return_index=True)[1]
    mychr = [intervals['chr'].values[index] for index in sorted(chr_indexes)]

    if len(mychr)==1:
        starts = genes_to_plot.iloc[:,1].values
        stops = genes_to_plot.iloc[:,2].values
        
        delta = 0
        if zero_offset:
            delta = intervals.iat[0,1]
            starts= starts - delta
            stops=stops- delta
        
        starts = starts / unit_size
        stops = stops / unit_size

    #ax.hlines(np.random(starts, stops)
    ngroups = np.max(genes_to_plot['group'])+1
    #ax.hlines(-1.0 * (genes_to_plot['group'].values +1), starts, stops, color=genes_to_plot['strand'].map(cmap).values, linewidth = 2)
    
    gps = genes_to_plot['group'].values
    clrs = genes_to_plot['strand'].map(cmap).values
    gnames = genes_to_plot['name'].values
    for i in range(genes_to_plot.shape[0]):
        rect=Rectangle((starts[i],-1.0 * (gps[i]+1)),stops[i]-starts[i],0.4,facecolor=clrs[i],edgecolor='black') #cmap(idx%20))
        ax.add_patch(rect)
        ax.text((starts[i]+stops[i])/2, -1.0 * (gps[i] +1)+0.6, gnames[i], fontsize='small', horizontalalignment='center')
    
    ax.set_ylim([-ngroups-0.5,0.1])

    #     start = starts[i]
    #     length = stops[i]-starts[i]
        
        #rect = Rectangle((start,0),length,1,facecolor=cmap(i%20),edgecolor=None) #cmap(idx%20))
        #ax.add_patch(rect)
    #ax.plot([0,n],[0,0],'k-',linewidth=1,zorder=-1)

def plt_tracks(signals, names, intervals, genes = None, figsize=None, nticks = 4, unit_size=1e6, zero_offset=True, ymax=None, link_y=False, frmt='{0:,}'):
    n = signals.shape[0]
    n0 = n
    if figsize is None:
        figsize = (6,n)
    
    fig = plt.figure(figsize = figsize)
    
    height_ratios = []

    
    for i in range(n0):
        height_ratios.append(0.5)

    if not genes is None:
        genes_to_plot = read_genes(genes, intervals)
        ngroups = np.max(genes_to_plot['group'])+1
        #print(ngroups)
        height_ratios.append(0.7 * ngroups)
        n+=1
    
    #print(n0, n)
    gs = gridspec.GridSpec(figure=fig,ncols=1,nrows=n,
        height_ratios=height_ratios)

    colors=sns.husl_palette(n0,l=0.55)
    

    if ymax is None:
        ymax_vec=np.nanmax(signals, axis=1)
    else:
        if type(ymax)==list:
            ymax_vec=np.array(ymax)
        elif type(ymax)==np.ndarray:
            ymax_vec=ymax.copy()
        else:
            ymax_vec = ymax * np.ones(n0)
    if type(link_y)==bool:
        if link_y:
            ymax_vec = np.max(ymax_vec) * np.ones(n0)
    else:
        for s in set(link_y):
            these_ixs = [i for i, k in enumerate(link_y) if k==s]
            this_max = np.max(ymax_vec[these_ixs])
            ymax_vec[these_ixs]=this_max
        
    axs = []

    chr_indexes = np.unique(intervals['chr'], return_index=True)[1]
    mychr = [intervals['chr'].values[index] for index in sorted(chr_indexes)]
    if len(mychr)>1:
        
        mychr_length=intervals[['chr','start']].groupby('chr').count().reset_index().sort_values('chr', key= lambda h: h.map({k: i for i, k in enumerate(mychr)}))['start'].values

        mychr_length_cum = np.cumsum(np.array(mychr_length))

    for i in range(n0):
        ax = fig.add_subplot(gs[i,0],frame_on=False)
        axs.append(ax)
        
        if len(mychr)>1:
            for ii in range(len(mychr)): 
                #ax.vlines(mychr_length_cum[ii], 0, ymax_vec[i], color='k')
                #ax.pac
                draw_chr_patches(ax, mychr_length, ymax_vec[i])

        drawTrack(ax, signals[i,:], intervals, colors[i], nticks = nticks, unit_size=unit_size, zero_offset=zero_offset)
        ax.set_title(names[i], loc='left', horizontalalignment='left', pad=-14)
        ax.set_ylim([0, ymax_vec[i]])
        ax.yaxis.set_ticks(np.linspace(0, ymax_vec[i], 2))

        
        ylim = ax.get_ylim()
        #print(ylim)

        
                #ax.set_ylim([-ymax_vec[i]/5, ymax_vec[i]])
    x_labels = list(ax.get_xticklabels())
    ticks = ax.get_xticks()
    if len(mychr)>1:
        tick_labels = mychr #+ [""]
        ax.tick_params(axis='x', which='major', labelbottom=False)
        ax.tick_params(axis='x', which='minor', labelbottom=True, bottom=False)
        ax.set_xticks(mychr_length_cum - mychr_length/2.0, minor=True)
        ax.set_xticklabels(tick_labels, ha = 'center',fontsize='small', minor=True, rotation=45)
    else:
        tick_labels = [frmt.format(x) for x in ticks]
#         print(tick_labels)
        ax.set_xticklabels(tick_labels, ha = 'center',fontsize='small')
        ax.tick_params(labelbottom=True)

    if not genes is None:
        #genes_to_plot = read_genes(genes, intervals)
        ax = fig.add_subplot(gs[n-1,0],frame_on=False)
        axs.append(ax)
        
        draw_genes(ax, genes_to_plot, intervals, zero_offset = zero_offset, unit_size = unit_size)
        #ax.set_ylim([-np.max(genes_to_plot['group'])-1,1])
        ax.set_xlim(axs[0].get_xlim())


    gs.tight_layout(fig)



def magic_plot(track_def, intervals, genes, plot_genes = None, plot_order=None, oper_bw="mean", normalize_tgt = True, normalize_src = False, 
figsize=None, nticks = 4, unit_size=1e6, zero_offset=True, ymax=None, link_y=False, frmt='{0:,}'):
    #normalization options: 
    #[chardata..., ['SOZ17'],'ES', sf]
    names = []
    signals = []
    for d in track_def:
        if type(d[0])==str:
            sf=1.0
            if len(d)>2:
                sf=d[2]
            signals+=[readBigWig(d[0], intervals, oper=oper_bw)/sf]
            names+=[d[1]]

        else:
            gene_ids = [genes.loc[genes['name']==n].index.to_list()[0] for n in d[1]]
            names += [k + d[2] for k in d[1]]
            sf=1.0
            if len(d)>3:
                sf=d[3]
            mygenes = genes.reindex(gene_ids)
            track, nreads, ntot = get_tracks(d[0],mygenes, intervals)
            norm_intervals = None
            
            if normalize_tgt:
                norm_intervals = intervals

            if normalize_src:
                track = normalize_tracks(track, nreads=nreads, intervals=norm_intervals, scale_factor=sf)
            else:
                track = normalize_tracks(track, nreads=None, intervals=norm_intervals, scale_factor=sf)
            
            signals+=[track]

    signals = np.vstack(signals)
    if not plot_order is None:
        signals = signals[plot_order, :]
        names = [names[i] for i in plot_order]
    
    plt_tracks(signals, names, intervals, genes = plot_genes, figsize=figsize, nticks = nticks, unit_size=unit_size, zero_offset=zero_offset, ymax=ymax, link_y=link_y, frmt=frmt)
    


# def make_meta_track(d, pr, tracks, beddf, outsuffix, normalize=True, normalize_reads=True, global_divide=1.0):
#     for t, g in tracks.items():
#         if g is None:
#             this_track=d.sum(axis=0).project(pr)
#         else:
#             this_track=d.select_annots(g).sum(axis=0).project(pr)
        
#         this_bw='data-raw/density.'+t+"."+outsuffix
#         this_track.to_bigwig(beddf, this_bw, normalize=normalize, normalize_reads=normalize_reads, global_divide=global_divide)


# def make_meta_track2(d, tracks, outsuffix, global_divide=1.0, comp=False):
#     for t, g in tracks.items():
#         if g is None:
#             this_track=d.sum(axis=0)
#     #this_track.compress()
#         else:
#             this_track=d.select_annots(g).sum(axis=0)
#         this_track.counts=this_track.counts/global_divide
        
#         if comp:
#             this_track.compress()
#         this_bw='data-raw/density.'+t+"."+outsuffix
#         this_track.to_bigwig(this_bw)

#     ids = np.zeros(len(track_def))
#     is_char = np.zeros(len(track_def), dtype=bool)
    
#     id_to_data = {}
#     for i, data in enumerate(track_def):
#         if type(data)==str:
#             ids[i] = id(data)
            
#         else:
#             ids[i] = id(data[0])
#             is_char[i] = True
#             if not ids[i] in id_to_data:

#                 tracks = get_tracks(data[0],genes.loc[genes['name']==data[1]], intervals)

# def magic_plot(track_def, intervals, plot_order, link_y=False, oper_bw="mean"):
    
#     ids = np.zeros(len(track_def))
#     is_char = np.zeros(len(track_def), dtype=bool)
    
#     id_to_data = {}
#     for i, data in enumerate(track_def):
#         if type(data)==str:
#             ids[i] = id(data)
            
#         else:
#             ids[i] = id(data[0])
#             is_char[i] = True
#             if not ids[i] in id_to_data:

#                 tracks = get_tracks(data[0],genes.loc[genes['name']==data[1]], intervals)

    
#     for i, data in enumerate(track_def):


#     ids_unq = {v:i for i, v in enumerate(set(ids))}
#     for this_in in set(ids)
#     gps = np.array([ids_unq[id0] for id0 in ids])

#     signals = []
    

    

#     for i, data in enumerate(track_def):
#         gps_char + [ids_char_unq.get(id(data[0]))  if type(data)==list else -1]


    
#     signals = {}


#     for grp, data in zip(ids, track_def):
#         if not grp in signals:
#             if type(data)==str:
#                 signals{grp} = [readBigWig(data, intervals, oper_bw="mean")]
#             else:

#         else:
            
    
#     #[path, data_ix, y_link_group]

#     [dES,g], [dDE, g]
#     trtools.get_tracks(chardata['dna']['DE']['exons'].copy(),mygenes, intervals)
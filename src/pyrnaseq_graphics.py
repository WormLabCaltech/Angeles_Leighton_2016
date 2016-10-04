# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 17:21:21 2016

A library of useful functions for RNA-seq analysis.

Accepts as input dataframes generated exclusively by sleuth.
May work with others but it certainly isn't optimized to do so.

@author: David Angeles Albores
@contact: dangeles@caltech.edu

A note on nomenclature. This library is meant to work with output from sleuth
I assume you  will have a single dataframe for each regression value in any
model you have.

Any such dataframe is referred to here as 'dfplot', since it has the values to
be plotted.

Crucially, 'dfplot' should have the following columns:
dfplot.b
dfplot.qval
dfplot.ens_gene  <--- must contain WBIDs only!

In order to select broad categories of genes, you will need dataframes that
contain the names of these genes. These dataframes are generically referred to
as dfindex dataframes. There are multiple possible formats for these, but they
all must have the following column:

dfindex.gene <--- WBIDs only!

Broadly speaking, there are three categories of dfindex so far
dfindex -- tissue data
dfindex -- gold standard comparison
dfindex -- genes with an expected effect

As such, there are 3 categories of functions:
Tissue-plotting functions
Gold-standard plotting
Effect plots

Tissue functions
To generate an appropriate tissue dfindex dataframe, simply use the
'organize' function to crawl through all tissues and assemble lists of genes
that are annotated in the tissue of interest.

Beware, 'organize' simply searches
for nodes that share substrings. i.e., if you said 'neuron', you'd get
genes annotated in either 'neuron' and 'dopaminergic neuron'. Likewise,
for 'tail' you may get 'tail neuron' or 'male tail'. So use organize carefully.
Longer 'names' may serve you better (i.e. 'dopaminergic neuron',
instead of just 'neuron').
Use the resulting dataframe to call on all the tissue related functions

In case you want your own dataframe, the df for tissue should have the
following columns:

dfindex.gene -- WBID
dfindex.tissue -- a string denoting a specific anatomic location
dfindex.expressed -- binary, denoting whether a gene is expressed there



Gold-standard functions
Gold-standard dataframes are dataframes that contain a list of genes that
we are interested in querying. In particular, we are mainly interested in
presence/absence of these genes, and the data have typically been generated
by previous publications. The dataframe is fairly small as a result

Columns:
dfindex.gene -- WBID
dfindex.origin -- a string, the name of the dataset this gene came from

Plot-by-value functions
These functions are designed to allow you to compare the value associated with
a list of genes BETWEEN datasets. for example, you may be interested in aging
With a two factor design, you may have an aging coefficient, and a genotype
coefficient. Your mutant of interest might be expected to have positive effects
on aging. As such, you might want to study how the distribution of genes
associated 'positively' with lifespan varies with age, and how it varies by
genotype on the same graph. These functions are designed to do exactly this.

In order to carry out these comparisons, your dataframe must have the
following columns:
dfindex.gene
dfindex.effect -- some categorical variable

"""

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


def wbid_extractor(tissue_df, main_tissue):
    """
    Given a string 'main_tissue', find all columns that
    have a substring equal to it in tissue_df. Then,
    extract all the wbids that are expressed in any of these
    columns and return a non-redundant list of these wbids
    """
    if type(main_tissue) != str:
        raise ValueError('please input a string in main tissue')

    matching = [s for s in tissue_df.columns if main_tissue in s]
    names = []
    for i in matching:
        if len(names) == 0:
            names = tissue_df.wbid[tissue_df[i] == 1].values
        else:
            names = np.append(names, tissue_df.wbid[tissue_df[i] == 1].values)
    names = list(set(names))

    return names


def organize(names, tissue_df):
    """
    Pick your favourite tissues, place them in a list
    provide the tissue dictionary and they will be assembled
    into a dataframe of all tissues that include that word

    INPUT:
    names - a list of tissues
    tissue_df - our standard tissue dictionary

    OUTPUT:
    df1 - a 'tidy' dataframe with three cols: 'wbid' 'tissue' 'expressed'
    """
    # guarantee its iterable
    if type(names) in [str]:
        names = [names]

    names = ['gene']+names

    df = pd.DataFrame(index=np.arange(0, len(tissue_df)), columns=names)
    df.gene = tissue_df.wbid
    for i, value in enumerate(names):
        genes_in_tissue = wbid_extractor(tissue_df, value)
        df[value][df.gene.isin(genes_in_tissue)] = 1

    df.fillna(0, inplace=True)
    df1 = pd.melt(df, id_vars='gene', value_vars=names[1:])
    df1.columns = ['gene', 'tissue', 'expressed']
    return df1


def fix_axes(**kwargs):
    """Makes modifications to axes to ensure proper plotting."""
    title = kwargs.pop('title',  '')
    savename = kwargs.pop('savename', '')
    xlab = kwargs.pop('xlab', r'$\beta$')
    ylab = kwargs.pop('ylab', r'log$_{10}Q$')
    yscale = kwargs.pop('yscale', 'log')
    xscale = kwargs.pop('xscale', 'symlog')
    xlim = kwargs.pop('xlim', [])
    ylim = kwargs.pop('xlim', [])
    loc = kwargs.pop('loc', 1)

    if len(xlim) != 0:
        plt.gca().set_xlim(xlim)
    if len(ylim) != 0:
        plt.gca().set_ylim(ylim)

    plt.gca().legend(fontsize=15.5, loc=loc, frameon=True)
    plt.title(title)
    plt.xlabel(xlab, fontsize=15)
    plt.ylabel(ylab, fontsize=15)
    plt.xscale(xscale)
    plt.yscale(yscale)

    if savename:
        plt.savefig(savename)


def volcano_plot_tissue(tissue, xname, yname, q, dfplot, dfindex, ax,
                        label, col='b', a=.8):
    """
    Plots all the tissue specific genes,i.e. all genes that appear in one
    and only one 'tissue'.


    tissue -- name of the tissue to be plotted
    xname, yname= name of column 'x' and 'y' in dfplot
    dfindex= must be exactly the tissue dataframe generated by organize
    """
    g = lambda x: ((dfindex.expressed == 1) & (dfindex.tissue == x))\
        & (~dfindex[dfindex.expressed == 1].duplicated('gene'))
    f = lambda x: (dfplot.ens_gene.isin(x)) & (dfplot.qval < q)

    gene_selection = g(tissue)
    genes_to_plot = dfindex[gene_selection].gene

    ind = f(genes_to_plot)
    x = dfplot[ind][xname]
    y = dfplot[ind][yname]
    plt.gca().plot(x, -np.log10(y), 'o', color=col, ms=6, alpha=a, label=label)


def explode(x, y, q, dfplot, dfindex, colors, **kwargs):
    """
    A function that makes a complete volcano plot of single-tissue specific
    genes. Warning: Genes annotated in multiple tissues are NOT plotted as
    tissue specific genes at the moment.
    It will plot:
    Genes not statistically significantly expressed
    Genes stat. sig. expressed but not in one of the dfindex lists
    Genes in dfindex (each tissue should have a different color)

    Parameters:
    x= name of the column containing fold-change or regression values for genes
    y= qvalue (NOT negative logged!)
    dfplot= df containing columns x, y
    dfindex= tissue dataframe output by 'organize' function
    colors= a vector of colors. Should have as many colors as 2 + no. tissues
    kwargs= all the kwargs that can be passed to fix-axes
    """
    a = kwargs.pop('a', .8)

    # provide enough colors
    nvals = len(dfindex.tissue.unique()) + 2
    if nvals > colors:
        raise ValueError('len(colors)  is < than len(tissues) +2!')

    # make sure x, y in dfplot
    if not [1 for s in dfplot.columns if s == x][0]:
        raise ValueError('x must be a dataframe column in dfplot')
    if not [1 for s in dfplot.columns if s == y][0]:
        raise ValueError('x must be a dataframe column in dfplot')

    # find all genes not expressed anywhere, or genes that are in multiple
    # 'tissues'
    #
    ind1 = (dfindex.expressed == 0) | \
           (dfindex[dfindex.expressed == 1].duplicated('gene'))
    ind2 = (dfplot.ens_gene.isin(dfindex[ind1].gene)) & (dfplot.qval < q)

    xnotisssig = dfplot[ind2][x]
    ynotisssig = dfplot[ind2][y]

    fig, ax = plt.subplots()
    plt.plot(xnotisssig, -np.log10(ynotisssig), 'o',
             color=colors[0], ms=6, alpha=a, label='all others')

    ind = (dfplot.qval > q)
    xnotsig = dfplot[ind].b
    ynotsig = dfplot[ind].qval

    fig, ax = plt.subplots()
    plt.plot(xnotsig, -np.log10(ynotsig), 'o',
             color=colors[0], ms=6, alpha=a, label='all others')

    values = dfindex.tissue.unique()
    # plot all the points not associated with a tissue
    for i, value in enumerate(values):
        volcano_plot_tissue(value, x, y, q, dfplot, dfindex, label=value,
                            col=colors[i+2], ax=ax, a=a)
    fix_axes(**kwargs)


def volcano_plot_cool_genes(q, x, y, gene_name, genes, dfplot, dfgenes,
                            ax, colors, a=1, dfplotgenes='ens_gene',
                            downsample=1, **kwargs):
    """
    Plots all the tissue specific genes, i.e. all genes that appear in one
    and only one 'tissue'
    """
    f = lambda x: (dfplot[dfplotgenes].isin(x))  # & (dfplot.qval < q)

    nvals = len(dfgenes[gene_name].unique())
    ncolors = len(colors)
    if nvals > ncolors:
        raise ValueError('Please provide as many colors as there are datasets.\
                          {0} {1}'.format(ncolors, nvals))

    for i, gname in enumerate(dfgenes[gene_name].unique()):
        # slice out desired genes
        selected_genes = (dfgenes[gene_name] == gname)
        ind = f(dfgenes[selected_genes][genes])
        sliced = dfplot[ind].dropna().copy()

        # no. genes showing up in assay
        ngsig = len(sliced[(sliced[y] < q)][dfplotgenes].unique())

        ds_select = np.floor(len(sliced)*downsample)
        ds = np.random.randint(0, len(sliced), int(ds_select))
        xcoord = sliced[x].values[ds]
        ycoord = sliced[y].values[ds]

        # no. of genes in dataset
        tg = len(dfgenes[dfgenes[gene_name] == gname][genes].unique())
        label = '{0} {1}= {2}, {3}= {4}'.format(gname, r'$n_{sig}$', 
                                                ngsig, r'$n_{tot}$', tg)
        plt.gca().plot(xcoord, -np.log10(ycoord), 'o',
                       color=colors[i], ms=6, alpha=a, label=label)


def explode_cool_genes(q, x, y, gene_name, genes, dfplot, dfgenes, colors,
                       dfplotgenes='ens_gene', **kwargs):
    """A function that generates all the relevant volcano plots."""
    a = kwargs.pop('a', .6)
    loc = kwargs.pop('loc', 'lower right')
    savename = kwargs.pop('savename', '')
    xlim = kwargs.pop('xlim', '')
    ylim = kwargs.pop('ylim', '')
    downsample = kwargs.pop('downsample', 1)

    # sig genes not in any given df
    ind1 = (~dfplot[dfplotgenes].isin(dfgenes[genes])) & (dfplot[y] < q)
    ind4 = (dfplot[y] < q) # for calculating total no. of sig genes
    # no sig genes
    # ind 2 is for plotting purposes
    # ind 3 is to calculate total number of non-sig genes
    ind2 = (~dfplot[dfplotgenes].isin(dfgenes[genes])) & (dfplot[y] > q)
    ind3 = (dfplot[y] > q)

    xnotsig = dfplot[ind2][x].values
    ynotsig = dfplot[ind2][y].values

    xsig = dfplot[ind1][x].values
    ysig = dfplot[ind1][y].values

    # down_sample selection:
    ds_select = np.floor(len(xnotsig)*downsample)
    ds = np.random.randint(0, len(xnotsig), int(ds_select))

    nnotsig = len(dfplot[ind3][dfplotgenes].unique())
    fig, ax = plt.subplots()
    plt.plot(xnotsig[ds], -np.log10(ynotsig[ds]), 'o',
             color=colors[0], ms=6, alpha=.15,
             label=r'not significant, $n$='+'{0}'.format(nnotsig))

    # down_sample selection:
    ds_select = np.floor(len(xsig)*downsample)
    ds = np.random.randint(0, len(xsig), int(ds_select))

    nsig = len(dfplot[ind4][dfplotgenes].unique())
    plt.plot(xsig[ds], -np.log10(ysig[ds]), 'o',
             color=colors[1], ms=6, alpha=.25,
             label=r'diff. exp. genes $n$=' '{0}'.format(nsig))

    # plot all the points not associated with a tissue
    volcano_plot_cool_genes(q, x, y, gene_name, genes, dfplot,
                            dfgenes, colors=colors[2:], ax=ax, a=a,
                            downsample=downsample)

    fix_axes(loc=loc, **kwargs)
    leg = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=True)
    leg.get_frame().set_facecolor('#EAEAF4')  # seaborn color

    plt.gca().set_xlim(xlim)
    plt.gca().set_ylim(ylim)
    if savename:
        fig.savefig(savename, dpi=300, format='png', bbox_extra_artists=(leg,),
                    bbox_inches='tight')
    
    return fig, ax


def kde_tissue(tissue, q, genes, x, y, dfplot, dfindex, ax, label, col='b'):
    """
    Plots all the tissue specific genes, i.e. all genes that appear in one
    and only one 'tissue'
    tissue -- tissue to plot
    q -- qvalue to slice on
    dfindex -- the dataframe generated by organizer
    dfplot -- the dataframe containing columns x, y and genes
    x -- the name of the column containing the values  to plot in the histogram
    y -- name of the column with which to slice the dataframe (q or p value)
    genes -- the name of the column containing the WBID names
    label -- name of the plot just made
    ax -- axis to plot in
    col -- color
    """
    g = lambda x: ((dfindex.expressed == 1) & (dfindex.tissue == x))
    # & (~dfindex[dfindex.expressed == 1].duplicated('gene'))
    f = lambda x: (dfplot[genes].isin(x)) & (dfplot[y] < q)

    gene_selection = g(tissue)
    genes_to_plot = dfindex[gene_selection].gene

    ind = f(genes_to_plot)
    to_plot = dfplot[ind][x]

    n = len(dfplot[ind][genes].unique())
    if len(to_plot) > 15:
        sns.kdeplot(to_plot, color=col, label=label+' $n$= {0}'.format(n),
                    ax=ax, lw=4, cut=0.5)
        if len(to_plot) <= 20:
            sns.rugplot(to_plot, color=col, ax=ax, height=.07, lw=2)


def tissue_kegg(q, genes, x, y, dfplot, dfindex, colors, **kwargs):
    """
    A function that given a dataframe dfplot and the tissue dataframe dfindex
    generates a KDE plot for each tissue that exists in dfindex

    q -- qvalue to slice on
    dfindex -- the dataframe as generated by organizer
    dfplot -- df containing columns x, y, and gene
    x -- the name of the column containing the values  to plot in the histogram
    y -- name of the column with which to slice the dataframe (q or p value)
    genes -- the name of the column containing the WBID names
    label -- name of the plot just made
    ax -- axis to plot in
    col -- color
    """
    # set scale parameters
    yscale = kwargs.pop('yscale', 'linear')
    xscale = kwargs.pop('xscale', 'linear')
    xlim = kwargs.pop('xlim', [-8, 8])
    ylim = kwargs.pop('ylim', [0, .5])

    ind1 = (dfindex.expressed == 0)
    # | (dftiss[dftiss.value==1].duplicated('gene'))
    ind2 = (dfplot[genes].isin(dfindex[ind1].gene)) & (dfplot[y] < q)

    xnotisssig = dfplot[ind2][x]

    n = len(dfplot[ind2][genes].unique())
    fig, ax = plt.subplots()
    sns.kdeplot(xnotisssig,
                color=colors[0], label='all others $n$= {0}'.format(n), ax=ax,
                lw=4, cut=0.5)
    plt.axvline(0, ls='--', color='black', lw=3)
    # plot all the points not associated with a tissue
    values = dfindex.tissue.unique()
    for i, value in enumerate(values):
        kde_tissue(value, q, genes, x, y, dfplot, dfindex, label=value,
                   col=colors[i+1], ax=ax)

    fix_axes(xscale=xscale, yscale=yscale, xlim=xlim, ylim=ylim, **kwargs)
    return fig, ax

def kde_value(value, q, dfplot, dfindex, ax, label,
              col='b', min_length=10, rug_length=20):
    """Plots all the value specific genes."""
    g = (dfindex.effect == value)
    f = lambda x: (dfplot.ens_gene.isin(x)) & (dfplot.qval < q)

    genes_to_plot = dfindex[g].gene

    ind = f(genes_to_plot)
    x = dfplot[ind].b
    n = len(dfplot[ind].ens_gene.unique())
    if len(x) > min_length:
        sns.kdeplot(x, color=col, label=label+' $n$= {0}'.format(n), ax=ax,
                    lw=5, cut=0.5)

        if len(x) < rug_length:
            sns.rugplot(x, color=col, ax=ax, height=.1, lw=2)
    else:
        print('too few values to plot {0}'.format(label+' n= {0}'.format(n)))


def kegg_compare_byval(value, q, Ldf, dfindex, colors, **kwargs):
    """
    Given a list of dataframes, Ldf, and a list of target genes dfindex,
    compare the distributions of genes within dfindex throughout every
    list for genes with trait 'value'.
    Ldf= a list of dataframes. must have df.ens_gene, df.b, df.qval exactly
    as written.

    dfindex= a list of genes to select, must have df.gene, df.effect
    value= a trait associated with genes in dfindex, (an entry in df.effect)
    colors= an array of len(Ldf) of colors
    """
    dfnames = kwargs.pop('dfnames', ['']*len(Ldf))
    xlim = kwargs.pop('xlim', [-10, 10])
    ylim = kwargs.pop('ylim', [0, 1])
    zeroline = kwargs.pop('zeroline', True)
    xscale = kwargs.pop('xscale', 'linear')
    yscale = kwargs.pop('yscale', 'linear')
    save = kwargs.pop('save', False)

    if len(Ldf) > len(colors):
        raise ValueError('Please provide as many colors as dataframes')

    if len(dfindex[dfindex.effect == value]) == 0:
        raise ValueError('Value \'{0}\' is not contained\
                          within dfindex'.format(value))

    if dfnames:
        if len(Ldf) != len(dfnames):
            raise ValueError('dfnames must be the same length as Ldf')

    fig, ax = plt.subplots()
    for i, df in enumerate(Ldf):
        kde_value(value, q, df, dfindex, ax, dfnames[i], colors[i])

    if zeroline:
        plt.axvline(0, ls='--', color='black', lw=2.5)

    if save:
        sv = '../output/Graphs/effect_'+value
        fix_axes(xlim=xlim, ylim=ylim, xscale=xscale, yscale=yscale,
                 savename=sv, **kwargs)

    else:
        fix_axes(xlim=xlim, ylim=ylim, xscale=xscale, yscale=yscale, **kwargs)


def kegg_compareall_byval(q, Ldf, dfindex, colors, **kwargs):
    """
    Given a list of dataframes Ldf, and a list of selection genes with criteria
    make all the plots of interest
    """
    vals = dfindex.effect.unique()
    titles = kwargs.pop('titles', ['']*len(vals))

    if len(titles) < len(vals):
        errormess = 'There are not enough titles for plots'.format(len(titles),
                                                                   len(vals))
        raise ValueError(errormess)

    for i, value in enumerate(vals):
        if titles[i]:
            title = titles[i]
        else:
            title = value
        kegg_compare_byval(value, q, Ldf, dfindex, colors,
                           title=title, **kwargs)


def line_prepender(filename, line):
    """
    Given a filename, opens it and prepends the line 'line'
    at the beginning of the file
    """
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)

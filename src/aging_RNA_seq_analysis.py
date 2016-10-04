#  -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 08:56:36 2016

@author: davidangeles
"""

#  -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import tissue_enrichment_analysis as tea
import os
import mpl_toolkits.mplot3d
import pyrnaseq_graphics as rsq
from sklearn.preprocessing import StandardScaler

#  Package to perform PCA
import sklearn.datasets
import sklearn.decomposition

sns.set_context("notebook")
mag = 2  # value of beta from regression
qval = .1  # qvalue from regression
qvalEn = 0.05  # q value for enrichment analysis (tissues)

dirLists = '../output/Gene_lists_for_analysis'
if not os.path.exists(dirLists):
    os.makedirs(dirLists)

dirGraphs = '../output/Graphs'
if not os.path.exists(dirLists):
    os.makedirs(dirGraphs)

os.chdir('./')
# gene_lists from sleuth
# tpm vals for PCA
dfTPM = pd.read_csv("../input/tpm_table.csv")
dfTPM.dropna(inplace=True)
# pos beta means high old adults
dfBetaA = pd.read_csv("../input/agebeta_wt.csv")
dfBetaA.dropna(inplace=True)
# pos beta means high in fog2
dfBetaG = pd.read_csv("../input/genotypebeta_wt.csv")
dfBetaG.dropna(inplace=True)
# pos beta means high in fog2-aged
dfBetaAG = pd.read_csv("../input/genotypecrossagebeta_wt.csv")
dfBetaAG.dropna(inplace=True)
# likelihood ratio test results
dfLRT = pd.read_csv("../input/lrt.csv")
dfLRT.dropna(inplace=True)

# sort by target_id
dfBetaA.sort_values('target_id', inplace=True)
dfBetaG.sort_values('target_id', inplace=True)
dfBetaAG.sort_values('target_id', inplace=True)
dfLRT.sort_values('target_id', inplace=True)

# gold standard datasets
dfDaf12 = pd.read_csv('../input/daf12genes.csv')
dfDaf16 = pd.read_csv('../input/daf16genes.csv')
dfLund = pd.read_csv('../input/lund_data.csv', header=None, names=['gene'])
dfEckley = pd.read_csv('../input/eckley_data.csv', header=None, names=['gene'])
dfMurphyUp = pd.read_csv('../input/murphy_data_lifespan_extension.csv')
dfMurphyDown = pd.read_csv('../input/murphy_data_lifespan_decrease.csv')
dfHalaschek = pd.read_csv('../input/Halaschek-Wiener_data.csv')

# gpcrs
dfGPCR = pd.read_csv('../input/all_gpcrs.csv')
dfICh = pd.read_csv('../input/select_ion_transport_genes.csv')
dfAxon = pd.read_csv('../input/axonogenesis_genes.csv')
dfNP = pd.read_csv('../input/neuropeptides.csv')

# gpcr is going to go into a gold standard fxn so add an 'origin' colmn
dfGPCR['origin'] = 'gpcrs'
dfICh['origin'] = 'select ion transport genes'
dfAxon['origin'] = 'axonogenesis genes'
dfNP['origin'] = 'neuropeptide genes'
frames = [dfGPCR, dfICh, dfAxon, dfNP]
dfTargets = pd.concat(frames)

# place all the gold standards in a single dataframe:
dfDaf12['origin'] = 'daf-12'
dfDaf16['origin'] = 'daf-16'
dfEckley['origin'] = 'Eckley'
dfLund['origin'] = 'Lund'
dfMurphyUp['origin'] = 'MurphyExt'
dfMurphyDown['origin'] = 'MurphyDec'
dfHalaschek['origin'] = 'Halaschek'
frames = [dfDaf12, dfDaf16, dfEckley, dfLund,
          dfMurphyDown, dfMurphyUp, dfHalaschek]
dfGoldStandard = pd.concat(frames)

# from wormbase
dfLifespanGenes = pd.read_csv('../input/lifespan gene list complete.csv')

# dfPAN = pd.read_csv()

# tissue dictionary--cite David Angeles et al TEA publication (forthcoming)
# if using the enrichment tool
tissue_df = pd.read_csv("../input/final_cutoff33_threshold0.95_methodany.csv")

# =============================================================================
#  KDE for select tissues
# =============================================================================
# color vector:
colors = ['# 696969', '# e41a1c', '# 377eb8',
          '# 4daf4a', '# 984ea3', '# ff7f00']
tissues = ['P3', 'extracellular', 'tail', 'dopaminergic']
# tissues = ['mu_int', 'sex organ', 'excretory', 'gonadal' ]
# tissues = ['gonad', 'sensillum', 'intestine', 'sex organ' ]
# tissues = ['head', 'tail', 'embryonic' ]
# tissues = ['muscle', 'coel', 'hyp' ]
# tissues = ['sperm', 'extracellular', 'tail', 'dopaminergic' ]
df_exp = rsq.organize(tissues, tissue_df)


colors2 = ['# ffff33', '# e41a1c', '# 377eb8',
           '# 4daf4a', '# 984ea3', '# ff7f00']
genes = 'ens_gene'
x = 'b'
y = 'qval'
rsq.tissue_kegg(qval, genes, x, y,  dfBetaA, df_exp, colors2, title='Aging',
                savename='../output/Graphs/aging_tissue_specific_kde.pdf',
                xlab=r'$\beta_{\mathrm{Aging}}$', ylab='Density')
rsq.tissue_kegg(qval, genes, x, y, dfBetaG, df_exp, colors2, title='Genotype',
                savename='../output/Graphs/genotype_tissue_specific_kde.pdf',
                xlab=r'$\beta_{\mathrm{Genotype}}$', ylab='Density')
rsq.tissue_kegg(qval, genes, x, y, dfBetaAG, df_exp, colors2,
                title='Aging::Genotype',
                savename='../output/Graphs/agingxgenotype_tissue_specific_kde.pdf',
                xlab=r'$\beta_{\mathrm{Aging::Genotype}}$', ylab='Density')
#  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Volcano plots for gold standards
#  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# colors
colors = ['# 999999', '# ffff33', '# f781bf', '# ff7f00',
          '# 984ea3', '# 4daf4a', '# 377eb8', '# e41a1c', '# a65628']

rsq.explode_cool_genes(qval, 'b', 'qval', 'origin', 'gene', dfBetaA,
                       dfGoldStandard, colors=colors, title='Aging',
                       savename='../output/Graphs/aging_with_goldstandards_volcplot.pdf',
                       xlab=r'$\beta_{\mathrm{Aging}}$', xlim=[-12, 12],
                       ylim=[10**-4, 10**2])
rsq.explode_cool_genes(qval, 'b', 'qval', 'origin', 'gene', dfBetaG,
                       dfGoldStandard, colors=colors, title='Genotype',
                       savename='../output/Graphs/genotype_with_goldstandards_volcplot.pdf',
                       xlab=r'$\beta_{\mathrm{Genotype}}$', xlim=[-12, 12],
                       ylim=[10**-4, 10**2])
rsq.explode_cool_genes(qval, 'b', 'qval', 'origin', 'gene', dfBetaAG,
                       dfGoldStandard, colors=colors, title='Aging::Genotype',
                       savename='../output/Graphs/agingcrossgenotype_with_goldstandards_volcplot.pdf',
                       xlab=r'$\beta_{\mathrm{Aging::Genotype}}$',
                       xlim=[-12, 12], ylim=[10**-4, 10**2])
# =============================================================================
#  Volcano plots for target genes
# =============================================================================
colors = ['# 999999', '# ffff33', '# f781bf', '# ff7f00',
          '# 984ea3', '# 4daf4a', '# 377eb8', '# e41a1c', '# a65628']
x = 'b'
y = 'qval'
targ_x = 'effect'
targ_y = 'gene'
rsq.explode_cool_genes(qval, x, y, targ_x, targ_y, dfBetaA, dfTargets,
                       colors=colors, title='Aging',
                       savename='../output/Graphs/aging_targets_volcplot.pdf',
                       xlab=r'$\beta_{\mathrm{Aging}}$',
                       xlim=[-12, 12], ylim=[10**-4, 10**2])
rsq.explode_cool_genes(qval, x, y, targ_x, targ_y, dfBetaG, dfTargets,
                       colors=colors, title='Genotype',
                       savename='../output/Graphs/genotype_targets_volcplot.pdf',
                       xlab=r'$\beta_{\mathrm{Genotype}}$', xlim=[-12, 12],
                       ylim=[10**-4, 10**2])
rsq.explode_cool_genes(qval, x, y, targ_x, targ_y, dfBetaAG, dfTargets,
                       colors=colors, title='Aging::Genotype',
                       savename='../output/Graphs/agingcrossgenotype_targets_volcplot.pdf',
                       xlab=r'$\beta_{\mathrm{Aging::Genotype}}$',
                       xlim=[-12, 12], ylim=[10**-4, 10**2])

# =============================================================================
#  Lifespan genes
# =============================================================================
colors = ['# 999999', '# ffff33', '# e41a1c',
          '# 377eb8', '# 4daf4a', '# 984ea3']
x = 'b'
y = 'qval'
targ_x = 'effect'
targ_y = 'gene'
rsq.explode_cool_genes(qval, x, y, targ_x, targ_y, dfBetaA, dfTargets,
                       colors=colors, title='Aging',
                       savename='../output/Graphs/aging_lifespan_volcplot.pdf',
                       xlab=r'$\beta_{\mathrm{Aging}}$', xlim=[-12, 12],
                       ylim=[10**-4, 10**2])
rsq.explode_cool_genes(qval, x, y, targ_x, targ_y, dfBetaG, dfTargets,
                       colors=colors, title='Genotype',
                       savename='../output/Graphs/genotype_lifespan_volcplot.pdf',
                       xlab=r'$\beta_{\mathrm{Genotype}}$', xlim=[-12, 12],
                       ylim=[10**-4, 10**2])
rsq.explode_cool_genes(qval, x, y, targ_x, targ_y, dfBetaAG, dfTargets,
                       colors=colors, title='Aging::Genotype',
                       savename='../output/Graphs/agingcrossgenotype_lifespan_volcplot.pdf',
                       xlab=r'$\beta_{\mathrm{Aging::Genotype}}$',
                       xlim=[-12, 12], ylim=[10**-4, 10**2])

# =============================================================================
# =============================================================================
#  #  Count how many lifespan genes show up
# =============================================================================
# =============================================================================
# figure out how many genes in dfLIfespan show up in this analysis
f = lambda x: (dfBetaA.ens_gene.isin(x)) & (dfBetaA.qval < qval)
ind = f(dfLifespanGenes.gene.values)
m = dfBetaA[ind].ens_gene.values

f = lambda x: (dfBetaG.ens_gene.isin(x)) & (dfBetaG.qval < qval)
ind = f(dfLifespanGenes.gene.values)
m = np.append(m, dfBetaG[ind].ens_gene.values)

f = lambda x: (dfBetaAG.ens_gene.isin(x)) & (dfBetaAG.qval < qval)
ind = f(dfLifespanGenes.gene.values)
m = np.append(m, dfBetaAG[ind].ens_gene.values)

m = list(set(m))
with open('../output/lifespan_genes_that_show_up.csv', 'w') as f:
    f.write('WBID\n')
    for gene in m:
        f.write(gene)
        f.write('\n')
    f.close()

# =============================================================================
# =============================================================================
#  #
# =============================================================================
# =============================================================================
Ldf = [dfBetaA, dfBetaG, dfBetaAG]  # list of dataframes
dfnames = ['Age', 'Genotype', 'Age::Genotype']
colors = ['# 377eb8', '# e41a1c', '# 4daf4a']
fnames = ['../output/Graphs/positive_aging.pdf',
          '../output/Graphs/negative_aging.pdf',
          '../output/Graphs/variable_aging.pdf',
          '../output/Graphs/unannotated_aging.pdf']
# =============================================================================
#  KDE for Lifespan Genes from Wormbase
# =============================================================================
rsq.kegg_compareall_byval(qval, Ldf, dfLifespanGenes, colors, savenames=fnames,
                          dfnames=dfnames, xscale='symlog',
                          ylab='Density', save=True)

# =============================================================================
#  KDE for Gene Targets
# =============================================================================
dfTargets.columns = ['gene', 'effect']
fnames = ['../output/Graphs/gpcrs.pdf',
          '../output/Graphs/ion_transporters.pdf',
          '../output/Graphs/axonogenesis_genes.pdf',
          '../output/Graphs/neuropeptide_genes.pdf']
colors = ['# 984ea3', '# 4daf4a', '# 377eb8',
          '# e41a1c', '# a65628']

dfTargets.columns = ['gene', 'effect']
rsq.kegg_compareall_byval(qval, Ldf, dfTargets, colors, savenames=fnames,
                          dfnames=dfnames, xscale='symlog', ylab='Density',
                          save=True)

f = lambda x, y: ((x.qval < qval) &
                  (x.ens_gene.isin(dfTargets[dfTargets.effect == y].gene)))

gpcrsA = dfBetaA[f(dfBetaA, 'gpcrs')]
gpcrsG = dfBetaG[f(dfBetaG, 'gpcrs')]
gpcrsAG = dfBetaAG[f(dfBetaAG, 'gpcrs')]
ionA = dfBetaA[f(dfBetaA, 'select ion transport genes')]
ionG = dfBetaG[f(dfBetaG, 'select ion transport genes')]
ionAG = dfBetaAG[f(dfBetaAG, 'select ion transport genes')]
axonA = dfBetaA[f(dfBetaA, 'axonogenesis genes')]
axonG = dfBetaG[f(dfBetaG, 'axonogenesis genes')]
axonAG = dfBetaAG[f(dfBetaAG, 'axonogenesis genes')]
# =============================================================================
#  Identify candidates for RNAi
#  Good targets: small qval, large positive b vals, not previously described
# =============================================================================


def exclude(df, excluded_genes, col):
    ind = (~df[col].isin(excluded_genes))
    return df[ind]


def find_molecular_targets(df, to_be_removed, cols='ens_gene', x='b', q=0.1):
    """
    Given a dataframe df, return a new dataframe that:
    Doesn't have WBIDs present in the exclude series
    Has only genes that have q value < qval
    """
    if cols in [str, int, float]:
        cols = [cols]

    df.sort_values(x, inplace=True)
    sig = (df.qval < q)  # take only sig genes

    temp = df[sig].copy()  # remove all non-sig genes and make a temp copy
    if isinstance(to_be_removed, list):
        for i, excluded_gene_list in enumerate(to_be_removed):
            temp = exclude(temp, excluded_gene_list, cols[i])
    return temp


def direction_specific_tissue_analysis(anames, fnames, df, Lind,
                                       genes='ens_gene'):
    """
    Given a single dataframe, perform a tissue analysis for genes
    using the selection indices in Lind
    """
    if not isinstance(anames, list):
        raise ValueError('anames must be a list!')
    if not isinstance(fnames, list):
        raise ValueError('fnames must be a list!')
    if len(anames) != len(fnames):
        raise ValueError('fnames and anames must match length')

    for i, aname in enumerate(anames):
        fname = fnames[i]
        ind = inds[i]
        df_results, unused = tea.enrichment_analysis(df[ind][genes], tissue_df,
                                                     qvalEn, show=False)

        # save results to csv file
        df_results.to_csv('../output/EnrichmentAnalysisResults/'+fname,
                          index=False)
        # reopen the file and add a comment with relevant info for the file
        line = '# ' + aname+'\n'
        rsq.line_prepender('../output/EnrichmentAnalysisResults/'+fname, line)
        # plot top fold change tissues
        tea.plot_enrichment_results(df_results, title=aname,
                                    dirGraphs=dirGraphs)
        plt.close()


# set the path to sve the rnai candidates
# and make sure to exclude known genes from analysis
path = '../output/RNAi Candidates/'
excluded1 = pd.Series(dfLifespanGenes.gene.append(dfGoldStandard.gene).unique())
# =============================================================================
# =============================================================================
# also exclude genes that have significant betas in any other list
x = dfBetaG[dfBetaG.qval < qval].target_id
y = dfBetaAG[dfBetaAG.qval < qval].target_id
excluded2 = pd.Series(x.append(y).unique())
excluded = [excluded1, excluded2]
cols = ['ens_gene', 'target_id']

aging_set = find_molecular_targets(dfBetaA, excluded, cols, q=qval)
aging_set.to_csv('../output/AgingGeneSet.csv')
aging_set.tail(55).to_csv(path + 'CandidatesAge_HighInOld.csv')
aging_set.head(55).to_csv(path + 'CandidatesAge_LowInOld.csv')

# enrichment on genes associated only with aging
aging_set = find_molecular_targets(dfBetaA, [excluded2],
                                   ['target_id'], q=qval)
aname1 = 'aging_up (genes assoc. with aging only)'
fname1 = 'aging_up.csv'
aname2 = 'aging_down (genes assoc. with aging only)'
fname2 = 'aging_down.csv'
anames = [aname1, aname2]
fnames = [fname1, fname2]
inds = [(aging_set.b > 0), (aging_set.b < 0)]
direction_specific_tissue_analysis(anames, fnames, aging_set, inds)
# =============================================================================
# =============================================================================
# now for genotype, same thing
x = dfBetaA[dfBetaA.qval < qval].target_id
y = dfBetaAG[dfBetaAG.qval < qval].target_id
excluded2 = pd.Series(x.append(y).unique())
excluded = [excluded1, excluded2]
genotype_set = find_molecular_targets(dfBetaG, excluded, cols, q=qval)
genotype_set.to_csv('../output/GenotypeGeneSet.csv')
genotype_set.tail(55).to_csv(path + 'CandidatesGenotype_HighInOld.csv')
genotype_set.head(55).to_csv(path + 'CandidatesGenotype_LowInOld.csv')


genotype_set = find_molecular_targets(dfBetaG, [excluded2],
                                      ['target_id'], q=qval)
aname1 = 'genotype up (genes assoc. with genotype only)'
fname1 = 'genotype up.csv'
aname2 = 'genotype down (genes assoc. with genotype only)'
fname2 = 'genotype down.csv'
anames = [aname1, aname2]
fnames = [fname1, fname2]
inds = [(genotype_set.b > 0), (genotype_set.b < 0)]
direction_specific_tissue_analysis(anames, fnames, genotype_set, inds)

# interactions
# more than linear increase in age, up in fog2 more than wt
# only take genes that go way up during age
ind = (dfBetaA.ix[dfBetaAG[dfBetaAG.qval < qval].index].b >
       dfBetaA.b.quantile(.9))
ind2 = (dfBetaAG.target_id.isin(dfLRT[dfLRT.qval < qval].target_id))
interaction_set = dfBetaAG[(dfBetaAG.qval < qval) & ind2]

interaction_set.tail(55).to_csv('../output/RNAi Candidates/CandidatesAgingXGenotypeMoreThanExp.csv')
dfBetaAG[dfBetaAG.qval < qval][ind & ind2].head(55).to_csv('../output/RNAi Candidates/CandidatesAgingXGenotypeLessThanExp.csv')

aname1 = 'genotypeXaging up (genes assoc. with a..g only)'
fname1 = 'agingCrossgenotype up.csv'
aname2 = 'genotypeXaging down (genes assoc. with a..g only)'
fname2 = 'agingCrossgenotype down.csv'
anames = [aname1, aname2]
fnames = [fname1, fname2]


inds = [(dfBetaAG[(dfBetaAG.qval < qval) & ind2].b > 0),
        (dfBetaAG[(dfBetaAG.qval < qval) & ind2].b < 0)]

sig = (dfBetaAG.qval < qval)
find_all_with_interaction_term = (dfBetaA.target_id.isin(
                                  dfBetaAG[sig & ind2].target_id))

direction_specific_tissue_analysis(anames, fnames,
                                   dfBetaAG[(dfBetaAG.qval < qval) & ind2],
                                   inds)
#
# s1 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[0]].est_counts.values
# s2 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[1]].est_counts.values
# s3 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[2]].est_counts.values
# s4 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[3]].est_counts.values
# s5 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[4]].est_counts.values
# s6 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[5]].est_counts.values
# s7 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[6]].est_counts.values
# s8 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[7]].est_counts.values
# s9 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[8]].est_counts.values
# s10 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[9]].est_counts.values
# s11 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[10]].est_counts.values
# s12 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[11]].est_counts.values
#
# s1 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[0]].tpm.values
# s2 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[1]].tpm.values
# s3 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[2]].tpm.values
# s4 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[3]].tpm.values
# s5 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[4]].tpm.values
# s6 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[5]].tpm.values
# s7 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[6]].tpm.values
# s8 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[7]].tpm.values
# s9 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[8]].tpm.values
# s10 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[9]].tpm.values
# s11 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[10]].tpm.values
# s12 = dfTPM[dfTPM['sample'] == dfTPM['sample'].unique()[11]].tpm.values
#
# mat = np.matrix([s1, s2, s3,
#                  s4, s5, s6,
#                  s7, s8, s9,
#                  s10, s11, s12]).transpose()
# sklearn_pca = sklearn.decomposition.PCA(n_components=3)
# sklearn_pca.fit(mat)
# x, y, z = sklearn_pca.components_
# #  Project the data into this 2D space and convert it back to a tidy dataframe
#
# colors = ['# e41a1c', '# 377eb8', '# 4daf4a', '# 984ea3']
# sample_label = ['mt a', 'mt a', 'mt a',
#                 'mt y', 'mt y', 'mt y',
#                 'wt a', 'wt a', 'wt a',
#                 'wt y', 'wt y', 'wt y']
# fig = plt.figure(1, figsize=(8, 6))
# ax = mpl_toolkits.mplot3d.Axes3D(fig)
# for i in np.arange(0, 12):
#     j = int(np.floor(float(i)/3))
#     print(i, j)
#     ax.plot([x[i], x[i]], [y[i], y[i]], [z[i], z[i]],
#             'o', color=colors[j], label=sample_label[i])
#
# ax.set_title("First three PC directions")
# ax.set_xlabel("PC 1")
# ax.set_ylabel("PC 2")
# ax.set_zlabel("PC 3")
# ax.legend(loc='upper left', fontsize=15)
#
#
# fig, ax = plt.subplots()
# for i in np.arange(0, 12):
#     j = int(np.floor(float(i)/3))
#     print(i, j)
#     ax.plot(x[i], y[i], 'o', color=colors[j], label=sample_label[i])
#
# ax.set_title("First three PC directions")
# ax.set_xlabel("PC 1")
# ax.set_ylabel("PC 2")
# ax.legend(loc='upper left', fontsize=15)
#
#
# # =============================================================================
# #  ''normalized'' pca whatever that means...
# # =============================================================================
# mat = np.matrix([s1, s2, s3,
#                  s4, s5, s6,
#                  s7, s8, s9,
#                  s10, s11, s12]).transpose()
# sklearn_pca = sklearn.decomposition.PCA(n_components=3)
# mat_std = StandardScaler().fit_transform(mat)
# sklearn_pca.fit_transform(mat_std)
# x, y, z = sklearn_pca.components_
# #  Project the data into this 2D space and convert it back to a tidy dataframe
#
# colors = ['# e41a1c', '# 377eb8', '# 4daf4a', '# 984ea3']
# sample_label = ['mt a', 'mt a', 'mt a',
#                 'mt y', 'mt y', 'mt y',
#                 'wt a', 'wt a', 'wt a',
#                 'wt y', 'wt y', 'wt y']
# fig = plt.figure(1, figsize=(8, 6))
# ax = mpl_toolkits.mplot3d.Axes3D(fig)
# for i in np.arange(0, 12):
#     j = int(np.floor(float(i)/3))
#     print(i, j)
#     ax.plot([x[i], x[i]], [y[i], y[i]], [z[i], z[i]],
#             'o', color=colors[j], label=sample_label[i])
#
# ax.set_title("First three PC directions")
# ax.set_xlabel("PC 1")
# ax.set_ylabel("PC 2")
# ax.set_zlabel("PC 3")
# ax.legend(loc='upper left', fontsize=15)
#
#
# fig, ax = plt.subplots()
# for i in np.arange(0, 12):
#     j = int(np.floor(float(i)/3))
#     if float(i) % 3 == 0:
#         ax.plot(x[i], y[i], 'o', color=colors[j], label=sample_label[i])
#     else:
#         ax.plot(x[i], y[i], 'o', color=colors[j])
#
# ax.set_title("First three PC directions")
# ax.set_xlabel("PC 1")
# ax.set_ylabel("PC 2")
# ax.legend(loc='upper left', fontsize=15)
#
# fig, ax = plt.subplots()
# for i in np.arange(0, 12):
#     j = int(np.floor(float(i)/3))
#     if float(i) % 3 == 0:
#         ax.plot(x[i], z[i], 'o', color=colors[j], label=sample_label[i])
#     else:
#         ax.plot(x[i], z[i], 'o', color=colors[j])
#
# ax.set_title("First three PC directions")
# ax.set_xlabel("PC 1")
# ax.set_ylabel("PC 2")
# ax.legend(loc='upper left', fontsize=15)
#
#
# fig, ax = plt.subplots()
# for i in np.arange(0, 12):
#     j = int(np.floor(float(i)/3))
#     if float(i) % 3 == 0:
#         ax.plot(z[i], y[i], 'o', color=colors[j], label=sample_label[i])
#     else:
#         ax.plot(z[i], y[i], 'o', color=colors[j])
#
# ax.set_title("First three PC directions")
# ax.set_xlabel("PC 1")
# ax.set_ylabel("PC 2")
# ax.legend(loc='upper left', fontsize=15)
#
#
# # =============================================================================
# # =============================================================================
# #  #  Regression examples
# # =============================================================================
# # =============================================================================
# # Make figures to show what regressions look like.....
# def example(B, x, y):
#     """
#     just a function to show a contour plot example for a linear regression with
#     interactions
#     """
#     b1, b2, b3 = B  # unpack
#     xx, yy = np.meshgrid(x, y)    # mesh
#     out = b1*xx + b2*yy + b3*xx*yy  # model
#     # output a normalized model that always has max absolute value 1 and
#     # min of 0
#     return (out+np.abs(np.min(out)))/np.max(np.abs(out))
#
# x = np.linspace(0, 10)
# y = np.linspace(0, 10)
# xx, yy = np.meshgrid(x, y)
# # set up the figure
# fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharey=True, figsize=(16, 4))
# cmap = sns.light_palette("navy", reverse=True, as_cmap=True)
#
# # positive interaction
# regress = example([1, 1, .5], x, y)
# cs = ax1.contourf(xx, yy, regress, cmap=cmap, alpha=0.7)
# ax1.set_title('Linear Model w +ive interactions')
# ax1.set_ylabel('y', fontsize=15)
#
# # no interaction
# regress = example([1, 1, 0], x, y)
# ax2.contourf(xx, yy, regress, cmap=cmap, alpha=0.7)
# ax2.set_title('Linear Model w/o interactions')
# ax2.set_xlabel('x', fontsize=15)
#
# # negative interaction
# regress = example([1, 1, -.5], x, y)
# ax3.contourf(xx, yy, regress, cmap=cmap, alpha=0.7)
# ax3.set_title('Linear Model w -ive interactions')
# fig.colorbar(cs, ax=ax3, shrink=0.9)
# # save
# plt.savefig('../output/ExamplesRegression.pdf')

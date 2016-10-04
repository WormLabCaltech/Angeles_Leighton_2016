# -*- coding: utf-8 -*-
"""
A script to analyze RNAi brood size assays.

author: David Angeles
contact: dangeles@caltech.edu
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

input_path = '../input/rnai_screen_results/'

df = pd.read_csv(input_path + 'brood_assay_screen.csv')
names = pd.read_csv(input_path + 'rnai_genes_dict.csv')

# Drop all useless columns that excel makes
# also drop any and all entries that have NaNs
names.drop('Unnamed: 4', 1, inplace=True)
names.dropna(0, 'any', inplace=True)

# make all codes upper or lower, not both
names.code = names.code.apply(str.lower)
df.rnai = df.rnai.apply(str.lower)

# extract the names that have been assayed so far
translate = lambda x: names[names.code == x].gene_name.values[0]
df['gene'] = df.rnai.apply(translate)

#fill all NANs with the mean value for that column
df.fillna(df.mean(), inplace=True)



df['total'] = df.d1 + df.d2 + df.d3
df.boxplot('total', by='rnai')
plt.show()

df.boxplot('total', by='gene')
plt.show()


sns.swarmplot(x='gene', y='total', data=df, size=10)
plt.show()

sns.swarmplot(x='rnai', y='total', data=df, size=10)
plt.show()


# trim outliers

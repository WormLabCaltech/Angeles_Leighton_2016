# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 13:28:07 2014
A script to re-format a file
@author: davidaalbores
"""
from __future__ import division, print_function, absolute_import
import pandas as pd

#Filenames!
fileInput= '../input/lifespan gene list complete.csv'; sepInput= ','; comInput= '#';
newFile= '../input/eckley_data.csv'

filename= '../input/c_elegans.PRJNA13758.WS241.livegeneIDs.unmaprm.txt'
separator= "\t"
comments= "#"

#Load the gene name dictionary
names= pd.read_csv(filename, sep= separator, comment= comments)

#Load the control database
df= pd.read_csv(fileInput, sep= sepInput, comment= comInput)


names[(names.HumanReadable.isin(df.gene)) | (names.GeneName.isin(df.gene))].WBID
df= df[df.gene.isin(names.HumanReadable)]
df['wbid']= names[(names.HumanReadable.isin(df.gene)) | (names.GeneName.isin(df.gene))].WBID#. to_csv(newFile, index= False)
#df[df.gene.isin(names.HumanReadable)].gene= names[(names.HumanReadable.isin(df.gene))].WBID.values


#df= df.dropna()


#Write to main file
#df.to_csv(newFile, sep=',', index= False)
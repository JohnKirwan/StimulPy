# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 15:54:26 2019

@author: John D. Kirwan
"""
# make an object containing every line
with open('test_data.dat','r') as specs:
    lines = specs.readlines()
specs.close() # immediately close connection

# Need to make a dataframe to
import pandas as pd
specdf = pd.DataFrame(columns=["material","dark","rep","wv","specrad"])
def num1dec (str):
    if "NAN" in str:
        str = "0"
    return round(float(str),1)

# run through the lines to extract the spec data
i = 0 # initialize counter at line 0
j = 0 # initialize row of spec df
while i < len(lines): # not <= because indexing begins at 0

    if  lines[i].__contains__("Comment"): # get spec description
        # assumes naming convention: material_dark_rep
        [material,dark,rep] = lines[i][21:-1].split('_')
        # subset string to get description
        i = i + 30 # skip down to start of spec data
        while lines[i].__contains__("[END]") == False: # while still spec data
            [wv, specrad] = lines[i].split(' ')[1:3]
            [wv, specrad] = map(num1dec, [wv, specrad]) # make numeric
            specdf.loc[j] = [material, dark, rep, wv, specrad]
            j = j + 1   # iterate j to move to next df row
            i = i + 1   # iterate i to move to next line of text

    i = i + 1  # iterate through lines in while loop

specdf.to_csv('spec_data.txt', sep='\t')
print('Dataframe is added to file')

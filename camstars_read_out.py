#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 13:17:04 2019

@author: Alex Hackett
An all-purpose out plotting package for Cambridge Stars
Reads in the Out file and allows the model to be read in model number by number
and for the meshpoint by meshpoint data to be utilized
"""
import numpy as np
import scipy as sp
import matplotlib.pylab as plt
import astropy as ap
import os
import pandas as pd
from tqdm import tqdm
import glob



class outReader:
    def __init__(self, file_dir):
        self.file_dir = file_dir
        
    def readOut(self):
        #split the out file
        if glob.glob(self.file_dir + '/src/xx*') == []:
            os.system("csplit -s -n 9 " + self.file_dir+'/out' + " '% K.*%' '/ K.*/-1' '{20}'")
        #Glob all of the xx split files
        self.xfiles = glob.glob(self.file_dir + '/src/xx*')
        self.model_data_array = []
        self.current_df1 = 0
        self.current_df2 = 0
        self.current_df3 = 0
        i = 0
        while i <= (len(self.xfiles)):
            print('Reading Model: ', i/3)
            print(np.genfromtxt(self.xfiles[i]))
            #self.current_df1 = pd.read_csv(self.xfiles[i], sep = ' ')
            #self.current_df2 = pd.read_table(self.xfiles[i+1], sep = ' ')
            #self.current_df3 = pd.read_table(self.xfiles[i+2], sep = ' ')
            #self.current_df3.drop(self.current_df3.tail(12).index, inplace=True)
            i += 3
            #self.model_data_array.append(pd.concat([self.current_df3,self.current_df2,self.current_df1], axis = 1))
            #self.model_data_array  =self.model_data_array self.current_df1
            
        


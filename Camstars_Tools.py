#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
import matplotlib.pylab as plt
import astropy as ap
import os
import pandas as pd

class CamRun:
    def __init__(self, run_dir):
        self.run_dir = str(run_dir)
        self.data_file_name = self.run_dir + '/data'
        self.modin_file_name = self.run_dir + '/modin'
        
    def readData(self):
        '''
        Opens the data file in the run dir and reads all the vars into a
        pandas dataframe
        '''
        df = pd.DataFrame()
        with open(self.data_file_name, 'r') as f:
            for i, line in enumerate(f):
                if i < 23:
                    df = pd.concat( [df, pd.DataFrame([tuple(line.strip().split(' '))])], ignore_index=True )
        self.datafile = df
        #Reading all the data vars in
        self.data_NH2 = self.datafile[0][0]
        self.data_ITER1 = self.datafile[1][0]
        self.data_ITER2 = self.datafile[2][0]
        self.data_JIN = self.datafile[3][0]
        self.data_JOUT = self.datafile[4][0]
        self.data_NCH = self.datafile[5][0]
        self.data_JP = self.datafile[6][0]
        self.data_ITH = self.datafile[7][0]
        self.data_IX = self.datafile[8][0]
        self.data_IY = self.datafile[9][0]
        self.data_IZ = self.datafile[10][0]
        self.data_IMODE = self.datafile[11][0]
        
        self.data_ICL = self.datafile[0][1]
        self.data_ION = self.datafile[1][1]
        self.data_IAM = self.datafile[2][1]
        self.data_IOP = self.datafile[3][1]
        self.data_INUC = self.datafile[4][1]
        self.data_IBC = self.datafile[5][1]
        self.data_ICN = self.datafile[6][1]
        self.data_IML1 = self.datafile[7][1]
        self.data_IML2 = self.datafile[8][1]
        self.data_ISGTH = self.datafile[9][1]
        self.data_IMO = self.datafile[10][1]
        self.data_IDIFF = self.datafile[11][1]
        
        self.data_NWRT1 = self.datafile[0][2]
        self.data_NWRT2 = self.datafile[1][2]
        self.data_NWRT3 = self.datafile[2][2]
        self.data_NWRT4 = self.datafile[3][2]
        self.data_NWRT5 = self.datafile[4][2]
        self.data_NSAVE = self.datafile[5][2]
        self.data_NMONT = self.datafile[6][2]
        
        self.data_EPS = self.datafile[0][3]
        self.data_DEL = self.datafile[1][3]
        self.data_DHO = self.datafile[2][3]
        self.data_DT3 = self.datafile[3][3]
        self.data_DDD = self.datafile[4][3]
        
        self.data_NE1 = self.datafile[0][4]
        self.data_NE2 = self.datafile[1][4]
        self.data_NE3 = self.datafile[2][4]
        self.data_NB = self.datafile[3][4]
        self.data_NEV = self.datafile[4][4]
        self.data_NF = self.datafile[5][4]
        self.data_J1 = self.datafile[6][4]
        self.data_J2 = self.datafile[7][4]
        self.data_IH = self.datafile[8][4]
        self.data_JH = self.datafile[9][4]
        
        self.data_DT1 = self.datafile[0][15]
        self.data_DT2 = self.datafile[1][15]
        
        self.data_ZS = self.datafile[0][16]
        self.data_ALPHA = self.datafile[1][16]
        self.data_CH = self.datafile[2][16]
        self.data_CC = self.datafile[3][16]
        self.data_CN = self.datafile[4][16]
        self.data_CO = self.datafile[5][16]
        self.data_CNE = self.datafile[6][16]
        self.data_CMG = self.datafile[7][16]
        self.data_CSI = self.datafile[8][16]
        self.data_CFE = self.datafile[9][16]
        
        self.data_RCD = self.datafile[0][17]
        self.data_OS = self.datafile[1][17]
        self.data_RML = self.datafile[2][17]
        self.data_RMG = self.datafile[3][17]
        self.data_ECA = self.datafile[4][17]
        self.data_XF = self.datafile[5][17]
        self.data_DR = self.datafile[6][17]
        
        self.data_RMT = self.datafile[0][18]
        self.data_RHL = self.datafile[1][18]
        self.data_AC = self.datafile[2][18]
        self.data_AK1 = self.datafile[3][18]
        self.data_AK2 = self.datafile[4][18]
        self.data_ECT = self.datafile[5][18]
        self.data_TRB = self.datafile[6][18]
        
        self.data_IRAM = self.datafile[0][19]
        self.data_IRS1 = self.datafile[1][19]
        self.data_VROT1 = self.datafile[2][19]
        self.data_IRS2 = self.datafile[3][19]
        self.data_VROT2 = self.datafile[4][19]
        self.data_FMAC = self.datafile[5][19]
        self.data_FAM = self.datafile[6][19]
        
        self.data_IVMC = self.datafile[0][20]
        self.data_TRC1 = self.datafile[1][20]
        self.data_IVMS = self.datafile[2][20]
        self.data_TRC2 = self.datafile[3][20]
        self.data_MWTS = self.datafile[4][20]
        self.data_IAGB = self.datafile[5][20]
        self.data_ISGFAC = self.datafile[6][20]
        self.data_FACSGIM = self.datafile[7][20]
        self.data_SGTHFAC = self.datafile[8][20]
        
        self.data_TENG = self.datafile[0][21]
        self.data_SMASS = self.datafile[1][21]
        self.data_FMASS = self.datafile[2][21]
        self.data_INJMD = self.datafile[3][21]
        self.data_STARTTIMEINJ = self.datafile[4][21]
        self.data_ENDTIMEINJ = self.datafile[5][21]
        
        self.data_ENDAGE = self.datafile[0][22]
        
        
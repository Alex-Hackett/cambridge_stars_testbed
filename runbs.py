#!/usr/bin/env python3

"""
Created on Wed Jan  8 16:15:49 2020

@author: ajh291
This is the script that runs the cambridge stars code, allowing you to select 
the modin, nucmodin, and data files that you wish to utilize, then asking you
where you would like to have the output files placed, and, most importantly
then performing the compuation, and running the liveplotter!
"""

import os
import sys
import subprocess
import numpy as np
import matplotlib.pylab as plt
import scipy as sp
#import easygui as eg
#from tkinter import tk
from tkinter.filedialog import askopenfilename, askdirectory


#System Calls
#set the stars working dir, just in case
stars_dir = os.path.dirname(os.path.abspath(__file__))
#Input files
data_filename = askopenfilename(initialdir=stars_dir, title='Select Datafile')
modin_filename = askopenfilename(initialdir=stars_dir, title = 'Select Modin File')
nucmodin_filename = askopenfilename(initialdir=stars_dir, title = 'Select Nuclear Synthesis File')
tzo_data_filename = askopenfilename(initialdir=stars_dir, title = 'Select TZO Control File')

#output files
out_dirname = askdirectory(initialdir=stars_dir, title = 'Select Output Directory')

subprocess.call('rm -f fort.*',shell=True)

#Link the code parameters for the run
os.symlink(data_filename, stars_dir + '/fort.1')
os.symlink(tzo_data_filename, stars_dir + '/fort.70')

#Link the physical data files, if COtables were chosen in the datafile,
#printa.f in the code handles this automatically
os.symlink(stars_dir + '/dat/phys02.dat', stars_dir + '/fort.11')
os.symlink(stars_dir + '/dat/nrate.dat', stars_dir + '/fort.13')
os.symlink(stars_dir + '/dat/splinecoefficients.dat', stars_dir + '/fort.14')

#Link inputs
os.symlink(modin_filename, stars_dir + '/fort.30')
os.symlink(nucmodin_filename, stars_dir + '/fort.31')

#Link outputs
os.symlink(out_dirname + '/out.out', stars_dir + '/fort.32')
os.symlink(out_dirname + '/plot.out', stars_dir + '/fort.33')

os.symlink(out_dirname + '/modout.out', stars_dir + '/fort.34')
os.symlink(out_dirname + '/nucmodout.out', stars_dir + '/fort.35')

os.symlink(out_dirname + '/syntha.out', stars_dir + '/fort.36')
os.symlink(out_dirname + '/synthb.out', stars_dir + '/fort.37')
os.symlink(out_dirname + '/synthc.out', stars_dir + '/fort.38')
os.symlink(out_dirname + '/surface.out', stars_dir + '/fort.39')
os.symlink(out_dirname + '/centre.out', stars_dir + '/fort.40')
os.symlink(out_dirname + '/sprocess.out', stars_dir + '/fort.41')
os.symlink(out_dirname + '/montage.out', stars_dir + '/fort.42')
os.symlink(out_dirname + '/adata.out', stars_dir + '/fort.45')
os.symlink(out_dirname + '/extra.out', stars_dir + '/fort.47')

#os.symlink(out_dirname + '/modin_last.in', stars_dir + '/fort.134')
#os.symlink(out_dirname + '/nucmodin_last.nucin', stars_dir + '/fort.135')


plotcaller = 'python ' + stars_dir + '/liveplotter.py ' + out_dirname + '/plot.out' + ' &'
subprocess.call(plotcaller,shell=True)
subprocess.call(stars_dir + '/bs',shell=True)

subprocess.call('mv ' + stars_dir + '/temp_last_modin ' + out_dirname + '/modin_last.in',shell=True)
subprocess.call('mv ' + stars_dir + '/temp_last_nucmodin ' + out_dirname + '/nucmodin_last.nucin',shell=True)

subprocess.call('rm -f fort.*', shell=True)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 17:46:14 2019

@author: ajh291

A live plotting utility for Cambridge Stars, to be bundled with it
Takes some of stars' live updating tools, namely plot, centre and surface
This will open up multiple, easy to manipulate windows, each of which will
display four plots, giving elemental information, ages, integrated values etc
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import os.path
import matplotlib.gridspec as gridspec
import pandas as pd
import sys

def animate(i):
    '''
    Main animation loop
    '''
    #Axis labels
    plot_axname = ["Model Number", "Age (years)", "Radius (log)",
                   "Effective Temperature (log)", "Luminosity (log)", "Mass",
                   "Mass of He Core", "Mass of C Core", "H Luminosity (log)", 
                   "He Luminosity (log)", "C Luminosity (log)",
                   "Mass Co-oridnate of Convective Boundary",
                   "Mass Co-oridnate of Convective Boundary",
                   "Mass Co-oridnate of Convective Boundary",
                   "Mass Co-oridnate of Convective Boundary",
                   "Mass Co-oridnate of Convective Boundary",
                   "Mass Co-oridnate of Convective Boundary",
                   "Mass Co-oridnate of Convective Boundary",
                   "Mass Co-oridnate of Convective Boundary",
                   "Mass Co-oridnate of Convective Boundary",
                   "Mass Co-oridnate of Convective Boundary",
                   "Mass Co-oridnate of Convective Boundary",
                   "Mass Co-oridnate of Convective Boundary",
                   "Mass Co-ordinate of Max H Energy Generation",
                   "Mass Co-ordinate of Max He Energy Generation",
                   "Opacity (log)", "Timestep", "Surface Abundance H",
                   "Surface Abundance He", "Surface Abundance C",
                   "Surface Abundance N", "Surface Abundance O",
                   "Surface Abundance 3He",
                   "Radius/Roche Lobe Radius","Spin Angular Momentum (Star 1)",
                   "Binar Period", "Binary Separation", "Binary Mass",
                   "Orbital Angular Momentum",
                   "Total Spin Angular Momentum (Binary)",
                   "Total Angular Momentum", "Angular Frequncy of Orbit",
                   "Angular Frequency of Star", "Moment of Inertia of Star",
                   "Moment of Intertia of Orbit", "Mass Loss Rate",
                   "Mass Location of Shell Boundary",
                   "Mass Location of Shell Boundary",
                   "Mass Location of Shell Boundary",
                   "Mass Location of Shell Boundary",
                   "Mass Location of Shell Boundary",
                   "Mass Location of Shell Boundary",
                   "Mass Location of Shell Boundary",
                   "Mass Location of Shell Boundary",
                   "Mass Location of Shell Boundary",
                   "Mass Location of Shell Boundary",
                   "Mass Location of Shell Boundary",
                   "Mass Location of Shell Boundary",
                   "Thermohaline Mixing Boundary Region",
                   "Thermohaline Mixing Boundary Region",
                   "Thermohaline Mixing Boundary Region",
                   "Thermohaline Mixing Boundary Region",
                   "Thermohaline Mixing Boundary Region",
                   "Thermohaline Mixing Boundary Region",
                   "Thermohaline Mixing Boundary Region",
                   "Thermohaline Mixing Boundary Region",
                   "Thermohaline Mixing Boundary Region",
                   "Thermohaline Mixing Boundary Region",
                   "Thermohaline Mixing Boundary Region",
                   "Thermohaline Mixing Boundary Region",
                   "Mass in Convective Envelope",
                   "Radius of Base of Convective Envelope",
                   "Central Density (log)", "Central Temperature (log)"]
    
    centre_axname = ["Model Number", "Age (Yrs)", "Gallinoes", "Neutrons",
                     r"$^{2}$H", r"$^{3}$He", r"$^{7}$Li", r"$^{7}$Be",
                     r"$^{11}$B", r"$^{13}$C", r"$^{14}$C", r"$^{15}$N",
                     r"$^{17}$O", r"$^{18}$O", r"$^{19}$F", r"$^{21}$Ne",
                     r"$^{22}$Ne", r"$^{22}$Na", r"$^{23}$Na", r"$^{24}$Mg",
                     r"$^{wtts/Lab4-2009.pdf25}$Mg", r"$^{26}$Mg", r"$^{26}$Al$^{m}$",
                     r"$^{26}$Al$^{g}$", r"$^{27}$Al", r"$^{28}$Si",
                     r"$^{29}$Si", r"$^{30}$Si", r"$^{31}$P", r"$^{32}$S",
                     r"$^{33}$S", r"$^{34}$S", r"$^{56}$Fe", r"$^{57}$Fe",
                     r"$^{58}$Fe", r"$^{59}$Fe", r"$^{60}$Fe", r"$^{59}$Co",
                     r"$^{58}$Ni", r"$^{59}$Ni", r"$^{60}$Ni", r"$^{61}$Ni",
                     r"$^{1}$H", r"$^{4}$He", r"$^{12}$C", r"$^{14}$N",
                     r"$^{16}$O", r"$^{20}$Ne"]
    #Entire loop depends on existence of files!
    #if (os.path.exists('plot') and os.path.exists('centre')):
    if True:
        #Load in both files
        plotdata = pd.read_table(sys.argv[1], header=None, sep='\s+').fillna(0).values
        #centredata = np.loadtxt('centre')/data/ajh291/ross_code/cont.36900b/plot.1
        
        '''
        Define all the plot vars
        '''
        
        #Panel 1
        #HRD
        x_1_1 = 4
        x_1_1 = x_1_1 - 1
        
        y_1_1 = 5
        y_1_1 = y_1_1 -1
        
        #Rhoc tc
        x_1_2 = 73
        x_1_2 = x_1_2 -1
        
        y_1_2 = 74
        y_1_2 = y_1_2 -1
        
        #Kipp diag
        x_1_3 = 2
        x_1_3 = x_1_3 -1
        
        #y_1_3
        
        #Central H frac and He Frac
        x_1_4 = 0
        x_1_4 = x_1_4 + 1
        
        y_1_4_a1 =1 +3 
        y_1_4_a2 = 1 +41
        
        y_1_4_b1 = 1+ 4
        y_1_4_b2 = 1+ 42
        
        '''
        #Panel 2
        x_2_1
        y_2_1
        
        x_2_2
        y_2_2
        
        x_2_3
        y_2_3
        
        x_2_4
        y_2_4
        '''
        #Set up HR diagram
        ax1_1.clear()
        ax1_1.plot(plotdata[:,x_1_1], plotdata[:,y_1_1])
        ax1_1.set_title('HR Diagram')
        ax1_1.set_xlim((np.max(plotdata[:,x_1_1]) + 0.01* np.max(plotdata[:,x_1_1]),
                        np.min(plotdata[:,x_1_1])- 0.01*np.min(plotdata[:,x_1_1])))
        ax1_1.set_xlabel(plot_axname[x_1_1])
        ax1_1.set_ylabel(plot_axname[y_1_1])
        
        #Set up Rhoc Tc
        ax1_2.clear()
        ax1_2.plot(plotdata[:,x_1_2], plotdata[:,y_1_2])
        ax1_2.set_xlabel(plot_axname[x_1_2])
        ax1_2.set_ylabel(plot_axname[y_1_2])
        
        #The very messy Kipp Diag
        ax1_3.clear()
        ax1_3.scatter(plotdata[:,1], plotdata[:,48], s = 0.03)
        ax1_3.scatter(plotdata[:,1], plotdata[:,49], s = 0.03)
        ax1_3.scatter(plotdata[:,1], plotdata[:,50], s = 0.03)
        ax1_3.scatter(plotdata[:,1], plotdata[:,51], s = 0.03)
        ax1_3.scatter(plotdata[:,1], plotdata[:,52], s = 0.03)
        ax1_3.scatter(plotdata[:,1], plotdata[:,53], s = 0.03)
        ax1_3.scatter(plotdata[:,1], plotdata[:,54], s = 0.03)
        ax1_3.scatter(plotdata[:,1], plotdata[:,55], s = 0.03)
        ax1_3.scatter(plotdata[:,1], plotdata[:,56], s = 0.03)
        ax1_3.scatter(plotdata[:,1], plotdata[:,57], s = 0.03)
        ax1_3.scatter(plotdata[:,1], plotdata[:,58], s = 0.03)
        ax1_3.scatter(plotdata[:,1], plotdata[:,59], s = 0.03)
        ax1_3.set_xlabel(plot_axname[1])
        ax1_3.set_ylabel(plot_axname[55])
        ax1_3.set_title("Kippenhahn Diagram (Burning Shells)")
        ax1_3.set_ylim(-0.1, np.max(plotdata[:,5]))
        
        #The h/he frac plot
        ax1_4.clear()
        
        #ax1_4.plot(plotdata[:,1], plotdata[:,8])
        
        
        
        ax1_4.plot(plotdata[:,1], plotdata[:,23], label = 'H', color = 'white')
        ax1_4.plot(plotdata[:,1], plotdata[:,24], label = 'He', color = 'red')
        ax1_4.plot(plotdata[:,1], plotdata[:,24], label = 'C', color = 'blue')
        ax1_4.legend()
        ax1_4.set_xlabel(plot_axname[1])
        ax1_4.set_ylabel('Mass Coordinate of Max Energy Generation')
        
        '''
        ax1_4.plot(centredata[:,x_1_4], np.log10(centredata[:,y_1_4_a1]+centredata[:,y_1_4_a2]), color='white', label=r'H Mass Frac')
        ax1_4.plot(centredata[:,x_1_4], np.log10(centredata[:,y_1_4_b1]+centredata[:,y_1_4_b2]), color='red', label=r'He Mass Frac')
        #ax1_4.set_ylim(0,1)
        ax1_4.set_xlabel(centre_axname[x_1_4])
        ax1_4.set_ylabel(r'Mass Fraction')
        ax1_4.set_title(r'Central Mass Fractions')
        ax1_4.legend()
        '''
        
        
        
        
        
        
        

'''
Set up the two intital plot windows, and the four subwindows each
'''



plt.style.use('dark_background')
fig1 = plt.figure(figsize = (20, 12), facecolor = 'black')

gs1 = gridspec.GridSpec(2,2)
ax1_1 = plt.subplot(gs1[0,0])
ax1_2 = plt.subplot(gs1[0,1])
ax1_3 = plt.subplot(gs1[1,0])
ax1_4 = plt.subplot(gs1[1,1])
'''
fig2 = plt.figure(figsize = (20, 12), facecolor = 'black')

gs2 = gridspec.GridSpec(2,2)
ax2_1 = plt.subplot(gs2[0,0])
ax2_2 = plt.subplot(gs2[0,1])
ax2_3 = plt.subplot(gs2[1,0])
ax2_4 = plt.subplot(gs2[1,1])
'''

ani = animation.FuncAnimation(fig1, animate, interval = 100)
plt.show()        
    
    

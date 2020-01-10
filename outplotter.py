#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 12:19:23 2019

@author: ajh291
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import os.path
import matplotlib.gridspec as gridspec
import pandas as pd
from tkinter.filedialog import askopenfilename, askdirectory
import numpy as n
import scipy.interpolate
import scipy.ndimage
from tqdm import tqdm

def congrid(a, newdims, method='linear', centre=False, minusone=False):
    import numpy as n
    import scipy.interpolate
    import scipy.ndimage
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).
    
    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''
    if not a.dtype in [n.float64, n.float32]:
        a = n.cast[float](a)

    m1 = n.cast[int](minusone)
    ofs = n.cast[int](centre) * 0.5
    old = n.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        print("[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions.")
        return None
    newdims = n.asarray( newdims, dtype=float )
    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = n.indices(newdims)[i]
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        cd = n.array( dimlist ).round().astype(int)
        newa = a[list( cd )]
        return newa

    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = n.arange( newdims[i] )
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        # specify old dims
        olddims = [n.arange(i, dtype = n.float) for i in list( a.shape )]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
        newa = mint( dimlist[-1] )

        trorder = [ndims - 1] + range( ndims - 1 )
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )

            mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
            newa = mint( dimlist[i] )

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )

        return newa
    elif method in ['spline']:
        oslices = [ slice(0,j) for j in old ]
        oldcoords = n.ogrid[oslices]
        nslices = [ slice(0,j) for j in list(newdims) ]
        newcoords = n.mgrid[nslices]

        newcoords_dims = range(n.rank(newcoords))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (n.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print("Congrid error: Unrecognized interpolation type.\n", \
              "Currently only \'neighbour\', \'nearest\',\'linear\',", \
              "and \'spline\' are supported.")
        return None

#Names
plot_names = ["Model Number", "Age (years)", "Radius (log)",
                   "Effective Temperature (log)", "Luminosity (log)", "Mass",
                   "Mass of H Core", "Mass of He Core", "H Luminosity (log)", 
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

out_names = ["K, Mesh Point", r"$\psi$ Degeneracy Parameter",
             r"P, Pressure (dyn/cm$^{3}$)", r"$\rho$, Density ($g/cm^{3}$)",
             r"T, Temperature (K)", r"$\kappa$, Opacity", r"$\nabla$",
             r"$\nabla _{ad}$", r"$\nabla_{rad} - \nabla_{ad}$",
             r"Mass, ($M/M_{\odot}$)", r"$^{1}$H", r"$^{4}$He", r"$^{12}$C",
             r"$^{14}$N", r"$^{16}$O", r"$^{20}$Ne", r"$^{24}$Mg",
             r"Radius ($R_{\odot}$)", r"Luminosity ($L_{\odot}$)",
             r"$E_{th}$, Thermodynamic Energy Generation ($erg/s/g$)",
             r"$E_{nuc}$, Nuclear Burning Energy Generation ($erg/s/g$)",
             r"$E_{\nu}$, Neutrino Energy Loss ($erg/s/g$)",
             r"$\delta$M, Mass Flux ($g/s$)",
             r"$K^{2}$, Squared Radius of Gyration",r"$n/(n+1)$",r"$U_{hom}$",
             r"$V_{hom}$", r"U, Internal Energy ($erg$)", r"S, Entropy",
             r"$L/L_{edd}$, Fraction of Eddington Luminosity", r"$\mu$ Mean Molecular Mass",
             r"$\mu_{ideal}$", r"$\sigma_{thermohaline}$", r"$E_{bind}$", r"$\beta$ ($P_{rad}/P_{tot}$)"]

#Read the out file
outfile = askopenfilename(initialdir='/data/ajh291/camstars_testbed/models/', title='Select Out File')
#outfile = '/data/ajh291/camstars_testbed/models/5_solar_mass/evo_output/final_few_models.out'
out_data = open(outfile, 'r').read()
out_data = out_data.split('\n')
out_data = [[float(x) for x in out_data[f].split()] for f in tqdm(range(len(out_data)))]

model_nums = []
mod_ages = []
mod_inter = []
i = 0
while True:
    try:
        modnum = out_data[i][0]
        modage = out_data[i][1]
        meshpnts = out_data[i][2]
        
        model_nums.append(modnum)
        mod_ages.append(modage)
        
        inter_details = []
        for j in range(int(meshpnts)):
            inter_details.append(out_data[j+1])
        
        mod_inter.append(inter_details)
        
        i = i + int(meshpnts) + 1
    except:
        break


mod_inter = np.swapaxes(mod_inter, 1, 2)





#Read the plot file
#plot_data = np.loadtxt('plot')

#Plotting Stuff
plt.style.use('dark_background')

#HRD
def HRD():
    fig, ax = plt.subplots()
    ax.plot(plot_data[:,3], plot_data[:,4])
    ax.set_title('HR Diagram')
    ax.set_xlabel(plot_names[3])
    ax.set_ylabel(plot_names[4])
    ax.set_xlim((np.max(plot_data[:,3]) + 0.01* np.max(plot_data[:,3]),
                        np.min(plot_data[:,3])-
                        0.01*np.min(plot_data[:,3])))

def KippDiag(var_track, lttc=False, yn=False):
    fig, ax = plt.subplots()
    age_of_mod = []
    mass_in_sh = []
    track_quant = []
    track_quant_yn = []
    for i in (range(len(model_nums))):
        for j in range(int(mod_inter[i][0][0])):
            age_of_mod.append(mod_ages[i])
            mass_in_sh.append(mod_inter[i][j][9])
            track_quant.append(mod_inter[i][j][var_track])
            if mod_inter[i][j][var_track] >= 0:
                track_quant_yn.append(1)
            else:
                track_quant_yn.append(0)
            
            
    if lttc:
        age_of_mod = np.log10(age_of_mod[-1] - np.asarray(age_of_mod))

    if yn:
        im = ax.scatter(age_of_mod, mass_in_sh, c=track_quant_yn)
    else:
        im = ax.scatter(age_of_mod, mass_in_sh, c=track_quant)
    if lttc:
        ax.set_xlim(np.nanmax(age_of_mod[np.isfinite(age_of_mod)]), np.nanmin(age_of_mod[np.isfinite(age_of_mod)]))
    ax.set_xlabel('Age (Yrs)')
    if lttc:
        ax.set_xlabel('log Time to Collapse (Yrs)')
    ax.set_ylabel(r'Mass Coordinate ($M/M_{\odot}$)')
    cbar = fig.colorbar(im, ax=ax)
    #cbar.set_label(r'$\nabla _{rad} - \nabla _{ad}$')
    


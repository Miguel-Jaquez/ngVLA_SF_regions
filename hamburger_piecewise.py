#!/usr/bin/env python
# coding: utf-8

#******************************
#Handling fits, pvs, and coords
#******************************
import sys
sys.path.append('/fs/posgrado30/other0/opt/star-forming-regions')
sys.path.append('/fs/posgrado30/other0/jesus/radmc3dPy/lib/python3.9/site-packages')

import sf3dmodels
import astropy.wcs as wcs
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
from spectral_cube import SpectralCube
from pvextractor import extract_pv_slice, Path, PathFromCenter
import radmc3dPy.image as image
from plot_helpers import plot_1d_props
#------------------
#Import the package
#------------------
from sf3dmodels import Model, Plot_model
from sf3dmodels import Resolution as Res
import sf3dmodels.utils.units as u
import sf3dmodels.utils.constants as ct
import sf3dmodels.rt as rt
#-----------------
#Extra libraries
#-----------------
from matplotlib import colors
import numpy as np
import os
import subprocess
import time
import astropy.units as asu
import astropy.constants as cts
import pandas as pd
sys.path.append("/home/jesus/Documents/paper2/ngVLA_SF_regions")
from utils_run_radmc import *

#------------------
#General Stellar Parameters
#------------------
MStar = 10.0           #[Msun]  stellar mas
MStar = MStar * u.MSun #past to SI units
LStar = u.LSun * ( MStar/u.MSun )**4   #this is need for ? for nothing 
RStar = u.RSun * ( MStar/u.MSun )**0.8
TStar = u.TSun * ( (LStar/u.LSun) / (RStar/u.RSun)**2 )**0.25
print ('STARs PARAMETERS   RStar:', RStar/u.RSun, ', LStar:', LStar/u.LSun, ', TStar:', TStar)
# Gravitational radius
## --- Gravitational radius
G = cts.G.si.value
k = cts.k_B.si
R = cts.R.si
#print(k)
T =9e3 * asu.K # 9.e3 * asu.K
mu = 1. #*asu.kg/asu.mol
mH = cts.m_p.si
gama = 5/3  # monoatomic gas
# for a polytropic gas, P \propto \rho^\gamma, cs = \sqrt(\gamma P/ \rho)
# for a isothermal gas P = \rho K_B T / (\mu m_H)
cs = np.sqrt(k * T/(mu*mH))
print (cs.to(asu.km/asu.s))
##### Roberto 
MStar_0 = MStar * asu.kg
#cs = 8650.704500032553 #m/s for a gas at T = 1e4
Rg = G * MStar /(cs.value**2) #m gravitational radius
print ("Gravitational radius [m]",Rg)
Rd = Rg #* u.au  # m ,Rg = gravitational radius #
"Gravitational Radius in AU = {:.2e}".format((Rg*asu.m).to(asu.au))

#-------------------------------
#Parameters for the disk
############# look for the denisty of reference rho(R=1 AU)
#-------------------------------

MRate = 5e-6 * u.MSun_yr #based on the ionizing photon rate                          ###########3.33e-5 * u.MSun_yr   #From Roberto Paper 2023
RhoE0 = Res.Rho0(MRate, Rd, MStar)    #density from Ulrich envelope
i = 60                #inclination
Rrot_pivot = Rd / 2 #150 #AU    #This value was the half of Rdics
pR_list_0 = -0.5
pR_list_1 = -0.5
# radial velocity parameters 
rrad_int = 150          # radius at what beggin the radial velocity
vr0 = 0                 # radial velocity
pr_list_0 = 0.0         # maybe a power law index

### from density_Env_Disc () #function
Arho = 1.0 #5.0 #5.29   # Roberto paper
rhoD0 = Arho * RhoE0      #Normalization factor based on the envelope
print ("density at R_g = {:.5e}".format(rhoD0))
#H = H0 * (RList / RStar)**p  #Scaleheight # equation Roberto paper , # = Z = 0
p_list_0 = -0.674
p_list = [p_list_0] 
#-------------------
#####      PHYSICAL PROPERTIES FOR THE MODEL
#-------------------
#DISC
#-------
Rmin = 1 * u.au #* H0_factor      #m,  Rmin dijo Roberto que lo voy a dejar en 1 au, 
### this is now our rho0 to create the disk density
rho_disk_Rmin = rhoD0 * (Rmin / Rd)**p_list_0 #* np.exp(-0.5 * zList**2 / H**2), 1.0)#1.0e9 / 2)
print ("density at r0(1 AU) = {:.5e}".format(rho_disk_Rmin))

Rdisc = 250 * u.au
R_list_rho = [Rmin, Rdisc]           #Interval limits for piecewise density profile
p_list = [p_list_0]                  #Powerlaws within intervals
rho0 = rho_disk_Rmin #* H0_factor**p_list_0   # in m^-3  # density pivot at Rmin
RH_list = R_list_rho
q_list = [1. for p in p_list]  # H = H0 * (R / R0)**(1 + 0.5*(1-q))

#from Izquierdo paper
#H = H0(R/Rstar)---> for q = 1  ---> H0 = H Rstar/Rdisc  where H = Rgrav 
print ("Rg =", Rg)
print ("Rmin" ,Rmin)
H0_0 = Rg * (Rmin / Rg)**1 #in m
print ("H0_0 = {:.2e}".format(H0_0))
H0_factor = 0.2                    # I'm going to leave it as is
H0 = H0_0 #* H0_factor
#### print properties
#print (RStar)
#print (R_list_rho)
#print (p_list)
print (rho0)

#---------------
#Density, Temperature and Velocity GRID Definition
#---------------
#box limits [-250,250]
sizex = sizey = sizez = 250 * u.au #halfsize
Nx = Ny = Nz = 201                   #Number of divisions for each axis
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], rt_code = 'radmc3d')#, include_zero=False)
NPoints = GRID.NPoints               # Final number of nodes in the grid

density = Model.density_Hamburgers_piecewise(RStar, H0, R_list_rho, p_list, rho0, GRID,
                                              RH_list = RH_list, q_list = q_list, rho_thres = 1e7,
                                             rho_min = 0.0, Rt = False)

#---------------------
# MODEL TEMPERATURE
#---------------------
# Ionized gas temperature
t_e = 9.0e3 #K
temperature = Model.temperature_Constant(density, GRID, discTemp=t_e, backTemp=2.725480)

#---------------------
# Calculation of recombination rate 
#---------------------
print('############################################')
print('The GRID cell lengths are: ({0},{1},{2}) meters'.format(GRID.step[0],GRID.step[1],GRID.step[2]))
dv = GRID.step[0]*GRID.step[1]*GRID.step[2] #m^3
mp = 1.6726e-27 #kg
Msun = 1.989e30 #kg
mass = mp*np.sum(density.total)*dv # [kg m ^-3]*[m^3] = kg
print('Total mass in model is {0} kg'.format(mass))
print('Total mass in model is {0} Msun'.format(mass/Msun))
alpha_t = 2.6e-19 #case B recombination coeff in m^3 s^-1
phi_ion = alpha_t*dv*np.sum(density.total**2)
print('The ionizing photon rate phi is {0:.3e} s^-1'.format(phi_ion))
print('The logarithmic ionizing photon rate phi is {0:.3f} s^-1'.format(np.log10(phi_ion)))
print('############################################')

#--------
#VELOCITY
#--------
# rotation
Rrot_pivot = Rrot_pivot * u.au
R_list_rot = [Rmin, Rrot_pivot, Rdisc] #Interval limits for piecewise rotation profile
pR_list = [pR_list_0, pR_list_1] #powerlaw exponents
vR0 = (ct.G * MStar / R_list_rot[0])**0.5	# Circular velocity at R_list_rot[0]
norm_vR = [0,0,vR0] #Components of rotation velocity. Only index [2] is taken.

rrad_int = rrad_int * u.au
r_list = [rrad_int, Rdisc]
pr_list = [pr_list_0] #powerlaw exponents
vr0 = vr0             #m/s
norm_vr = [vr0,0,0]   #Components of infall. Only [0] is taken. sign determines direction, '-' for infall

vel = Model.velocity_piecewise(density, GRID,
                               R_list=R_list_rot, pR_list=pR_list, v0R=norm_vR, #polar piecewise velocity
                               r_list=r_list, pr_list=pr_list, v0r=norm_vr,     #radial piecewise velocity
                               )
Model.PrintProperties(density, temperature, GRID, species='dens_ion')
Model.PrintProperties(density, temperature, GRID, species='dens_e')
#***************
#prop DICTIONARY
#***************
prop = {'vel_x' : vel.x, 'vel_y' : vel.y, 'vel_z' : vel.z,
    'dens_e' : density.total, 'dens_ion' : density.total,
    'temp_gas' : temperature.total}
#-----------------------------------------------
#3D Points Distribution (weighting with density)
#-----------------------------------------------
tag = 'Main'
dens_plot = density.total / 1e6
weight = 10*RhoE0
r = GRID.rRTP[0] / u.au #GRID.rRTP hosts [r, R, Theta, Phi] --> Polar GRID
Plot_model.scatter3D(GRID, density.total, weight, NRand = 4000, 
                     colordim = r, axisunit = u.au, cmap = 'jet', colorscale = 'log', 
                     colorlabel = r'${\rm log}_{10}(r [au])$',output = '3Dpoints%s.png'%tag, show = False)
#---------------------
#2D PLOTTING (Density)
#---------------------

vmin, vmax = np.array([2e13, 1e19]) / 1e6
norm = colors.LogNorm(vmin=vmin, vmax=vmax)

Plot_model.plane2D(GRID, dens_plot, axisunit = u.au,
                   cmap = 'jet', plane = {'z': 0*u.au},
                   norm = "log", colorlabel = r'$[\rm cm^{-3}]$',
                   output = 'DensMidplane_%s.png'%tag, show = True)

Plot_model.plane2D(GRID, dens_plot, axisunit = u.au,
                   cmap = 'jet', plane = {'y': 0*u.au},
                   norm = "log", colorlabel = r'$[\rm cm^{-3}]$',
                   output = 'DensMidplane_%s.png'%tag, show = False)

plot_1d_props(GRID, prop, density, tag='single_model')


#############----------------------------- Radiative transfer

################
#read the lines properties to run RADMC
#lines_df = pd.read_csv("/fs/posgrado30/other0/jesus/paper2/Hydrogen_recom_lines_in_use.csv")
lines_df =  pd.read_csv("/home/jesus/Documents/paper2/Hydrogen_recom_lines_in_use.csv")

#run radmc
i = 0   #this is the index of the line that you want to simulate.
run_radmc(lines_df.iloc[i],GRID,prop,beam=(2e-3,2e-3))



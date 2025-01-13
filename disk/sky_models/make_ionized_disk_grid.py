#!/usr/bin/env python
# coding: utf-8

import sys
import os
myhost = os.uname()[1]
if myhost == 'posgrado30.crya.privado':
    ### import sf3dmodels path
    sys.path.append('/fs/posgrado30/other0/opt/star-forming-regions')
    ### import radmc3dPy
    sys.path.append('/fs/posgrado30/other0/opt/star-forming-regions')
    sys.path.append('/fs/posgrado30/other0/jesus/radmc3dPy/lib/python3.9/site-packages')
    ### import utils folder to run radmc
    sys.path.append("/fs/posgrado30/other0/jesus/respaldo_paper2/utils")
    lines_properties_file = "/fs/posgrado30/other0/jesus/respaldo_paper2/ngVLA_SF_regions/Hydrogen_recom_lines_in_use.csv"
else:
    ### import sf3dmodels path
    sys.path.append('/fs/posgrado30/other0/opt/star-forming-regions')
    ### import radmc3dPy
    sys.path.append('/fs/posgrado30/other0/opt/star-forming-regions')
    sys.path.append('/fs/posgrado30/other0/jesus/radmc3dPy/lib/python3.9/site-packages')
    ### import utils folder to run radmc
    sys.path.append("/fs/posgrado30/other0/jesus/respaldo_paper2/utils")
    lines_properties_file = "/fs/posgrado30/other0/jesus/respaldo_paper2/ngVLA_SF_regions/Hydrogen_recom_lines_in_use.csv"
from utils_run_radmc import * # this run radmc
##------
# import modules
#----------------
import astropy.wcs as wcs
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
#from spectral_cube import SpectralCube
#from pvextractor import extract_pv_slice, Path, PathFromCenter
import radmc3dPy.image as image
from plot_helpers import plot_1d_props, plot_1d_props_jaquez_dens_vs_radius
import matplotlib.pyplot as plt
#------------------
#Import sf3d models package
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
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import os
import subprocess
import time
import astropy.units as asu
import astropy.constants as cts
import pandas as pd
radmc_path = "/fs/posgrado30/other0/jesus/radmc-3d/version_0.40/examples/run_recomblines_userdef"

######## radio velocity to frequency definition
def freq2vel(freq, freq_0):
    return (freq_0 - freq) * (cts.c.cgs.to(asu.km/asu.s)) / freq_0
def vel2freq(vel,freq_0):
    return freq_0 * (1.-(vel/(cts.c.cgs.to(asu.km/asu.s))))

# compare with the old model
#------------------
#General Stellar Parameters
#------------------
MStar = 20.0           #[Msun]  stellar mas
MStar = MStar * u.MSun #past to SI units
LStar = u.LSun * ( MStar/u.MSun )**4   #this is need for ? for nothing 
RStar = u.RSun * ( MStar/u.MSun )**0.8
TStar = u.TSun * ( (LStar/u.LSun) / (RStar/u.RSun)**2 )**0.25
print ('STARs PARAMETERS   RStar:', RStar/u.RSun, ', LStar:', LStar/u.LSun, ', TStar:', TStar)

## ------------ Gravitational radius
G = cts.G.si.value
k = cts.k_B.si
R = cts.R.si
#print(k)
T =9.e3 * asu.K 
mu = 1. #*asu.kg/asu.mol
mH = cts.m_p.si
gama = 5./3  # monoatomic gas
# for a polytropic gas, P \propto \rho^\gamma, cs = \sqrt(\gamma P/ \rho)
# for a isothermal gas P = \rho K_B T / (\mu m_H)
#cs = np.sqrt(k * T/(mu*mH))
cs = np.sqrt(k * T/(mu*mH))
print ('velocidad del sonido en el gas ionizado = ', cs.to(asu.km/asu.s)) ##### Roberto  paper #cs = 8650.704500032553 #m/s for a gas at T = 1e4
MStar_0 = MStar * asu.kg
Rg = 0.5 * G * MStar /(cs.value**2) #m Sonic Point (Vel scape = Vel sound)
print ("Gravitational radius [au]", Rg/u.au)
Rd = Rg #* u.au  # m ,Rg = gravitational radius # 

#-------------------------------
#Parameters for the disk density
############# look for the denisty of reference rho(R=1 AU)
#-------------------------------
MRate = 3.2e-5 * u.MSun_yr #based on the ionizing photon rate  ########### 3.33e-5 * u.MSun_yr   #From Roberto Paper 2023
RhoE0 = Res.Rho0(MRate, Rd, MStar)    #density from Ulrich Envelope at Rd
print ('The density envelope at Rd is {:.2e}'.format(RhoE0))
i = 30                #inclination

### from density_Env_Disc () #function
Arho = 1.0 #5.0 #5.29   # Roberto paper
rhoD0 = Arho * RhoE0      #Normalization factor based on the envelope
print ("density at R_g in the mid plane = {:.5e}".format(rhoD0))
#H = H0 * (RList / RStar)**p  #Scaleheight # equation Roberto paper , # = Z = 0
p_list_0 = -0.674  
p_list = [p_list_0] 
Rmin = 1 * u.au #* H0_factor      #m,
### this is now our rho0 to create the disk density
rho_disk_Rmin = rhoD0 * (Rmin / Rd)**p_list_0 # * np.exp(-0.5 * zList**2 / H**2), 1.0)#1.0e9 / 2)
print ("density at R(R=1 AU) = {:.5e}".format(rho_disk_Rmin))

#-------------------
#####      PHYSICAL PROPERTIES FOR THE MODEL
Rdisc = Rd  # Disk size at mid plane
R_list_rho = [Rmin, Rdisc]           #Interval limits for piecewise density profile
## why we leave the inner radius at 1 au
p_list = [p_list_0]                  #Powerlaws within intervals
rho0 = rho_disk_Rmin #* H0_factor**p_list_0   # in m^-3  # density pivot at Rmin
RH_list = R_list_rho
q_list = [1. for p in p_list]  # H = H0 * (R / R0)**(1 + 0.5*(1-q))

#from Izquierdo paper
#H = H0(R/Rstar)---> for q = 1  ---> H0 = H Rstar/Rdisc  where H = Rgrav 
print ("Rg = {:.2e}".format((Rg*asu.m).to(asu.au)))
print ("Rmin = {:.2e}".format((Rmin*asu.m).to(asu.au)))
H0_0_old = Rg * (Rmin / Rdisc) #**1 #in m
print ("H0_0 = {:.2e}".format(H0_0_old))
H0_factor = 0.2                    # I'm going to leave it as is
H0 = H0_0_old #* H0_factor
print ((H0*asu.m*Rmin/Rd).to(asu.au))
print ((np.sqrt(2)*Rg*asu.m).to(asu.au))

###### radial velocity parameters (not used for the model)
Rrot_pivot = 0.5 * Rd    # this is for calculate the velocities
pR_list_0 = -0.5
pR_list_1 = -0.5
rrad_int = 150          # radius at what beggin the radial velocity
vr0 = 0                 # radial velocity
pr_list_0 = 0.0         # maybe a power law index

#---------------
#Density, Temperature and Velocity GRID Definition
#---------------
#box limits [-250,250]
sizex = sizey = sizez = 200 * u.au #halfsize
Nx = Ny = Nz = 375     #cdelt = 0.2 au      #Number of divisions for each axis
GRID_old = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], rt_code = 'radmc3d')#, include_zero=False)
NPoints = GRID_old.NPoints               # Final number of nodes in the grid

density_old = Model.density_Hamburgers_piecewise_bound(RStar, H0, R_list_rho, p_list, rho0, GRID_old,
                                              RH_list = RH_list, q_list = q_list, rho_thres = 1e7,
                                             rho_min = 0.0, Rt = False)

#---------------------
# MODEL TEMPERATURE
#---------------------
# Ionized gas temperature
t_e = 9.0e3 #K
temperature_old = Model.temperature_Constant(density_old, GRID_old, discTemp=t_e, backTemp=2.725480)

#---------------------
# Calculation of recombination rate 
#---------------------
print('############################################')
print('The GRID cell lengths are: ({0},{1},{2}) meters'.format(GRID_old.step[0],GRID_old.step[1],GRID_old.step[2]))
dv_old = GRID_old.step[0]*GRID_old.step[1]*GRID_old.step[2] #m^3
mp = 1.6726e-27 #kg
Msun = 1.989e30 #kg
mass = mp*np.sum(density_old.total)*dv_old # [kg m ^-3]*[m^3] = kg
print('Total mass in model is {0} kg'.format(mass))
print('Total mass in model is {0} Msun'.format(mass/Msun))
alpha_t = 2.6e-19 #case B recombination coeff in m^3 s^-1
phi_ion = alpha_t*dv_old*np.sum(density_old.total**2)
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

vel_old = Model.velocity_piecewise(density_old, GRID_old,
                               R_list=R_list_rot, pR_list=pR_list, v0R=norm_vR, #polar piecewise velocity
                               r_list=r_list, pr_list=pr_list, v0r=norm_vr,     #radial piecewise velocity
                               )

Model.PrintProperties(density_old, temperature_old, GRID_old, species='dens_ion')
Model.PrintProperties(density_old, temperature_old, GRID_old, species='dens_e')
#***************
#prop DICTIONARY
#***************
prop_old = {'vel_x' : vel_old.x, 'vel_y' : vel_old.y, 'vel_z' : vel_old.z,
    'dens_e' : density_old.total, 'dens_ion' : density_old.total,
    'temp_gas' : temperature_old.total}


#############
# Radiative transfer
#############
#read the lines properties to run RADMC
lines_df = pd.read_csv(lines_properties_file)

# RADMC PARAMETERS}
BEAM = beam = (3.7e-3,2.6e-3)   #arcsec
TIME_INT= time_int = 3600*5.0 #1hr #second 
DISTANCE= distance = 4000.0    # pc
CHANNEL_WIDTH= channel_width_vel = 5   #km/s
NUMBER_CHANNELS= number_channels = 60 # INT
INCLINATION= inclination = 30     #DEG
SIZE_REGION=size_region = 2400    #AU  that the image will be cover
NPIX=npix = 750*4 # to have cdelt = 0.19 mas, if the beam was 1mas the 4 pixels in the beam.

line = 1 #H38a
line_prop = lines_df.iloc[line]
n_sup = line_prop["nsup"]
n_inf = line_prop["ninf"]
line_freq = line_prop["Frequency_radmc"]
channel_width_velocity = channel_width_vel * asu.km/asu.s
channel_width_freq_0 = vel2freq(channel_width_velocity,line_freq*asu.Hz)
channel_width_freq = line_freq*asu.Hz - channel_width_freq_0 
print ("channel width {} km/s or {} MHz".format(channel_width_velocity, channel_width_freq.to(asu.MHz)))
freq_ini = line_freq*asu.Hz - ((channel_width_freq)*(number_channels/2))
freq_fin = line_freq*asu.Hz + ((channel_width_freq)*(number_channels/2))
print ("{:.7e}, {:.7e}".format(freq_ini,freq_fin))
wave_ini = freq_fin.to(asu.micron, equivalencies=asu.spectral())
wave_fin = freq_ini.to(asu.micron, equivalencies=asu.spectral())
wave_range=[abs(wave_ini.value),abs(wave_fin.value)]
print ("the wave range is {} to {} microns".format(wave_ini,wave_fin))
#--------------------------------
i = inclination #30

#********************
#WRITING FOR RADMC-3D
#********************
print ("################ Writing for RADMC-3D ####################")
A = rt.Radmc3dDefaults(GRID_old)
A.recomblines(prop_old, [n_sup,n_inf], kwargs_control = {'userdef_nonlte': 1}) #Writing prop with radmc3d format, non-LTE



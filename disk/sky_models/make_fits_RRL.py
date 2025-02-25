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
    radmc_path = "/fs/posgrado30/other0/jesus/radmc-3d/version_0.40/examples/run_recomblines_userdef"
    sys.path.append("/fs/posgrado30/other0/jesus/respaldo_paper2/utils")
    lines_properties_file = "/fs/posgrado30/other0/jesus/respaldo_paper2/ngVLA_SF_regions/Hydrogen_recom_lines_in_use.csv"
    # import ngVLA calculator
    sys.path.append("../../ngvla_sensitivity_calculator_master")
else:
    ### import sf3dmodels path
    sys.path.append('/fs/posgrado30/other0/opt/star-forming-regions')
    ### import radmc3dPy
    sys.path.append('/fs/posgrado30/other0/opt/star-forming-regions')
    sys.path.append('/fs/posgrado30/other0/jesus/radmc3dPy/lib/python3.9/site-packages')
    ### import utils folder to run radmc
    sys.path.append("/fs/posgrado30/other0/jesus/respaldo_paper2/utils")
    lines_properties_file = "/fs/posgrado30/other0/jesus/respaldo_paper2/ngVLA_SF_regions/Hydrogen_recom_lines_in_use.csv"

### import ngVLA noise calculator
from ngVLA_sensitivity_calculator import *

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



#read the lines properties to run RADMC
lines_df = pd.read_csv(lines_properties_file)

# RADMC PARAMETERS}
BEAM = beam = (4.8e-3,2.6e-3)   #arcsec
TIME_INT= time_int = 3600*5.0 #1hr #second 
DISTANCE= distance = 4000    # pc
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

#make fits
im = image.readImage('model_disk_line_incl30.out')    
model = 'disk_hipotenuse_corrected'
output = 'img_rl_{}_H{}_{}pc.fits'.format('disk','H38',distance)
im.writeFits(fname=output, dpc=distance, coord='18h00m00.0s +25d00m00.00s') #writting 3d fits cube
#fix header of fits file written by radmc3dpy
hdul = fits.open(output, mode='update')
hdul[0].header['CTYPE3'] = 'FREQ'
hdul[0].header['CUNIT3'] = 'Hz'
hdul[0].header['RESTFRQ'] = line_freq       #for H38                 #original #2.319935e11 # 2.3190093e+11
hdul[0].header['SPECSYS'] = 'LSRK'
hdul[0].header['ALTRVAL'] = 0.0
hdul[0].header['ALTRPIX'] = 25.0
hdul[0].header['VELREF'] =  257             # Local Standar of Rest: mean motion of material in the Milky Way in the neighborhood of the Sun
hdul.flush()   #### save the file resulting from the radiative transfer.  
############ separate the line cube from the continuum cube
print ("separating the continuum data cube and the line data cube")
#print('Separatig Starting continuum_subtraction (convl_csub.py)')
hdul = fits.open(output)
############################### fin add noise ####################
data = hdul[0].data  #units of Jy/pix  #output from radmc
data_line = data.copy()
data_cont = data.copy()
min_values_per_pixel = np.min(data,axis=0) #Ya prove que esto si funcione bien
for slide in range(data.shape[0]):
    data_line[slide,y,x] = data[slide,y,x] - min_values_per_pixel
    data_cont[slide,y,x] = min_values_per_pixel
csub_image = 'img_rl_{}_line_H{}_{}pc.fits'.format('disk','H38',distance)
cont_image = 'img_rl_{}_cont_H{}_{}pc.fits'.format('disk','H38',distance)
fits.writeto(csub_image,data_line,hdul[0].header,overwrite=True)
fits.writeto(cont_image,data_cont,hdul[0].header,overwrite=True)	#QZ
############--------------------------##############
#Noise calculation
frequency = line_freq/1e9
#continuum for now the noise is in the entire bandwidth 20GHz
#a_cont,b_cont,c_cont = sigma_ps_fn("main",  freq=frequency,type_cal = 'continuum', theta=beam[0],t_int=time_int, verbose=True )
#noise_input_continuum = a_cont[1]*1e-6   #Jy/beam  #original in uJy/beam
#print ("the continuum noise = {}".format(noise_input_continuum))
#line
a_line,b_line,c_line = sigma_ps_fn("main",  freq=frequency, type_cal = 'line', theta=beam[0],t_int=time_int, delta_v=channel_width_vel*1e3, verbose=True )
noise_input_line = 0.46e-6 #a_line[1]*1e-6     #Jy/beam  #original in uJy/beam
print ("the line noise = {}".format(noise_input_line))
#a[0] = ""#frequency
#a[1] is the simga_rms
#a[2] is the T_rms
#### end: Noise calculation
#--------------------------------------------------------------------------------
# line convolution and add noise
#------------------------------------------------------------------------------
print('\n')  
print('Starting line cube convolution and adding noise (convl_csub.py)')
convl_csub_image = 'img_rl_{}_line_convl_H{}_{}pc.fits'.format('disk','H38',distance)
hdul = fits.open(csub_image)
new_beam = beam[0] #arcsec = 2 mas  
npixsigma = (1/2.3548)*(new_beam/(3600*hdul[0].header['CDELT2']))
print ("are {} pixels per beam".format(npixsigma))
gauss_kernel = Gaussian2DKernel(npixsigma)
pixperbeam = 1.442*(np.pi/4)*(new_beam**2)/((3600*hdul[0].header['CDELT2'])**2)
convl_data = hdul[0].data.copy()
convl_data_noise = hdul[0].data.copy()

# PSF convolved
for element in range(hdul[0].data.shape[0]):
        fft_mp = lambda a: scipy.fft.fftn(a, workers=-1)  # use all available cores
        ifft_mp = lambda a: scipy.fft.ifftn(a, workers=-1)
        convl_data[element,:,:] = pixperbeam*convolve_fft(data_line[element,:,:],gauss_kernel, fftn=fft_mp, ifftn=ifft_mp)
        #add the noise after convolutio
#adding noise
for element in range(hdul[0].data.shape[0]):
    convl_data_noise[element,:,:] = convl_data[element,:,:] + np.random.normal(0,scale=noise_input_line, size=(hdul[0].data.shape[1],hdul[0].data.shape[2])) #Jy/beam

    #for x in range(hdul[0].data.shape[2]):
     #   for y in range(hdul[0].data.shape[1]):
      #      convl_data_noise[element,y,x] = convl_data[element,y,x] + np.random.normal(0,scale=noise_input_line) #Jy/beam

hdul[0].header['BMIN']= new_beam/3600
hdul[0].header['BMAJ']= new_beam/3600
hdul[0].header['BPA']= 0
hdul[0].header['BUNIT']= 'JY/BEAM'
fits.writeto(convl_csub_image,convl_data_noise,hdul[0].header,overwrite=True)
hdul.close()
#------------------------------------------------------

print('\n')
print('Finished running convolution and continuum subtraction')
print ('Ellapsed time for run RADMC: %.3fs' % (time.time() - t0))
print ('-------------------------------------------------\n-------------------------------------------------\n')
###########################################################################################################

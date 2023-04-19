#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
from astropy.io import fits
import matplotlib.pyplot as plt
import aplpy
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
# run radmc


# In[10]:


#------ Parameters
MStar = 10. * u.MSun
MRate = 3.33e-5 * u.MSun_yr   #Mass accretion rate From Roberto Paper 2022  4e-4 * u.MSun_yr #Mass accretion rate
RStar = 26 * u.RSun * ( MStar/u.MSun )**0.27 * ( MRate / (1e-3*u.MSun_yr) )**0.41
LStar = 3.2e4 * u.LSun
TStar = u.TSun * ( (LStar/u.LSun) / (RStar/u.RSun)**2 )**0.25
shFactor = 0.01 #Izquierdo 2018
Rd = 250 * u.au
RhoE0 =  Res.Rho0(MRate, Rd, MStar) #Base density for Ulrich's model
Arho =  24.1 #Disc-envelope density factor
p = 2.25  #from single disk de Roberto 2.25 is the original value 

#---------------
#GRID Definition
#---------------
# Cubic grid, each edge ranges [-sizex, sizex] au.
sizex = sizey = sizez = 300 * u.au #halfsize
Nx = Ny = Nz = 241 #Number of divisions for each axis
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], rt_code = 'radmc3d')
NPoints = GRID.NPoints #Number of nodes in the grid

density = Model.density_Hamburgers(RStar, shFactor, Rd, RhoE0, Arho, GRID, 
                       p = 2.25, q = 0.5, rho_thres = 10.0, rho_min = 1.0, 
                       Rt = False, discFlag=True, rdisc_max = False)

#-----------------------------------------------
#3D Points Distribution (weighting with density)
#-----------------------------------------------
tag = 'Main'
dens_plot = density.total / 1e6

weight = 10*RhoE0
r = GRID.rRTP[0] / u.au #GRID.rRTP hosts [r, R, Theta, Phi] --> Polar GRID
Plot_model.scatter3D(GRID, density.total, weight,
                     NRand = 4000, colordim = r, axisunit = u.au,
                     cmap = 'jet', colorscale = 'log',
                     colorlabel = r'${\rm log}_{10}(r [au])$',
                     output = '3Dpoints%s.png'%tag, show = True)
#---------------------
#2D PLOTTING (Density)
#---------------------

vmin, vmax = np.array([2e13, 1e19]) / 1e6
norm = colors.LogNorm(vmin=vmin, vmax=vmax)

Plot_model.plane2D(GRID, dens_plot, axisunit = u.au,
                   cmap = 'jet', plane = {'z': 0*u.au},
                   norm = norm, colorlabel = r'$[\rm cm^{-3}]$',
                   output = 'DensMidplane_%s.png'%tag, show = False)


Plot_model.plane2D(GRID, dens_plot, axisunit = u.au,
                   cmap = 'jet', plane = {'y': 0*u.au},
                   norm = norm, colorlabel = r'$[\rm cm^{-3}]$',
                   output = 'DensVertical_%s.png'%tag, show = False)

#---------------------
# MODEL TEMPERATURE
#---------------------
# Ionized gas temperature
t_e = 9.0e3 #K
temperature = Model.temperature_Constant(density, GRID, discTemp=t_e, backTemp=2.725480)

#--------
#VELOCITY
#--------
#Keplerian
vel = Model.velocity(RStar, MStar, Rd, density, GRID)
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

Model.PrintProperties(density, temperature, GRID, species='dens_ion')
Model.PrintProperties(density, temperature, GRID, species='dens_e')
#***************
#prop DICTIONARY
#***************
prop = {'vel_x' : vel.x, 'vel_y' : vel.y, 'vel_z' : vel.z,
    'dens_e' : density.total, 'dens_ion' : density.total,
    'temp_gas' : temperature.total}


plot_1d_props(GRID, prop, density, tag='single_model')

def freq2vel(freq, freq_0):
    return (freq_0 - freq) * (cts.c.cgs.to(asu.km/asu.s)) / freq_0
def vel2freq(vel,freq_0):
    return freq_0 * (1.-(vel/(cts.c.cgs.to(asu.km/asu.s))))

def run_radmc(line_prop,GRID,prop,beam=(2e-3,2e-3),noise_input = 5.3481272e-5,inclination=53,channel_width_0=4):
    '''
    noise in Jy/beam
    '''
    t0 = time.time()
    radmc_path = "/home/jesus/radmc-3d/version_0.40/examples/run_recomblines_userdef"
    #radmc_path = "/fs/posgrado30/other0/jesus/radmc-3d/version_0.40/examples/run_recomblines_userdef"
    #-------------------------------
    #channel_width_0 = 4 #km/s
    #ancho de banda
    n_sup = line_prop["nsup"]
    n_inf = line_prop["ninf"]
    line_freq = line_prop["Frequency_radmc"]
    #print ("line frequency",line_freq)
    #print ("control: velocity line", freq2vel(line_freq,line_freq))
    #print (n_sup,n_inf)
    channel_width = channel_width_0 * asu.km/asu.s
    print ("ancho del canal en km/s", channel_width)
    a = vel2freq(channel_width,line_freq*asu.Hz)
    #print ("ancho del canal en HZ = {:.2e} ".format(a))
    #cdelt
    print ("ancho del canal en HZ = {:.2e} ".format(a - line_freq*asu.Hz))
    #40 divitions
    freq_ini = (a - line_freq*asu.Hz)*20 - line_freq*asu.Hz
    freq_fin = (a - line_freq*asu.Hz)*20 + line_freq*asu.Hz
    #print ("{:.2e}, {:.2e}".format(freq_ini,freq_fin))
    line_wave = (line_prop["Frequency_radmc"]*asu.Hz).to(asu.micron, equivalencies=asu.spectral())
    #print (line_wave)
    wave_ini = freq_ini.to(asu.micron, equivalencies=asu.spectral())
    wave_fin = freq_fin.to(asu.micron, equivalencies=asu.spectral())
    wave_range=[abs(wave_ini.value),abs(wave_fin.value)]
    #print (wave_ini,wave_fin)
    #--------------------------------
    i = inclination #53
    #print (wave_range)  #control
    #********************
    #WRITING FOR RADMC-3D
    #********************
    print ("################ Writing for RADMC-3D ####################")
    A = rt.Radmc3dDefaults(GRID)
    A.recomblines(prop, [n_sup,n_inf], kwargs_control = {'userdef_nonlte': 1}) #Writing prop with radmc3d format, non-LTE

    ########################################## RADIATIVE TRANSFER WITH RADMC-3D ####################################
    radmc_rl = os.path.join(radmc_path, 'radmc3d')
    ##### JAQUEZ H38_\alpha -> 2600.6855 microns
    subprocess.run(radmc_rl+' image lambdarange {} {} nlam 40 incl {} npix 221 sizeau 650'.format(wave_range[0],wave_range[1],str(i)), shell=True, executable='/bin/bash') 
     
    ################################### CONVOLUTION AND CONTINUUM SUBTRACTION ######################################
    # Use radmc3dpy to write fits file
    output = 'img_rl_disk_H{}.fits'.format(n_inf)
    #subprocess.call(['rm',output], shell=True)
    os.system('rm -f '+output)
    dist = 5000.   #49000 = nube de magallanes
    im = image.readImage()    
    #data = im.image #Accessing image data
    #plt.imshow(data[:,:,29], origin='lower', cmap='cool') #plotting a single channel
    #plt.show()
    im.writeFits(fname=output, dpc=dist, coord='18h10m28.652s -19d55m49.66s') #writting 3d fits cube
        
    #fix header of fits file written by radmc3dpy
    hdul = fits.open(output, mode='update')
    hdul[0].header['CTYPE3'] = 'FREQ'
    hdul[0].header['CUNIT3'] = 'Hz'
    hdul[0].header['RESTFRQ'] = line_freq       #for H38                 #original #2.319935e11 # 2.3190093e+11
    hdul[0].header['SPECSYS'] = 'LSRK'
    hdul[0].header['ALTRVAL'] = 0.0
    hdul[0].header['ALTRPIX'] = 25.0
    hdul[0].header['VELREF'] = 257
    hdul.flush()   #### save the file resulting from the radiative transfer.     
    print('\n')
    ############ No convolutions
    print('Starting continuum_subtraction (convl_csub.py)')
    hdul = fits.open(output)
    ################################ add noise #######################
    #beam = beam[0] #arcsec 
    cdelt = hdul[0].header["CDELT1"]*3600   #this is in arcsec
    ## change the noise from Jy/beam to Jy/Pix to adding to the data
    pixperbeam = 1.442*(np.pi/4)*(beam[0]*beam[1])/((3600*hdul[0].header['CDELT2'])**2)
    noiseJyPix = (noise_input)*(1/pixperbeam)   #(1/(beam[0]*beam[1]))*((cdelt)**2)
    #check how past from beam to pixels
    ############################### fin add noise ####################
    data = hdul[0].data  #units of Jy/pix  #output from radmc
    data_csub = data.copy()
    data_cont = data.copy()
    for slide in range(data.shape[0]):
        for x in range(data.shape[2]):
            for y in range(data.shape[1]):
                data_csub[slide,y,x] = data[slide,y,x] - data[:,y,x].min()
                data_cont[slide,y,x] = data[:,y,x].min()
    csub_image = 'img_rl_disk_csub_H{}.fits'.format(n_inf)
    cont_image = 'img_rl_disk_cont_H{}.fits'.format(n_inf)
    fits.writeto(csub_image,data_csub,hdul[0].header,overwrite=True)
    fits.writeto(cont_image,data_cont,hdul[0].header,overwrite=True)	#QZ
    

    ########################    convolution
    print('\n')  
    print('Starting full image convolution (convl_csub.py)')
    hdul = fits.open(output)
    new_beam = beam[0] #arcsec = 2 mas  
    npixsigma = (1/2.3548)*(new_beam/(3600*hdul[0].header['CDELT2']))
    print ("are {} pixels per beam".format(npixsigma))
    gauss_kernel = Gaussian2DKernel(npixsigma)
    pixperbeam = 1.442*(np.pi/4)*(new_beam**2)/((3600*hdul[0].header['CDELT2'])**2)
    convl_data = hdul[0].data.copy()
    convl_data_noise = hdul[0].data.copy()
    for element in range(hdul[0].data.shape[0]):
        convl_data[element,:,:] = pixperbeam*convolve(hdul[0].data[element,:,:],gauss_kernel)
        #add the noise after convolution
    #adding noise
    for element in range(hdul[0].data.shape[0]):
        for x in range(hdul[0].data.shape[2]):
            for y in range(hdul[0].data.shape[1]):
                convl_data_noise[element,y,x] = convl_data[element,y,x] + np.random.normal(0,scale=noise_input) #Jy/beam
    hdul[0].header['BMIN']= new_beam/3600
    hdul[0].header['BMAJ']= new_beam/3600
    hdul[0].header['BPA']= 0
    hdul[0].header['BUNIT']= 'JY/BEAM'
    
    convl_image = 'img_rl_disk_convl_H{}.fits'.format(n_inf)
    fits.writeto(convl_image,convl_data_noise,hdul[0].header,overwrite=True)
    
    
    #--------------------------------------------------------------------------------
    # csub convolution
    #------------------------------------------------------------------------------
    print('\n')  
    print('Starting continum substracted convolution (convl_csub.py)')
    convl_csub_image = 'img_rl_disk_csub_convl_H{}.fits'.format(n_inf)
    hdul = fits.open(csub_image)
    new_beam = beam[0] #arcsec = 2 mas  
    npixsigma = (1/2.3548)*(new_beam/(3600*hdul[0].header['CDELT2']))
    print ("are {} pixels per beam".format(npixsigma))
    gauss_kernel = Gaussian2DKernel(npixsigma)
    pixperbeam = 1.442*(np.pi/4)*(new_beam**2)/((3600*hdul[0].header['CDELT2'])**2)
    convl_data = hdul[0].data.copy()
    convl_data_noise = hdul[0].data.copy()

    for element in range(hdul[0].data.shape[0]):
        convl_data[element,:,:] = pixperbeam*convolve(hdul[0].data[element,:,:],gauss_kernel)
        #add the noise after convolution
    #adding noise
    for element in range(hdul[0].data.shape[0]):
        for x in range(hdul[0].data.shape[2]):
            for y in range(hdul[0].data.shape[1]):
                convl_data_noise[element,y,x] = convl_data[element,y,x] + np.random.normal(0,scale=noise_input) #Jy/beam
        
    hdul[0].header['BMIN']= new_beam/3600
    hdul[0].header['BMAJ']= new_beam/3600
    hdul[0].header['BPA']= 0
    hdul[0].header['BUNIT']= 'JY/BEAM'
    fits.writeto(convl_csub_image,convl_data_noise,hdul[0].header,overwrite=True)
    hdul.close()
    #------------------------------------------------------
    ################ continuum image + noise
    print('\n')
    print('Continuum +noise image (convl_csub.py)')

    convl_c_image = 'img_rl_disk_cont_convl_H{}.fits'.format(n_inf)	#QZ
    hdul = fits.open(cont_image)
    new_beam = beam[0] #arcsec
    npixsigma = (1/2.3548)*(new_beam/(3600*hdul[0].header['CDELT2']))
    gauss_kernel = Gaussian2DKernel(npixsigma)
    pixperbeam = 1.442*(np.pi/4)*(new_beam**2)/((3600*hdul[0].header['CDELT2'])**2)
    convl_data = hdul[0].data.copy()
    convl_data_noise = hdul[0].data.copy()
    for element in range(hdul[0].data.shape[0]):
        convl_data[element,:,:] = pixperbeam*convolve(hdul[0].data[element,:,:],gauss_kernel)
    #adding noise
    for element in range(hdul[0].data.shape[0]):
        for x in range(hdul[0].data.shape[2]):
            for y in range(hdul[0].data.shape[1]):
                convl_data_noise[element,y,x] = convl_data[element,y,x] + np.random.normal(0,scale=noise_input) #Jy/beam

        
    hdul[0].header['BMIN']= new_beam/3600
    hdul[0].header['BMAJ']= new_beam/3600
    hdul[0].header['BPA']= 0
    hdul[0].header['BUNIT']= 'JY/BEAM'
    fits.writeto(convl_c_image,convl_data,hdul[0].header,overwrite=True)
    hdul.close()


    print('\n')
    print('Finished running convolution and continuum subtraction')
    print ('Ellapsed time for run RADMC: %.3fs' % (time.time() - t0))
    print ('-------------------------------------------------\n-------------------------------------------------\n')
    ###########################################################################################################
### agregar si queremos que saque las imagenes sin ruido


#lines_df = pd.read_csv("/fs/posgrado30/other0/jesus/paper2/Hydrogen_recom_lines_in_use.csv")
lines_df =  pd.read_csv("/home/jesus/Documents/paper2/Hydrogen_recom_lines_in_use.csv")
lines_df



from joblib import Parallel, delayed
## run all lines 
Parallel(n_jobs=1)(delayed(run_radmc)(lines_df.iloc[i],GRID,prop) for i in range(len(lines_df)))


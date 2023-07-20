import time
from matplotlib import colors
import numpy as np
import os
import subprocess
import time
import astropy.units as asu
import astropy.constants as cts
import pandas as pd
from sf3dmodels import Model, Plot_model
from sf3dmodels import Resolution as Res
import sf3dmodels.utils.units as u
import sf3dmodels.utils.constants as ct
import sf3dmodels.rt as rt
import sys
sys.path.append('/fs/posgrado30/other0/opt/star-forming-regions')
import sf3dmodels
import astropy.wcs as wcs
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
from spectral_cube import SpectralCube
from pvextractor import extract_pv_slice, Path, PathFromCenter
import radmc3dPy.image as image
from plot_helpers import plot_1d_props
#for the noise 
from ngVLA_sensitivity_calculator import *

def freq2vel(freq, freq_0):
    return (freq_0 - freq) * (cts.c.cgs.to(asu.km/asu.s)) / freq_0
def vel2freq(vel,freq_0):
    return freq_0 * (1.-(vel/(cts.c.cgs.to(asu.km/asu.s))))

#this is a copy of create full model
#def run_radmc(line_prop,GRID,prop,beam=(2e-3,2e-3),noise_input_line = 17.29055e-6, noise_input_continuum = 0.13844e-6,inclination=60,channel_width_0=4,number_channels=20):
def run_radmc(line_prop,GRID,prop,beam=(2e-3,2e-3),frequency = 93, time_int = 36000, inclination=60,channel_width=4,number_channels=20):
    '''
    vel_region in km/s
    noise in Jy/beam
    '''
    t0 = time.time()
    radmc_path = "/home/jesus/radmc-3d/version_0.40/examples/run_recomblines_userdef"
    #radmc_path = "/fs/posgrado30/other0/jesus/radmc-3d/version_0.40/examples/run_recomblines_userdef"
    #-------------------------------
    n_sup = line_prop["nsup"]
    n_inf = line_prop["ninf"]
    line_freq = line_prop["Frequency_radmc"]
    #print ("control: velocity line", freq2vel(line_freq,line_freq))
    channel_width_vel = channel_width * asu.km/asu.s
    channel_width_freq_0 = vel2freq(channel_width_vel,line_freq*asu.Hz)
    channel_width_freq = line_freq*asu.Hz - channel_width_freq_0 
    print ("channel width {} or {}".format(channel_width_vel, channel_width_freq.to(asu.MHz)))
    freq_ini = line_freq*asu.Hz - (channel_width_freq)*(number_channels/2)
    freq_fin = line_freq*asu.Hz + (channel_width_freq)*(number_channels/2)
    print ("{:.7e}, {:.7e}".format(freq_ini,freq_fin))
    line_wave = (line_prop["Frequency_radmc"]*asu.Hz).to(asu.micron, equivalencies=asu.spectral())
    #print (line_wave)
    wave_ini = freq_fin.to(asu.micron, equivalencies=asu.spectral())
    wave_fin = freq_ini.to(asu.micron, equivalencies=asu.spectral())
    wave_range=[abs(wave_ini.value),abs(wave_fin.value)]
    print (wave_ini,wave_fin)
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
    subprocess.run(radmc_rl+' image lambdarange {} {} nlam {} incl {} npix 221 sizeau 550'.format(wave_range[0],wave_range[1],number_channels,str(i)), shell=True, executable='/bin/bash') 
     
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
    hdul[0].header['VELREF'] =  257             # Local Standar of Rest: mean motion of material in the Milky Way in the neighborhood of the Sun
    hdul.flush()   #### save the file resulting from the radiative transfer.     
    print('\n')
    ############ separate the line cube from the continuum cube
    print ("separating the continuum data cube and the line data cube")
    #print('Separatig Starting continuum_subtraction (convl_csub.py)')
    hdul = fits.open(output)
    ############################### fin add noise ####################
    data = hdul[0].data  #units of Jy/pix  #output from radmc
    data_line = data.copy()
    data_cont = data.copy()
    for slide in range(data.shape[0]):
        for x in range(data.shape[2]):
            for y in range(data.shape[1]):
                data_line[slide,y,x] = data[slide,y,x] - data[:,y,x].min()
                data_cont[slide,y,x] = data[:,y,x].min()
    csub_image = 'img_rl_disk_line_H{}.fits'.format(n_inf)
    cont_image = 'img_rl_disk_cont_H{}.fits'.format(n_inf)
    fits.writeto(csub_image,data_line,hdul[0].header,overwrite=True)
    fits.writeto(cont_image,data_cont,hdul[0].header,overwrite=True)	#QZ
    ############--------------------------##############
    #Noise calculation
    #------  noises
    #continuum 
    a_cont,b_cont,c_cont = sigma_ps_fn("main+lba",  freq=frequency,type_cal = 'continuum', theta=beam[0],t_int=time_int, delta_v=channel_width*1e3, verbose=True )
    noise_input_continuum = a_cont[1]*1e-6   #Jy/beam  #original in uJy/beam
    print ("the continuum noise = {}".format(noise_input_continuum))
    #line
    a_line,b_line,c_line = sigma_ps_fn("main+lba",  freq=frequency, type_cal = 'line', theta=beam[0],t_int=time_int, delta_v=channel_width*1e3, verbose=True )
    noise_input_line = a_line[1]*1e-6     #Jy/beam  #original in uJy/beam
    print ("the line noise = {}".format(noise_input_line))
    #a[0] = ""#frequency
    #a[1] is the simga_rms
    #a[2] is the T_rms
    #### end: Noise calculation
    
    
    ########################    convolution and add noise
    print('\n')  
    print('Starting full image convolution (convl_csub.py)')
    hdul = fits.open(output)
    new_beam = beam[0] #arcsec = 2 mas  
    npixsigma = (1/2.3548)*(new_beam/(3600*hdul[0].header['CDELT2']))
    print ("are {} pixels per beam".format(npixsigma))
    gauss_kernel = Gaussian2DKernel(npixsigma)
    pixperbeam = 1.442*(np.pi/4)*(new_beam**2)/((3600*hdul[0].header['CDELT2'])**2)
    convl_data = hdul[0].data.copy()
    for element in range(hdul[0].data.shape[0]):
        convl_data[element,:,:] = pixperbeam*convolve(hdul[0].data[element,:,:],gauss_kernel)
    hdul[0].header['BMIN']= new_beam/3600
    hdul[0].header['BMAJ']= new_beam/3600
    hdul[0].header['BPA']= 0
    hdul[0].header['BUNIT']= 'JY/BEAM'
    convl_image = 'img_rl_disk_convl_H{}_control.fits'.format(n_inf)
    fits.writeto(convl_image,convl_data,hdul[0].header,overwrite=True)
    
    
    #--------------------------------------------------------------------------------
    # line convolution and add noise
    #------------------------------------------------------------------------------
    print('\n')  
    print('Starting line cube convolution and adding noise (convl_csub.py)')
    convl_csub_image = 'img_rl_disk_line_convl_H{}.fits'.format(n_inf)
    hdul = fits.open(csub_image)
    new_beam = beam[0] #arcsec = 2 mas  
    npixsigma = (1/2.3548)*(new_beam/(3600*hdul[0].header['CDELT2']))
    print ("are {} pixels per beam".format(npixsigma))
    gauss_kernel = Gaussian2DKernel(npixsigma)
    pixperbeam = 1.442*(np.pi/4)*(new_beam**2)/((3600*hdul[0].header['CDELT2'])**2)
    convl_data = hdul[0].data.copy()
    convl_data_noise = hdul[0].data.copy()
    for element in range(hdul[0].data.shape[0]):
        convl_data[element,:,:] = pixperbeam*convolve(data_line[element,:,:],gauss_kernel)
        #add the noise after convolution
    #adding noise
    for element in range(hdul[0].data.shape[0]):
        for x in range(hdul[0].data.shape[2]):
            for y in range(hdul[0].data.shape[1]):
                convl_data_noise[element,y,x] = convl_data[element,y,x] + np.random.normal(0,scale=noise_input_line) #Jy/beam
        
    hdul[0].header['BMIN']= new_beam/3600
    hdul[0].header['BMAJ']= new_beam/3600
    hdul[0].header['BPA']= 0
    hdul[0].header['BUNIT']= 'JY/BEAM'
    fits.writeto(convl_csub_image,convl_data_noise,hdul[0].header,overwrite=True)
    hdul.close()
 
    #------------------------------------------------------
    ################ continuum cube + noise
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
        convl_data[element,:,:] = pixperbeam*convolve(data_cont[element,:,:],gauss_kernel)
    #adding noise
    for element in range(hdul[0].data.shape[0]):
        for x in range(hdul[0].data.shape[2]):
            for y in range(hdul[0].data.shape[1]):
                convl_data_noise[element,y,x] = convl_data[element,y,x] + np.random.normal(0,scale=noise_input_continuum) #Jy/beam

        
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

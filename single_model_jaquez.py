'''
Run radmc Hamburger disk, optimized for many lines
'''
#******************************
#Handling fits, pvs, and coords
#******************************
import astropy.wcs as wcs
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
from spectral_cube import SpectralCube
from pvextractor import extract_pv_slice, Path, PathFromCenter
import radmc3dPy.image as image

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


#################################### run radmc ##########################
def run_radmc(line_prop,GRID,prop):
    t0 = time.time()
    radmc_path = "/home/jesus/radmc-3d/version_0.40/examples/run_recomblines_userdef"
    n_sup = line_prop["nsup"]
    n_inf = line_prop["ninf"]
    line_freq = line_prop["Frequency_radmc"]
    print (n_sup,n_inf)
    line_wave = (line_prop["Frequency_radmc"]*asu.Hz).to(asu.micron, equivalencies=asu.spectral())
    wave_range = [line_wave.value-1, line_wave.value+1]
    i = 53
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
    subprocess.run(radmc_rl+' image lambdarange {} {} nlam 60 incl '.format(wave_range[0],wave_range[1])+str(i)+' npix 201 sizeau 500', shell=True, executable='/bin/bash') 
     
    ################################### CONVOLUTION AND CONTINUUM SUBTRACTION ######################################
    # Use radmc3dpy to write fits file
    output = 'img_rl_disk.fits'
    #subprocess.call(['rm',output], shell=True)
    os.system('rm -f '+output)
    dist = 5000.
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
    hdul.flush()        
    print('\n')
    print('Starting continuum_subtraction (convl_csub.py)')
    csub_image = 'img_rl_disk_csub.fits'
    cont_image = 'img_rl_disk_cont.fits'
    hdul = fits.open(output)
    data = hdul[0].data
    data2 = data.copy()
    datac = data.copy()	#QZ
    for x in range(data.shape[2]):
    	for y in range(data.shape[1]):
    		data2[:,y,x] = data[:,y,x] - data[:,y,x].min()
    		datac[:,y,x] = data[:,y,x].min()

    fits.writeto(csub_image,data2,hdul[0].header,overwrite=True)
    fits.writeto(cont_image,datac,hdul[0].header,overwrite=True)	#QZ
    

    #convolve
    print('\n')
    print('Starting convolution (convl_csub.py)')

    convl_image = 'img_rl_disk_csub_convl.fits'
    hdul = fits.open(csub_image)
    new_beam = 2e-3 #arcsec = 2 mas  
    npixsigma = (1/2.3548)*(new_beam/(3600*hdul[0].header['CDELT2']))
    gauss_kernel = Gaussian2DKernel(npixsigma)
    pixperbeam = 1.442*(np.pi/4)*(new_beam**2)/((3600*hdul[0].header['CDELT2'])**2)
    convl_data = hdul[0].data.copy()
    for element in range(hdul[0].data.shape[0]):
    	convl_data[element,:,:] = pixperbeam*convolve(hdul[0].data[element,:,:],gauss_kernel)

    hdul[0].header['BMIN']= new_beam/3600
    hdul[0].header['BMAJ']= new_beam/3600
    hdul[0].header['BPA']= 0
    hdul[0].header['BUNIT']= 'JY/BEAM'
    fits.writeto(convl_image,convl_data,hdul[0].header,overwrite=True)
    hdul.close()
    #
    # QZ added, extract continuum image
    print('\n')
    print('Extract continuum image (convl_csub.py)')

    convl_c_image = 'img_rl_disk_cont_convl.fits'	#QZ
    hdul = fits.open(cont_image)
    new_beam = 2e-3 #arcsec
    npixsigma = (1/2.3548)*(new_beam/(3600*hdul[0].header['CDELT2']))
    gauss_kernel = Gaussian2DKernel(npixsigma)
    pixperbeam = 1.442*(np.pi/4)*(new_beam**2)/((3600*hdul[0].header['CDELT2'])**2)
    convl_data = hdul[0].data.copy()
    for element in range(hdul[0].data.shape[0]):
    	convl_data[element,:,:] = pixperbeam*convolve(hdul[0].data[element,:,:],gauss_kernel)

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
    
#####
# sf3dmodels
####

# from single_model_piecewise.py
# --------------
##### Parameters
# --------------
#### from piecewise_pv
#model constructions
#------------------
#General Parameters
#------------------
MStar = 10.0           #[Msun]  stellar mas
MStar = MStar * u.MSun
LStar = u.LSun * ( MStar/u.MSun )**4   #this is need for ?

#-------------------------------
#Parameters for the Pringle disc
#-------------------------------
RStar = u.RSun * ( MStar/u.MSun )**0.8
TStar = u.TSun * ( (LStar/u.LSun) / (RStar/u.RSun)**2 )**0.25
Rdisc = 250             #
Rd = Rdisc * u.au
print ('STARs PARAMETERS   RStar:', RStar/u.RSun, ', LStar:', LStar/u.LSun, ', TStar:', TStar)
############# rho(R=1 AU)
MRate = 3.33e-5 * u.MSun_yr   #From Roberto Paper 2022
RhoE0 = Res.Rho0(MRate, Rd, MStar)    #density from Ulrich envelope
i = 51.2                #inclination
Rrot_pivot = 150 #AU    #This value was the half of Rdics
pR_list_0 = -0.5
pR_list_1 = -0.5
# radial velocity parameters 
rrad_int = 150          # radius at what beggin the radial velocity
vr0 = 0                 # radial velocity
pr_list_0 = 0.0         # maybe a power law index

### from density_Env_Disc () #function
Arho = 5.29   #
rhoD0 = Arho * RhoE0      #Normalization factor based on the envelope
print ("density at Rdisc = {:.5e}".format(rhoD0))
#H0 = 0.01 * RStar         #Scaleheight at RStar
#H = H0 * (RList / RStar)**1.25  #Scaleheight
# equation Roberto paper , # = Z = 0
#Estamos normalizando a R=Rdisc # donde topa el disco con la envvoltura.
p_list_0 = -0.674
Rmin = 1 * u.au #* H0_factor      # Rmin dijo Roberto que lo voy a dejar en 1 au, 
### this is now our rho0 to create the disk density
rho_disk_Rmin = rhoD0 * (Rmin / Rd)**p_list_0 #* np.exp(-0.5 * zList**2 / H**2), 1.0)#1.0e9 / 2)
print ("density at r0(1 AU) = {:.5e}".format(rho_disk_Rmin))


print('Running model with:')
print('rho0: {}'.format(rho_disk_Rmin))
print('Rdisc: {}'.format(Rdisc))
print('p_list_0: {}'.format(p_list_0))
print('i: {}'.format(i))
print('MStar: {}'.format(MStar))
print('Rrot_pivot: {}'.format(Rrot_pivot))
print('pR_list_0: {}'.format(pR_list_0))
print('pR_list_1: {}'.format(pR_list_1))
print('rrad_int: {}'.format(rrad_int))
print('vr0: {}'.format(vr0))
print('pr_list_0: {}'.format(pr_list_0))


#---------------
#GRID Definition
#---------------
#box limits [-250,250]
sizex = sizey = sizez = 250 * u.au #halfsize
Nx = Ny = Nz = 201                   #Number of divisions for each axis
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], rt_code = 'radmc3d')#, include_zero=False)
NPoints = GRID.NPoints               # Final number of nodes in the grid
#-------------------
#PHYSICAL PROPERTIES
#-------------------
#-------
#DISC
#-------
H0_factor = 0.2                    # I'm going to leave it as is
#Rmin = 1 * sfu.au #* H0_factor      # Rmin dijo Roberto que lo voy a dejar en 1 au, 
#R_rho_1 = R_rho_1 * sfu.au
Rdisc = Rd
R_list_rho = [Rmin, Rdisc]           #Interval limits for piecewise density profile
p_list = [p_list_0]                  #Powerlaws within intervals
rho0 = rho_disk_Rmin*H0_factor**p_list_0     # in m^-3  #check why this.

RH_list = R_list_rho
q_list = [1. for p in p_list]  # H = H0 * (R / R0)**(1 + 0.5*(1-q))
#H0 = 1.058e11 #m at Rmin = 10au, from Hollenbach+94, such that H(R=2000au) = 2000 au, valid for q=1.5
#H0 = 1.496e12 * H0_factor #m at Rmin = 10au, such that H(R=2000au) = 2000 au, valid for q = 1.0
#1.496e12 = 10AU
H0 = 1.496e11 * H0_factor

density = Model.density_Hamburgers_piecewise(RStar, H0, R_list_rho, p_list, rho0, GRID,
                                              RH_list = RH_list, q_list = q_list,
                                              rho_thres = 1e7, rho_min = 0.0, Rt = False)


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
    
    
########## plots 

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
                   output = 'DensMidplane_%s.png'%tag, show = True)

######################
# read the lines data from Hydrogen_recom_lines_in_use.csv
######################
lines_df = pd.read_csv("/home/jesus/Documents/paper2/Hydrogen_recom_lines_in_use.csv")
line = lines_df.iloc[0]
lines_df.head()      
run_radmc(line,GRID,prop)
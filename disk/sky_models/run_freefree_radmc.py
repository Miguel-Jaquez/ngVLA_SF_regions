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



########################################## RADIATIVE TRANSFER WITH RADMC-3D #################################### 
lambda_ini = (116 * asu.GHz).to(asu.micron, equivalencies=asu.spectral())
lambda_fin = (96 * asu.GHz).to(asu.micron, equivalencies=asu.spectral())
#lambda_ini, lambda_fin

########################### RRL H38##############################################  
radmc_rl = os.path.join(radmc_path, 'radmc3d')
##### JAQUEZ H38_\alpha -> 2600.6855 microns
#subprocess.run(radmc_rl+' image lambdarange {} {} nlam 60 incl 30 npix 3000 sizeau 2400 setthreads 14 '.format(lambda_ini.value, lambda_fin.value), shell=True, executable='/bin/bash')

#os.system('mv image.out model_disk_continuum_incl30.out')


subprocess.run(radmc_rl+' image lambdarange 2584.4177 3944.6376 nlam 2 incl 30 npix 3000 sizeau 2400 setthreads 15', shell=True, executable='/bin/bash')

os.system('mv image.out model_disk_continuum_incl30_1channel.out')
 

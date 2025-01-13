import sys
#change this path for the correct ones in your computer.
import numpy as np
import os
sys.path.append("/share/Part1/jesus/SF_regions/utils/")
from casatasks.private import simutil
## Set the name of the configuration file for the ngVLA Main subarray:
conf_file = 'ngvla-revD.main.cfg' # main = core + spiral+ mid subarrays #Which anntena configuration we want 
#conf_dir = os.getenv('CASAPATH').split(' ')[0] + '/data/alma/simmos/' #'/home/jesus/CASA/casa-6.5.3-28-py3.8/data/alma/simmos/'
conf_dir = '/share/Part1/jesus/SF_regions/utils/simmos/'
conf_path = os.path.join(conf_dir, conf_file)
## Use simutil to read the .cfg file:
u = simutil.simutil()
xx,yy,zz,diam,padnames,_,telescope,posobs = u.readantenna(conf_dir+conf_file)

line = 'H38'
#for the noise 
from ngVLA_sensitivity_calculator import *
tsource = 18000

#a_cont,b_cont,c_cont = sigma_ps_fn("main",  freq=93,type_cal = 'continuum', theta=-1, t_int=tsource, verbose=True )
a_cont, b_cont, c_cont = sigma_ps_fn("main", freq = 93, type_cal = 'line',delta_v = 814653.42, theta = -1, t_int = tsource, verbose=True) #delta_v  = 250 MHz
sigma_NA = a_cont[1]*1e-6   #in Jy/beam  #original in uJy/beam 
print ("continuum noise per channel = {} Jy/beam".format(sigma_NA))

nantennas = len(xx)
nbaselines = nantennas*(nantennas-1)/2
print ('number of baselines: ', nbaselines)
npol = 2.0
nchan = 1.0
tint = 10
nintegrations = tsource/tint
noise_input_continuum = sigma_NA * np.sqrt(nchan* npol * nbaselines * nintegrations ) # jy

sigma_simple ='{}Jy'.format(noise_input_continuum)
print ("the input noise is {}".format(sigma_simple))


# In CASA

## Create a copy of the noise-free MS:
os.system('cp -r jet_93GHz_continuum/jet_93GHz_continuum.ngvla-revD.main.ms jet_93GHz_continuum_noisy.ms')

## Open the MS we want to add noise to with the sm tool:
sm.openfromms('jet_93GHz_continuum_noisy.ms')

## Set the noise level using the simplenoise parameter estimated in the section on Estimating the Scaling Parameter for Adding Thermal Noise:


sm.setnoise(mode = 'simplenoise', simplenoise = sigma_simple)

## Add noise to the 'DATA' column (and the 'CORRECTED_DATA' column if present):
sm.corrupt()

## Close the sm tool:
sm.done()

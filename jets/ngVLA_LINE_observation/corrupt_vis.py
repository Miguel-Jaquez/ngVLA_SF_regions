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

line = 'H41'
line_freq = 92.07 
#for the noise 
from ngVLA_sensitivity_calculator import *
tsource = 3600*5
a_line,b_line,c_line = sigma_ps_fn("main",  freq=line_freq,type_cal = 'line', theta=-1, t_int=tsource, delta_v=30.e3, verbose=True )
sigma_NA = a_line[1]*1e-6   #in Jy/beam  #original in uJy/beam 
print ("line noise = {} Jy/beam".format(sigma_NA))

nantennas = len(xx)
nbaselines = (nantennas*(nantennas-1))/2
print ('number of baselines: ', nbaselines)
npol = 2.0
nchan = 1.0
tint = 10
nintegrations = tsource/tint
noise_input = sigma_NA * np.sqrt(nchan* npol * nbaselines * nintegrations ) # jy

sigma_simple ='{}Jy'.format(noise_input)
print ("the input noise is {}".format(sigma_simple))


# In CASA

## Create a copy of the noise-free MS:
os.system('cp -r jet_93GHz_line_{0}_timeint_5/jet_93GHz_line_{0}_timeint_5.ngvla-revD.main.ms jet_93GHz_line_{0}_noisy.ms'.format(line))

## Open the MS we want to add noise to with the sm tool:
sm.openfromms('jet_93GHz_line_{}_noisy.ms'.format(line))

## Set the noise level using the simplenoise parameter estimated in the section on Estimating the Scaling Parameter for Adding Thermal Noise:


sm.setnoise(mode = 'simplenoise', simplenoise = sigma_simple)

## Add noise to the 'DATA' column (and the 'CORRECTED_DATA' column if present):
sm.corrupt()

## Close the sm tool:
sm.done()

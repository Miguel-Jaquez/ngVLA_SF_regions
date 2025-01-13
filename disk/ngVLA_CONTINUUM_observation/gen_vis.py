import os

## Set the name of the configuration file for the ngVLA Main subarray:
conf_file = 'ngvla-revD.main.cfg' # main = core + spiral+ mid subarrays #'ngvla-main-revC.cfg'   #Which anntena configuration we want 
conf_dir = '/share/Part1/jesus/SF_regions/utils/simmos/' #  '/home/jesus/CASA/casa-6.5.3-28-py3.8/data/alma/simmos/'
conf_path = os.path.join(conf_dir,conf_file)

## Position of the source that we want to observe: 
#'18h10m28.652s -19d55m49.66s'
# tuto 'J2000 00:00:00.0 +24.00.00.0'
direction = 'J2000 18:00:00.0 +25.00.00.0'

line = 'H38'

#model = 'disk_cont_BW20GHz_{}'.format(line)
model = '../model/img_rl_disk_cont_H38'
simobserve(project = 'disk_93GHz_cont_{}'.format(line), 
                    skymodel = model+'.fits', 
                    setpointings = True,
                    indirection = direction,
                    incenter = '106GHz',
                    inwidth = '333.33MHz',
                    integration = '10s',  
                    obsmode = 'int', 
                    antennalist = conf_path, 
                    hourangle = 'transit', 
                    totaltime = '18000s',  
                    thermalnoise = '',
                    overwrite=True,
                    graphics = 'none',
                    verbose = True)


importfits( fitsimage = model+'.fits', imagename = model+'.image',overwrite=True)
# exportfits ( '')

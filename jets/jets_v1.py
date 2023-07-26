#### create a jet
"""
Basic docstring explaining example
"""
from __future__ import print_function
import sys
#change this path for the correct ones in your computer.
sys.path.append('/fs/posgrado30/other0/opt/star-forming-regions')
sys.path.append('/fs/posgrado30/other0/jesus/radmc3dPy/lib/python3.9/site-packages')

#********************
#sf3dmodels libraries
#********************
import sf3dmodels.Plot_model as Pm            #Plotting model
from sf3dmodels.outflow import OutflowModel   #Model functions
import sf3dmodels.utils.units as u            #Units
import sf3dmodels.rt as rt                    #Writing functions for radiative transfer
import sf3dmodels.Model as Model              #Grid
from sf3dmodels.grid import Overlap           #Overlap submodels
from radmc3dPy.image import *
from matplotlib import cm
from matplotlib import pyplot as plt
sys.path.append("/home/jesus/Documents/paper2/ngVLA_SF_regions")
from utils_run_radmc import *
#********************
#Extra libraries
#********************
import numpy as np
import time

########### control
GRID_construction = True
RADMC_run = True


t0 = time.time()
if GRID_construction == True:
    #--------------------
    ### at which resolutions we need to put the jets grid.
    dx_grid = 2*u.au   # esto es la resolucion
    #*********************
    #OUTFLOW 1: standar collimated
    #*********************
    tag = '_outflow_collimated'
    #---------------------
    #GEOMETRIC PARAMETERS
    #---------------------
    pos_c = np.array([0, 0, 0]) #Origin of coordinates on the axis.
    axis = np.array([0,1,0])    #Long axis direction, of arbitrary length. For example, an axis pointing to the z-direction: axis = [0,0,1] 
    z_min = 25*u.au      #Lower and upper limits to compute the physical properties  # convierte au a m
    z_max = 150*u.au 
    #---------------------
    #PHYSICAL PARAMETERS
    #---------------------
    ### I (Jaquez) put the same M_dot at in the disc
    mdot = 1e-6*u.MSun_yr    # convierte masa solar por año a kg/s   ##va 1.5; 2 es para la propuesta   
    mu = 1.6735e-27          # Kg ; masa H
    cte = 4.9975e-7
    #---------------------
    #PHYSICAL PROPERTIES OF THE REYNOLDS 86 JET MODEL
    #---------------------
    r0=z_min               #inicio del jet
    v0 = 500 * 1e3         #  Gas speed at `z0`. Speed normalization factor. 500 es km/s, pero convierto a metros porque el codigo trabaja en mks
    #####
    # The formula is w = w0 (r/r0)^\epsilon
    #where w is the widht of the jet
    #####
    w0 = 5*u.au            # parameters to compute the Jet half-width at each z. poner en au y se transforma a m. 
    eps = 2./3.            # eps = 1 == corresponds to a conical (constant opening angle) jet
    w = [w0, eps]

    s = np.pi*w0**(2)                                   # calculado con ec 12 paper. w0 en m
    n =cte*mdot*mu**(-1)*v0**(-1)*s**(-1)               # en cm-3 # Reynolds equaation 16.5
    print('Densidad jet en la base = ', n, ' cm-3')
    
    # power law (?) for the temperature
    T0 = 10000.
    qT = 0
    temp = [T0, qT]
    
    #power law (?) for the density 
    dens0 = n * 1e6                                     # paso la densidad a m-3 porque el codigo trabaja en mks
    qn = -2*eps
    dens = [dens0, qn]

    ionfrac = [1, 0]           #all the gas is ionized

    abund = [1e-4, 0]          # parameters to compute the molecular abundance
    gtd = 100                  # gas to dust ratio
    #---------------------
    #COMPUTING MODEL PROPERTIES
    #---------------------
    Outf1 = OutflowModel(pos_c, axis, z_min, z_max, dx_grid) #Initializing Class with grid parameters
    Outf1.reynolds86(w, dens, ionfrac, temp, v0, abund, gtd, [0,0,0], r0) #Invoking the outflow model from Reynolds et al. 1986

    #---------------------
    #WRITING FOR RADMC3D
    #---------------------
    #Using the Radmc3d class
    prop1 = {'vel_x' : Outf1.x, 'vel_y' : Outf1.y, 'vel_z' : Outf1.z,
             'dens_e': Outf1.density_ion,
             'dens_ion': Outf1.density_ion,
             'temp_gas': Outf1.temperature}
    radmc1 = rt.Radmc3d(Outf1.GRID)
    radmc1.submodel(prop1, output='datatab'+tag+'.dat')
    print('Output columns', radmc1.columns)


    #*******************************************************************
    #OUTFLOW 2: conical recombining
    #*********************
    tag = '_outflow_conical'
    #---------------------
    #GEOMETRY
    #---------------------
    pos_c = np.array([0, 0, 0])   #[200*u.au, -50*u.au, 0]
    axis = np.array([0,1,0])      #[1,1,-1]
    z_min = 3*u.au
    z_max = 150*u.au

    #---------------------
    #PHYSICAL PARAMETERS
    #---------------------
    #also put the same accretion rate
    mdot = 1e-6*u.MSun_yr  # convierte masa solar por año a kg/s
    #---------------------
    #PHYSICAL PROPERTIES
    #---------------------
    r0=z_min
    v0 = 100 * 1e3 #100 es km/s, pero convierto a metros porque el codigo trabaja en mks
    w0 = 5*u.au  #  poner en au. La unidad es m
    eps = 1
    w = [w0, eps]

    T0 = 10000.
    qT = -0.5   
    temp = [T0, qT]

    s = np.pi*w0**(2)   # calculado con ec 12 paper. w0 en m
    n =cte*mdot*mu**(-1)*v0**(-1)*s**(-1)   # en cm-3
    print('Densidad viento en la base = ',n ,' cm-3')

    dens0 = n* 1e6   # paso la densidad a m-3 porque el codigo trabaja en mks
    qn = -2*eps   # -2
    dens = [dens0, qn]

    ionfrac = [0.1,-0.5]      #qx = -0.5 es como varia # 10% de ionizacion inicialmente

    abund = [1e-4, 0]
    gtd = 100

    #---------------------
    #COMPUTING MODEL
    #---------------------
    Outf2 = OutflowModel(pos_c, axis, z_min, z_max, dx_grid) #Initializing Class with grid parameters
    Outf2.reynolds86(w, dens, ionfrac, temp, v0, abund, gtd, [0,0,0], r0) #Invoking the outflow model from Reynolds et al. 1986

    #---------------------
    #WRITING FOR RADMC3D
    #---------------------
    #Using the Radmc3d class
    prop2 = {'vel_x' : Outf2.x, 'vel_y' : Outf2.y, 'vel_z' : Outf2.z,
            'dens_e': Outf2.density_ion,
             'dens_ion': Outf2.density_ion,
             'temp_gas': Outf2.temperature}
    radmc2 = rt.Radmc3d(Outf2.GRID)
    radmc2.submodel(prop2, output='datatab'+tag+'.dat')
    print('Output columns', radmc2.columns)

    #-------
    #TIMING
    #-------
    print ('Ellapsed time to create the grids properties for the two (separated) jets: %.3fs' % (time.time() - t0))
    print ('-------------------------------------------------\n-------------------------------------------------\n')
else:
    print ("The grids were not created")

if RADMC_run == True:
    #####################################-------------------------
    #JOIN the JETS 
    ########################## -----------------------------------
    #********
    #GRIDDING
    #********
    sizex = 150 * u.au       #500
    sizey = sizez = 150 * u.au    #500
    Nx = Ny = Nz = 300
    GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], rt_code='radmc3d')

    files = ['datatab_outflow_collimated.dat', 'datatab_outflow_conical.dat'] #

    #### I add the columns for the vel_x, vel_y, vel_z in order to have the same columns that in the model disk
    #### and to run recomblines instead that freefree in RADMC
    columns = ['id', 'x', 'y', 'z', 'dens_ion', 'dens_e', 'temp_gas', 'vel_x', 'vel_y', 'vel_z']
    outflows = Overlap(GRID)
    finalprop = outflows.fromfiles(columns, submodels = files, rt_code = 'radmc3d')
    
    
    ###########--------------------------------------------
    # Run RADMC-3D
    
    lines_df =  pd.read_csv("/home/jesus/Documents/paper2/Hydrogen_recom_lines_in_use.csv")
    #run radmc
    i = 0   #this is the index of the line that you want to simulate.
    
    run_radmc(lines_df.iloc[i],GRID,finalprop,beam=(2e-3,2e-3),inclination=0,freefree=True)
    
####################################### plots 

#********
#PLOTTING --------------------cuidado que no se esten Mezclando las definiciones de densidad
#******** 
density = finalprop['dens_e']/ 1e6 #dens. in cm^-3  
temperature = finalprop['temp_gas']

weight = 100 * np.mean(density)

Pm.plane2D(GRID, density, axisunit = u.au,
                   cmap = 'jet', plane = {'z': 0*u.au},
                   norm = "log", colorlabel = r'$[\rm cm^{-3}]$',
                   output = 'DensMidplane_%s.png'%tag, show = False)


#-----------------
#Plot for DENSITY
#-----------------
Pm.scatter3D(GRID, density, weight, NRand = 4000, axisunit = u.au, colorscale = 'log', cmap = 'cool',
             colorlabel = r'${\rm log}_{10}(n [cm^{-3}])$', output = 'global_grid_dens.png', vmin = 5, show=False)




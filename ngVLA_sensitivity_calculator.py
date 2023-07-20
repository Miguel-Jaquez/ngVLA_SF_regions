#!/usr/bin/env python
'''ngVLA subarray sensitivity and key performance metric calculator 
V. Rosero '''



import sys
import pickle
import argparse        
import numpy as np
from scipy.interpolate import PPoly,BSpline,CubicSpline

master_path = '/home/jesus/Documents/paper2/ngVLA_SF_regions/'

## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&    
## 
## Function to calculate sensitivity of a subarray
##
## input parameters:
## 
##     subarray ## name of array (string): sba, core, spiral, mid, lba, main, 
##                 main+lba, spiral+mid, spiral+core mid+lba
## 
##     freq     ## frequency in GHz (float): e.g., 17
##         
##     theta    ## resolution in arcsec (float): e.g., 0.5
##                 theta=-1 will calculate native resolution (natural, no taper)
## 
##     t_obs    ## on-source time in hours (float): e.g., 1.  [optional: default=1.0]
## 
##     delta_v  ## channel width in m/s (float): e.g., 2e3.  [optional: default=10e3]
##
## 
## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


def calculate_sensitivity( subarray, freq=-1, theta=-1, t_obs=1., delta_v=10e3, verbose=True ):
    
    if subarray == 'sba':
        D = 6.            
        if theta != -1:
            print('\n***WARNING*** The calculator does not support tapering of the sba. Values will be reported for natural resolution.')
            theta=-1
    else:
        D = 18.

        
    a, b, c = sigma_ps_fn( subarray=subarray, freq = freq, D = D, theta= theta, delta_v=delta_v,  verbose=False)        

    type_cal = 'continuum'
    print('\n\n {0} calculation:'.format(type_cal))
    a, b, c = sigma_ps_fn( subarray=subarray, freq = freq, D = D, type_cal = type_cal, delta_v=delta_v, t_int= t_obs*3600., theta= theta, verbose=verbose)
    band_key, theta2, delta_freq, eta_w, freq2 = c[0], c[1], c[2], c[3], a[0]
    str1 = '{0} array, {4:.1f} hours at frequency {1} GHz ({2}) with resolution= {3} mas  (eta_w = {5:.2f})'.format(subarray, freq2, band_key, theta2*1e3, t_obs, eta_w)
    print('\n'+str1+'\n'+'-'*len(str1))
    print('continuum point source sensitivity: {0} uJy/beam'.format(a[1] ))
    print('continuum brightness sensitivity: {0} K'.format(a[2] ))

    type_cal = 'line'
    print('\n\n {0} calculation:'.format(type_cal))
    a, b, c = sigma_ps_fn( subarray=subarray, freq = freq, D = D, type_cal = type_cal, delta_v=delta_v, t_int= t_obs*3600., theta= theta, verbose=verbose)
    band_key, theta2, delta_freq, eta_w, freq2 = c[0], c[1], c[2], c[3], a[0]
    str1 = '{0} array, {4:.1f} hours at frequency {1} GHz ({2}) with resolution= {3} mas  (eta_w = {5:.2f})'.format(subarray, freq2, band_key, theta2*1e3, t_obs, eta_w)
    print('\n'+str1+'\n'+'-'*len(str1))
    print('line point source sensitivity in a {0} km/s ({1:.3f} kHz) channel: {2}uJy/beam'.format(delta_v/1e3, delta_freq/1e3, a[1] ))
    print('line brightness sensitivity in a {0} km/s ({1:.3f} kHz) channel: {2}K'.format(delta_v/1e3, delta_freq/1e3, a[2] ))



                       
## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&    
##                 
## Function to calculate the Point Source  Sensitivity  (sigma_rms)
## and the Surface Brightness Sensitivity  (sigma_T)
##
## using the equations and values from ngVLA memo 17     
##              
## c = 2.9979246e8 ## [m/s] speed of light
## k_B = 1.38e-23  ## [J/K]  Boltzmann Constant
## eta_c = 0.98    ## correlator efficiency 
## n_pol = 2.      ## number of polarizations
## eta_Q =  0.9625 ## Digitizer quantization efficiency
## D = 18.         ## [m] Diameter of the antenna 
## N_ant           ## Total number of antennas in subarray (from table)
## T_s             ## System Temperature (from table)
## t_int           ## [s] integration (or observation) time
## A               ## [m^2] Geometric area of single antenna
## eta_A           ## aperture (antenna) efficiency (from table)
## SEFD            ## [Jy] System Equivalent flux Density  of an antenna
## delta_v         ## [m/s] velocity resolution
## delta_nu        ## [Hz] correlated bandwidth  (continuum: from table, line: calculated from delta_v)
## theta           ## [arcsec] resolution (clean beam size)
## b_max           ## [m] physical length of longest baseline in subarray
##
## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


def sigma_ps_fn( subarray = 'core', freq = -1, type_cal = 'continuum', \
                t_int = 3600., N_ant = 214., k_B = 1.38e-23, delta_v = 10e3, theta = -1, verbose=False, \
                eta_c = 0.98, eta_Q =  0.9625, n_pol = 2., D = 18., c = 2.9979246e8 ):

    if subarray == 'sba':
        D = 6.            
        if theta != -1:
            print('\n***WARNING*** The calculator does not support tapering of the sba. Values will be reported for natural resolution.')
            theta=-1
    else:
        D = 18.
 
    
    if isinstance(freq,str):  nu = float( freq.rstrip('GHz') )
    else:  nu = freq

    if not ((theta == -1) or (theta > 0)):
        print ('theta {0} must be either -1 or >0'.format(theta))
        sys.exit(0)
        
    try:
        with open( master_path + 'receiver_data.pkl','rb') as in1:
            receiver_data = pickle.load(in1)                           
    except:
        with open( master_path + 'receiver_data.pkl','rb') as in1:
            receiver_data = pickle.load(in1,encoding='latin1')                           


    try:
        with open( master_path + 'subarray_data.pkl','rb') as in1:
            subarray_parameter_data = pickle.load(in1)                           
    except:
        with open( master_path + 'subarray_data.pkl','rb') as in1:
            subarray_parameter_data = pickle.load(in1,encoding='latin1')                                   
    
    try:
        subarray_parameters = subarray_parameter_data[subarray]
    except:
        print ('subarray {0} must be one of:\nsba, core, spiral, mid, lba, main, main+lba, spiral+mid, spiral+core, \nmid+lba'.format(subarray))
        sys.exit(0)



    # find band from frequency
    band_keys = receiver_data.keys()
    valid_band_keys = []
    for key in band_keys:
        band_freqs = receiver_data[key]['freq']             
        if (nu < np.max(band_freqs)) and (nu >= np.min(band_freqs)):
            valid_band_keys.append(key)


    if len(valid_band_keys) >= 1:
        if verbose: print('frequency {0} GHz found in band(s): {1}'.format(nu, valid_band_keys))
    if len(valid_band_keys) > 1:
        band_key = valid_band_keys[-1]
        if verbose: print('selecting band {0} based on lower Tsys'.format( band_key))
    elif len(valid_band_keys) == 0:
        print('frequency {0} GHz not found within any receiver bands'.format(nu))
        sys.exit(0)
    else:
        band_key = valid_band_keys[0]        
            
        
    max_bw = receiver_data[band_key]['max_bw']     
    band_freqs = receiver_data[band_key]['freq'] 
    fC = receiver_data[band_key]['freq_center'] 
    all_T_s = receiver_data[band_key]['tSys'] 
    all_eta_A = receiver_data[band_key]['antEff'] 

    fL, fH = np.min(band_freqs), np.max(band_freqs)
    cs_T_s = CubicSpline(band_freqs, all_T_s)
    cs_eta_A = CubicSpline(band_freqs, all_eta_A) 
    
    b_max = subarray_parameters['b_max']
    b_min = subarray_parameters['b_min']
    N_ant = subarray_parameters['N_ant']
    
    if subarray_parameters['spline_type'] == 'cubic':
        cs = PPoly(*subarray_parameters['spline_params'])
    elif subarray_parameters['spline_type'] == 'univariate':
        cs = BSpline(*subarray_parameters['spline_params'])
    


    if type_cal == 'line':
        T_s = cs_T_s(nu)
        eta_A = cs_eta_A(nu)
    else:
        if band_key != 'BAND_6':
            x_freqs = np.linspace(fL,fH,100)
            nu = np.mean([fL,fH])
            if verbose: print('using entire continuum bandwidth ({1}-{2} GHz) at center frequency {0} GHz'.format(nu,fL,fH))

        else:
            if (nu-10 < fL): 
                if verbose: print('continuum bandwidth extends beyond receiver edges, shifting center frequency from {0} to {1}'.format(nu, fL+10))
                nu = fL + 10
                if verbose: print('using continuum bandwidth ({1}-{2} GHz) at center frequency {0} GHz'.format(nu,nu-10,nu+10))
            elif (nu+10 > fH): 
                if verbose: print('continuum bandwidth extends beyond receiver edges, shifting center frequency from {0} to {1}'.format(nu, fH-10))
                nu = fH - 10
                if verbose: print('using entire continuum bandwidth ({1}-{2} GHz) at center frequency {0} GHz'.format(nu,nu-10,nu+10))
            x_freqs = np.linspace(nu-10,nu+10,100) 

        T_s = np.mean(cs_T_s(x_freqs))
        eta_A = np.mean(cs_eta_A(x_freqs))
    if verbose: print('using interpolated  Tsys: {0} K at frequency {2} GHz (band average: {1})'.format(T_s,receiver_data[band_key]['Tsys'],nu))
    if verbose: print('using interpolated eta_A: {0} at frequency {2} GHz (band average: {1})'.format(eta_A,receiver_data[band_key]['eta_A'],nu))        

    
    if theta == -1: 
        eta_w = 1.
        theta = subarray_parameters['theta_nat_30_GHz'] * 30./nu /1e3 #arcsec
        if verbose: print('using native resolution {0} mas at frequency {1} GHz'.format(theta*1e3,nu))        
    else:
        eta_w = float( cs( np.log10( theta*1e3*nu/30. ) ) )
        if verbose: print('using inefficiency factor eta_w: {0} at frequency {1} GHz and resolution {2} mas)'.format(eta_w,nu,theta*1e3))

    
    delta_nu = nu * delta_v/c * 1e9  ## line width in Hz    
    A = np.pi * (D/2.)**2  ## [m**2]
    Field_view = 1.02*(206265/60.)*c/(D*nu*1e9)   ## arcmin  Uniform illumination
    total_eff_A= eta_A * N_ant * A
    Res_max_base = (206265.)*c*1e3/(nu*1e9*b_max)  ## marcsec
    LAS = (206265.)*c/(nu*1e9*b_min)  ## arcsec
    
    
    SEFD = (2* k_B * T_s/ (eta_Q * eta_A * A))/1e-26   ## [Jy]

    avg_SEFD = SEFD * (receiver_data[band_key]['Tsys'] / T_s) * (receiver_data[band_key]['eta_A'] / eta_A)

    if verbose: print('calculated SEFD: {0} Jy at frequency {1} GHz (band average: {2})'.format(SEFD,nu,avg_SEFD))        


    data_to_print_performance = (nu, fL, fH,   max_bw, \
                                 Field_view, eta_A, \
                                 total_eff_A/1e3, T_s,\
                                 SEFD, Res_max_base, LAS)

    band_info = (band_key, theta, delta_nu, eta_w)
    
    if type_cal == 'continuum':
        delta_nu = max_bw * 1e9  ## GHz
    elif type_cal == 'line':        
        delta_nu = nu * delta_v/c * 1e9  ## line width in GHz


    sigma_ps = (SEFD / (eta_c * np.sqrt(n_pol* delta_nu * t_int * N_ant*(N_ant-1))))/1e-6  ##uJy        
    sigma_rms = eta_w * sigma_ps
        
    sigma_T = 1.216 * (sigma_rms/(nu**2 * theta**2))   ## K when sigma_ps in uJy
    data_to_print = (nu,   sigma_rms,  sigma_T)

    
    return data_to_print, data_to_print_performance, band_info




## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&    
##                 
## Code to support command line execution
##    
## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


subarray_help = '''name of array (string): sba, core, spiral, mid, lba, main, 
                main+lba, spiral+mid, spiral+core, mid+lba'''

freq_help = ''' frequency in GHz (float): e.g., 24.5'''

theta_help = '''resolution in arcsec (float or -1): e.g., 0.5. 
                theta=-1 will calculate native resolution (natural, no taper)'''

t_obs_help = '''on-source time in hours (float): e.g., 4.  Default: 1'''

delta_v_help = '''channel width in m/s (float): e.g., 2e3.  Default: 10e3'''


    
if __name__ == '__main__':     

    parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('subarray', type=str, help=subarray_help)
    parser.add_argument('frequency', type=float, help=freq_help)
    parser.add_argument('theta', type=float, help=theta_help)
    parser.add_argument('--t_obs', type=float, default=1.0, help=t_obs_help)
    parser.add_argument('--delta_v',type=float, default=10e3, help=delta_v_help)
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    args = parser.parse_args()    
    print(__doc__)
    calculate_sensitivity( args.subarray, args.frequency, args.theta, args.t_obs, args.delta_v, args.verbose )    
    


#

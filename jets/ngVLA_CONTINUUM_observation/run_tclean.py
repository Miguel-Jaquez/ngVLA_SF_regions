import time
#### first  without noise
#### lo desmodifique, revisarlo si se vuelve a correr
## noisy
t0 = time.time()
line = 'H41'
def run_tclean(msfile,n,r):
    if n == 0 :
       imagename = 'psf_model_n0_r{}'.format(r)
    else:
        imagename = 'jet_continuum_r{}_n{}'.format(r,1)
        interactive = True
    #os.system('rm -rf '+imagename+'.*')
    tclean(vis=msfile, 
            datacolumn='data', 
            imagename=imagename, 
            imsize = [2400], 
            startmodel='', 
            cell='0.355mas', 
            specmode='mfs',  
            gridder='standard', 
            deconvolver='mtmfs', scales=[0,6,18,54], 
            weighting='briggs',robust=r, uvtaper='2.0mas',
            threshold = '1.2uJy',niter=n, interactive=True)

    return n


msfile = 'jet_93GHz_continuum_noisy.ms'
n = 10000000
r = -0.5 

run_tclean(msfile,n,r)

print ('the tclean finish in {} seconds or {} minuts'.format(time.time()-t0, (time.time()-t0)/60))

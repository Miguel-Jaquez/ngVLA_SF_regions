import time
#### first  without noise
#### lo desmodifique, revisarlo si se vuelve a correr
## noisy
t0 = time.time()
line = 'H41'
def run_tclean(msfile,n,r):
    if n == 0 :
       imagename = 'psf_model_n0_r{}_uvtaper1.5mas'.format(r)
    else:
        imagename = 'jet_{}_line_r{}_n{}_uvtaper1.5mas'.format(line,r,1)
        interactive = True
    #os.system('rm -rf '+imagename+'.*')
    tclean(vis=msfile, 
            datacolumn='data', 
            imagename=imagename, 
            imsize =[1080], startmodel='', 
            cell='0.8mas', 
            specmode='cube',  
            gridder='standard', 
            deconvolver='multiscale', scales=[0,6,18], 
            weighting='briggs',robust=r, uvtaper = '1.5mas', 
            threshold = '10.0uJy',niter=n, interactive=True)
    return n


msfile = 'jet_93GHz_line_{}_noisy.ms'.format(line)
n = 1000000
r = 0.5


run_tclean(msfile,n,r)

print ('the tclean finish in {} seconds or {} minuts'.format(time.time()-t0, (time.time()-t0)/60))


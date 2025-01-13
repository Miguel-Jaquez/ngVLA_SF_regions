import time
#### first  without noise
#### lo desmodifique, revisarlo si se vuelve a correr


t0 = time.time()
line = 'H38'
def run_tclean(msfile,n,r):
    if n == 0 :
       imagename = 'psf_model_n0_r{}'.format(r)
    else:
        imagename = 'disk_{}_line_r{}_n{}'.format(line,r,1)
        interactive = True
    #os.system('rm -rf '+imagename+'.*')
    tclean(vis=msfile, 
            datacolumn='data', imagename=imagename, 
            imsize = [1200],startmodel='', 
            cell='0.5mas', specmode='cube',  
            gridder='standard', 
            deconvolver='multiscale',scales=[0,7,21], 
            weighting='briggs',robust=r, uvtaper='2.0mas', 
            threshold = '50.0uJy',niter=n, interactive=True)
    return n


msfile = 'disk_93GHz_line_{}_noisy.ms'.format(line)
n =  100000
r = - 0.5

run_tclean(msfile,n,r)

print ('the tclean finish in {} seconds or {} minuts'.format(time.time()-t0, (time.time()-t0)/60))


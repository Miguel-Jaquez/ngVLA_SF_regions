import time
#### first  without noise
#### lo desmodifique, revisarlo si se vuelve a correr
## noisy
t0 = time.time()
line = 'H38'
def run_tclean(msfile,n,r):
    if n == 0 :
       imagename = 'psf_model_n0_r{}'.format(r)
    else:
        imagename = 'disk_{}_cont_r{}_n{}'.format(line,r,1)
        interactive = True
    #os.system('rm -rf '+imagename+'.*')
    tclean(vis=msfile, datacolumn='data', imagename=imagename, 
            imsize = [2400],startmodel='', 
            cell='0.25mas', 
            specmode='mfs',  gridder='standard', 
            deconvolver='mtmfs', scales=[0,8,24,60], 
            weighting='briggs',robust=r, uvtaper='2.0mas', 
            threshold = '0.40uJy',niter=n,interactive=True)
    return n


msfile = 'disk_93GHz_cont_{}_noisy.ms'.format(line)
n = 1000000
r = -0.5

run_tclean(msfile,n,r)

print ('the tclean finish in {} seconds or {} minuts'.format(time.time()-t0, (time.time()-t0)/60))


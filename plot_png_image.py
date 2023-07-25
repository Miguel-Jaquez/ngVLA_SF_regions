from radmc3dPy.image import *
from matplotlib import cm
from matplotlib import pyplot as plt

a=readImage()
plotImage(a,log=True,maxlog=4,cmap=cm.hot,bunit='snu',dpc=725.,au=True)
a.writeFits('Jet_wind.fits', dpc=725., coord='22h56m17.985s 62d01m49.55s')


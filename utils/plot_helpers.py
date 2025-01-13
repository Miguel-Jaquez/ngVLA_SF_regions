"""
Plotting functions, should be used with posterior distribution of parameters.
Making plt plots during emcee execution must be avoided.
"""

import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.ticker import FormatStrFormatter
import sf3dmodels.utils.units as sfu
from myplt_style_andizq import *

#SMALL_SIZE = 10
#MEDIUM_SIZE = 15

class PlotTools: #Extracted from sf3dmodels.model.disc2d
    @staticmethod
    def mod_nticks_cbars(cbars, nbins=5):
        for cb in cbars:
            cb.locator = ticker.MaxNLocator(nbins=nbins)
            cb.update_ticks()

    @staticmethod
    def mod_major_ticks(ax, axis='both', nbins=6):
        ax.locator_params(axis=axis, nbins=nbins)

    @staticmethod
    def mod_minor_ticks(ax):
        ax.minorticks_on()
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2)) #1 minor tick per major interval                                                                                                             
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))

    @classmethod
    def make_up_ax(cls, ax, xlims=(None, None), ylims=(None, None),
                   mod_minor=True, mod_major=True, **kwargs_tick_params):
        '''
        Set plot lims and tick params all at once.
        If mod_minor and/or mod_major call mod_minor_ticks and/or mod_major_ticks functions to control tick freq.
        '''
        kwargs_t = dict(labeltop=True, labelbottom=False, top=True, right=True, which='both', direction='in')
        kwargs_t.update(kwargs_tick_params)
        if mod_major: cls.mod_major_ticks(ax)
        if mod_minor: cls.mod_minor_ticks(ax)
        ax.set_xlim(*xlims)
        ax.set_ylim(*ylims)
        ax.tick_params(**kwargs_t)

make_up_ax = PlotTools.make_up_ax
mod_nticks_cbars = PlotTools.mod_nticks_cbars
cmap_im = 'gnuplot2_r'

#******************************
#SPECIFIC PLOTS FOR G10 PROJECT
#******************************
def plot_corner(samples, labels=None, quantiles=None):
    """Plot corner plot to check parameter correlations. Requires the 'corner' module"""
    import corner
    quantiles = [0.16, 0.5, 0.84] if quantiles is None else quantiles
    corner.corner(samples, labels=labels, title_fmt='.4f', bins=30,
                  quantiles=quantiles, show_titles=False)

def plot_walkers(samples, best_params, nstats=None, header=None, kind=None, tag=''):
    """function from sf3dmodels.model.disc2d, copy-pasted to avoid unnecessary imports"""
    npars, nsteps, nwalkers = samples.shape
    if kind is not None:
        ukind, neach = np.unique(kind, return_counts=True)
        ncols = len(ukind)
        nrows = np.max(neach)
    else:
        ukind = [''] 
        ncols = 1
        nrows = npars
        kind = ['' for i in range(nrows)] 

    if header is not None:
        if len(header) != npars: raise InputError(header, 'Number of headers must be equal to number of parameters')

    kind_col = {ukind[i]: i for i in range(ncols)}
    col_count = np.zeros(ncols).astype('int')
    print (kind_col)
    
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(3*ncols, 3*nrows))
    if nrows==1 and ncols==1: ax = [ax] 
    x0_hline = 0
    for k, key in enumerate(kind):
        j = kind_col[key]
        i = col_count[j] 
        for walker in samples[k].T:
            if ncols==1: axij = ax[i]
            elif nrows==1: axij = ax[j]
            else: axij = ax[i][j]
            axij.plot(walker, alpha=0.1, lw=1.0, color='k')
            if header is not None: 
                axij.text(0.1, 0.1, header[k], va='center', ha='left', fontsize=MEDIUM_SIZE+1, transform=axij.transAxes, rotation=0) #0.06, 0.95, va top, rot 90
        if i==0: axij.set_title(key, pad=10, fontsize=MEDIUM_SIZE+2)
        if nstats is not None: 
            axij.axvline(nsteps-nstats, ls=':', lw=2, color='r')
            x0_hline = nsteps-nstats

        axij.plot([x0_hline, nsteps], [best_params[k]]*2, ls='-', lw=3, color='dodgerblue')
        axij.text((nsteps-1)+0.03*nsteps, best_params[k], '%.3f'%best_params[k], va='center', color='dodgerblue', fontsize=MEDIUM_SIZE, rotation=90) 
        axij.tick_params(axis='y', which='major', labelsize=SMALL_SIZE, rotation=45)
        axij.set_xlim(None, nsteps-1 + 0.01*nsteps)
        col_count[j]+=1

    for j in range(ncols):
        i_last = col_count[j]-1
        if ncols==1: ax[i_last].set_xlabel('Steps')
        elif nrows==1: ax[j].set_xlabel('Steps')
        else: ax[i_last][j].set_xlabel('Steps')
        if nrows>1 and i_last<nrows-1: #Remove empty axes
            for k in range((nrows-1)-i_last): ax[nrows-1-k][j].axis('off')



def plot_pv_comparison(obs_data, model_data, pos_pixscale=100, velo_pixscale=2.0, nametag='Keplerian'):
    # Definir objetos de figura y ejes con su forma y tamaño
    fig, ax = plt.subplots(1,3,figsize=(16,4)) #,sharey='row'
    #axpos = ax[1].get_position()
    vmin, vmax = obs_data.min(), obs_data.max()

    # Llenar los arreglos con datos sobre los respectivos objetos eje con imshow    
    im0 = ax[0].imshow(obs_data, cmap='viridis', origin='lower', norm=colors.PowerNorm(gamma=0.5, vmin=vmin, vmax=vmax))
    im1 = ax[1].imshow(model_data, cmap='viridis', origin='lower', norm=colors.PowerNorm(gamma=0.5, vmin=vmin, vmax=vmax))

    # Proporciones de los ejes
    ax[0].set_aspect(0.5)
    ax[1].set_aspect(0.5)

    # Etiquetas de los ejes
    ax[0].set_title('Data')
    ax[1].set_title(nametag+' model')
    ax[0].set(xlabel='Position [au]', ylabel=r'Velocity [km s$^{-1}$]')
    ax[1].set(xlabel='Position [au]')

    # Agregar contornos a los objetos eje
    ax[0].contour(obs_data, colors='black', vmin=vmin, vmax=vmax)
    ax[1].contour(model_data, colors='black', vmin=vmin, vmax=vmax)

    #plt.colorbar(im1, ax=ax[1], shrink=0.8, pad=0.01)
    plt.colorbar(im1, ax=ax[0], shrink=1.0, pad=0.02)
    
    # Make residuals plot
    #axr = fig.add_axes([axpos.x1+0.05, axpos.y0, axpos.width, axpos.height])
    
    residuals = obs_data-model_data
    maxr = np.max(np.abs(residuals))
    #imr = axr.imshow(residuals, cmap='RdBu_r', origin='lower', vmin=-maxr, vmax=maxr)
    imr = ax[2].imshow(residuals, cmap='RdBu_r', origin='lower', vmin=-maxr, vmax=maxr)
    ax[2].set_aspect(0.5)

    cbar = plt.colorbar(imr, ax=ax[2], shrink=1.0, pad=0.02)
    cbar.set_label(r'Jy beam$^{-1}$')
    ax[2].set_title('Residuals')
    ax[2].set(xlabel='Position [au]')
    ax[1].tick_params(labelleft=False)
        
    #Custom ticks
    v0 = (obs_data.shape[0]-1)/2
    v_ticks = [v0-30, v0-20, v0-10, v0, v0+10, v0+20, v0+30]
    v_ticklabels = []
    for i in range(len(v_ticks)):
        v_label = '{:.1f}'.format((v_ticks[i]-v0)*velo_pixscale)
        v_ticklabels.append(v_label)     
    ax[0].set_yticks(v_ticks)
    ax[0].set_yticklabels(v_ticklabels)
    ax[1].set_yticks(v_ticks)
    ax[1].set_yticklabels(v_ticklabels)
    ax[2].set_yticks(v_ticks)
    ax[2].set_yticklabels(v_ticklabels)
    
    x0 = (obs_data.shape[1]-1)/2 #14
    x_ticks = [x0-14, x0-7, x0, x0+7, x0+14]
    x_ticklabels = []
    for i in range(len(x_ticks)):
        x_label = str(int((x_ticks[i]-x0)*pos_pixscale))
        x_ticklabels.append(x_label)
    ax[0].set_xticks(x_ticks)
    ax[0].set_xticklabels(x_ticklabels)
    ax[1].set_xticks(x_ticks)
    ax[1].set_xticklabels(x_ticklabels)
    ax[2].set_xticks(x_ticks)
    ax[2].set_xticklabels(x_ticklabels)    
 
    #set spaces between subplots   
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.05)

    #Save and show figure
    figname= 'pv_comparison_'+nametag+'.png'
    plt.savefig(figname, bbox_inches='tight',dpi=150)
    plt.show()
    plt.close()

def plot_cont_comparison(obs_data, model_data, pos_pixscale=100, figname='cont_comparison.png'):
    # Definir objetos de figura y ejes con su forma y tamaño
    fig, ax = plt.subplots(2,1,sharex=True,figsize=(10,5))
    vmin, vmax = obs_data.min(), obs_data.max()

    # Llenar los arreglos con datos sobre los respectivos objetos eje con imshow
    im0 = ax[0].imshow(obs_data, cmap='viridis', origin='lower', norm=colors.PowerNorm(gamma=0.5, vmin=vmin, vmax=vmax))
    im1 = ax[1].imshow(model_data, cmap='viridis', origin='lower', norm=colors.PowerNorm(gamma=0.5, vmin=vmin, vmax=vmax))

    # Etiquetas de los ejes
    ax[0].set_title('Continuum Data', fontdict={'fontsize': 25})
    ax[1].set_title('Model', fontdict={'fontsize': 25})
    ax[1].set_xlabel('Position [au]', fontdict={'fontsize': 22})
    
    #Custom ticks
    x0 = (obs_data.shape[1]-1)/2
    x_ticks = [x0-4,x0-3,x0-2,x0-1,x0,x0+1,x0+2,x0+3,x0+4]
    x_ticklabels = []
    for i in range(len(x_ticks)):
        x_label = str(int((x_ticks[i]- x0)*pos_pixscale))
        x_ticklabels.append(x_label)
    ax[1].set_xticks(x_ticks)
    ax[1].set_xticklabels(x_ticklabels)

    for axi in ax: axi.tick_params(labelleft=False)
    
    cbar = plt.colorbar(im1, ax=ax, shrink=0.8, pad=0.01)
    cbar.set_label(r'Jy beam$^{-1}$', fontdict={'fontsize': 22})
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(16)
    plt.savefig(figname, bbox_inches='tight',dpi=150)
    plt.show()
    plt.close()

def plot_1d_props(GRID, prop, density_obj, tag=''):
    #*********************
    #1D PLOTTING
    #*********************
    kwargs_sc_filled = dict(s=50, facecolor='none', linewidth=1.5)
    kwargs_sc_unfill = dict(s=80, linewidth=2.0)
    sizex = np.max(GRID.XYZgrid[0])
    #planez = np.unique(abs(GRID.XYZ[2]))[5] #plane value on z to be analyzed
    #print ('Considering the followingn Z to plot 1D densities and velocities:', planez/sfu.au)
    #indz = GRID.XYZ[2] == planez
    indz = np.abs(GRID.XYZ[2]) < 1.0 #indices where z==0 (not exactly zero due to small machine errors, that's why < 1.0)
    #indz = (np.abs(GRID.XYZ[2]) < np.min(np.abs(GRID.XYZ[2]))+0.1*sfu.au) & (GRID.XYZ[2]>0) #avoid midplane
    indphi = GRID.rRTP[3][indz] == np.pi/4 #indices where phi, on z=0, equals pi/4 from the grid

    fig, ax = plt.subplots(ncols=2,nrows=2, figsize=(14,11))
    prop2D = prop['dens_ion'][indz].reshape(GRID.Nodes[:2])/1e6 #Reshaping to Nx,Ny matrix
    prop1D = np.ma.array(np.ones(len(indphi)), mask=~indphi).reshape(GRID.Nodes[:2]) #Array masked where ind != indphi
    im0 = ax[0,0].imshow(prop2D, norm=colors.LogNorm(), cmap=cmap_im, extent = np.array([-sizex,sizex,-sizex,sizex])/sfu.au)
    ax[0,0].imshow(prop1D, cmap='binary_r', origin='lower', extent = np.array([-sizex,sizex,-sizex,sizex])/sfu.au) #Showing where the 1D profile comes from
    #ax[0,0].set_title('Ionized density on z=%.1f au'%(planez/sfu.au))
    ax[0,0].set_title('$n_{ion}$ [cm$^{-3}$] on z=0')
    ax[0,0].set_xlabel('x [au]')
    ax[0,0].set_ylabel('y [au]')
    cbar0 = plt.colorbar(im0, ax = ax[0,0])
    cbar0.ax.tick_params(which='major', direction='in', width=2.3, size=5, pad=7, labelsize=MEDIUM_SIZE)

    #vel x,y,z along 1D profile on ax[0,0]
    phi_line = GRID.rRTP[0][indz][indphi]/sfu.au
    ax[0,1].scatter(phi_line, prop['vel_x'][indz][indphi]/1e3, marker='p', edgecolor='tomato', label=r'$\upsilon_x$', **kwargs_sc_filled)
    ax[0,1].scatter(phi_line, prop['vel_y'][indz][indphi]/1e3, marker='o', edgecolor='forestgreen', label=r'$\upsilon_y$', **kwargs_sc_filled)
    ax[0,1].scatter(phi_line, prop['vel_z'][indz][indphi]/1e3, marker='X', edgecolor='dodgerblue', label=r'$\upsilon_z$', zorder=-1, **kwargs_sc_filled)
    ax[0,1].set_xlim(0, sizex/sfu.au)
    ax[0,1].set_xlabel('r [au]')
    ax[0,1].set_ylabel('Velocity [km s$^{-1}$]')
    ax[0,1].set_title(tag)
    ax[0,1].legend()

    xplot, yplot = phi_line, prop['dens_ion'][indz][indphi]*1e-6
    ax[1,0].scatter(xplot, yplot, marker='+', facecolor='k', **kwargs_sc_unfill) #prop along the 1D profile
    ax[1,0].plot(xplot, yplot, lw=4, color='k', alpha=0.2)
    ax[1,0].set_xlim(0, sizex/sfu.au)
    ax[1,0].set_xlabel('r [au]')
    ax[1,0].set_ylabel(r'$n_{ion}$ [cm$^{-3}$]')
    ax[1,0].yaxis.set_major_formatter(FormatStrFormatter('%.0e'))

    xplot, yplot = phi_line, density_obj.H[indz][indphi]/sfu.au
    ax[1,1].scatter(xplot, yplot, marker='+', facecolor='k', **kwargs_sc_unfill) #prop along the 1D profile
    ax[1,1].plot(xplot, yplot, lw=4, color='k', alpha=0.2)
    ax[1,1].set_xlim(0, sizex/sfu.au)
    ax[1,1].set_xlabel('r [au]')
    ax[1,1].set_ylabel('$H$ [au]')

    for axi in ax:
        for j,axj in enumerate(axi):
            make_up_ax(axj)
            axj.tick_params(labelbottom=True, labeltop=False)
            if j==1:
                axj.yaxis.set_label_position('right')
                axj.tick_params(labelright=True, labelleft=False)
            
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=0.3)
    plt.savefig('3Dvely_%s.png'%tag, bbox_inches='tight')
    plt.show()
    plt.close()

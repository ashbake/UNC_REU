#!/usr/bin/python
#execfile('plot_outliers_racz.py')

path = '/srv/two/ashbake/REU/eco_data/'
import numpy as np
from pylab import *
import math
from matplotlib import rc, rcParams
rcParams.update({'figure.autolayout': True,'figure.subplot.left'  : 0.15
                 ,'figure.subplot.bottom'  : 0.2})
import matplotlib.pyplot as plt
from numpy import *
import random
import matplotlib.ticker as mticker

Ho = 100. #units of h
execfile(path+'eco_wrapper.py')

#select & plot galaxies incorrectly grouped
bad = np.where(gal_dmass - mass_grp > 1)
plot(gal_dmass,mass_grp,'k.',label=bz + '_' + bp,markeredgecolor='b',markerfacecolor='b')
plot(gal_dmass[bad],mass_grp[bad],'r.',label=bz + '_' + bp,markeredgecolor='r',markerfacecolor='r')
savefig('plotbad_masses.png')

number        = 10    #min number of galaxies in a group
parameter     = 'meanssfr'    #key name (dextrmagnew, goodmhi, meanssfr, vrot_tf)

############## RA/CZ PLOTS ##############
howmanygroups = len(np.where(num[ind] > 10)[0])
clf()
mr = sn['dextrmagnew']
cz = sn['goodz']*3e5

#col = sn[parameter]                    # Magnitude
#xxx= sn[parameter] - sn['dextumagnew']
#col = xxx
#col = log10(sn['goodmhi'])  #GasMass

ssfr = sn[parameter]
ssfr[np.where(ssfr > 2)] = 2
col = ssfr                   #Spec Star Form Rate

#col = sn['vrot_tf']
j=0
save_col = array([[0.0]*howmanygroups]*4)        #for plotting mass vs mean group 'parameter'
for i in range(0,len(ind)-1):
    n=ind[i+1] - ind[i]
    if n > number:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        # SIZE based on MSTARS
        sizes = 50*(6**((sn['vrot'][ind[i]:ind[i+1]]-
                         min(sn['vrot'][ind[i]:ind[i+1]]))/
                        (max(sn['vrot'][ind[i]:ind[i+1]] - 
                             min(sn['vrot'][ind[i]:ind[i+1]])))))
        # COLOR based on R BAND MAGNITUDE
        #col = mr
        #col = 
        # GAS MASS or STELLAR MASS dominates?
        gasmass = sn['goodmstarsnew'][ind[i]:ind[i+1]]/log10(sn['goodmhi'][ind[i]:ind[i+1]])
        gas = np.where(gasmass < 1)
        star = np.where(gasmass > 1)

        # LOAD colormap & PLOT
        cm = plt.cm.get_cmap('RdYlBu')
        im = plt.scatter(ra[ind[i]:ind[i+1]],cz[ind[i]:ind[i+1]],
                         c=col[ind[i]:ind[i+1]],
                         cmap=cm,s=sizes,edgecolor='none',vmin = min(col),
                         label='Stellar Mass',vmax = max(col))
        cbar = plt.colorbar(im)
        cbar.set_label(parameter,size=12)
        #legend(loc=0,scatterpoints = 3)
        plt.scatter(ra[ind[i]:ind[i+1]][gas],cz[ind[i]:ind[i+1]][gas],marker='o',
                    s=sizes[gas]+35,edgecolor='k',facecolors='none') 

        # MARK BAD galaxies NOT in AM group
        bad = np.where(gal_dmass[ind[i]:ind[i+1]] - mass_grp[ind[i]:ind[i+1]] > 1)
        scatter(ra[ind[i]:ind[i+1]][bad],cz[ind[i]:ind[i+1]][bad],
                marker='x',c='k',lw=1,s=sizes[bad])
        
        # LABEL plot w/ mass and stuff
        #text(0.2,0.9,'log ' + r'$M_{halo}$' + ' = ' + str(round(gal_dmass[ind[i]],2)),ha='center', va='center',transform=ax.transAxes,fontsize=12)
        #text(0.2,0.85,'Skewness =',ha='center', va='center',transform=ax.transAxes,fontsize=12)
        font={'family' : 'Courier',
              'weight' : 'normal',
              'size'   :  20}

        suptitle(r'log $\mathsf{M_{halo}}$ = ' + str(round(gal_dmass[ind[i]],2)),fontsize=20)
        xlabel('Right Ascension (degrees)',fontsize=20)
        ylabel('CZ (km/s)',fontsize=20)

        #myLocator = mticker.MultipleLocator(4)
        #ax.xaxis.set_major_locator(myLocator)
        #gcf().subplots_adjust(left=0.15)
        plt.xticks(rotation=45)
        plt.matplotlib.rc('font', **font)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        ax.get_xaxis().get_major_formatter().set_scientific(False)

        savefig('racz_' +str(i) + '.png')
        clf()
        
        save_col[0,j] = gal_dmass[ind[i]]
        save_col[1,j] = median(col[ind[i]:ind[i+1]])
        save_col[2,j] = max(col[ind[i]:ind[i+1]])
        save_col[3,j] = min(col[ind[i]:ind[i+1]])
        j+=1
        
        
xlabel('RA')
ylabel('CZ')
plt.savefig('plotbad_racz.png')

clf()
plot(save_col[0,:],save_col[1,:],'k.')#,save_col[0,:],save_col[2,:],'b.',save_col[0,:],save_col[3,:],'r.')
xlabel('Dynamical Mass of Group')
ylabel(parameter)
savefig('test.png')


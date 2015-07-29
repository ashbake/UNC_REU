#!/usr/bin/python
#execfile('plot_outliers.py')

path = '/srv/two/ashbake/REU/eco_data/'
import numpy as np
from pylab import *
import math
from matplotlib import rc, rcParams
import matplotlib.pyplot as plt
from numpy import *
import random



Ho = 100. #units of h
execfile(path+'eco_wrapper.py')

#select galaxies incorrectly grouped
bad = np.where(gal_dmass - mass_grp > 1)

plot(gal_dmass,mass_grp,'k.',label=bz + '_' + bp,markeredgecolor='b',markerfacecolor='b')
plot(gal_dmass[bad],mass_grp[bad],'r.',label=bz + '_' + bp,markeredgecolor='r',markerfacecolor='r')
savefig('plotbad_masses.png')


############## RA/DEC PLOTS ##############
mr = sn['dextrmagnew']       # Magnitude
#col = log10(sn['goodmhi'])  #GasMass
ssfr = sn['meanssfr']
ssfr[np.where(ssfr > 4)] = 4
col = ssfr                   #Spec Star Form Rate

clf()
for i in range(0,len(ind)-1):
    n=ind[i+1] - ind[i]
    if n > 20:
        fig = plt.figure()
        ax = fig.add_subplot(111)

        # SIZE based on MSTARS
        sizes = 20*(5**((sn['goodmstarsnew'][ind[i]:ind[i+1]]-min(sn['goodmstarsnew'][ind[i]:ind[i+1]]))/(max(sn['goodmstars'][ind[i]:ind[i+1]] - min(sn['goodmstars'][ind[i]:ind[i+1]])))))
        # COLOR based on R BAND MAGNITUDE
        #colors = (mr-max(mr))/(min(mr)-max(mr))
        # GAS MASS or STELLAR MASS dominates?
        gasmass = sn['goodmstarsnew'][ind[i]:ind[i+1]]/log10(sn['goodmhi'][ind[i]:ind[i+1]])
        gas = np.where(gasmass < 1)
        star = np.where(gasmass > 1)

        # LOAD colormap & PLOT
        cm = plt.cm.get_cmap('RdYlBu')
        im = plt.scatter(ra[ind[i]:ind[i+1]],dec[ind[i]:ind[i+1]],c=col[ind[i]:ind[i+1]],cmap=cm,s=sizes,edgecolor='none',vmin = min(col), vmax = max(col))
        cbar = plt.colorbar(im)
        cbar.set_label(parameter,size=18)
        cm = plt.cm.get_cmap('winter')
        plt.scatter(ra[ind[i]:ind[i+1]][gas],dec[ind[i]:ind[i+1]][gas],marker='o',s=sizes[gas]+35,edgecolor='k',facecolors='none') 

        # MARK BAD galaxies NOT in AM group
        bad = np.where(gal_dmass[ind[i]:ind[i+1]] - mass_grp[ind[i]:ind[i+1]] > 1)
        scatter(ra[ind[i]:ind[i+1]][bad],dec[ind[i]:ind[i+1]][bad],marker='x',c='k',lw=1.2,s=sizes[bad])
        
        # LABEL plot w/ mass and stuff
        text(0.2,0.9,'log /M_{halo} = ' + str(round(gal_dmass[ind[i]],2)),ha='center', va='center',transform=ax.transAxes,fontsize=12)
        text(0.2,0.85,'HEYY sk ku =',ha='center', va='center',transform=ax.transAxes,fontsize=12)
        
        xlabel('Right Ascension')
        ylabel('Declination')
        plt.savefig('radec_' +str(i) + '.png')
        clf()

xlabel('RA')
ylabel('DEC')
plt.savefig('plotbad_radec.png')


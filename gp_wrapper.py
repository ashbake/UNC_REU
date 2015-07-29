#!/usr/bin/python
#execfile('gp_wrapper.py')
#run this before all other gp files

path = '/srv/two/ashbake/REU/'
mm=0
pp=3  # 3 and 5 
method = ['avg','cen'][mm]
llens = ['new1','new','old','12_13','11_14','14_75','12_14','12_17'][pp]
shift = [0,0,0,0.7797,0.8252,0.8611,0.7636,0.704]
Ho = 100.0 #h
print 'method = ', method, llens
import numpy as np
from pylab import *
import math 
from matplotlib import rc, rcParams
import matplotlib.pyplot as plt
from numpy import *
import cPickle
from numpy import arange,array,ones,linalg
from pylab import plot,show
from numpy import arange,array,ones#,random,linalg
import pylab
pylab.ion()

#import scipy as scipy
#from scipy import stats

filename = path + 'gp_dynm_out_' + method + '_' + llens + '.dat' #avg method & good linking lengths
file=np.loadtxt(filename)

fof_dmass = log10((1.27*3*3.14/(2*.7))*file[:,0])# + shift[pp] 
fof_vel	  = file[:,1]
fof_rad   = file[:,2]
num       = file[:,3]
inew      = file[:,4]
nold      = file[:,5]
iold      = file[:,6]
rmag      = file[:,7]
halom     = file[:,8] + 0.155
ofind	  = file[:,9]
ra     	  = file[:,10]
dec    	  = file[:,11]
dmass	  = log10((1.27*3*3.14/(2*.7))*(file[:,12])) #+ shift[pp]
vel       = file[:,13] #from dynmass.py
rad       = file[:,14]
galr      = file[:,15] #each galaxies position/vel in halo
galv      = file[:,16]
gd_galr   = file[:,17]
gd_galv   = file[:,18]
cz        = file[:,19]

#rho0      = 1.477 * 10**11 #h^2 Msun/Mpc^3
rho0      = .75*2.775*10**11 #h^2 Msun/Mpc^3 from Andreas's code

r200       = ((10**(fof_dmass))*(3./(rho0*200*4.*np.pi)))**(1./3.)
r200gd      = ((10**(dmass))*(3./(rho0*200*4.*np.pi)))**(1./3.)



#match_mass = array(cPickle.load(open(path + 'match_mass_'+ method + '_' + llens + '.p','rb')))
match_mass = fof_dmass

######### new index ###########
ind=[0]    #make index
i=1
while i < (len(inew)-1):
  	tempind=i
  	ind.extend([tempind])
  	while (i < (len(inew)-2)) & (inew[i+1] == inew[i]):
 	    i = i + 1
 	i=i + 1

ind.extend([tempind+1]) #boundary condition

#HAM MASSES
filename = path + 'AM_code/out_Resolve_5001_hod1_mock_zspace.dat'
file=np.loadtxt(filename)
hamass = file[:,1]
hamr200 = file[:,2]
gal_hamass = [0] * (len(ind) + 1)
for i in range(0,len(ind)-1):
    n = ind[i+1] - ind[i]
    gal_hamass[ind[i]:ind[i+1]] = [hamass[i]] * n

gal_hamass.extend([hamass[-1]])
gal_hamass = array(gal_hamass) + .155


############## old index #############
filename = path + 'dynm_out_'+ method + '.dat' #avg method & good linking lengths
file=np.loadtxt(filename)
czold = file[:,6]
oldid = file[:,7]
oldnum = file[:,8]
gd_galr   = file[:,4]
gd_galv   = file[:,5]
indold=[0]
i=1
while i < (len(czold)-1):
    tempind = i
    indold.extend([tempind])
    while (oldid[i+1] == oldid[i]):
        i = i + 1
    i = i + 1
indold.extend([tempind+1]) #boundary condition

class BinFun:
    'given x and y, bin along x and find dispersion in y'
    def __init__(self,xin,yin,nbin):
        binsize   = (max(xin) - min(xin))/nbin
        xout      = [0] * nbin      #store middle of bin
        yout      = [0] * nbin      #store avg halo_m for bin
        errs      = [0] * nbin
        for pizza in range(0,nbin):
            temp_y       =  yin[(xin > (binsize*pizza + min(xin))) & (xin < (binsize*(pizza + 1)+min(xin)))] 
            errs[pizza]  =  sqrt(sum(abs(temp_y - median(temp_y))**2)/(len(temp_y)-1))
            xout[pizza]  =  (binsize*pizza + binsize*(pizza+1))/2 + min(xin)
            yout[pizza]  =  median(temp_y)
        self.xout = xout
        self.yout = yout
        self.errs = errs

class MAD:
    'Median absolute deviation'
    def __init__(self,xin,yin,nbin):
        binsize   = (max(xin) - min(xin))/nbin
        xout      = [0] * nbin      #store middle of bin
        yout      = [0] * nbin      #store avg halo_m for bin
        errs      = [0] * nbin
        for pizza in range(0,nbin):
            temp_y       =  yin[(xin > (binsize*pizza + min(xin))) & (xin < (binsize*(pizza + 1)+min(xin)))] 
            errs[pizza]  =  median(abs(temp_y - median(temp_y)))
            xout[pizza]  =  (binsize*pizza + binsize*(pizza+1))/2 + min(xin)
            yout[pizza]  =  median(temp_y)
        self.xout = xout
        self.yout = yout
        self.errs = errs

#!/usr/bin/python
#execfile('eco_makepickle.py')

path = '/srv/two/ashbake/REU/eco_data/'
import numpy as np
from pylab import *
import math
from matplotlib import rc, rcParams
import matplotlib.pyplot as plt
from numpy import *
import random
from scipy.io.idl import readsav
import scipy as scipy
import cPickle as pickle
import pylab
pylab.ion()

bzs = ['.12','.13','.14','.12','.11']
bps = ['2.2','1.5','.75','1.3','1.4']
shifts = [0,0,0.8611,.7797,0.8252]
############
llpick = 3 #2 and 3 only options
############

bz = bzs[llpick]
bp = bps[llpick]
shift = shifts[llpick]

#LOAD MY DICTIONARY & add things
eco_dict = pickle.load( open( "eco_dict" + bz + bp+ ".p", "rb" ) )
locals().update(eco_dict) #restore all variables
rho0       = .7*.7*.75*2.775* 10**11 #h^2 Msun/Mpc^3
eco_dict['r200']       = ((10.**shift)*eco_dict['gal_dmass']*
                          (3./(rho0*200*4.*np.pi)))**(1./3.)
eco_dict['rvir']       = 1.27 * gal_rad       #1.27 from plotting Mhalo & Mdyn
eco_dict['gal_dmass'] = log10((1.27*3*np.pi/2)*gal_dmass)  #redefine w/ shift

## ADD HALO ABUNDANCE MATCHING MASSES & R200
filename = path + 'out_ecospring_hammass.txt'
file = np.loadtxt(filename)

hamass = file[:,1]
hamr200 = file[:,2]

gal_hamass = [0]   #define per galaxy array
gal_hamr200 = [0]
for i in range(0,len(ind)-1):
    n = ind[i+1] - ind[i]
    gal_hamass[ind[i]:ind[i+1]] = [hamass[i]] * n
    gal_hamr200[ind[i]:ind[i+1]] = [hamass[i]] * n
gal_hamass.extend([hamass[-1]])  #gotta add the last value
gal_hamass = array(gal_hamass) + 0.15 #convert by subtracting log.7 for little h
gal_hamr200.extend([hamr200[-1]])  #gotta add the last value

eco_dict['gal_hamass'] = gal_hamass
eco_dict['hamass'] = hamass
eco_dict['hamr200'] = hamr200
eco_dict['gal_hamr200'] = hamr200


#LOAD AMANDA'S STUFF
path2 = '/srv/two/ashbake/SENIOR/RESOLVE_data/'
s = readsav(path2 + 'jushr_alfalfa_clean_102113.dat',python_dict=True)
xx = (np.where(s['goodnewabsr'] < -17.33))[0]
sn = {}
for k in s.keys():
    eco_dict[k] = s[k][xx][list(of)]


pickle.dump(eco_dict, open( 'eco_dict_all' + bz + bp+ '.p', "wb" ) )

#######   Pick In Sample #######
#xx = (np.where((s['dextrmagnew'] < -17.33) & (s['groupcentcz'] > 3000.) & (s['groupcentcz'] < 7000.) & (radiir > 1) & (radiid > 1)))[0]

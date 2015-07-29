#!/usr/bin/python
#execfile('eco_wrapper.py')

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
h = 0.7
Ho=100*h

eco_dict = pickle.load( open( path + "eco_dict_all" + bz + bp+ ".p", "rb" ) )

##########   FOR APPLYING BOUNDARY CUTS / TO DEFINE SAMPLE  ############

#inbound = (np.where((array(eco_dict_all['avg_cz']) > 3000.) & (array(eco_dict_al#l['avg_cz']) < 7000.)))[0]
#indall = eco_dict_all['ind']
#inbound2 = (np.where((array(eco_dict_all['avg_cz'])[indall] > 3000.) & (array(ec#o_dict_all['avg_cz'])[indall] < 7000.)))[0]
#
#eco_dict = {}
#for k in eco_dict_all.keys():
#    if len(eco_dict_all[k]) > 10000:
#        eco_dict[k] = array(eco_dict_all[k])
#    else:
#        eco_dict[k] = array(eco_dict_all[k])[inbound2]
#        print k, 'should be only key to change'


locals().update(eco_dict) #restore all variables




class BinFun:
    'given x and y, bin along x and find dispersion in y'
    def __init__(self,xin,yin,nbin):
        binsize   = (max(xin) - min(xin))/nbin
        xout      = [0] * nbin      #store middle of bin
        yout      = [0] * nbin      #store avg halo_m for bin
        errs      = [0] * nbin
        for pizza in range(0,nbin):
            temp_y       =  array(yin)[(xin > (binsize*pizza + min(xin))) & (xin < (binsize*(pizza + 1)+min(xin)))] 
            errs[pizza]  =  sqrt(sum(abs(temp_y - median(temp_y))**2)/(len(temp_y)-1))
            xout[pizza]  =  (binsize*pizza + binsize*(pizza+1))/2 + min(xin)
            yout[pizza]  =  median(temp_y)
        self.xout = xout
        self.yout = yout
        self.errs = errs

class MAD:
    'given x and y, bin along x and find dispersion in y'
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


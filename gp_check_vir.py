#!/usr/bin/python
#execfile('gp_check_vir.py')

path = '/srv/two/ashbake/REU/'
import numpy as np
from pylab import *
import math
from matplotlib import rc, rcParams
import matplotlib.pyplot as plt
from numpy import *
import random
import matplotlib.ticker as mticker
from scipy.special import erf
from scipy.stats import anderson

#llpick=3
#execfile(path+'gp_wrapper.py')

number        = 5            #min number of galaxies in a group


howmanygroups = len(np.where(num[ind] > number)[0])
virline = 0 #0 if defining virialized line to be rvir, neg if want rvir bigger

A2 = [0]*len(cz)
alpha = [0]*len(cz)
count = 0
for i in range(0,len(ind)-1):
    n=ind[i+1] - ind[i]
    if n > number:
        czs_grp = cz[ind[i]:ind[i]+n]
        sort = np.argsort(czs_grp)
        czs_grp = czs_grp[sort]
        
        andtest = anderson(czs_grp,dist='norm')
        A2notstar = andtest[0]
        A2[ind[i]:ind[i+1]] = [A2notstar*(1 + 0.75/n + 2.25/(n**2))]*n
        
        a = 3.6789468
        b = 0.1749916  #from Hou et al 2009
        alpha[ind[i]:ind[i+1]] = [a * exp(-A2[ind[i]]/b)]*n  #alpha < 5% not gaussian
        
        count += 1

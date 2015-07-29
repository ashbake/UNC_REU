#!/usr/bin/python
#execfile('makefile.py')

path = '/srv/two/ashbake/REU/'
import numpy as np
from pylab import *
import math
from matplotlib import rc, rcParams
import matplotlib.pyplot as plt
from numpy import *
import random
import matplotlib.ticker as mticker


filename = path + 'AM_code/Resolve_5001_hod1_mock_zspace.dat'
file=np.loadtxt(filename)
ra = file[:,0]
dec = file[:,1]
cz = file[:,2]
id = file[:,4]
mr = file[:,3]
nu = file[:,6]

ind=[0]
i=1
while i < (len(id)-1):
    tempind = i
    ind.extend([tempind])
    while (id[i+1] == id[i]):
        i = i + 1
    i = i + 1

ind.extend([tempind+1]) #boundary condition

mrtot=[0.0]*len(id)
num = [0]*len(id)

for i in range(0,len(ind)-1):
    if (ind[i+1]-ind[i]) != nu[ind[i]]: #remove halos w/ missing centrals
        num[ind[i]:ind[i+1]] = [ind[i+1] - ind[i]]*(ind[i+1] - ind[i])
        mrtot[ind[i]:ind[i+1]] = [-99]*(ind[i+1] - ind[i])
    else:
        ntemp = ind[i+1] - ind[i]
        num[ind[i]:ind[i+1]] = [ntemp]*ntemp
        mrtot[ind[i]:ind[i+1]] = [-2.5*log10(sum(10**(-0.4*mr[ind[i]:ind[i+1]])))] * ntemp
 

mrtot[-1] = mr[-1]
num[-1]=1

num = array(num)
mrtot=array(mrtot)
Gs = array(['G']*len(num))

filename2 = 'G_Resolve_5001_hod1_mock_zspace.txt'
file2=np.loadtxt(filename2,dtype = str)
newid = arange(1,len(id[ind])+1)

X = [Gs[ind],newid,ra[ind],dec[ind],cz[ind],num[ind],cz[ind],cz[ind],mrtot[ind]]

np.savetxt('G_goodResolvegroups2.txt',np.c_[ra[ind],dec[ind],
                                                cz[ind],num[ind],cz[ind],
                                                cz[ind],mrtot[ind]],fmt='%1.3f')

np.savetxt('numberstemp.txt',np.c_[num[ind]],fmt='%1i')
np.savetxt('mrtottemp.txt',np.c_[mrtot[ind]],fmt='%1.4f')

##################load
filename3 = 'out_goodResolvegroups2_mass.txt'
filename2='Resolve_5001_hod1_mock_zspace.dat'

file3=np.loadtxt(filename3)
newm  = file3[:,1]

file2 = np.loadtxt(filename2)
truem = file2[:,5]

plot(newm,truem[ind],'r.')

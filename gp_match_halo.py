#!/usr/bin/python
#gp_match_halo.py
#execfile('gp_match_halo.py')
#Use gp_dynmass.py to get dynmass, vel, sig, etc. pick best halo mass to
#associate with the new groups the group finder found

import numpy as np
from pylab import *
import math
from matplotlib import rc, rcParams
import matplotlib.pyplot as plt
from numpy import *
from scipy.stats import mode
import cPickle

filename="./gp_dynm_out_avg_new1.dat"
file=np.loadtxt(filename)

fof_dmass = log10(file[:,0])
fof_vel	  = file[:,1]
fof_rad   = file[:,2]
num       = file[:,3]
inew      = file[:,4]
nold      = file[:,5]
iold      = file[:,6]
rmag      = file[:,7]
halom     = file[:,8]
ofind	  = file[:,9]
ra     	  = file[:,10]
dec    	  = file[:,11]
dmass	  = log10(file[:,12])
vel       = file[:,13]
rad       = file[:,14]

ind=[0]    #initial condition
i=1
while i < (len(inew)-1):
  	tempind=i
  	ind.extend([tempind])
  	while (i < (len(inew)-2)) & (inew[i+1] == inew[i]):
 	    i = i + 1
 	i=i+1
 	
ind.extend([tempind+1]) #boundary condition

numbers = np.arange(0, max(iold), 1)

class GetMode:
	'for a set of galaxies in one group, gives the max recurring old index'
	'select first based on max contribution from original halo, then on max luminosity'
	def __init__(self,i):
		mode = {}
		flag = {}
	        alist = iold[ind[i]:ind[i+1]]
		tempbins = np.bincount(list(alist))
		if size(tempbins[tempbins == max(tempbins)]) == 1:
			mode = tempbins.argmax()
			flag = 0
		else:
			#select based on luminosity
			maxes = numbers[[tempbins == max(tempbins)]]
			lums  = [0] * len(maxes)
			for j in range(0,len(maxes)):
				lums[j] = sum(rmag[iold == maxes[j]])
			mode=array(maxes)[lums == min(lums)][0]
			flag = 1
		self.mode = mode
		self.flag = flag
		
#make index showing if mode was unique or not so can color code plot based on that

flag = []
match_mass = []
for i in range(0,len(ind)-1):
	n = ind[i+1] - ind[i]
	if ind[i+1] - ind[i] < 3:
		match_mass[ind[i]:ind[i+1]] = [0] * n #[halom[ind[i]]]
		flag[ind[i]:ind[i+1]] = [0] * n
	#elif mode == 0:
	#	match_mass[ind[i]:ind[i+1]] = [0] * n
	else:
		imode = GetMode(i)
		mode = imode.mode
		match_mass[ind[i]:ind[i+1]] = [halom[iold == mode][0]] * n
		flag[ind[i]:ind[i+1]] = [imode.flag] * n
		print i

match_mass.append(0)
flag.append(0)


filename='match_mass_avg_new1.p'
file=open(filename,'w')
cPickle.dump(match_mass, file, protocol=2)
file.close()

match_mass = array(match_mass)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#	#	#	#	#	#	PLOT	#	#	#	#	#	#	 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
plot(halom,match_mass,'k.',label="Matched Mass vs. Galaxies' originial Halo Mass") #should be 1:1
plot(fof_dmass,match_mass,'y.',label='Matched Halo Mass vs calculated Group Dyn Mass',zorder = -100)
plot(fof_dmass[num > 5],match_mass[num > 5],'r.',label='num galaxies > 5',zorder = -100)
out = BinFun(fof_dmass[fof_dmass > 0],match_mass[fof_dmass > 0],10)
plt.errorbar(out.xout,out.yout,yerr=out.errs,barsabove=True,ls='None', capsize=4.5, color='black',zorder = 300,elinewidth=2.5)
plot(out.xout,out.yout,'ko')
ylim(10,16)
xlim(8,15)
xlabel('Dynamical Mass post fof9')
ylabel('Matched Halo Mass')
title('Halo Matching using mode halo ID & max lum for tie breaking')
legend(numpoints=1)
savefig('match_mass1.png')
show()

out = BinFun(fof_dmass[fof_dmass > 0],match_mass[fof_dmass > 0],10)
plt.errorbar(out.xout,out.yout,yerr=out.errs,barsabove=True,ls='None', capsize=4.5, color='black',zorder = 300,elinewidth=2.5)
plot(out.xout,out.yout,'ko',label='match_mass')

out = BinFun(fof_dmass[fof_dmass > 0],halom[fof_dmass > 0],10)
plt.errorbar(out.xout,out.yout,yerr=out.errs,barsabove=True,ls='None', capsize=4.5, color='blue',zorder = 300,elinewidth=2.5)
plot(out.xout,out.yout,'bo',label='after fof')

out = BinFun(dmass[dmass > 0],halom[dmass > 0],10)
plt.errorbar(out.xout,out.yout,yerr=out.errs,barsabove=True,ls='None', capsize=4.5, color='red',zorder = 300,elinewidth=2.5)
plot(out.xout,out.yout,'ro',label='true dmass')
legend(numpoints=1,loc=2)
savefig('test.png')
show()


plot(fof_dmass,match_mass,'r.', label='unique mode')
plot(fof_dmass[flag == 1],match_mass[flag == 1],'k.',label='tie breakers')
ylabel('matched mass')
xlabel('dynamical mass')
savefig('match_mass3.png')
legend(numpoints=1)
show()

#!/usr/bin/python
#execfile('eco_check_gaus.py')
#this program determines how gaussian a group's vel dispersion and then checks
#to see if it correlates w/ halo mass dispersion
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata

execfile('eco_wrapper.py') #choose which linking lengths

class SkewKurtReal():
	'find the skewness and kurtosis of real velocity data'
	def __init__(self,galv,ind,num): 
		skew = [0] * len(ind)
		kurt = [0] * len(ind)
		suck = []
		for i in range(0,len(ind)-1):
			n = ind[i+1] - ind[i]
			if n <=  3:
				skew[ind[i]:ind[i+1]] = [0] * n
				kurt[ind[i]:ind[i+1]] = [0] * n
				suck[ind[i]:ind[i+1]] = [0] * n
			elif n != num[ind[i]]:
				skew[ind[i]:ind[i+1]] = [0] * n
				kurt[ind[i]:ind[i+1]] = [0] * n
				suck[ind[i]:ind[i+1]] = [0] * n
			elif n > 3:
				subv = array(galv[ind[i]:ind[i+1]])
				m3  =  (1./n) * sum((subv - subv.mean())**3.)
				m2  =  ((1./n) * sum((subv - subv.mean())**2.))**(3./2)
				skew[ind[i]:ind[i+1]]  =  [(sqrt(n*(n-1.))/(n-2.))*m3/m2] * n
				g1  =  (1./n) * sum((subv - subv.mean())**4.)
				g2  =  ((1./n) * sum((subv - subv.mean())**2.))**2.
				k   =  g1/g2 - 3
				kurt[ind[i]:ind[i+1]] = [((n-1.)/((n-2.)*(n-3.))) * ((n+1.) * k + 6.)] * n
				#suck[ind[i]:ind[i+1]] = [mean(abs(fof_dmass[ind[i]:ind[i+1]] - halom[ind[i]:ind[i+1]]))] * n
		kurt.append(0)
		skew.append(0)
		suck.append(0)
		self.skew = skew
		self.kurt = kurt
		self.suck = suck

import random

class SkewKurtFake():
	'CHECK RANDOM SKEWNESS & KURTOSIS from random gaussian distribution'
       	def __init__(self,num): 
		size = 500
		skew = [0] * size
		kurt = [0] * size
		numb = [0] * size
		for it in range(0,size):
			#n = random.randint(1,200)
			n = random.sample(num[num > 2],1)[0]
			#velocities = array(range(size)) + 500.  #flat distribution
			velocities = np.random.normal(random.randint(472,523),random.uniform(2,400),10000)  #gaussian
			cz = random.sample(velocities,int(n))
			if n <= 3:
				skew[it] = 0
				kurt[it] = 0
				numb[it] = float(n)
			if n > 3:
				subv = array(cz)
				m3  =  (1./n) * sum((subv - subv.mean())**3.)
				m2  =  ((1./n) * sum((subv - subv.mean())**2.))**(3./2)
				skew[it]  =  (sqrt(n*(n-1.))/(n-2.))*m3/m2
				g1  =  (1./n) * sum((subv - subv.mean())**4.)
				g2  =  ((1./n) * sum((subv - subv.mean())**2.))**2.
				k   =  g1/g2 - 3
				kurt[it] = ((n-1.)/((n-2.)*(n-3.))) * ((n+1.) * k + 6.)
				numb[it] = n
		self.skew = skew
		self.kurt = kurt
		self.numb = numb

real = SkewKurtReal(galv,ind,num)
fake = SkewKurtFake(num)
cols=num
cols[cols > 30] = 30

filename = '../gp_dynm_out_avg_12_13.dat' #avg method & good linking lengths
file=np.loadtxt(filename)
inew      = file[:,4]
galv      = file[:,16]
num       = file[:,3]
ind=[0]    #make index
i=1
while i < (len(inew)-1):
  	tempind=i
  	ind.extend([tempind])
  	while (i < (len(inew)-2)) & (inew[i+1] == inew[i]):
 	    i = i + 1
 	i=i + 1

ind.extend([tempind+1]) #boundary condition

real2= SkewKurtReal(galv,ind,num)

cols2=num
cols2[cols2>30]=30

fig=plt.figure()
#plot in one window
plt.subplot(2,3,1)
scatter(real.skew,real.kurt,marker='.',c=cols,edgecolors='none')
xlabel('skewness')
ylabel('kurtosis')
ylim(-6,6)
xlim(-2,2)
title('ECO data')
colorbar()

plt.subplot(2,3,2)
scatter(fake.skew,fake.kurt,marker='.',c=fake.numb,edgecolors='none')
xlabel('skewness')
ylabel('kurtosis')
title('From Gaussian Sample')
ylim(-6,6)
xlim(-2,2)
colorbar()

plt.subplot(2,3,3)
scatter(real2.skew,real2.kurt,marker='.',c=cols2,edgecolors='none')
xlabel('skewness')
ylabel('kurtosis')
ylim(-6,6)
xlim(-2,2)
title('Mock Data post FOF')
colorbar()

#contours
plt.subplot(2,3,4)
ngrid = 30
xi = np.linspace(min(real.skew),max(real.skew),ngrid)
yi = np.linspace(min(real.kurt),max(real.kurt),ngrid)
zi = griddata(real.skew,real.kurt,cols,xi,yi)
plt.contour(xi,yi,zi,20,linewidths=1)
#plt.scatter(real.skew,real.kurt,c=num,s=20)
plt.xlim(min(real.skew),max(real.skew))
plt.ylim(min(real.kurt),max(real.kurt))

plt.subplot(2,3,5)
ngrid = 30
xi = np.linspace(min(fake.skew),max(fake.skew),ngrid)
yi = np.linspace(min(fake.kurt),max(fake.kurt),ngrid)
zi = griddata(fake.skew,fake.kurt,fake.numb,xi,yi)
plt.contour(xi,yi,zi,20,linewidths=1)
#plt.scatter(fake.skew,fake.kurt,c=fake.numb,s=20)
plt.xlim(min(fake.skew),max(fake.skew))
plt.ylim(min(fake.kurt),max(fake.kurt))

plt.subplot(2,3,6)
ngrid = 20
xi = np.linspace(min(real2.skew),max(real2.skew),ngrid)
yi = np.linspace(min(real2.kurt),max(real2.kurt),ngrid)
zi = griddata(real2.skew,real2.kurt,cols2,xi,yi)
plt.contour(xi,yi,zi,20,linewidths=1)
#plt.scatter(real.skew,real.kurt,c=num,s=20)
plt.xlim(min(real2.skew),max(real2.skew))
plt.ylim(min(real2.kurt),max(real2.kurt))

savefig('figures/plotgaus.png')
plt.show()

clf()

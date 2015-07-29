#!/usr/bin/python
#execfile('dynmass.py')


#### CHOOSE MASS ESTIMATOR ####
mass_est = 'pow(velarr,2) * avg_rad / G'
#mass_est = 'pow(velarr,2)  / G'
#mass_est = '(pow(velarr,2) / avg_rad) / G' #(used with angles3)

option = 'avg'
          #'centrals' or 'averages'
inverse = 'no' 
          #only pick yes if using 3rd mass_est where need it


import numpy as np
from pylab import *
import math
from matplotlib import rc, rcParams
import matplotlib.pyplot as plt
from numpy import *
import pylab
pylab.ion()

filename="./Resolve_5001_hod1_mock_zspace.dat"
file=np.loadtxt(filename)

ra = file[:,0]
dec= file[:,1]
cz = file[:,2] #real cz, no peculiars for rspace
mr = file[:,3] #absolute R band mag
id = file[:,4] #halo id number
lm = file[:,5] #log halo mass
nu = file[:,6] #number in halo
ch = file[:,7] #central 0 satellite 1

G  = 4.302E-9   #Mpc(km/s)^2 Msun^-1
c  = 3.0*10**5  #km/sec
Ho = 100.0      #+/-0.80 (km/s)/Mpc  was using 69.3 but use 100 and then mass will be in terms of inverse little h

### DEFINE INDEX MATRIX ###
N   = (len(cz)-1)
ch = ch[0:N+1]
ind2 = [0]
xx  = [i+1 for i in xrange(N)]
ind2.extend(xx)
ind2 = array(ind2)
ind2 = ind2[ch == 0]                     #define index as start of new group. sometimes a group doesnt have a central...? causing scatter
#fixing missing central problem (redefine better 'ind')
ind=[0]
i=1
while i < (len(cz)-1):
    tempind = i
    ind.extend([tempind])
    while (id[i+1] == id[i]):
        i = i + 1
    i = i + 1
ind.extend([tempind+1]) #boundary condition

### DEFINE VEL DISPERSION && AVG RADIUS ###
velarr   = [0] * len(ind)
avg_rad  = [0] * len(ind)
dyn_mass = [0] * len(ind)
galr     = [0] * len(ind)
galv     = [0] * len(ind)
for i in range(0,len(ind)-1):
    if (ind[i+1]-ind[i]) != nu[ind[i]]: #remove halos w/ missing centrals
        dyn_mass[i] = 0
        galr[ind[i]:ind[i+1]] = [0]*(ind[i+1] - ind[i])
        galv[ind[i]:ind[i+1]] = [0]*(ind[i+1] - ind[i])
    elif ind[i+1] == ind[i] + 1:     #rid of halos w/ 1 galaxy
        dyn_mass[i] = 0.
        galr[ind[i]:ind[i+1]] = [0]
        galv[ind[i]:ind[i+1]] = [0]
    elif ind[i+1] == ind[i] + 2:     #rid of halos w/ 2 galaxies
        dyn_mass[i] = 0.
        galr[ind[i]:ind[i+1]] = [0]*2
        galv[ind[i]:ind[i+1]] = [0]*2
    else:
        if option == 'cen':
            #vel dispersion
            czs  = cz[ind[i]:ind[i+1]]
            num  = (ind[i+1] - ind[i])
            ceni = [mr[ind[i]:ind[i+1]] == min(mr[ind[i]:ind[i+1]])]
            if len(czs[ceni]) > 1:
                mrs = mr[ind[i]:ind[i+1]] + .001*array([random.randint(1,100) for _ in range(num)])
                ceni = [mrs == min(mrs)]
            sig  = sqrt(sum((czs - czs[ceni])**2)/(num - 1))
            galv[ind[i]:ind[i+1]] = czs - median(czs)
            velarr[i] = sig
            #proj radii
            ras     = ra[ind[i]:ind[i+1]]
            decs    = dec[ind[i]:ind[i+1]]
            angles  = sqrt(((ras - ras[ceni])*cos((math.pi/180.)*decs[ceni]))**2 + (decs - decs[ceni])**2)
            radii   = (math.pi/180.)* angles  * czs[ceni] / Ho
            galr[ind[i]:ind[i+1]] = radii
            if inverse == 'no':
                avg_rad[i] = median(radii)
            elif inverse == 'yes':
                dists=array([[0.]*len(ras)]*len(ras))
                for ii in range(0,len(ras)):
                    for j in range(0,len(ras)):
                        if ii != j:
                            dists[ii,j] = (sqrt(((ras[ii] - ras[j])*cos((math.pi/180.)*decs[j]))**2 + (decs[ii] - decs[j])**2))
                radii  = (math.pi/180.)* dists * avg_cz / Ho
                avg_rad[i] = sum(1/radii[radii > 0])/(size(radii)-len(ras))  
            else:
                print 'do you want the inverse of r or not??'
        elif option == 'avg':
            #vel dispersion
            czs = cz[ind[i]:ind[i+1]]
            num = ind[i+1] - ind[i]
            avg_cz  = mean(czs) #sum(czs)/num
            sig = sqrt(sum((czs - avg_cz)**2)/(num - 1))
            velarr[i] = sig
            galv[ind[i]:ind[i+1]] = czs - mean(czs)
            #proj radii
            ras     = ra[ind[i]:ind[i+1]]
            decs    = dec[ind[i]:ind[i+1]]
            avg_ra  = sum(ras)/num
            avg_dec = sum(decs)/num
            angles  = sqrt(((ras - avg_ra)*cos((math.pi/180.)*avg_dec))**2 + (decs - avg_dec)**2)
            radii   = (math.pi/180.)* angles * avg_cz / Ho
            galr[ind[i]:ind[i+1]] = radii
            if inverse == 'no':
                avg_rad[i] = mean(radii)
            elif inverse == 'yes':
                dists=array([[0.]*len(ras)]*len(ras))
                for ii in range(0,len(ras)):
                    for j in range(0,len(ras)):
                        if ii != j:
                            dists[ii,j] = (sqrt(((ras[ii] - ras[j])*cos((math.pi/180.)*decs[j]))**2 + (decs[ii] - decs[j])**2))
                radii  = (math.pi/180.)* dists * avg_cz / Ho
                avg_rad[i] = sum(1/radii[radii > 0])/(size(radii)-len(ras))
            else:
                print 'do you want the inverse of r or not??'

velarr   = array(velarr)
radarr  = array(avg_rad)
dyn_mass = eval(mass_est)    #choose mass estimator above
halo_m   = lm[ind]           #actual mass of halo
galr.append(0)
galv.append(0)


##~~~~~~~~~define per galaxy arrays ~~~~~~~~~~~##
gal_dmass = [0] #* (len(ind) + 1)
gal_vel   = [0] #* (len(ind) + 1)
gal_rad   = [0] #* (len(ind) + 1)
numgood = [0] #* (len(ind) + 1)
for i in range(0,len(ind)-1):
    n = ind[i+1] - ind[i]
    n2 = nu[ind[i]]
    numgood[ind[i]:ind[i+1]] = [n] * n2
    gal_dmass[ind[i]:ind[i+1]] = [dyn_mass[i]] * n2
    gal_vel[ind[i]:ind[i+1]] = [velarr[i]] * n2
    gal_rad[ind[i]:ind[i+1]] = [radarr[i]] * n2

gal_vel.extend([0])
gal_rad.extend([0])
gal_dmass.extend([0])
numgood.extend([0])
gal_vel   = array(gal_vel)
gal_rad   = array(gal_rad)
gal_dmass = array(gal_dmass)
num = array(nu)

DataOut = column_stack((gal_dmass,gal_vel,gal_rad,num,galr,galv,cz,id,nu))
savetxt('dynm_out_' + option + '.dat', DataOut)




### BIN for error bars - error in halo mass for given dynmass###
big_dm  = dyn_mass[(dyn_mass > 0) & (dyn_mass < Infinity) & (nu[ind] > 3)]
big_hm  = halo_m[(dyn_mass > 0) & (dyn_mass < Infinity)   & (nu[ind] > 3)]
bins    = [0] * 20      #store errors
xx      = [0] * 20      #store middle of bin
yy      = [0] * 20      #store avg halo_m for bin
start   = min(log10(big_dm[big_dm > 0])) 
end     = max(log10(big_dm[big_dm > 0]))
binsize = (end - start)/len(bins)
for pizza in range(0,len(bins)):
    temp_hm      =  big_hm[(log10(big_dm[big_dm > 0]) > (binsize*pizza + start)) &(log10(big_dm) < (binsize*(pizza + 1)+start))] 
    bins[pizza]  =  sqrt(sum(abs(temp_hm - temp_hm.mean())**2)/(len(temp_hm)-1))
    xx[pizza]    =  (binsize*pizza + binsize*(pizza+1))/2 + start
    yy[pizza]    =  mean(temp_hm)

norm_err  =  array(bins)/array(yy)     #normalized error

DataOut = column_stack((xx,norm_err))
savetxt('./norm_err/' + option + '.dat', DataOut)






# -~-~-~-~-~-~-~PLOT~-~-~-~-~-~-~#
#dyn_mass = 4.466*dyn_mass
plot(log10(dyn_mass[dyn_mass != 0]),halo_m[dyn_mass != 0],'y.',zorder=-100,label='$N_\mathrm{galaxies}$ > 2')
#plot(log10(dyn_mass[nu[ind] > 3]),halo_m[nu[ind] > 3],'g.', zorder=-100,label='$N_\mathrm{galaxies}$ > 3')                            #overplot halos w/ 4 members & up
plot(log10(dyn_mass[nu[ind] > 4]),halo_m[nu[ind] > 4],'b.',zorder=-100,label='$N_\mathrm{galaxies}$ $>=$ 5')                            #overplot halos w/ 4 members & up
#plt.errorbar(array(xx),yy,yerr=bins,barsabove=True,ls='None',label='halo mass dispersion', capsize=4.5, color='black',zorder = 300,elinewidth=2.5)
#plot(xx,yy,'ko')
x = np.arange(0, 15, 0.5);
y = x #+ .7
plot(x,y,color='r',label='1 to 1 line')        #plot 1:1 line
ylim(10,15)
xlim(9,14.5)
xlabel('log Dynamical Mass')
ylabel('log Halo Mass')
title(mass_est)
#title('')
legend(loc=2,numpoints=1)
#savefig('bigmasses.png')


# - - - - - - make error plot - - - - - - #
      #norm_err_files - lists xx & norm_err 
          #for several modes/options
fig=plt.figure()
norm_err_files = loadtxt('norm_err/err_files2',dtype='str',unpack=True)
for i in range(0,len(norm_err_files)):
    xx1,err1 = loadtxt('norm_err/' + norm_err_files[i] + '.dat', unpack=True)
    plot(xx1,err1,label=norm_err_files[i])

title('Normalized Errors for various mass estimators')
xlabel('log Dyn Mass')
ylabel('Normalized Scatter')
legend(loc=8,numpoints=1)
savefig('errors.png')



#!/usr/bin/python
#gp_dynmass.py
#execfile('gp_dynmass.py')
#find dynamical mass after applying group finder
################## CHOOSE MASS ESTIMATOR ###################
mass_est = 'pow(velarr,2) * radarr / G'
#mass_est = 'pow(velarr,2)  / G'
#mass_est = '(3. * pow(velarr,2) / radarr) / G' #(used with angles3)

option = 'avg'
          #'centrals' or 'averages'
inverse = 'no' 
          #only pick yes if using 3rd mass_est where need it
############################################################
  
  #import things
import numpy as np
from pylab import *
import math
from matplotlib import rc, rcParams
import matplotlib.pyplot as plt
from numpy import *
import random

#linking length combos
bzs = ['14','13','12','12','11','12','12'] #parallel, vel
bps = ['75','15','22','13','14','14','17'] #perpend, rad
#       0    1    2     3    4    5    6
llpick = 0  #0 to 4

filename = "fof9_noG_" + bzs[llpick] + "_" + bps[llpick] + ".dat"
#filename="fof9_noG_14_75.dat"   #linking lengths .14 and .75
#name = 'old'
#filename="fof9_noG_13_15.dat"  #linking lengths .13 and 1.5
#name = 'new'
#filename="fof9_noG_12_22.dat"  #linking lengths .13 and 1.5
#name = 'new1'


file=np.loadtxt(filename)

id = file[:,0] #id number
of = file[:,1] #location in original catalog
ra = file[:,2] #real cz, no peculiars for rspace
dec= file[:,3] #absolute R band mag
cz = file[:,4] #halo id number

filename="./Resolve_5001_hod1_mock_zspace.dat"

file=np.loadtxt(filename)

mr = file[:,3][list(of)] #absolute R band mag
hm = file[:,5][list(of)] #log halo mass
n2 = file[:,6][list(of)] #number in halo
ch = file[:,7][list(of)] #central 0 satellite 1
iold = file[:,4][list(of)]

#define constants
G  = 4.302E-9   #Mpc(km/s)^2 Msun^-1
c  = 3.0*10**5  #km/sec
h = .69
Ho = 100.0      #+/-0.80 (km/s)/Mpc was using 69.3 but use 100 and then mass will be in terms of inverse little h


#~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ DEFINE INDEX MATRIX ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#based on the halo ID, not whether a central, assuming IDs in order
ind=[0]    #initial condition
i=1
while i < (len(cz)-1):
    tempind = i
    ind.extend([tempind])
    while (i < (len(cz)-2)) & (id[i+1] == id[i]):
        i = i + 1
    i=i+1

ind.extend([tempind+1]) #boundary condition

#define nu (number in each group)
nu = array(ind[1::]) - array(ind[0:len(ind)-1])
nu = list(nu)
nu.append(1)            #fill in last group
nu = array(nu)




#~ ~ ~ ~ ~ ~ ~ DEFINE VEL DISPERSION && AVG RADIUS~ ~ ~ ~ ~ ~ ~ ~#
velarr   = [0] * len(ind)
radarr   = [0] * len(ind)
dyn_mass = [0] * len(ind)
galr     = [0] * len(ind)
galv     = [0] * len(ind)
for i in range(0,len(ind)-1):
    if ind[i+1] == ind[i] + 1:     #rid of halos w/ 1 galaxy
        dyn_mass[i] = 0.
        galr[ind[i]:ind[i+1]] = [0]
        galv[ind[i]:ind[i+1]] = [0]
    elif ind[i+1] == ind[i] + 2:   #rid of halos w/ 2 galaxies
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
                mrs = mr[ind[i]:ind[i+1]] + .0001*array([random.randint(1,1000) for _ in range(num)])
                ceni = [mrs == min(mrs)]
            sig  = sqrt(sum((czs - czs[ceni])**2)/(num - 1))
            galv[ind[i]:ind[i+1]] = czs - mean(czs)
            velarr[i] = sig
            #proj radii
            ras     = ra[ind[i]:ind[i+1]]
            decs    = dec[ind[i]:ind[i+1]]
            angles  = sqrt(((ras - ras[ceni])*cos((math.pi/180.)*decs[ceni]))**2 + (decs - decs[ceni])**2)
            radii   = (math.pi/180.)* angles  * czs[ceni] / Ho
            galr[ind[i]:ind[i+1]] = radii
            if inverse == 'no':
                radarr[i] = mean(radii)
            elif inverse == 'yes':
                radarr[i] = mean(1/radii)
            else:
                print 'do you want the inverse of r or not??'
        elif option == 'avg':
            #vel dispersion
            czs = cz[ind[i]:ind[i+1]]
            num = ind[i+1] - ind[i]
            avg_cz  = sum(czs)/num
            sig = sqrt(sum((czs - avg_cz)**2)/(num - 1))
            galv[ind[i]:ind[i+1]] = czs - avg_cz
            velarr[i] = sig
            #proj radii
            ras     = ra[ind[i]:ind[i+1]]
            decs    = dec[ind[i]:ind[i+1]]
            avg_ra  = sum(ras)/num
            avg_dec = sum(decs)/num
            angles  = sqrt(((ras - avg_ra)*cos((math.pi/180.)*avg_dec))**2 + (decs - avg_dec)**2)
            radii   = (math.pi/180.)* angles * avg_cz / Ho
            galr[ind[i]:ind[i+1]] = radii
            if inverse == 'no':
                radarr[i] = mean(radii)
            elif inverse == 'yes':
                radarr[i] = mean(1/radii)
            else:
                print 'do you want the inverse of r or not??'

velarr   = array(velarr)            #velocity dispersion for each group
radarr   = array(radarr)            #average projected radius for e/ group
dyn_mass = eval(mass_est)           #choose mass estimator above
galr.append(0)
galv.append(0)


# ------------ save & define per galaxy arrays -------------- #

gal_dmass = [0] * (len(ind) + 1)
num       = [0] * len(ind)       #no plus one bc last index isnt 0
gal_vel   = [0] * (len(ind) + 1)
gal_rad   = [0] * (len(ind) + 1)
for i in range(0,len(ind)-1):
    n = ind[i+1] - ind[i]
    gal_dmass[ind[i]:ind[i+1]] = [dyn_mass[i]] * n
    num[ind[i]:ind[i+1]] = [nu[i]] * n
    gal_vel[ind[i]:ind[i+1]] = [velarr[i]] * n
    gal_rad[ind[i]:ind[i+1]] = [radarr[i]] * n

num.extend([1])
gal_vel.extend([0])
gal_rad.extend([0])
gal_dmass.extend([0])
gal_halom = hm
gal_vel   = array(gal_vel)
gal_rad   = array(gal_rad)
gal_dmass = array(gal_dmass)
num = array(num)

filename="./dynm_out_" + option + ".dat"
file=np.loadtxt(filename)
gd_mass = file[:,0][list(of)]
gd_vels = file[:,1][list(of)] 
gd_rads = file[:,2][list(of)]
gd_galr = file[:,4][list(of)]
gd_galv = file[:,5][list(of)]

DataOut = column_stack((gal_dmass,gal_vel,gal_rad,num,id,
	n2,iold,mr,hm,of,ra,dec,gd_mass,gd_vels,gd_rads,
        galr,galv,gd_galr,gd_galv,cz))  #save file of info
savetxt('gp_dynm_out_'+ option + '_' + bzs[llpick] + '_' + bps[llpick] + '.dat', DataOut)

#f.write("/n".join(map(lambda x: str(x), match_mass)))


### ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ BIN for error bars ~ ~ ~ ~ ~ ~ ~ ~ ~ ###

#big_dm  = gal_dmass[(gal_dmass > 0) & (gal_dmass < Infinity) & (num > 4)]
#big_hm  = gal_halom[(gal_dmass > 0) & (gal_dmass < Infinity) & (num > 4)]
#bins    = [0] * 20      #store errors here
#xx      = [0] * 20      #store middle of bin
#yy      = [0] * 20      #store avg halo_m for bin
#start   = min(log10(big_dm[big_dm > 0]))
#end     = max(log10(big_dm[big_dm > 0]))
#binsize = (end - start)/len(bins)
#for pizza in range(0,len(bins)):
#    temp_hm      =  big_hm[(log10(big_dm[big_dm > 0]) > (binsize*p#izza + start)) &(log10(big_dm) < (binsize*(pizza + 1)+start))] 
#    bins[pizza]  =  sqrt(sum(abs(temp_hm - temp_hm.mean())**2)/(le#n(temp_hm)-1))
#    xx[pizza]    =  (binsize*pizza + binsize*(pizza+1))/2 + start
#    yy[pizza]    =  mean(temp_hm)

#norm_err  =  array(bins)/array(yy)     #normalized error
#DataOut = column_stack((xx,norm_err))  #save file of info
#savetxt('./norm_err/' + 'gp_' + option + '.dat', DataOut)

#PICKLINGG
#filename='galr.p'
#file=open(filename,'w')
#cPickle.dump(galr, file, protocol=2)
#file.close()

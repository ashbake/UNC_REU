#!/usr/bin/python
#execfile('eco_dynmass.py')
#find dynamical mass after applying group finder


mass_est = 'pow(velarr,2) * radarr / G'
          #mass_est = 'pow(velarr,2)  / G'
          #mass_est = '(3. * pow(velarr,2) / radarr) / G' #(used with angles3)
option = 'avg'
          #'centrals' or 'averages'
inverse = 'no' 
          #only pick yes if using 3rd mass_est where need it
bzs = ['.12','.13','.14','.12','.11']
bps = ['2.2','1.5','.75','1.3','1.4']
llpick = 2
##########################################################################

import numpy as np
from pylab import *
import math
from matplotlib import rc, rcParams
import matplotlib.pyplot as plt
from numpy import *
import random
import re 
import cPickle as pickle

bz = bzs[llpick]
bp = bps[llpick]
bbz = re.sub('[.]','',bz)
bbp = re.sub('[.]','',bp)
path = "/srv/two/ashbake/REU/"
filename=path + "eco_data/noG_" + bz + "_" + bp + "_infof_ecospring.txt"       #real data .14 & .75


file=np.loadtxt(filename)

id = file[:,0] #id number
of = file[:,1] #location in original catalog
ra = file[:,2] #real cz, no peculiars for rspace
dec= file[:,3] #absolute R band mag
cz = file[:,4] #halo id number

file=np.loadtxt(path + 'eco_data/in_ecofall.txt')
mr = file[:,3][list(of)]

#constants
G  = 4.302E-9   #Mpc(km/s)^2 Msun^-1
c  = 3.0*10**5  #km/sec
Ho = 100.*.7      # (km/s)/Mpc was using 69.3 but use 100 and then mass will be in terms of inverse little h. amanda uses Ho=70 so use keep it as 70

#based on the halo ID, not whether a central, assuming IDs in order
ind=[0]    #initial condition
i=1
while i < (len(cz)-1):
    while (i < (len(cz)-2)) & (id[i+1] == id[i]):
        i = i + 1
    ind.extend([i+1])
    i=i+1

#define num (number in each group)
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
avg_cz   = [0] * len(ind)
avg_dec  = [0] * len(ind)
avg_ra   = [0] * len(ind)

for i in range(0,len(ind)-1):
    if ind[i+1] == ind[i] + 1:     #rid of halos w/ 1 galaxy
        dyn_mass[i] = 0.
        galr[ind[i]:ind[i+1]] = [0]
        galv[ind[i]:ind[i+1]] = [0]
        avg_cz[ind[i]:ind[i+1]] = [0]
        avg_ra[ind[i]:ind[i+1]] = [0]
        avg_dec[ind[i]:ind[i+1]] = [0]
    elif ind[i+1] == ind[i] + 2:   #rid of halos w/ 2 galaxies
        dyn_mass[i] = 0.
        galr[ind[i]:ind[i+1]] = [0]*2
        galv[ind[i]:ind[i+1]] = [0]*2
        avg_cz[ind[i]:ind[i+1]] = [0]*2
        avg_ra[ind[i]:ind[i+1]] = [0]*2
        avg_dec[ind[i]:ind[i+1]] = [0]*2
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
            avg_cz[ind[i]:ind[i+1]]  = [sum(czs)/num]*num
            sig = sqrt(sum((czs - avg_cz[ind[i]])**2)/(num - 1))
            galv[ind[i]:ind[i+1]] = czs - avg_cz[ind[i]]
            velarr[i] = sig
            #proj radii
            ras     = ra[ind[i]:ind[i+1]]
            decs    = dec[ind[i]:ind[i+1]]
            avg_ra[ind[i]:ind[i+1]]  = [sum(ras)/num]*num
            avg_dec[ind[i]:ind[i+1]] = [sum(decs)/num]*num
            angles  = sqrt(((ras - avg_ra[ind[i]])*cos((math.pi/180.)*avg_dec[ind[i]]))**2 + (decs - avg_dec[ind[i]])**2)
            radii   = (math.pi/180.)* angles * avg_cz[ind[i]] / Ho
            galr[ind[i]:ind[i+1]] = radii
            radarr[i] = mean(radii)



velarr   = array(velarr)            #velocity dispersion for each group
radarr   = array(radarr)            #average projected radius for e/ group
dyn_mass = eval(mass_est)           #choose mass estimator above
galr.append(0)
galv.append(0)
avg_cz.append(0)
avg_ra.append(0)
avg_dec.append(0)

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
gal_vel   = array(gal_vel)
gal_rad   = array(gal_rad)
gal_dmass = array(gal_dmass)
num = array(num)

DataOut = column_stack((gal_dmass,gal_vel,gal_rad,num,id,mr,of,ra,dec,galr,galv,cz))  #save file of info
savetxt('eco_dynm_out_'+ bz + '_' + bp + '_spr021114.dat', DataOut)

eco_dict = {}
eco_dict['gal_dmass']  =  gal_dmass
eco_dict['gal_vel']    =  gal_vel
eco_dict['gal_rad']    =  gal_rad
eco_dict['num']        =  num
eco_dict['id']         =  id
eco_dict['mr']         =  mr
eco_dict['of']         =  of
eco_dict['ra']         =  ra
eco_dict['dec']        =  dec
eco_dict['galr']       =  galr
eco_dict['galv']       =  galv
eco_dict['cz']         =  cz
eco_dict['ind']        =  ind
eco_dict['avg_cz']     = avg_cz
eco_dict['avg_ra']     = avg_ra
eco_dict['avg_dec']    = avg_dec

pickle.dump(eco_dict, open( 'eco_dict' + bz + bp+ '.p', "wb" ) )
#eco_dict = pickle.load( open( "eco_dict" + bz + bp+ ".p", "rb" ) )

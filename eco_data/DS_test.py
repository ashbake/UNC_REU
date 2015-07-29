#!/usr/bin/python
#execfile('AS_test.py')
#see cohen 2014 and anderson shectman 1988

#execfile('../eco_wrapper.py')
path = '/srv/two/ashbake/REU/eco_data/'

def scrambled(orig):
    dest = orig[:]
    random.shuffle(dest)
    return dest

del_sim = array([[0.]*100]*len(cz))
delta = [0.] * len(cz)
DELTA = [0.] * len(cz)
DELSIM= [0.] * len(cz)
p_value = [0.] * len(cz)
for i in range(0,len(ind)-1):
    n=int(num[ind[i]])
    #size = round(n/4)
    size=11
    if n > size:
        czs = cz[ind[i]:ind[i]+n]
        ras = ra[ind[i]:ind[i]+n]
        decs = dec[ind[i]:ind[i]+n]
        gals_near = [0] * n
        for j in range(0,n):
            angles  = sqrt(((ras - ras[j])*cos((math.pi/180.)*dec[j]))**2 + (decs - decs[j])**2)
            dists = (math.pi/180.)* angles * czs[j] / Ho
            close = argsort(dists)[1:size+1] #indices of closest galaxies
            sig = gal_vel[ind[i]]
            czs_loc = czs[close]    #czs of closest gals
            sig_loc = sqrt(sum((czs_loc - mean(czs_loc))**2)/(n - 1))
            delta[ind[i]+j] = (size/gal_vel[ind[i]]**2) * ((mean(czs_loc) - mean(czs))**2 + (sig_loc - sig)**2)
            #simulated values
            for k in range(0,100):
                cz_scram = array(scrambled(list(czs)))
                angles  = sqrt(((ras - ras[j])*cos((math.pi/180.)*dec[j]))**2 + (decs - decs[j])**2)
                dists = (math.pi/180.)* angles * cz_scram[j] / Ho
                close = argsort(dists)[1:size] #indices of closest galaxies
                sig = gal_vel[ind[i]]
                czs_loc = cz_scram[close]    #czs of closest gals
                sig_loc = sqrt(sum((czs_loc - mean(czs_loc))**2)/(n - 1))
                del_sim[ind[i]+j][k] = (size/gal_vel[ind[i]]**2) * ((mean(czs_loc) - mean(czs))**2 + (sig_loc - sig)**2)
            
        DELTA[ind[i]:ind[i+1]] = [sum(sqrt(delta[ind[i]:ind[i]+n]))]*n
        DELSIM[ind[i]:ind[i+1]] = [sum(sqrt(del_sim[ind[i]:ind[i]+n]),axis=0)]*n
        p_value[ind[i]:ind[i+1]] = [len(where(DELSIM[ind[i]] > DELTA[ind[i]])[0])/100.]*n


execfile(path + 'calc_grps/check_vir.py')   #do anderson darling test to get alpha values
alpha=array(alpha)
DELTA = array(DELTA)
p_value = array(p_value)

#PLOT ALLLLLLL THE GROUPS###
import numpy as np
import matplotlib.pyplot as plt

big = where(num[ind] > size)[0]
for hi in range(0,len(big)):
    i = big[hi]
    n=int(num[ind[i]])
    czs = cz[ind[i]:ind[i]+n]
    ras = ra[ind[i]:ind[i]+n]
    decs = dec[ind[i]:ind[i]+n]
    data=[ras,decs,array(delta[ind[i]:ind[i+1]]),czs,goodgroupn[ind[i]:ind[i+1]]]
    labels = ['{0}'.format(k) for k in data[4]]
    plt.subplots_adjust(bottom = 0.1)
    plt.scatter(
        data[0], data[1], marker = 'o', c = data[2],s=50+abs((czs - mean(czs))))
    #for label, x, y in zip(labels, data[0], data[1]):
    #    plt.annotate(
    #        label, 
    #        xy = (x, y), xytext = (-20, 20),
    #        textcoords = 'offset points', ha = 'right', va = 'bottom',
    #        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.3),
    #        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    title('P val = ' + str(p_value[ind[i]]) + ', alpha = ' + str(alpha[ind[i]]))
    colorbar()
    plt.show()

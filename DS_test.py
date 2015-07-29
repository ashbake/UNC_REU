#!/usr/bin/python
#execfile('AS_test.py')
#see cohen et al 2014 and dressler shectman 1988 & Hou et al 2014

execfile('gp_wrapper.py')

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
    size=11
    if n >= size:
        czs = cz[ind[i]:ind[i]+n]
        ras = ra[ind[i]:ind[i]+n]
        decs = dec[ind[i]:ind[i]+n]
        gals_near = [0] * n
        for j in range(0,n):
            angles  = sqrt(((ras - ras[j])*cos((math.pi/180.)*dec[j]))**2 + (decs - decs[j])**2)
            dists = (math.pi/180.)* angles * czs[j] / Ho
            close = argsort(dists)[1:size] #indices of closest galaxies
            sig = fof_vel[ind[i]]
            czs_loc = czs[close]    #czs of closest gals
            sig_loc = sqrt(sum((czs_loc - mean(czs_loc))**2)/(n - 1))
            delta[ind[i]+j] = (size/fof_vel[ind[i]]**2) * ((mean(czs_loc) - mean(czs))**2 + (sig_loc - sig)**2)
            #simulated values
            for k in range(0,100):
                cz_scram = array(scrambled(list(czs)))
                angles  = sqrt(((ras - ras[j])*cos((math.pi/180.)*dec[j]))**2 + (decs - decs[j])**2)
                dists = (math.pi/180.)* angles * cz_scram[j] / Ho
                close = argsort(dists)[1:size] #indices of closest galaxies
                sig = fof_vel[ind[i]]
                czs_loc = cz_scram[close]    #czs of closest gals
                sig_loc = sqrt(sum((czs_loc - mean(czs_loc))**2)/(n - 1))
                del_sim[ind[i]+j][k] = (size/fof_vel[ind[i]]**2) * ((mean(czs_loc) - mean(czs))**2 + (sig_loc - sig)**2)
            
        DELTA[ind[i]:ind[i+1]] = [sum(sqrt(delta[ind[i]:ind[i]+n]))]*n
        DELSIM[ind[i]:ind[i+1]] = [sum(sqrt(del_sim[ind[i]:ind[i]+n]),axis=0)]*n
        p_value[ind[i]:ind[i+1]] = [len(where(DELSIM[ind[i]] > DELTA[ind[i]])[0])/100.]*n


execfile('gp_check_vir.py')   #do anderson darling test to get alpha values
alpha=array(alpha)
DELTA = array(DELTA)

#PLOT ALLLLLLL THE GROUPS###
import numpy as np
import matplotlib.pyplot as plt

big = where(num[ind] > size)[0]
for hi in range(0,len(big)):
    plt.figure()
    i = big[hi]
    n=int(num[ind[i]])
    czs = cz[ind[i]:ind[i]+n]
    ras = ra[ind[i]:ind[i]+n]
    decs = dec[ind[i]:ind[i]+n]
    data=[ras,decs,array(delta[ind[i]:ind[i+1]]),czs,nold[ind[i]:ind[i+1]]]
    labels = ['{0}'.format(k) for k in data[4]]
    plt.subplots_adjust(bottom = 0.1)
    plt.scatter(
        data[0], data[1], marker = 'o', c = data[2],s=60+abs(data[3]-mean(data[3])))
    for label, x, y in zip(labels, data[0], data[1]):
        plt.annotate(
            label, 
            xy = (x, y), xytext = (-20, 20),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.3),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    title('P val = ' + str(p_value[ind[i]]) + ', alpha = ' + str(alpha[ind[i]]))
    colorbar()
    


#PLOT TWO METHODS
scatter(array(alpha)[num > 10],array(delta)[num > 10],c=log10(gal_dmass)[num > 10])
xlabel('alpha')
ylabel('delta')
text(1.3,23,'MOCK',fontsize=20)
xlim(-.05,2)
ylim(-.03,25)
cbar = colorbar()
cbar.set_label("Dynamical Mass",rotation=270)
savefig('alphadeltaMOCK')


x= array(alpha)
y= array(delta)
xx = MAD(x[where((num > 10) & (num < 25))],y[where((num > 10) & (num < 25))],3)
yy = MAD(x[where((num >= 25) & (num < 50 ))],y[where((num >= 25) & (num < 50))],3)
zz = MAD(x[where(num >= 50)],y[where(num >= 50)],3)

#plot(array(xx.xout),xx.yout,'ro',zorder = 450,ms = 10)
plot(xx.xout,xx.yout,'b-',xx.xout,xx.yout,'bo',zorder = 400,ms = 5)
plt.errorbar(array(xx.xout),xx.yout,yerr=xx.errs,xerr=None,barsabove=True,
             ls='None', capsize=4.5, color='blue',zorder = 300,
             elinewidth=3)
plot(yy.xout,yy.yout,'-',yy.xout,yy.yout,'o',c='orange',zorder = 400,ms = 5)
plt.errorbar(array(yy.xout),yy.yout,yerr=yy.errs,xerr=None,barsabove=True,
             ls='None', capsize=4.5, color='orange',zorder = 300,
             elinewidth=3)
plot(zz.xout,zz.yout,'r-',zz.xout,zz.yout,'ro',zorder = 400,ms = 5)
plt.errorbar(array(zz.xout),zz.yout,yerr=zz.errs,xerr=None,barsabove=True,
             ls='None', capsize=4.5, color='red',zorder = 300,
             elinewidth=3)


matplotlib.rcParams.update({'font.size': 16})
text(1.3,9,'Mock',fontsize=20)

ylim(0,10)
xlabel('alpha')
ylabel('delta')

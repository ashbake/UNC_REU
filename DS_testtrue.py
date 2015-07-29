#!/usr/bin/python
#execfile('AS2_test.py')
#see cohen 2014 and anderson shectman 1988
#computes AS test for mock for true groups

execfile('gp_dynmass.py')
execfile('gp_check_vir.py')   #do anderson darling test to get alpha values

def scrambled(orig):
    dest = orig[:]
    random.shuffle(dest)
    return dest
gal_vel=array(gal_vel)

del_sim = array([[0.]*100]*len(cz))
delta = [0.] * len(cz)
DELTA = [0.] * len(cz)
DELSIM= [0.] * len(cz)
p_value = [0.] * len(cz)
for i in range(0,len(ind)-1):
    n=int(ind[i+1] - ind[i])
    #size = round(n/4)
    size=11
    if (n >= size) & (gal_vel[ind[i]] !=0):
        czs = cz[ind[i]:ind[i]+n]
        ras = ra[ind[i]:ind[i]+n]
        decs = dec[ind[i]:ind[i]+n]
        gals_near = [0] * n
        for j in range(0,n):
            angles  = sqrt(((ras - ras[j])*cos((math.pi/180.)*dec[j]))**2 + (decs - decs[j])**2)
            dists = (math.pi/180.)* angles * czs[j] / Ho
            close = argsort(dists)[1:size] #indices of closest galaxies
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


#execfile('gp_check_vir.py')   #do anderson darling test to get alpI hha values
#alpha=array(alpha)
DELTA = array(DELTA)

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
    plt.show()


#PLOT TWO METHODS
scatter(array(alpha)[num > 10],array(delta)[num > 10],c=log10(gal_dmass)[num > 10])
scatter(array(alpha)[num > 75],array(delta)[num > 75],c='k')
xlabel('alpha')
ylabel('delta')
text(1.1,30,'MOCK',fontsize=20)
xlim(-.05,1.6)
ylim(-.03,38)
cbar = colorbar()
cbar.set_label("Dynamical Mass",rotation=270)
savefig('alphadeltaMOCK.pdf')



x= array(alpha)
y= array(delta)
x1 = x[where((num > 10) & (num < 25))]
y1 = y[where((num > 10) & (num < 25))]
x2 = x[where((num >= 25) & (num < 50))]
y2 =y[where((num >= 25) & (num < 50))]
x3 = x[where((num >= 50) & (y > 0))]
y3 = y[where((num >= 50) & (y > 0))]

#plot(array(xx.xout),xx.yout,'ro',zorder = 450,ms = 10)
fig=plt.figure()
ax = fig.add_subplot(111)

ax1 = subplot(313)
plot(x1,y1,'bo',markeredgecolor='b',label=r'$11 \leq N_{grp} < 25$')
ax1.set_ylim(0,28)
ax1.set_xlim(-0.1,2)
ax1.set_xlabel('alpha')

ax2 = subplot(312)
plot(x2,y2,'o',c='orange',markeredgecolor='orange',label=r'$25 \leq N_{grp} < 50$')
ax2.set_ylim(0,28)
ax2.set_xlim(-0.1,2)
setp( ax2.get_xticklabels(), visible=False)
ax2.set_ylabel('delta')

ax3 = subplot(311)
plot(x3,y3,'ro',markeredgecolor='red',label=r'$50 \leq N_{grp} \leq 100$')
ax3.set_ylim(0,28)
ax3.set_xlim(-0.1,2)
setp( ax3.get_xticklabels(), visible=False)

matplotlib.rcParams.update({'font.size': 16})

plt.setp(ax1, yticks=[0,5,10,15,20,25])
plt.setp(ax2, yticks=[5,10,15,20,25])
plt.setp(ax3, yticks=[5,10,15,20,25])
fig.subplots_adjust(hspace=0)
ax.tick_params(direction='in', length=2, width=6, labeltop = 'off')

text(1.3,20,'MOCK',fontsize=20)
#text(1.2,22,'ECO',fontsize=20)









plot(dd.xout,dd.yout,'m-',dd.xout,dd.yout,'mo',zorder = 400,ms = 5)
plt.errorbar(array(dd.xout),dd.yout,yerr=dd.errs,xerr=None,barsabove=True,
             ls='None', capsize=4.5, color='m',zorder = 300,
             elinewidth=3)

#!/usr/bin/python
#execfile('gp_plots.py')
execfile('gp_wrapper.py')
#execfile('gp_check_vir.py')
from scipy.stats import spearmanr


# for making plots of dyn mass/vel/radii before and after
fig = plt.figure()
ax = fig.add_subplot(111)


plt.scatter(fof_dmass,dmass,marker='o',edgecolor='None',
	    facecolors='gray',label=r'$\mathsf{3 < N_{gal} < 9}$')

plt.scatter(fof_dmass[where(num > 9)],dmass[where(num > 9)],
	    marker='D',edgecolor='None',facecolors='skyblue',
	    s=20,label=r'$\mathsf{N_{gal} > 9}$')
font={'family' : 'normal',
      'weight' : 'normal',
      'size'   :  20}

sr_mass = round((spearmanr(fof_dmass[where((fof_dmass > 0) & (dmass > 0))],
          dmass[where((fof_dmass > 0) & (dmass > 0))]))[0],2)

xx = fof_dmass[where((fof_dmass > 0) & (dmass > 0))]
yy = dmass[where((fof_dmass > 0) & (dmass > 0))]
SS = sum(abs(yy - xx))/len(xx)

text(10,13.8,r'$b_{\perp} = 0.12$',fontsize=18,horizontalalignment='center',
     verticalalignment='center')
text(10.5,14.6,'a)  ' + r'$\rho = $ ' + str(sr_mass)+ ',  S =  ' + str(round(SS,2)),fontsize=16,horizontalalignment='center',
     verticalalignment='center')

text(10,13.2,r'$b_{z} = 1.3 $',fontsize=18,horizontalalignment='center',
     verticalalignment='center')


xlabel('Dyn Mass Post Group Finder',fontsize=20)
ylabel('Dyn Mass of True Groups',fontsize=20)
xlim(9,15)
ylim(9,15)
x=y=[9,15]
plot(x,y,'k-')
matplotlib.rc('font', **font)
legend(loc=4,scatterpoints=1,prop={'size':16,'weight':'bold'})
savefig('plotmass_' + llens + '.pdf')
show()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ vel ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clf()
plt.scatter(fof_vel,vel,marker='o',edgecolor='None',
	    facecolors='gray',label=r'$\mathsf{3 < N_{gal} < 9}$')

plt.scatter(fof_vel[where(num > 9)],vel[where(num > 9)],#
	    marker='D',edgecolor='None',facecolors='lightcoral',
	    s=20,label=r'$\mathsf{N_{gal} > 9}$')
#font={'family' : 'normal',
#      'weight' : 'normal',
#      'size'   :  20}
xx = fof_vel[where((fof_dmass > 0) & (dmass > 0))]
yy = vel[where((fof_dmass > 0) & (dmass > 0))]
SS = sum(abs(yy - xx))/len(xx)

sr_vel = round((spearmanr(fof_vel[where((fof_dmass > 0) & (dmass > 0))],
          vel[where((fof_dmass > 0) & (dmass > 0))]))[0],2)
text(150,450,'e)    ' + r'$\rho = $ ' + str(sr_vel)+ ',    S =  ' + str(round(SS,0)),fontsize=16,horizontalalignment='center',
     verticalalignment='center')

xlabel('Velocity Dispersion Post Group Finder',fontsize=20)
ylabel('Velocity Dispersion of True Groups',fontsize=20)
xlim(0,500)
ylim(0,500)
x=y=[0,500]
plot(x,y,'k-')
#matplotlib.rc('font', **font)
legend(loc=4,scatterpoints=1,prop={'size':16,'weight':'bold'})
savefig('plotvel_' + llens + '.pdf')
show()


#~~~~~~~~~~~~~~~~~~~~~~~~ radius ~~~~~~~~~~~~~~~~~~~~~

clf()
plt.scatter(fof_rad,rad,marker='o',edgecolor='None',
	    facecolors='gray',label=r'$\mathsf{3 < N_{gal} < 9}$')

plt.scatter(fof_rad[where(num > 9)],rad[where(num > 9)],#
	    marker='D',edgecolor='None',facecolors='forestgreen',
	    s=20,label=r'$\mathsf{N_{gal} > 9}$')
#font={'family' : 'normal',
#      'weight' : 'normal',
#      'size'   :  20}

sr_rad = round((spearmanr(fof_rad[where((fof_dmass > 0) & (dmass > 0))],
          rad[where((fof_dmass > 0) & (dmass > 0))]))[0],2)
xx = fof_rad[where((fof_dmass > 0) & (dmass > 0))]
yy = rad[where((fof_dmass > 0) & (dmass > 0))]
SS = sum(abs(yy - xx))/len(xx)

text(.18,.55,'f)   ' + r'$\rho = $ ' + str(sr_rad) + ',  S = ' + str(round(SS,3)),fontsize=16,horizontalalignment='center',
     verticalalignment='center')

xlabel('Group Radii Post Group Finder',fontsize=20)
ylabel('True Group Radii',fontsize=20)
xlim(0,.8)
ylim(0,.6)
x=y=[0,.8]
plot(x,y,'k-')
#matplotlib.rc('font', **font)
legend(loc=4,scatterpoints=1,prop={'size':16,'weight':'bold'})
savefig('plotrad_' + llens + '.pdf')
show()


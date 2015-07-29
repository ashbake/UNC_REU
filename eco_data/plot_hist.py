#!/usr/bin/python
#execfile('plottest.py')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

path = '/srv/two/ashbake/REU/eco_data/'
execfile(path+'eco_wrapper.py')
#execfile('AS_test.py')
execfile('calc_grps/check_vir.py')
clf()
fig = plt.figure()

plt.scatter(gal_hamass[where(gal_dmass > 0)],
            gal_dmass[where(gal_dmass > 0)],label=r'$3 \ \leq  \ N_{grp} \ \leq \ 5$',
            edgecolor='None',c='salmon')

xx = BinFun(gal_hamass[where(gal_dmass > 0)],
            gal_dmass[where(gal_dmass > 0)],5)

dm = gal_dmass[where(num > 5)]
hm = gal_hamass[where(num > 5)]

plt.scatter(hm,dm,edgecolor='None',c='skyblue',label=r'$N_{grp} \ > \  5$')

dm = gal_dmass[where(num > 5)]
hm = gal_hamass[where(num > 5)]
alphas = array(alpha)[num > 5]

plt.scatter(hm[alphas < .05],dm[alphas < .05],marker='*',
               edgecolor='None',c='navy',s=80,label=r'$\alpha \ < \ 5 \%$')

plot(xx.xout,xx.yout,'ko')
plt.errorbar(array(xx.xout),xx.yout,yerr=xx.errs,xerr=None,barsabove=True,
             ls='None', capsize=4.5, color='black',zorder = 300,
             elinewidth=1.5)

x = y = [11,15.5]
plot(x,y,'k-',ms=6)
xlabel('HAM Masses')
ylabel('Dynamical Masses')
ylim(9,15.2)
xlim(11,15)

plt.legend(loc = 4,scatterpoints=1,prop={'size':16})
text(11.3,14.5,'(a) ECO Catalog')
plt.savefig('figures/plot_masses_'+bz + '_' + bp+'.pdf')
show()


####### HISTOGRAM PLOT ########
clf()
fig = plt.figure(figsize=(5,5)) 
plt.subplots(2, sharex=True)
ax=fig.add_subplot(1,1,1)
gcf().subplots_adjust(bottom=0.20)
#plt.gca().xaxis.set_major_locator(MaxNLocator(prune='lower'))
#ax.hist(gal_dmass[ind][gal_dmass[ind] > 0],np.linspace(min(mass_grp),
#         max(mass_grp),15),facecolor='red',edgecolor='None',alpha=.3,lw=2,
#         label=r'$\mathsf{M_{dynamical}}$')
#ax.hist(mass_grp[ind][gal_dmass[ind] > 0],np.linspace(min(mass_grp),
#         max(mass_grp),15),alpha=1,lw=2,hatch="///",histtype='step',
#         label=r'$\mathsf{M_{luminosity}}$',edgecolor='navy')
ax.hist(gal_dmass[ind][num[ind] > 3],np.linspace(min(gal_hamass),
         max(gal_hamass),15),facecolor='red',edgecolor='None',alpha=.3,lw=2,
         label=r'$\mathsf{M_{dynamical}}$')
ax.hist(gal_hamass[ind][num[ind] > 3],
        np.linspace(min(gal_hamass),max(gal_hamass),15),alpha=1,lw=2,
        hatch="///",histtype='step',
        label=r'$\mathsf{M_{luminosity}}$',edgecolor='navy')
font={'family' : 'normal',
      'weight' : 'normal',
      'size'   :  20}
ax.legend(loc=1,prop={'size':24,'weight':'bold'})
#axarr[0].xlabel(r'$\mathsf{log(M_{group})}$',fontsize=30)
ax.set_ylabel('Number of Groups',fontsize=25)
ax.set_xlim(11.0,15)
ax.set_title('ECO Sample')
matplotlib.rc('font', **font)
savefig('figures/ECOmasshists.pdf')
show()

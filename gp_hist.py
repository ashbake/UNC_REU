#!/usr/bin/python
#execfile('gp_hist.py')

#plot histograms for galaxies' masses & velocities

#execfile('eco_data/eco_wrapper.py')
execfile('gp_wrapper.py')

blah = 20  #number of bins

plt.subplot(2,1,1)

fof_dmass = fof_dmass - log10(0.7)

plt.hist(fof_dmass[ind][num[ind] >= 3],np.linspace(min(halom),
     max(halom),blah),facecolor='red',edgecolor='None',alpha=.6,lw=2,
     label=r'$\mathsf{M_{dynamical}}$',normed=True)
plt.hist(halom[ind][nold[ind] >= 3],np.linspace(min(halom),
     max(halom),blah),alpha=1,lw=2,hatch="///",histtype='step',
     label=r'$\mathsf{M_{halo}}$',edgecolor='navy',normed=True)
plt.hist(dmass[ind][nold[ind] >= 3],np.linspace(min(halom),
     max(halom),blah),alpha=1,lw=2,hatch="|",histtype='step',
     label=r'$\mathsf{M_{dyn,true}}$',edgecolor='green',normed=True)


xlim(min(fof_dmass[fof_dmass > 0]),max(fof_dmass[fof_dmass > 0]))
ylim(0,1)
ylabel('number')
legend(loc=0)
#title('Mock, Linking Lengths = ' + bz+' , '+bp)

plt.subplot(2,1,2)
plt.hist(fof_dmass[ind][num[ind] >= 3],np.linspace(min(halom),
     max(halom),blah),facecolor='red',edgecolor='None',alpha=.6,lw=2,
     label=r'$\mathsf{M_{dynamical}}$',normed=True)
plt.hist(halom[ind][nold[ind] >= 3],np.linspace(min(halom),
     max(halom),blah),alpha=1,lw=2,hatch="///",histtype='step',
     label=r'$\mathsf{M_{halo}}$',edgecolor='green',normed=True)
execfile('eco_data/eco_wrapper.py')
plt.hist(gal_dmass[ind][gal_dmass[ind] > 0],np.linspace(min(halom),
     max(halom),blah),facecolor='orange',edgecolor='None',alpha=.6,lw=2,
     label=r'$\mathsf{M_{dyn,eco}}$',normed=True)
plt.hist(gal_hamass[ind][gal_dmass[ind] > 0],np.linspace(min(halom),
     max(halom),blah),lw=2,hatch="///",histtype='step',edgecolor='navy',
     label=r'$\mathsf{M_{lum,eco}}$',normed=True)

xlim(min(fof_dmass[fof_dmass > 0]),max(fof_dmass[fof_dmass > 0]))
ylim(0,1)
xlabel('log dynamical masses')
ylabel('number')
legend(loc=0)
title('ECO')


#savefig('plothist_' + bz + '_' + bp + '.pdf')
plt.show()
clf()


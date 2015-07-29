#!/usr/bin/python
#execfile('plot_sample.py')
import matplotlib.gridspec as gridspec
#plot histograms for galaxies' masses & velocities

#execfile('eco_data/eco_wrapper.py')
execfile('gp_wrapper.py')

blah = 20  #number of bins
fig, (ax1, ax2)  = plt.subplots(2, 1,sharex=True)
halom = halom + .15
ax3 = plt.axes([.65, .6, .2, .2])

ax1.hist(halom[ind],np.linspace(min(halom),
     max(halom),blah),alpha=.3,lw=2,histtype='stepfilled',
     label=r'$\mathsf{Mock \ M_{halo}}$',facecolor='red',normed=True)
ax1.hist(array(gal_hamass)[ind],np.linspace(min(halom),
     max(halom),blah),lw=3,histtype='stepfilled',edgecolor='green',
     label=r'$\mathsf{Mock \ M_{HAM}}$',normed=True,facecolor='None')
ax2.hist(dmass[ind][nold[ind] >= 5],np.linspace(min(gd),
     max(gd),blah),alpha=1,lw=3,histtype='step',edgecolor='green',
     label=r'$\mathsf{Mock \ M_{dyn,true}}$',facecolor='None',normed=True)
ax2.hist(fof_dmass[ind][num[ind] >= 5],np.linspace(min(gd),
     max(gd),blah),edgecolor='red',histtype='step',
         alpha=.8,lw=3,label=r'$\mathsf{Mock \  M_{dyn, FoF}}$',normed=True)
ax3.hist(halom[ind],np.linspace(min(halom),
     max(halom),blah),alpha=.3,lw=2,histtype='stepfilled',
     label=r'$\mathsf{Mock \ M_{halo}}$',facecolor='red',normed=True)
ax3.hist(array(gal_hamass)[ind],np.linspace(min(halom),
     max(halom),blah),lw=3,histtype='stepfilled',edgecolor='green',
     label=r'$\mathsf{Mock \ M_{HAM}}$',normed=True,facecolor='None')
ax1.text(8.5,.2,'(a)',fontsize=20)
ax2.set_ylabel('Normalized Frequency')
ax1.set_ylabel('Normalized Frequency')
ax1.set_ylim(0,1.7)

execfile('eco_data/eco_wrapper.py')
gd = gal_dmass[ind][gal_dmass[ind] > 0]
hm = gal_hamass[ind]

ax2.hist(gal_dmass[ind][num[ind] >= 5],np.linspace(min(gd),
     max(gd),blah),histtype='step',facecolor='None',
         edgecolor='navy',alpha=1,lw=1,hatch='/',
         label=r'$\mathsf{ECO \ M_{dyn}}$',normed=True)
ax1.hist(gal_hamass[ind],np.linspace(min(halom),
     max(halom),blah),lw=1,hatch="/",histtype='step',edgecolor='navy',
     label=r'$\mathsf{ECO \ M_{HAM}}$',normed=True)
ax2.hist(gal_hamass[ind][num[ind] >= 5],np.linspace(min(halom),
     max(halom),blah),lw=0,alpha=.6,histtype='stepfilled',facecolor='skyblue',
     label=r'$\mathsf{ECO \ M_{HAM}}$',normed=True,zorder=-100)
ax1.legend(loc=2,prop={'size':14})
ax2.legend(loc=2,prop={'size':14})
ax2.text(8.5,.1,'(b)',fontsize=20)
ax2.set_xlim(8,16)
ax2.set_ylim(0,1)
ax2.set_xlabel('Mass (log $M_\odot$)')
fig.subplots_adjust(hspace=0)
#savefig('plot_sample.pdf')

plt.setp(ax2, yticks=[0,.2,.4,.6,.8])

ax3.hist(hm,np.linspace(min(halom),
     max(halom),blah),lw=1,hatch="/",histtype='step',edgecolor='navy',
     label=r'$\mathsf{ECO \ M_{HAM}}$',normed=True)
plt.setp(ax3, xticks=[12,13,14,15], yticks=[0,.1,.2,.3])
ax3.set_xlim(12.5,15)
ax3.set_ylim(0,.11)

ax1.set_ylabel('Normalized Frequency')
ax1.set_ylim(0,1.7)


#font = {'size'   : 14}
#matplotlib.rc('font', **font)
#
savefig('plot_sample.pdf')
show()

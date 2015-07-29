
execfile('eco_wrapper.py')

#plot group  number versus mass
fig = plt.figure()
ax = plt.gca()
cm = ax.scatter(gal_dmass,num,c=gal_hamass,edgecolor='none')
ax.plot([8,15.5],[20,20])
ax.set_yscale('log')
xlim(8,15.5)
xlabel('Dynamical Mass')
ylabel(r'$\mathsf{N_{grp}}$')
font={'size'   :  20}
matplotlib.rc('font', **font)
cb = plt.colorbar(cm)
cb.set_label("HAM Masses",rotation=270)
text(8.5,450,'ECO Groups')
savefig('numvsmass.pdf')
show()

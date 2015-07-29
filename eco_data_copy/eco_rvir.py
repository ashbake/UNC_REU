#!/usr/bin/python
#execfile('eco_rvir.py')

execfile('eco_wrapper.py')

rho0            = .7*.7*.75*2.775* 10**11 #h^2 Msun/Mpc^3
halo_rvir       = ((10**gal_hamass)*(3./(rho0*200*4.*np.pi)))**(1./3.)
dm_rvir         = (((10.**(gal_dmass))*(3./(rho0*200*4.*np.pi))))**(1./3.)
groupn = 5
big = where(num[ind] > groupn)

x=[0,5]
scatter(rvir,r200,c='gray',marker='o')
plot(x,x)

ylim(0,1.2)
xlim(0,1.2)
xlabel('$\mathsf{R_{vir}}$')
ylabel('$\mathsf{R_{200}}$')
title('Comparison of $\mathsf{R_{200}}$ and $\mathsf{R_{vir}}$')


matplotlib.rcParams.update({'font.size': 16})

savefig('figures/plotvir.pdf')

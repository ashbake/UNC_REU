#to run AM code:
#run_AM_code,filename,bxy,bz,volume
#14011.0 Mpc^3
#4805.773 mpc3/h3
#massmatch3 MassFunction.dat 207034.45 < G_infof_ecospring.txt  > out_hamass_ecospring_.12_1.3.txt 
#execfile('AM_plots.py')

execfile('../gp_wrapper.py')

filename = 'out_Resolve_5001_hod1_mock_zspace.dat'
file=np.loadtxt(filename)
hamass = file[:,1]
hamr200 = file[:,2]
gal_hamass = [0] * (len(ind) + 1)
for i in range(0,len(ind)-1):
    n = ind[i+1] - ind[i]
    gal_hamass[ind[i]:ind[i+1]] = [hamass[i]] * n

gal_hamass.extend([hamass[-1]])

G = 4.302*10**-9 #;pc Msun-1 (km/s)**2
rho0 = .25*(3*(100**2)/ (8*np.pi * G))  #;H0 in km/s/Mpc
rho0      = 1.477 * 10**11 #h^2 Msun/Mpc^3
r200ham      = ((10**(hamass))*(3./(rho0*200*4.*np.pi)))**(1./3.)

f, axarr  = plt.subplots(2, 2)
x = [7,14.5]
y = [10,15]
onetoone = [10,16]
newcol = nold
newcol[newcol > 50] = 50
#A~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
axarr[0,0].scatter(dmass[num > 2],halom[num > 2],c = log10(newcol[num > 2]),
                   edgecolor='none', marker='.',zorder = -200)#,cmap=mpl.cm.gray)
axarr[0,0].plot(onetoone,onetoone,'g-',zorder=-10)
axarr[0,0].set_xlabel('Dyn Mass "True"')
axarr[0,0].set_ylabel('Original Halo Mass')
axarr[0,0].set_ylim(10,15)
axarr[0,0].set_xlim(8,15)
out = BinFun(dmass[dmass > 0],halom[dmass > 0],7)
axarr[0,0].errorbar(out.xout,out.yout,yerr=out.errs,
                    label='halo mass dispersion',barsabove=True,
                    ls='None', capsize=4.5, color='m',zorder = 300,
                    elinewidth=2)
axarr[0,0].plot(out.xout,out.yout,'ro',mec='w')
axarr[0,0].text(9,14.5,'a)',fontsize=16,horizontalalignment='center',
     verticalalignment='center')
#C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
axarr[1,0].scatter(fof_dmass,halom,c=log10(newcol),edgecolor='none',
                   marker='.',zorder=-200)# cmap=mpl.cm.gray)
axarr[1,0].plot(onetoone,onetoone,'g-',zorder=-10)
axarr[1,0].set_xlabel('Dynamical Mass')
axarr[1,0].set_ylabel('Original Halo Mass')
axarr[1,0].set_ylim(10,15)
axarr[1,0].set_xlim(8,15)
out = BinFun(fof_dmass[where((fof_dmass > 12) & (nold > 3) & (halom > 12))],
             halom[where((fof_dmass > 12) & (nold > 3) & (halom > 12))],5)
axarr[1,0].errorbar(out.xout,out.yout,yerr=out.errs,
                    label='halo mass dispersion',barsabove=True,
                    ls='None', capsize=4.5, color='m',zorder = 300,
                    elinewidth=2)
axarr[1,0].plot(out.xout,out.yout,'ro',mec='w')
axarr[1,0].text(9,14.5,'c)',fontsize=16,horizontalalignment='center',
     verticalalignment='center')
#D~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
img = axarr[1,1].scatter(array(gal_hamass)[num > 2],halom[num > 2], marker='.',
                   c=log10(newcol[num > 2]),edgecolor='none',zorder=-200)#, cmap=mpl.cm.gray)
axarr[1,1].plot(onetoone,onetoone,'g-',zorder=-10)
axarr[1,1].set_xlabel('HAM Mass')
axarr[1,1].set_ylabel('Original Halo Mass')
axarr[1,1].set_ylim(10,15)
axarr[1,1].set_xlim(8,15)
out = BinFun(array(gal_hamass)[where((array(gal_hamass) > 12.0) & (nold > 2) & 
                                     (halom > 12))],
             halom[where((array(gal_hamass) > 12.0) & (nold > 2) & 
                         (halom > 12))],5)

axarr[1,1].arrow(9,11,0.0,.2, head_width=0.1, head_length=0.1, fc='k', ec='k')
axarr[1,1].arrow(9,11.2,0.0,-.2, head_width=0.1, head_length=0.1, fc='k', ec='k')

axarr[1,1].errorbar(out.xout,out.yout,yerr=out.errs,
                    label='halo mass dispersion',barsabove=True,
                    ls='None', capsize=4.5, color='m',zorder = 300,
                    elinewidth=2)
axarr[1,1].plot(out.xout,out.yout,'ro',mec='w',zorder=400)
axarr[1,1].set_xticks([8,9,10,11, 12, 13,14,15])
axarr[1,1].text(9,14.5,'d)',fontsize=16,horizontalalignment='center',
     verticalalignment='center')
#B fake cartoon~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
axarr[0,1].plot(onetoone,onetoone,'g-')
#axarr[0,1].plot(onetoone,onetoone,'w--')
#axarr[0,1].scatter(array(gal_hamass)[num > 3],array(gal_hamass)[num > 3], marker='.',
#                   c=log10(newcol[num > 3]),edgecolor='none',zorder=-200)#, cmap=mpl.cm.gray)

file1 = 'out_forashley.dat'
file2 = 'noG_forashley.dat'
data1 = np.loadtxt(file1)
data2 = np.loadtxt(file2)
hamass=data1[:,1][35:269488]-log10(0.70)
halomm=data2[:,1]-log10(0.70)

axarr[0,1].plot(hamass[halomm > 11],halomm[halomm > 11],'.',c='gray',markeredgecolor='gray',zorder=-100,alpha=.1)

out = BinFun(hamass,halomm,15)

msun = [12,12.5,13,13.5,14,14.5] - log10(.7) #corresponding errors
masserrslow= [.2,.22,.22,.25,.32,.4]
masserrshi=[.15,.1,.09,.08,.08,.15]
#axarr[0,1].fill_between(msun, msun+masserrshi, msun-masserrslow,color='navy',zorder=-10,alpha=.2)

axarr[0,1].errorbar(out.xout[0::3],out.yout[0::3],yerr=out.errs[0::3],
                    barsabove=True,
                    ls='None', capsize=4.5, color='k',zorder = 300,
                    elinewidth=2)

axarr[0,1].arrow(9,11,0.0,.2, head_width=0.1, head_length=0.1, fc='k', ec='k')
axarr[0,1].arrow(9,11.2,0.0,-.2, head_width=0.1, head_length=0.1, fc='k', ec='k')

axarr[0,1].set_xlim(8,15)
axarr[0,1].set_ylim(10,15)
axarr[0,1].plot(out.xout[0::3],out.yout[0::3],'ro',mec='w',zorder=400)
axarr[0,1].set_xticks([8,9,10,11, 12, 13,14,15])
axarr[0,1].text(9,14.5,'b)',fontsize=16,horizontalalignment='center',
     verticalalignment='center')
axarr[0,1].set_xlabel('HAM Masses "True"')
axarr[0,1].set_ylabel('Original Halo Mass')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

f.subplots_adjust(right=0.8)
cbaxes = f.add_axes([0.85, 0.15, 0.05, 0.7])
#cb=colorbar()
#cb.set_label("True ($\mathsf{N_{grp}}$)")
cbar = f.colorbar(img,cax = cbaxes)
cbar.set_label('True $\mathsf{N_{grp}}$', rotation=270)
cbar.ax.set_yticklabels([' 1',' 2',' 3',' 4',' 6','10','16','25','>40'])
#cbar.ax.get_yaxis().set_ticks([])
#for j, lab in enumerate(['$1$','$4$','$12$','$>50$']):
#    cbar.ax.text(.5, (2 * j + 1) / 8.0, lab, ha='right', va='center')

plt.savefig('compare_masses_new.pdf')
plt.show()
#f.gca().xaxis.set_major_locator(MaxNLocator(prune='lower'))
#plt.show()

#~~~~~~~~~~~~~~~~~~~~~ECO Comp~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plt.figure()
plt.scatter(array(gal_hamass),fof_dmass,
                       c=log10(newcol), marker='.',edgecolors='none',
                       zorder=-200,label=r'$\mathsf{3 \ \leq N_{grp} \ \leq 5}$')#, cmap=mpl.cm.gray)
plt.scatter(array(gal_hamass)[num > 5],fof_dmass[num > 5],
                       c=log10(newcol)[num > 5], marker='^',edgecolors='none',
                       zorder=-200,s=30,label=r'$\mathsf{N_{grp} \ > \ 5}$')#, cmap=mpl.cm.gray)
plot(onetoone,onetoone,'m-',zorder = -10)
xlabel('HAM Masses')
ylabel('Dynamical Mass')
ylim(9,15)
xlim(11,15)
out = BinFun(array(gal_hamass)[num > 3],fof_dmass[num > 3],7)
plt.errorbar(out.xout,out.yout,yerr=out.errs,barsabove=True,
                    ls='None', capsize=4.5, color='r',zorder = 300,
                    elinewidth=3)
plt.xticks([8,9,10,11, 12, 13,14,15])
plt.plot(out.xout,out.yout,'ro',zorder=400,mec='w',ms=10)
text(12,
14.2,'(b) Mock Catalog',fontsize=18,horizontalalignment='center',
     verticalalignment='center')
xlim(11,15)
legend(loc = 4)
savefig('hamdmass_comp.pdf')
show()

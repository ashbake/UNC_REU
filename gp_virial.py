#!/usr/bin/python
#execfile('gp_virial.py')

execfile('gp_wrapper.py')

rho0            = 1.477 * 10**11 #h^2 Msun/Mpc^3
halo_rvir       = ((10**halom)*(3./(rho0*200*4.*np.pi)))**(1./3.)
dm_rvir         = (((10.**(0.8611+fof_dmass))*(3./(rho0*200*4.*np.pi))))**(1./3.)
groupn = 5
big = where(num[ind] > groupn)
#plot(halo_rvir[ind][num[ind] > groupn],dm_rvir[ind][num[ind] > groupn],'k.')
plot(rad[ind][num[ind] > groupn],(halo_rvir[ind][num[ind] > groupn]),'k.',mec='k')
plot(fof_rad[ind][num[ind] > groupn],(dm_rvir[ind][num[ind] > groupn]),'r.',mec='r')
#plot(rad[ind][num[ind] > groupn],(halo_rvir[ind][num[ind] > groupn])/(rad[ind][num[ind] > groupn]),'r.',mec='r')
#ylim(0,.7)
xlim(0,.4)
x=[0,.1,.2,.3,.4,.5,.6,.7]
plot(x,array(x)+.2,'g')
ylabel('R_vir according to True Halo Mass')
xlabel('Average radius for true group')
title('Constant to add to  avg radius')
savefig('compare_rvir')





plot(rad,rvirgd,'k.',mec='k')
#plot(rad,rvirgd,'r.',mec='r')
plot(fof_rad[ind],rvir[ind],'m.',mec='m')
#plot(fof_rad,rvirgd,'r.',mec='r')
#y_intgd = GetIntercept2(rad[rad > 0],rvirgd[rad > 0])
#bbgd = y_intgd.b

#y_int = GetIntercept2(fof_rad[fof_rad[ind] > 0],rvir[ind][fof_rad[ind] > 0])
#bb = y_int.b

x=[0,4]
y = [0,4*2]#1.27]
ygd= [0,4*1.34]
plot(x,y)
plot(x,ygd)
xlim(0,.8)
ylim(0,.8)
xlabel('Projected Group Average Radii')
ylabel('Group Virial Radius')
show()

ymass = GetIntercept(dmass[where((nold > 8) & (dmass > 0))],halom[where((nold > 8) & (dmass > 0))])
yfofmass = GetIntercept(fof_dmass[num > 8],halom[num > 8])

ymass.b
yfofmass.b

class GetIntercept2:
	def __init__(self,xi,yi):
		j = 1
		while j < 10:
			guess = np.arange(0,2,.01)
			chisquare=[0] * len(guess)
			for i in range(0,len(guess)):
				expected = xi*guess[i]
				chisquare[i] = sum(((yi - expected)**2))
			b = guess[chisquare == min(chisquare)]
			j = j*2
		self.b = b
		self.chisquare = min(chisquare)


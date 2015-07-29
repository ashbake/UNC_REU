#!/usr/bin/python
#execfile('gp_find_shift.py')

#find shift of dynmass vs halo mass/matched mass plots before
#and after group finder.

#execfile('gp_wrapper.py')

n=12

#NEED OLD IND BC CANT HAVE REPEATS WHEN DO X^2	

#-----------------------------------------------------
#               SHIFT OF NO FOF9/TRUE DMASSES
#-----------------------------------------------------

#fix slope to be 1
#minimize chi^2 iterating over guesses
class GetIntercept:
	def __init__(self,xi,yi):
		j = 1
		b = mean(yi - xi)
		while j < 10:
			guess = np.arange(b-1./(j**2),b+1./(j**2),.1/(j**4))
			chisquare=[0] * len(guess)
			for i in range(0,len(guess)):
				expected = xi + guess[i]
				chisquare[i] = sum(((yi - expected)**2))
			b = guess[chisquare == min(chisquare)]
			j = j*2
		self.b = b
		self.chisquare = chisquare


# - 1 - find shift for original dyn mass before group finder	
xi = dmass[(nold > n) & (dmass > 0)]  #some are -inf bc dmass=0 if group was cut off
yi = halom[(nold > n) & (dmass > 0)]

y_int = GetIntercept(xi,yi)
b_dmass = y_int.b
print 'b_dmass = ', b_dmass[0]
print 'min chi^2 =' , min(y_int.chisquare)

#plot
y=x = np.arange(min(xi), max(xi), 0.5);
plot(dmass,halom-b_dmass,'y.')
plot(xi,yi-b_dmass,'r.',x,y,'g-',label='1 to 1 line')
xlabel('perfect dyn mass')
ylabel('halo mass')
legend()
title('shift =' + str(b_dmass))
savefig('test_shift')
show()


# ~ 2 ~ find shift for halo mass matched dmass after group finder
xi = fof_dmass[(num > n) & (fof_dmass > 0) & (halom > 11.9)]  #some are -inf bc dmass=0 if group was cut off
yi = halom[(num > n) & (fof_dmass > 0) & (halom > 11.9)]

y_int = GetIntercept(xi,yi)
b_fof_dmass = y_int.b
print 'b_fof_dmass = ', b_fof_dmass[0]
print 'min chi^2 =' , min(y_int.chisquare)

#plot
y=x = np.arange(min(xi), max(xi), 0.5);
#plot(fof_dmass,halom-b_fof_dmass,'y.')
#plot(xi,yi-b_fof_dmass,'r.',x,y,'g-',label='1 to 1 line')
#xlabel('dyn mass after fof')
#ylabel('halo mass')
#title('shift =' + str(b_fof_dmass))
#savefig('test_shift2')
#show()


# ~ 3 ~ find shift for original dyn mass before group finder	
xi = fof_dmass[ind][(num[ind] > n) & (fof_dmass[ind] > 0)]  #some are -inf bc dmass=0 if group was cut off
yi = match_mass[ind][(num[ind] > n) & (fof_dmass[ind] > 0)]

y_int = GetIntercept(xi,yi)
b_matchmass = y_int.b
#print 'b_matchmass = ', b_matchmass[0]
#print 'min chi^2 =' , min(y_int.chisquare)

#plot
y = x = np.arange(min(xi), max(xi), 0.5);
#plot(fof_dmass,match_mass-b_matchmass,'y.')
#plot(xi,yi-b_matchmass,'r.',x,y,'g-',label='1 to 1 line')
#xlabel('dyn mass after fof')
#ylabel('halo mass')
#title('shift =' + str(b_matchmass))

#show()


#find best line
#dont fix slope at 1
A = array([xi, ones(len(xi))])
w = linalg.lstsq(A.T,yi)[0] # obtaining the parameters

line = w[0]*xi+w[1] # regression line
x=y= np.arange(12.5, 15, 0.5);
#plot(xi,line,'r-',xi,yi,'o',x,y,'g')
x=y= np.arange(12.5, 15, 0.5);
#show()

#plot(x,y+b_matchmass,label='matched mass shift')
#plot(x,y+b_fof_dmass,label='fof')




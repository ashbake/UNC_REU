#!/usr/bin/python
#execfile('gp_findcuts.py')
#pretty sure this plot only is meaningful if all points included bc if only select
#bad points, i'm shading based on how far points are from 1:1 line anyways, so 
#it won't be evident in final plot

execfile('gp_wrapper.py')
pylab.ioff()
xx = (galr/fof_rad)[num > 3]
yy = abs(galv/fof_vel)[num > 3]
fof_dmass = fof_dmass[num > 3]
halom = halom[num > 3]


#outliers = (halom > fof_dmass + 0.3)
#outliers1 = (halom < 11.5)
outliers2 = halom < fof_dmass - 0.7
i_out = outliers2

#select bad points
bad = nold[xx > 0] == 1
good = nold > 1


class BinXY:
    'given x and y, bin along x and find dispersion in y'
    def __init__(self,xin,yin,nbin):
    	binsizex  = (max(xin) - min(xin))/nbin
    	binsizey  = (max(yin) - min(yin))/nbin
    	xout	  = zeros(nbin**2).reshape((nbin,nbin))
    	yout	  = zeros(nbin**2).reshape((nbin,nbin))
        shade     = zeros(nbin**2).reshape((nbin,nbin))
        number    = zeros(nbin**2).reshape((nbin,nbin))
        for cake  in range(0,nbin):
        	bin_x       =  xin[(xin > (binsizex*cake + min(xin))) & (xin <= (binsizex*(cake + 1)+min(xin)))]
        	xout[cake,:]  =  binsizex*cake + min(xin)
        	ind_xbin    =  (xin > (binsizex*cake + min(xin))) & (xin <= (binsizex*(cake + 1)+min(xin)))
        	temp_y      =  yin[ind_xbin]
        	for pizza in range(0,nbin):
        		ind_ybin    = (temp_y > (binsizey*pizza + min(yin))) & (temp_y <= (binsizey*(pizza + 1)+min(yin)))
        		bin_y       =  temp_y[ind_ybin]
        		yout[cake,pizza] =  binsizey*pizza + min(yin)
        		hm_bin		=  halom[fof_dmass > 0][ind_xbin][ind_ybin]
        		dm_bin		=  fof_dmass[fof_dmass > 0][ind_xbin][ind_ybin]
        		num_bin		=  num[fof_dmass > 0][ind_xbin][ind_ybin]
        		if len(hm_bin) == 0:
        			shade[cake,pizza] = 0
        			number[cake,pizza] = 0
        		else:
        			number[cake,pizza] = mean(num_bin[np.isfinite(num_bin)])
        			shade[cake,pizza] = mean(abs(hm_bin[np.isfinite(dm_bin)] - dm_bin[np.isfinite(dm_bin)]))
        self.yout = yout
        self.xout = xout
        self.shade = shade
        self.number = number
        self.binsizey = binsizey
        self.binsizex = binsizex

xx = xx[xx > 0]
yy = yy[yy > 0]
nbin = 10
out = BinXY(xx,yy,nbin)
yout = out.yout
xout = out.xout
shade = out.shade
number = out.number
binsizex = out.binsizex
binsizey = out.binsizey

shadenew=[0]*nbin**2
xfg = [0]*nbin**2
yfg = [0]*nbin**2
nfg = [0] * nbin**2

for i in range(0,nbin):
	shadenew[i*nbin:nbin*i+nbin] = shade[i,:]
	xfg[i*nbin:nbin*i+nbin] = xout[i,:]
	yfg[i*nbin:nbin*i+nbin]= yout[i,:]
	nfg[i*nbin:nbin*i+nbin]= number[i,:]
	
#DataOut = column_stack((xfg,yfg,shadenew,nfg))  #save file of info for filtergraph
#savetxt('gridplot.dat', DataOut)

colorscale = (shade - min(shade[np.isfinite(shade)]))/(max(shade[np.isfinite(shade)]) - min(shade[np.isfinite(shade)]))
fig, ax = plt.subplots(1)

plt.subplot(2,1,1)

for i in range(0,nbin):
	for j in range(0,nbin):
		plot(xout[i,j],yout[i,j],alpha=colorscale[i,j],ms=27,marker='s',c='b')

#plot(xout,yout,'k.')
xlim(0-binsizex/2,xout[nbin-1,nbin-1]+binsizex/2)
ylim(0-binsizey/2,yout[nbin-1,nbin-1]+binsizey/2)
xlabel('galaxy radius / group radius')
ylabel('abs(cz gal)/ group vel dispersion')
savefig('sigrad_dispersion.png')

plt.subplot(2,1,2)
plot(xx,yy,'r.',markeredgecolor='r',alpha=.3)
#plot(xx[bad],yy[bad],'b.')
savefig('figures/sigrad_dispersion.png')
show()




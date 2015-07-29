#!/usr/bin/python
#execfile('gp_PDF.py')

#get probability distribution for galaxies' masses

execfile('gp_wrapper.py')

class BinFun2:
    'given x and y, bin along x and find dispersion in y'
    def __init__(self,yin,nbin):
        binsize   = (max(yin) - min(yin))/nbin
        ymid      = [0] * nbin
        size      = [0] * nbin
        yin.sort()
        for pizza in range(0,nbin):
            temp_y       =  yin[(yin > (binsize*pizza + min(yin))) 
            					& (yin <= (binsize*(pizza + 1)+min(yin)))] 
            size[pizza]  =  len(temp_y)
            ymid[pizza]  =  (binsize*pizza + binsize*(pizza+1))/2 + min(yin)
            print 'pizza binfun=' ,pizza
        self.ymid = ymid
        self.size = size


class YDist:
    'given x range and corresponding y-array, outputs histogram of y-distribution calling BinFun2'
    def __init__(self,xin,yin,nbin):
        binsize   = (max(xin) - min(xin))/(nbin)
        xout      = [0] * nbin      #store middle of bin
        yout      = [0] * nbin      #store avg halo_m for bin
        size      = [0] * nbin
        for pizza in range(0,nbin):
            temp_y       =  yin[(xin > (binsize*pizza + min(xin))) & (xin <= (binsize*(pizza + 1)+min(xin)))] 
            xout[pizza]  =  (binsize*pizza + binsize*(pizza+1))/2 + min(xin)
            output       =  BinFun2(temp_y,15)
            yout[pizza]  =  [output.ymid]
            size[pizza]  =  [output.size]
            print 'pizza ydist =',pizza
        self.xout = xout
        self.yout = yout
        self.size = size


#PDF for matched mass plot
xin = fof_dmass[ind][fof_dmass[ind] > 0]
yin = match_mass[ind][fof_dmass[ind] > 0]

xin = dmass[ind][dmass[ind] > 0]
yin = halom[ind][dmass[ind] > 0]

xin = fof_dmass[num > 3][fof_dmass[num > 3] > 0]
yin = halom[num > 3][fof_dmass[num > 3] > 0]

plot(fof_dmass[fof_dmass > 0],halom[fof_dmass > 0],'k.')
plot(xin,yin,'y.')
show()

out = YDist(xin,yin,4)

xout = out.xout
size = out.size
yout = out.yout

#for i in range(0,3): #len(xout)):#
	#x = array(yout[9-2*i][0]#)
	#y = array(size[9-2*i][0])
	#width = x[4] - x[3]
#	colors = ['r','g','b','y','m']
#	plt.bar(x,y,width,color=colors[i],alpha=.7,label='xslice ='+str(round(xout[9-2*i],3)))
#	xlim(min(yin),max(yin))
#	legend()
#	print width
#show()

plt.subplot(4,1,1)
plt.bar(yout[0][0],size[0][0],width=(yout[0][0][2]-yout[0][0][1]),color='r',alpha=.7,label='xslice ='+str(round(xout[0],3)))
title('xslice ='+str(round(xout[0],3)))
xlim(10,16)
plt.subplot(4,1,2)
plt.bar(yout[1][0],size[1][0],width=(yout[1][0][2]-yout[1][0][1]),color='g',alpha=.7,label='xslice ='+str(round(xout[1],3)))
title('xslice ='+str(round(xout[1],3)))
xlim(10,16)
plt.subplot(4,1,3)
plt.bar(yout[2][0],size[2][0],width=(yout[2][0][2]-yout[2][0][1]),color='y',alpha=.7,label='xslice ='+str(round(xout[2],3)))
title('xslice ='+str(round(xout[2],3)))
xlim(10,16)
plt.subplot(4,1,4)
plt.bar(yout[3][0],size[3][0],width=(yout[3][0][2]-yout[3][0][1]),color='b',alpha=.7,label='xslice ='+str(round(xout[3],3)))
title('xslice ='+str(round(xout[3],3)))
xlim(10,16)
savefig('histogramplot.png')
show()
#find standard deviations

class Deviation:
    'find sigmas'
    def __init__(self,halo):
        halo.sort()
        N = len(halo)
        if N > 16:
        	self.ymed  = np.median(halo)
        	self.y1slo = halo[round(0.16*N)]
        	self.y1shi = halo[round(0.84*N)]
        	self.y2slo = halo[round(0.025*N)]
        	self.y2shi = halo[round(0.975*N)]
        else:
        	self.ymed  = 0
        	self.y1slo = 0
        	self.y1shi = 0
        	self.y2slo = 0
        	self.y2shi = 0
        	
        	
class YDist2:
    'given x range and corresponding y-array, outputs histogram of y-distribution calling BinFun2'
    def __init__(self,xin,yin,nbin):
        binsize   = (max(xin) - min(xin))/(nbin)
        xout = [0] * nbin 
        ymed = [0] * nbin 
        y1slo = [0] * nbin 
        y1shi = [0] * nbin 
        y2slo = [0] * nbin 
        y2shi = [0] * nbin      #store avg halo_m for bin
        for pizza in range(0,nbin):
            temp_y       =  yin[(xin > (binsize*pizza + min(xin))) & (xin <= (binsize*(pizza + 1)+min(xin)))] 
            xout[pizza]  =  (binsize*pizza + binsize*(pizza+1))/2 + min(xin)
            output       =  Deviation(temp_y)
            ymed[pizza]  =  output.ymed
            y1slo[pizza] =output.y1slo
            y1shi[pizza] =output.y1shi
            y2slo[pizza] =output.y2slo
            y2shi[pizza] =output.y2shi
            print 'pizza ydist =',pizza
        self.xout  = array(xout)
        self.ymed  = array(ymed)
        self.y1slo = array(y1slo)
        self.y1shi = array(y1shi)
        self.y2slo = array(y2slo)
        self.y2shi = array(y2shi)
      
xin = fof_dmass[ind][fof_dmass[ind] > 0]
yin = match_mass[ind][fof_dmass[ind] > 0]

xin = dmass[ind][dmass[ind] > 0]
yin = halom[ind][dmass[ind] > 0]

xin = fof_dmass[ind][fof_dmass[ind] > 0]
yin = halom[ind][fof_dmass[ind] > 0]



outer = YDist2(xin,yin,8)
xout = outer.xout
ymed = array(outer.ymed)

fig, ax = plt.subplots(1)#,2,1)
   
ylim(10.5,15)
xlim(9,15)
plot(array(xout)[ymed > 0],ymed[ymed > 0],'k.',alpha=.3)
plot(xout[ymed > 0],outer.y1slo[ymed > 0],'b-',label='1 sigma',alpha=.5)
plot(xout[ymed > 0],outer.y1shi[ymed > 0],'b-',alpha=.5)
plot(xout[ymed > 0],outer.y2slo[ymed > 0],'b-',label='2 sigma',alpha=.3)
plot(xout[ymed > 0],outer.y2shi[ymed > 0],'b-',alpha=.3)
legend(numpoints = 1,loc=2)
xlabel('dynamical mass')
ylabel('matched halo mass')
ax.plot(xin,yin,'k.',alpha=.5)
ax.fill_between(xout[ymed > 0], outer.y1shi[ymed > 0], outer.y1slo[ymed > 0], facecolor='blue',alpha=.5)
ax.fill_between(xout[ymed > 0], outer.y2shi[ymed > 0], outer.y2slo[ymed > 0], facecolor='blue',alpha=0.3,zorder=-1)
savefig('fig3.png')

show()

#plt.subplot(1,2,1)

plot(xin,yin,'y.')    
ylim(10,15)
xlim(8,14)
plot(array(xout)[ymed > 0],ymed[ymed > 0],'k.')

fig, ax = plt.subplots(1)
ax.plot(xin,yin,'k.',alpha=.5)
ax.fill_between(xout[ymed > 0], outer.y1shi[ymed > 0], outer.y1slo[ymed > 0], facecolor='blue',alpha=.5)
ax.fill_between(xout[ymed > 0], outer.y2shi[ymed > 0], outer.y2slo[ymed > 0], facecolor='blue',alpha=0.3,zorder=-1)

show()


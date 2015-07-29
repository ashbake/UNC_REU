#!/usr/bin/python
#execfile('AS_test.py')
#fit line, rotate, take ratios of sigmas of histograms

execfile('gp_wrapper.py')

for i in range(0,len(ind)-1):
    n=int(num[ind[i]])
    #size = round(n/4)
    size=11
    if n > size:
        czs = cz[ind[i]:ind[i]+n]
        ras = ra[ind[i]:ind[i]+n]
        decs = dec[ind[i]:ind[i]+n]
        coords_rot=[0.0] * len(ras)
        #fit line
        fit = np.polyfit(ras, decs, 1, rcond=None, full=False)
        x = [min(ras),max(ras)]
        x = array(x)
        y = fit[0]*x +fit[1]
        #rotate points around line
        angle=(math.atan((y[1] - y[0])/(x[1] - x[0])))
        rotMatrix = array([[cos(angle), -sin(angle)], 
                   [sin(angle),  cos(angle)]])
        coords = zip(*[ras,decs])
        for j in range(0,len(ras)):
            coords_rot[j] = np.dot(coords[j],rotMatrix)
        ras_rot = zip(*coords_rot)[0]
        decs_rot = zip(*coords_rot)[1]
        #plot
        plot(ras,decs,'ro')
        plot(x,y,'k-')
        plot(ras_rot,decs_rot,'bo')
        xr = [min(ras_rot),max(ras_rot)]
        xr = array(xr)
        fit = np.polyfit(ras_rot,decs_rot, 1, rcond=None, full=False)
        yr = fit[0]*xr +fit[1]
        plot(xr,yr,'k-')
        show()

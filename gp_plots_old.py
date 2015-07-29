#!/usr/bin/python
#run gp_dynmass.py first
#plots shit
#execfile('gp_dynmass.py')
#execfile('gp_plots.py')
truem = gal_halom[gal_dmass != 0]
dynm  = log10(gal_dmass[gal_dmass != 0])
vel   = gal_vel[gal_dmass != 0]
rad   = gal_rad[gal_dmass != 0]
radarr= array(radarr)
nus   = num[gal_dmass != 0]


#LOAD ORIGINAL DYNMASS, VEL, & RAD#
filename="/home/bakerad/REU/dynm_out_avg.dat"
file=np.loadtxt(filename)
gd_mass = file[:,0][list(of)] #absolute R band mag
gd_vels = file[:,1][list(of)] 
gd_rads = file[:,2][list(of)]
true_mass = log10(gd_mass[gal_dmass != 0])
true_vel  = gd_vels[gal_dmass != 0]
true_rad  = gd_rads[gal_dmass != 0]

#~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
#~	~	~	~	~	~	PLOTTING	~	~	~	~	~	~
#~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
class BinFun:
    'given x and y, bin along x and find dispersion in y'
    def __init__(self,xin,yin,nbin):
        binsize   = (max(xin) - min(xin))/nbin
        xout      = [0] * nbin      #store middle of bin
        yout      = [0] * nbin      #store avg halo_m for bin
        errs      = [0] * nbin
        for pizza in range(0,nbin):
            temp_y       =  yin[(xin > (binsize*pizza + min(xin))) & (xin < (binsize*(pizza + 1)+min(xin)))] 
            errs[pizza]  =  sqrt(sum(abs(temp_y - temp_y.mean())**2)/(len(temp_y)-1))
            xout[pizza]  =  (binsize*pizza + binsize*(pizza+1))/2 + min(xin)
            yout[pizza]  =  mean(temp_y)
        self.xout = xout
        self.yout = yout
        self.errs = errs


def plotmass():
    "plot halo mass vs dynamical mass"
    plot(fof_dmass[fof_dmass > 0],halom[fof_dmass > 0],'y.',zorder=-100,label='$N_\mathrm{galaxies}$ = 3')
    plot(fof_dmass[num > 5],halom[num > 5],'b.',zorder=-100,label='$N_\mathrm{galaxies}$ = 3')
    #plot(log10(gal_dmass[gal_dmass != 0]),gal_halom[gal_dmass != 0],'y.',zorder=-100,label='$N_\mathrm{galaxies}$ = 3')
    #plot(log10(gal_dmass[num > 3]),gal_halom[num > 3],'g.', zorder=-100,label='$N_\mathrm{galaxies}$ = 4')                            #overplot halos w/ 4 members & up
    #plot(log10(gal_dmass[num > 5]),gal_halom[num > 5],'r.', zorder=-100,label='$N_\mathrm{galaxies}$ >= 5')                            #overplot halos w/ 4 members & up
    out = BinFun(dynm,truem,15)
    #plt.errorbar(out.xout,out.yout,yerr=out.errs,label='halo mass dispersion',barsabove=True,ls='None', capsize=4.5, color='black',zorder = 300,elinewidth=2.5)
    #plot(out.xout,out.yout,'ko')
    y = x = np.arange(0, 15, 0.5);
    plot(x,y,color='r',label='1 to 1 line',zorder=100)        #plot 1:1 line
    ylim(10,15)
    xlim(9,14.5)
    xlabel('log Dynamical Mass')
    ylabel('log Halo Mass')
    title('avg / better lnking lengths')
    legend(loc=2,numpoints=1)
    savefig('gp_galmass_old.png')
    show()



def newplot():
    'plot each galaxy with old dynmass found'
    
    x = np.arange(7, 15, 0.5);
    y = x #+ .7
    
    #plot(x,y,color='r',label='1 to 1 line')        #plot 1:1 line
    plot(fof_dmass,dmass,'y.',zorder = -100,label='3 < $N_\mathrm{galaxies}$ <= 8')
    plot(fof_dmass[num > 8],dmass[num > 8],'b.',label='$N_\mathrm{galaxies}$ > 8')
    #plot(dynm[nus > 15],true_mass[nus > 15],'r.',label='$N_\mathrm{galaxies}$ > 15')
    #oddities = ((dynm/truem) > 1) & (truem < 11.4)
    #plot(dynm[oddities],true_mass[oddities],'ko')
    #xx = BinFun(dynm[true_mass > 0],true_mass[true_mass > 0],15)
    #plt.errorbar(xx.xout,xx.yout,yerr=xx.errs,xerr=None,barsabove=True,ls='None', capsize=4.5, color='black',zorder = 300,elinewidth=2.5)
    #plot(xx.xout,xx.yout,'ko')
    xlabel('dynamical mass post group finder')
    ylabel('dynamical mass')
    ylim(8,15)
    xlim(8,15)
    title('Change in Dyn Mass Due to Group Finder')
    legend(loc=2,numpoints=1)
    savefig('newmassplot_old.png')
    show()
    
    plot(rad,true_rad,'y.',zorder=-100)
    plot(rad[nus > 10],true_rad[nus > 10],'b.',label='$N_\mathrm{galaxies}$ > 10')
    x = np.arange(min(rad), max(rad), 0.5);
    y = x
    plot(x,y,color='r',label='1 to 1 line')
    xx = BinFun(rad[true_rad > 0],true_rad[true_rad>0],15)
    #plt.errorbar(xx.xout,xx.yout,yerr=xx.errs,xerr=None,barsabove=True,ls='None', capsize=4.5, color='black',zorder = 300,elinewidth=2.5)
    xlabel('avg radius after group finder')
    ylabel('true avg radius')
    title('central method/ worse linking lengths/ radius')
    ylim(0,1.7)
    xlim(0,0.8)
    legend(loc=2,numpoints=1)
    savefig('radplot_old.png')
    show()

    x = np.arange(min(vel), max(vel), 0.5);
    y = x
    plot(vel,true_vel,'y.')
    plot(vel[nus > 10],true_vel[nus > 10],'g.',label='$N_\mathrm{galaxies}$ > 10')
    plot(x,y,color='r',label='1 to 1 line')
    xx = BinFun(rad[true_rad > 0],true_rad[true_rad>0],15)
    plt.errorbar(xx.xout,xx.yout,yerr=xx.errs,xerr=None,barsabove=True,ls='None', capsize=4.5, color='black',zorder = 300,elinewidth=2.5)
    xlabel('velocity dispersion after group finder')
    ylabel('original group velocity dispersion')
    title('central method / worse linking lengths 1 / vel dispersion')
    ylim(0,600)
    xlim(0,600)
    savefig('velplot_old.png')
    show()

    x = np.arange(min(nus), max(nus), 0.5);
    y = x
    plot(nus,n2[list(of)][gal_dmass != 0],'y.')
    plot(x,y,color='k',label='1 to 1 line')
    xlim(0,80)
    ylim(0,200)
    xlabel('size of group after group finder')
    ylabel('original galaxies group size')
    title('central method / worse linking lengths/ num in group')
    savefig('numplot_old.png')
    show()

#~ ~ ~ ~ ~ ~ ~ ~ ~ ~ plot errors ~ ~ ~ ~ ~ ~ ~ ~ ~#
#norm_err_files - lists xx & norm_err for several modes/options
def ploterrs():
    "plot normalized errors for each mass determination method"
    norm_err_files = loadtxt('norm_err/err_files2',dtype='str',unpack=True)
    for i in range(0,len(norm_err_files)):
        xx1,err1 = loadtxt('norm_err/' + norm_err_files[i] + '.dat', unpack=True)
        plot(xx1,err1,label=norm_err_files[i])
        
    title('Normalized Errors for various mass estimators')
    xlabel('log Dyn Mass')
    ylabel('Normalized Scatter')
    legend(loc=8,numpoints=1)
    savefig('gp_errors.png')
    show()

# - - - - - analyze sigma vs r - - - - - #
#dynmass matches ind. pick out bad mass estimates groups
#redefine
def plotsigrad(radarr,ind):
    'plot velocity dispersion and radius and analyze outlyers'
    truem = gal_halom[gal_dmass != 0]
    dynm  = log10(gal_dmass[gal_dmass != 0])
    vel   = gal_vel[gal_dmass != 0]
    rad   = gal_rad[gal_dmass != 0]
    radarr = array(radarr)
    
    ind = array(ind)
    weirdos = radarr[gal_dmass[ind] != 0] > 0.1
    crazies = velarr[gal_dmass[ind] != 0] > 100
    plot(dynm,truem,'y.',zorder=-100,label='$N_\mathrm{galaxies}$ > 2')
    plot(dynm[weirdos],truem[weirdos],'bo')
    plot(dynm[crazies],truem[crazies],'go')
    x = np.arange(10, 15.5, 0.5)
    y = x #+ .7
    plot(x,y,color='g',label='1 to 1 line')        #plot 1:1 line
    show()



#dudes below 1:1 line
def plotodds():
    'plot things/outliers'
    oddities = ((dynm/truem) > 1) & (truem < 11.4)
    plot(dynm[oddities],truem[oddities],'bo')
    xlabel('log Dynamical Mass')
    ylabel('log Halo Mass')
    xlim(8,14)
    ylim(10,15)
    title('oddities 1')
    savefig('gp_odd_m1.png')
    show()

    scatter(vel,rad,c='y',marker='o',edgecolor='None')
    scatter(vel[oddities],rad[oddities],c='b',marker='o')
    xlim(0,400)
    ylim(0,.6)
    xlabel('sigma - vel dispersion')
    ylabel('projected radius')
    title('oddities radius vs sigma')
    savefig('gp_odd_vr1.png')
    show()

def plotrainbow():
#plot sigma rad plot for various dyn masses
    'plot sigma v radius color coding on dmass, halom, then # galaxies'
    mind1 = dynm < 12.5
    mind2 = dynm < 12
    mind3 = dynm < 11.5
    mind4 = dynm < 11
    mind5 = dynm < 10.5
    plot(vel[mind1],rad[mind1],'r.')
    plot(vel[mind2],rad[mind2],'y.')
    plot(vel[mind3],rad[mind3],'g.')
    plot(vel[mind4],rad[mind4],'b.')
    plot(vel[mind5],rad[mind5],'m.')
    show()
#plot sigma rad plot for various true masses
    mind1 = truem < 15
    mind2 = truem < 14
    mind3 = truem < 13
    mind4 = truem < 12
    mind5 = truem < 11
    plot(vel[mind1],rad[mind1],'r.')
    plot(vel[mind2],rad[mind2],'y.')
    plot(vel[mind3],rad[mind3],'g.')
    plot(vel[mind4],rad[mind4],'b.')
    plot(vel[mind5],rad[mind5],'m.')
    show()
#plot sigma rad plot for various halo size by # galaxies
    mind1 = nus > 3
    mind3 = nus > 5
    mind5 = nus > 7
    plot(vel[mind1],rad[mind1],'r.')
    plot(vel[mind3],rad[mind3],'y.')
    plot(vel[mind5],rad[mind5],'b.')
    show()

def grp_hist():
    'plot histogram of # galaxies in groups (true and after group finder)'

    freq = [0] * max(num)
    for i in range(1,int(max(num))+1):      #where num is true numbers        
        freq[i-1] = size(num[num == i])
    x = np.arange(1, max(num)+1, 1)
    
    pos = np.arange(len(x))
    width = 1.0
    
    ax = plt.axes()
    ax.set_xticks(pos + (width / 2))
    ax.set_xticklabels(x)
    
    plt.bar(pos[0:50], np.log2(freq[0:50]),width,color='b',alpha=0.5,label='after fof')
    #plt.bar(pos[0:50], np.log2(freqs[0:50]),width,color='r',alpha=0.5,label='true grouping')
    legend()
    #savefig('histplot.png')
    show()

def sigrad2():
    'plot true sig and radius for each galaxy vs the new one from set group finder ran to create'
    
plotmass()
newplot()
#ploterrs()
#plotsigrad(radarr,ind)
#plotodds()
#plotrainbow()
#grp_hist()


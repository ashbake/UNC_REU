#!/usr/bin/python
#execfile('plotECO.py')


execfile('eco_wrapper.py')

fig = plt.figure(figsize=(15,5))

plot(ra,dec,'k,')
scatter(ra[gal_dmass > 10],dec[gal_dmass > 10],marker='o',c=gal_dmass[gal_dmass > 10],edgecolor='None',s=3,zorder=300)
cb=colorbar()
cb.set_label("Log halo mass ($\mathsf{M_\odot}$)",rotation=270)
xlabel('RA (deg)')
ylabel('Dec (deg)')
xlim(min(ra),max(ra))
ylim(min(dec),max(dec))

alfatop = where((goodinalfa == 1) & (dec > 20))
alfabot = where((goodinalfa == 1) & (dec < 20))

maxratop, maxdectop = max(ra[alfatop]), max(dec[alfatop])
minratop, mindectop = min(ra[alfatop]), min(dec[alfatop])
maxrabot, maxdecbot = max(ra[alfabot]), max(dec[alfabot]) 
minrabot, mindecbot = min(ra[alfabot]), min(dec[alfabot])

plot([minratop,maxratop],[maxdectop,maxdectop],'m',ls='dashed',lw=3)
plot([minratop,maxratop],[mindectop,mindectop],'m',ls='dashed',lw=3)
plot([minrabot,maxrabot],[maxdecbot,maxdecbot],'m',ls='dashed',lw=3)
plot([minrabot,maxrabot],[mindecbot,mindecbot],'m',ls='dashed',lw=3)
plt.fill_between([minratop,maxratop],maxdectop, mindectop, color='m',alpha=.1)
plt.fill_between([minratop,maxratop],maxdecbot, mindecbot, color='m',alpha=.1)

plot([131.25,236.25],[0,0],'darkolivegreen',ls='dashed',lw=4)
plot([131.25,236.25],[5,5],'darkolivegreen',ls='dashed',lw=4)
plot([131.25,131.25],[0,5],'darkolivegreen',ls='dashed',lw=4)
plot([236.25,236.25],[0,5],'darkolivegreen',ls='dashed',lw=4)

savefig('figures/ECOmap.pdf',dpi=300, bbox_inches='tight')
show()

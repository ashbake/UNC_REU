#!/usr/bin/python
#execfile('gp_velcor.py')
#this program determines error in velocity dispersion if use
#linking lengths that don't maximize for it

execfile('gp_wrapper.py')

#galaxy vel vs vel after fof
cols=nold
cols[cols > 100]=100
pick = all([num < 100,num > 4],axis=0)
scatter(list(vel[pick]),list(fof_vel[pick]),marker='.',c=sqrt(cols[pick]),edgecolors='none')
xlabel('old velocity dispersion')
ylabel('v disp after fof')
colorbar()
ylim(0,500)
xlim(0,600)
title(llens+'color = nold')
savefig('veltest'+llens)
clf()

pick = all([num < 100,num > 2],axis=0)
scatter(list(num[pick]),list(fof_vel[pick]),marker='.',c=sqrt(vel[pick]),edgecolor='none')
ylabel('fof vel disp')
xlabel('number in group')
colorbar()
ylim(0,400)
xlim(0,100)
title(llens + '  color = sqrt(old vel disp)')
savefig('numtest'+llens)
clf()

ncut = [3,4,5,6]
cuts =[100,140,160,190]

badpts0 = all([num == ncut[0],fof_vel > cuts[0]],axis=0)
badpts1 = all([num == ncut[1],fof_vel > cuts[1]],axis=0)
badpts2 = all([num == ncut[2],fof_vel > cuts[2]],axis=0)
badpts3 = all([num == ncut[3],fof_vel > cuts[3]],axis=0)

bad = any([badpts0,badpts1,badpts2,badpts3],axis=0)
scatter(list(vel[pick]),list(fof_vel[pick]),marker='.',c=sqrt(cols[pick]),edgecolors='none')
scatter(list(vel[bad]),list(fof_vel[bad]),marker='.',c='k',edgecolors='none')





pick = all([num < 100,num > 3],axis=0)
scatter(list(num[pick]),list(fof_vel[pick]),marker='.',c=sqrt(num[pick]),edgecolor='none')
ylabel('fof vel disp')
xlabel('number in group')
colorbar()
savefig('numtest2')
clf()

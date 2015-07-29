file = 'MassFunction.dat'
data =np.loadtxt(file)
num = data[:,1]

nums = [num[400],num[350],num[300],num[250],num[200],num[150],num[100]]  #abundance/yaxis
msun = [12,12.5,13,13.5,14,14.5,15] #corresponding errors
errs = [.015,.0175,.02,.022,.031,.09,.15]    #from graph
nums = array(nums)
errs1 = 1.0 - 17.5*array(errs)
errs2 = 1.0+17.5*array(errs)
msun=array(msun)-log10(.7)

plot(msun,nums)
plot(msun,nums*errs1)
plot(msun,nums*errs2)
plot(msun[0]+.2,nums[0]*errs1[0],'ro')
plot(msun[1]+.22,nums[1]*errs1[1],'ro')
plot(msun[2]+.22,nums[2]*errs1[2],'ro')
plot(msun[3]+.15,nums[3]*errs1[3],'ro')
plot(msun[4]+.32,nums[4]*errs1[4],'ro')
plot(msun[5]+.4,nums[5]*errs1[5],'ro')
plot(msun[6]+.5,nums[6]*errs1[6],'ro')

plot(msun[0]-.15,nums[0]*errs2[0],'ro')
plot(msun[1]-.1,nums[1]*errs2[1],'ro')
plot(msun[2]-.09,nums[2]*errs2[2],'ro')
plot(msun[3]-.08,nums[3]*errs2[3],'ro')
plot(msun[4]-.08,nums[4]*errs2[4],'ro')
plot(msun[5]-.15,nums[5]*errs2[5],'ro')
plot(msun[6]-.1,nums[6]*errs2[6],'ro')
xlabel('mass log msun')
ylabel('N(m)')

ax3 = plt.axes([.5, .5, .3, .3])
ax3.plot(msun,nums)
ax3.plot(msun,nums*errs1)
ax3.plot(msun,nums*errs2)
ax3.set_xlim(14,15.2)
ax3.set_ylim(0,.0002)

xlabel('mass log msun')
ylabel('N(m)')

masserrslow= [.2,.22,.22,.15,.32,.4,.5]
masserrshi=[.15,.1,.09,.08,.08,.15,.1]

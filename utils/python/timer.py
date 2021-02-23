from common_modules import *
tweak_rcParams(fs=12, home='/home1/07918/fyli')
mpl.use('TKAgg')

exec(open('./para.out').read())
data=np.loadtxt('../../data/time.dat')

plt.figure()
time=data[:,0]*dt
ind=np.argsort(data[-1,1:])
print(ind)
for i in range(7):
  j=ind[-(i+1)]+1
  plt.plot(time,data[:,j],label='timer %d'%j)

plt.legend(loc='best')
#plt.ylim(0,2000)
plt.grid()
plt.ylabel('timers (s)')
plt.xlabel(r'$t\Omega_i$')
# plt.savefig('timer.png')
plt.show()


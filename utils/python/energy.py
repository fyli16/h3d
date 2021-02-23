import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('GTK3Agg')
#matplotlib.get_backend()

exec(open('./para.out').read())
data=np.loadtxt('../../data/energy.dat')

time=data[:,0]*dt
energy_e=data[:,1]
energy_b=data[:,2]
#energy_b=data[:,2]
energy_f=data[:,3]
energy_t=data[:,4]
energy_p=data[:,5]
plt.plot(time,energy_e,label=r'$\delta E^2$')
plt.plot(time,energy_b,label=r'$\delta B^2$')
plt.plot(time,energy_f,label=r'flow energy')
plt.plot(time,energy_t,label=r'thermal energy')
plt.plot(time,energy_p,label=r'particle energy')
plt.plot(time,energy_e+energy_b+energy_p,label='total energy')
plt.grid()
plt.legend(loc='best')
plt.xlabel(r'$t\Omega_i$')
plt.savefig('energy.png')
plt.show()

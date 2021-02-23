import numpy as np

import matplotlib.pyplot as plt

runid=30
den=np.fromfile('den_fluc_h3d-%03d.dat'%runid)
dt=0.025
nt=20
exec(open('./para.out').read())
time=np.arange(nt)*dt*1000
plt.cla()
plt.plot(time,den[:nt],label=r'$\delta \rho^2$')

b=np.fromfile('magnetic_energy_h3d-%03d.dat'%runid)
e=np.fromfile('electric_energy_h3d-%03d.dat'%runid)
plt.plot(time,b[:nt],label=r'$\delta B^2$')
plt.plot(time,e[:nt],label=r'$\delta E^2$')
plt.grid()
plt.legend()
plt.xlabel(r'$t\Omega_i$')
plt.savefig('history.png')
plt.show()

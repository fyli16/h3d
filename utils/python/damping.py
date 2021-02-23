import numpy as np

import matplotlib.pyplot as plt

run_name='h3d-106'
dt=0.025
nt=20
ifshow=True
exec(open('./para.out').read())
#den=np.fromfile('den_fluc_h3d-%s.dat'%run_name)
e=np.fromfile('electric_energy_%s.dat'%run_name)
b=np.fromfile('magnetic_energy_%s.dat'%run_name)
time=np.arange(nt)*dt*tinterval
plt.cla()
#plt.semilogy(time,den[:nt],label=r'$\delta \rho^2$')
plt.semilogy(time,e[:nt],label=r'$\delta E^2$')
plt.semilogy(time,b[:nt],label=r'$\delta B^2$')
#plt.ylim(ymin=1e-3)
plt.legend()
plt.grid()
plt.xlabel(r'$t\Omega_i$')
plt.savefig('damping.png')
if ifshow:
  plt.show()

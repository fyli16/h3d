import numpy as np

import matplotlib.pyplot as plt

run_name='h3d-106'
dt=0.025
nt=20
ifshow=True
exec(open('./para.out').read())
data=np.loadtxt('decomp_%s.dat'%run_name)
ax=plt.subplot(111)
plt.plot(data[0,:]*dt,data[1,:],label='compressible flow')
plt.plot(data[0,:]*dt,data[2,:],label='solenoidal flow')
plt.legend()
plt.ylim(bottom=0)
plt.grid()
ax2=ax.twinx()
ax2.plot(data[0,:]*dt,data[1,:]/(data[1,:]+data[2,:]),'k')
ax.set_xlabel(r'$t\Omega_i$')
ax.set_ylabel('flow energy')
plt.xlim(data[0,0]*dt, data[0,-1]*dt)
plt.ylim(bottom=0)
plt.savefig('decomp.png')
if ifshow:
  plt.show()

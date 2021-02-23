import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

binary=False
binary=True
tracking_width=14
nop=500
nop=15
cova=100
if binary:
  f=FortranFile('../../data/tracking_b.dat')
  #f=FortranFile('../../data/tracking_mpi.dat')
  d=[]
  while True:
    try:
      tmp=f.read_reals(np.float64)
    except:
      break
    d.append(tmp)
  data=np.array(d).reshape(-1,nop,tracking_width) 
  # note, the order of dimensions is the reverse of that in Fortran code
  print(data.shape)
  print(data[0,0,:])
else:
  data=np.loadtxt('../../data/tracking.dat')

plt.figure(figsize=(12,8))
ax1=plt.subplot(2,2,1)
ax2=plt.subplot(2,2,2)
ax3=plt.subplot(2,2,3)
ax4=plt.subplot(2,2,4)
tags=[1]
names=['x','y','z']
if binary:
  for i in tags:
    # i starts from 0 for binary output
    i=i-1
    for j in range(3):
      x=data[:,i,j].flatten()
      ax1.plot(x,'-',label=names[j])
    for j in range(3):
      v=data[:,i,j+3].flatten()
      ax2.plot(v,'-',label='v'+names[j])
    for j in range(3):
      v=data[:,i,j+8].flatten()
      ax3.plot(v*cova**2,'-',label='E'+names[j])
    for j in range(3):
      v=data[:,i,j+11].flatten()
      ax4.plot(v*cova,'-',label='B'+names[j])
  ax1.legend()
  ax1.grid()
  ax2.legend()
  ax2.grid()
  ax3.legend()
  ax3.grid()
  ax3.set_xlabel('t')
  ax4.legend()
  ax4.grid()
  ax4.set_xlabel('t')
  plt.savefig('tracking_b.png')
else:
  for i in tags:
    ind=np.where(abs(data[:,7]-i)<1e-6)
    for j in range(3):
      x=data[ind,j].flatten()
      ax1.plot(x,'-',label=names[j])
    for j in range(3):
      v=data[ind,j+3].flatten()
      ax2.plot(v,'-',label='v'+names[j])
  ax1.legend()
  ax1.grid()
  ax2.set_xlabel('t')
  ax2.legend()
  ax2.grid()
  plt.savefig('tracking.png')
plt.show()

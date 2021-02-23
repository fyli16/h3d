import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

binary=False
binary=True
tracking_width=14
nop=500
nop=20
nop=15
nop=1000
cova=100
f=FortranFile('../particles.dat')
d=[]
while True:
  try:
    tmp=f.read_reals(np.float64)
  except:
    break
  d.append(tmp)
print(len(d))
data=np.array(d).reshape(nop,-1,tracking_width) 
# note, the order of dimensions is the reverse of that in Fortran code
print(data.shape)
print(data[0,0,:])

tags=[1,2,3,4,5]
for i in tags:
  print(i, data[i-1,0,6])
  plt.figure(figsize=(12,8))
  ax1=plt.subplot(2,2,1)
  ax2=plt.subplot(2,2,2)
  ax3=plt.subplot(2,2,3)
  ax4=plt.subplot(2,2,4)
  names=['x','y','z']
  # i starts from 0 for binary output
  for j in range(3):
    x=data[i-1,:,j].flatten()
    ax1.plot(x,'-',label=names[j])
  for j in range(3):
    v=data[i-1,:,j+3].flatten()
    ax2.plot(v,'-',label='v'+names[j])
  for j in range(3):
    v=data[i-1,:,j+8].flatten()
    ax3.plot(v*cova**2,'-',label='E'+names[j])
  for j in range(3):
    v=data[i-1,:,j+11].flatten()
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
  plt.suptitle('Particle #%d'%i)
  plt.savefig('particle_%04d.png'%i)
plt.show()

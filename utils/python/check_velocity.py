import numpy as np
import matplotlib.pyplot as plt
# plot flow velocities of different ion species in H3D

nx=256
ny=256
nz=256
dt=0.025

exec(open('./para.out').read())
nspec=len(wspec)

it=1
ix=2
iy=5
iz=19
cut='x'
cut='z'
cut='y'

plt.figure()
ax1=plt.subplot(311)
ax2=plt.subplot(312)
ax3=plt.subplot(313)
for s in range(nspec):
  vxs=np.fromfile('../../data/vxs_%d_%d.gda'%(s+1,it),dtype=np.float32)
  vxs=np.reshape(vxs,(nx,ny,nz),order='F')
  if cut=='z':
    ax1.plot(vxs[ix,iy,:],label='vx_%d'%(s+1))
  elif cut=='x':
    ax1.plot(vxs[:,iy,iz],label='vx_%d'%(s+1))
  else:
    ax1.plot(vxs[ix,:,iz],label='vx_%d'%(s+1))

  vys=np.fromfile('../../data/vys_%d_%d.gda'%(s+1,it),dtype=np.float32)
  vys=np.reshape(vxs,(nx,ny,nz),order='F')
  if cut=='z':
    ax2.plot(vys[ix,iy,:],label='vy_%d'%(s+1))
  elif cut=='x':
    ax2.plot(vys[:,iy,iz],label='vy_%d'%(s+1))
  else:
    ax2.plot(vys[ix,:,iz],label='vy_%d'%(s+1))

  vzs=np.fromfile('../../data/vzs_%d_%d.gda'%(s+1,it),dtype=np.float32)
  vzs=np.reshape(vzs,(nx,ny,nz),order='F')
  if cut=='z':
    ax3.plot(vzs[ix,iy,:],label='vz_%d'%(s+1))
  elif cut=='x':
    ax3.plot(vzs[:,iy,iz],label='vz_%d'%(s+1))
  else:
    ax3.plot(vzs[ix,:,iz],label='vz_%d'%(s+1))


ax1.legend()
ax2.legend()
ax3.legend()
if cut=='z':
  ax3.set_xlabel('z')
elif cut=='x':
  ax3.set_xlabel('x')
else:
  ax3.set_xlabel('y')
plt.suptitle(r'$t\Omega_i=%d$'%(it*dt))
plt.show()

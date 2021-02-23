import numpy as np
import matplotlib.pyplot as plt

nx=64
ny=64
nz=64
dt=0.025
exec(open('./para.out').read())

it=1
it=4000
l0e=l0b=l0v=0.50
iz=32
iy=9
ix=8
cut='y'
cut='x'
cut='z'

fixed_level=True
fixed_level=False

ex=np.fromfile('../../data/ex_%d.gda'%it,dtype=np.float32)
ex=np.reshape(ex,(nx,ny,nz),order='F')

ey=np.fromfile('../../data/ey_%d.gda'%it,dtype=np.float32)
ey=np.reshape(ey,(nx,ny,nz),order='F')

ez=np.fromfile('../../data/ez_%d.gda'%it,dtype=np.float32)
ez=np.reshape(ez,(nx,ny,nz),order='F')

bx=np.fromfile('../../data/bx_%d.gda'%it,dtype=np.float32)
bx=np.reshape(bx,(nx,ny,nz),order='F')

by=np.fromfile('../../data/by_%d.gda'%it,dtype=np.float32)
by=np.reshape(by,(nx,ny,nz),order='F')

bz=np.fromfile('../../data/bz_%d.gda'%it,dtype=np.float32)
bz=np.reshape(bz,(nx,ny,nz),order='F')

vix=np.fromfile('../../data/vix_%d.gda'%it,dtype=np.float32)
vix=np.reshape(vix,(nx,ny,nz),order='F')

viy=np.fromfile('../../data/viy_%d.gda'%it,dtype=np.float32)
viy=np.reshape(viy,(nx,ny,nz),order='F')

viz=np.fromfile('../../data/viz_%d.gda'%it,dtype=np.float32)
viz=np.reshape(viz,(nx,ny,nz),order='F')

plt.figure(figsize=(12,10))
ax=plt.subplot(3,3,1)
if fixed_level:
  levels=np.linspace(-l0e,l0e,20)
else:
  levels=None
if cut=='z':
  field=ex[:,:,iz].transpose()
elif cut=='y':
  field=ex[:,iy,:].transpose()
else:
  field=ex[ix,:,:].transpose()
m=ax.contourf(field,levels=levels)
plt.colorbar(m)
plt.title('Ex')
ax.set_xticklabels([])
if cut=='z':
  plt.ylabel('y')
else:
  plt.ylabel('z')
plt.grid()
ax=plt.subplot(3,3,2)
if cut=='z':
  field=ey[:,:,iz].transpose()
elif cut=='y':
  field=ey[:,iy,:].transpose()
else:
  field=ey[ix,:,:].transpose()
m=ax.contourf(field,levels=levels)
plt.colorbar(m)
plt.title('Ey')
ax.set_xticklabels([])
ax.set_yticklabels([])
plt.grid()
ax=plt.subplot(3,3,3)
if cut=='z':
  field=ez[:,:,iz].transpose()
elif cut=='y':
  field=ez[:,iy,:].transpose()
else:
  field=ez[ix,:,:].transpose()
m=ax.contourf(field,levels=levels)
plt.colorbar(m)
plt.title('Ez')
plt.grid()
ax.set_xticklabels([])
ax.set_yticklabels([])


if fixed_level:
  levels=np.linspace(-l0b,l0b,20)
else:
  levels=None
ax=plt.subplot(3,3,4)
if cut=='z':
  field=bx[:,:,iz].transpose()
elif cut=='y':
  field=bx[:,iy,:].transpose()
else:
  field=bx[ix,:,:].transpose()
m=ax.contourf(field,levels=levels)
plt.colorbar(m)
plt.title('Bx')
ax.set_xticklabels([])
if cut=='z':
  plt.ylabel('y')
else:
  plt.ylabel('z')
plt.grid()
ax=plt.subplot(3,3,5)
if cut=='z':
  field=by[:,:,iz].transpose()
elif cut=='y':
  field=by[:,iy,:].transpose()
else:
  field=by[ix,:,:].transpose()
m=ax.contourf(field,levels=levels)
plt.colorbar(m)
plt.title('By')
ax.set_xticklabels([])
ax.set_yticklabels([])
plt.grid()
ax=plt.subplot(3,3,6)
if cut=='z':
  field=bz[:,:,iz].transpose()
elif cut=='y':
  field=bz[:,iy,:].transpose()
else:
  field=bz[ix,:,:].transpose()
if fixed_level:
  m=ax.contourf(field,levels=levels+1)
else:
  m=ax.contourf(field)
plt.colorbar(m)
ax.set_xticklabels([])
ax.set_yticklabels([])
plt.title('Bz')
plt.grid()

if fixed_level:
  levels=np.linspace(-l0v,l0v,20)
else:
  levels=None
ax=plt.subplot(3,3,7)
if cut=='z':
  field=vix[:,:,iz].transpose()
elif cut=='y':
  field=vix[:,iy,:].transpose()
else:
  field=vix[ix,:,:].transpose()
m=ax.contourf(field,levels=levels)
plt.colorbar(m)
plt.title('vix')
if cut=='z':
  plt.xlabel('x')
  plt.ylabel('y')
elif cut=='y':
  plt.xlabel('x')
  plt.ylabel('z')
elif cut=='x':
  plt.xlabel('y')
  plt.ylabel('z')
plt.grid()
ax=plt.subplot(3,3,8)
if cut=='z':
  field=viy[:,:,iz].transpose()
elif cut=='y':
  field=viy[:,iy,:].transpose()
else:
  field=viy[ix,:,:].transpose()
m=ax.contourf(field,levels=levels)
plt.colorbar(m)
plt.title('viy')
if cut=='x':
  plt.xlabel('y')
else:
  plt.xlabel('x')
ax.set_yticklabels([])
plt.grid()
ax=plt.subplot(3,3,9)
if cut=='z':
  field=viz[:,:,iz].transpose()
elif cut=='y':
  field=viz[:,iy,:].transpose()
else:
  field=viz[ix,:,:].transpose()
m=ax.contourf(field,levels=levels)
plt.colorbar(m)
if cut=='x':
  plt.xlabel('y')
else:
  plt.xlabel('x')
ax.set_yticklabels([])
plt.title('viz')
plt.grid()
if cut=='z':
  plt.suptitle(r'$t\Omega_i=%d, z=%d$'%(it*dt,iz))
elif cut=='y':
  plt.suptitle(r'$t\Omega_i=%d, y=%d$'%(it*dt,iy))
else:
  plt.suptitle(r'$t\Omega_i=%d, x=%d$'%(it*dt,ix))
plt.savefig('2d_cut_%s_%d.png'%(cut,it*dt))
plt.show()

import numpy as np
import matplotlib.pyplot as plt

nx=256
ny=256
nz=256
L=256.

dt=0.025

# for particle 933
it=11000
iz=77
it=13000
iz=77
#for particle 193 and 195
it=1000
iz=120

ex=np.fromfile('../../data/ex_%d.gda'%it,dtype=np.float32)
ex=np.reshape(ex,(nx,ny,nz),order='F')

ey=np.fromfile('../../data/ey_%d.gda'%it,dtype=np.float32)
ey=np.reshape(ey,(nx,ny,nz),order='F')

bx=np.fromfile('../../data/bx_%d.gda'%it,dtype=np.float32)
bx=np.reshape(bx,(nx,ny,nz),order='F')

by=np.fromfile('../../data/by_%d.gda'%it,dtype=np.float32)
by=np.reshape(by,(nx,ny,nz),order='F')

bz=np.fromfile('../../data/bz_%d.gda'%it,dtype=np.float32)
bz=np.reshape(bz,(nx,ny,nz),order='F')

den=np.fromfile('../../data/den_%d.gda'%it,dtype=np.float32)
den=np.reshape(den,(nx,ny,nz),order='F')

data=np.load('j_%d.npz'%it)
jx=data['jx'].real
jy=data['jy'].real
jz=data['jz'].real
j=np.sqrt(jx**2+jy**2+jz**2)

plt.figure(figsize=(12,10))
ax=plt.subplot(2,2,1)
levels=np.linspace(-0.55,0.55,20)
field=ex[:,:,iz].transpose()
m=ax.contourf(field,levels=levels)
plt.colorbar(m)
plt.title('Ex')
plt.ylabel('y')
plt.grid()
ax=plt.subplot(2,2,2)
field=ey[:,:,iz].transpose()
m=ax.contourf(field,levels=levels)
plt.colorbar(m)
plt.title('Ey')
plt.grid()


#levels=np.linspace(-0.5,0.5,20) # for j_x
#levels=np.linspace(0.0,0.5,20) #for j
levels=np.linspace(0.0,0.7,20) #for j
#levels=np.linspace(0.0,4.5,20) #for density
ax=plt.subplot(2,2,3)
#field=bx[:,:,iz].transpose()
#field=jx[:,:,iz].transpose()
field=j[:,:,iz].transpose()
#field=den[:,:,iz].transpose()
m=ax.contourf(field,levels=levels)
plt.colorbar(m)
plt.title('|j|')
plt.xlabel('x')
plt.ylabel('y')
plt.grid()


levels=np.linspace(-0.45,0.45,20)
ax=plt.subplot(2,2,4)
field=bz[:,:,iz].transpose()
m=ax.contourf(field,levels=levels+1)
plt.colorbar(m)
plt.xlabel('x')
plt.title('Bz')
plt.grid()

plt.suptitle(r'$t\Omega_i=%d$'%(it*dt))
plt.savefig('2d_cut_%d.png'%(it*dt))
plt.show()

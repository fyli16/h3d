import numpy as np
import matplotlib.pyplot as plt

def read_binary(filename,shape):
  fd=open(filename,'rb')
  # read Fortran binary array
  data=np.fromfile(filename,dtype='f4').reshape((nx,ny,ny),order='F')
  fd.close()
  return data

nx=32
ny=32
nz=32
it=100
ndimy=2
ndimz=2
npy=ny/ndimy
npz=nz/ndimz
path='../../data/'
titles=np.array([['ex','ey','ez'],['bx','by','bz']])

data=np.loadtxt('../../data/probes.dat')
for myid in range(1):
  time=data[:,0]
  fig,axes=plt.subplots(2,3,figsize=(12,6))
  for i in range(2):
    for j in range(3):
      axes[i,j].plot(time,data[:,i*3+j+1+myid*6])
      axes[i,j].grid()
      axes[i,j].set_title(titles[i,j])
  for it in [1,200,400]:
    jb=npy*(int(myid/2))+1
    kb=npy*(myid%2)+1
    print(jb,kb)
    f=read_binary(path+'ex_%d.gda'%it,(nx,ny,nz))
    axes[0,0].plot(it,f[0,jb-1,kb-1],'*r') # note index starts at 0 in python vs
                                           # 1 in fortran
    f=read_binary(path+'ey_%d.gda'%it,(nx,ny,nz))
    axes[0,1].plot(it,f[0,jb-1,kb-1],'*r')
    f=read_binary(path+'ez_%d.gda'%it,(nx,ny,nz))
    axes[0,2].plot(it,f[0,jb-1,kb-1],'*r')
    #bz=read_binary(path+'bz_%d.gda'%it,(nx,ny,nz))
    f=read_binary(path+'bx_%d.gda'%it,(nx,ny,nz))
    axes[1,0].plot(it,f[0,jb-1,kb-1],'*g')
    f=read_binary(path+'by_%d.gda'%it,(nx,ny,nz))
    axes[1,1].plot(it,f[0,jb-1,kb-1],'*g')
    f=read_binary(path+'bz_%d.gda'%it,(nx,ny,nz))
    axes[1,2].plot(it,f[0,jb-1,kb-1],'*g')
  for j in range(3):
    axes[1,j].set_xlabel('t')
  plt.tight_layout()
  plt.savefig('probes_%d.png'%myid)
plt.show()

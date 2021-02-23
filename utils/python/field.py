import numpy as np
import matplotlib.pyplot as plt
#from scipy.io.numpyio import fread, fwrite

def read_binary(filename,shape):
  fd=open(filename,'rb')
  # read Fortran binary array
  data=np.fromfile(filename,dtype='f4').reshape((nx,ny,nz),order='F')
  fd.close()
  return data

nx=32
ny=32
nz=1
it=7500
it=3000
it=1000
it=1
exec(open('./para.out').read())
path='../../data/'
#for it in [1, 1000, 3000, 7000]:
xaxis='x'
xaxis='y'
#for it in [1000,2000, 3000]:
for it in range(1000,11000,2000):
#for it in [20000,40000,60000,80000]:
  print('it=',it)
  plt.figure(figsize=(8,8))
  i=1
  #for fname in ['fox','foy','foz']:
  #for fname in ['bx','by','bz']:
  for fname in ['ex','ey','ez']:
    ax=plt.subplot(3,1,i)
    print(fname)
    f=read_binary(path+fname+'_%d.gda'%it,(nx,ny,nz))
    print('max=',f.max())
    print('min=',f.min())
    print('avg=',f.mean())

    if xaxis=='x':
      for iy in [0,1,2, 3]:
        plt.plot(f[0,iy,:],label='iy=%d'%iy)
    elif xaxis=='y':
      #for ix in [0,1,2, 3]:
      for ix in [0]:
        plt.plot(f[ix,0,:],label='ix=%d'%ix)
    plt.grid()
    plt.ylabel(fname)
    plt.legend(loc='best')
    if i==1:
      plt.title('t=%d'%(it*dt))
    i+=1
  if xaxis=='x':
    #plt.savefig('friction_y_%d'%(it))
    #plt.savefig('B_y_%d'%(it))
    plt.savefig('E_y_%d'%(it))
  elif xaxis=='y':
    #plt.savefig('friction_x_%d'%(it))
    #plt.savefig('B_x_%d'%(it))
    plt.savefig('E_x_%d'%(it))
  plt.show()


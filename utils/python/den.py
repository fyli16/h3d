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
it=1
it=7500
it=3000
it=20000
it=15000
it=10000
it=1
it=1000
exec(open('./para.out').read())
path='../../data/'

fig=plt.figure(figsize=(8,6))
#for it in [1,1000,2000,3000,4000,5000]:
for it in range(1000,27000,1000):
  print(it)
  fig.clf()
  plt.subplot(111)
  den=read_binary(path+'den_%d.gda'%it,(nx,ny,nz))
  print('total den=%.8f'%np.sum(den))
  #plt.plot(den[:,16,0],label='n vs x')
  #plt.plot(den[16,:,0],label='n vs y')
  plt.contourf(den[:,:,0].transpose(),20)
  plt.colorbar()
  plt.grid()
  plt.title(r'$t\Omega_i=%d$'%(it*dt))
  plt.savefig('density_%d.png'%(it*dt))


#plt.figure(figsize=(8,6))
#plt.subplot(111)
#for iy in [1,5, 16, 31]:
#  plt.plot(den[:,iy,0],label='iy=%d'%iy)
#plt.grid()
#plt.legend(loc='best')

#plt.show()


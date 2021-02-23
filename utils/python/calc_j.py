import numpy as np
import matplotlib.pyplot as plt

nx=240
ny=240
nz=240
L=240.
it=1000
it=16000
it=13000
it=11000
twopi=np.pi*2


if False:
  # 1D test
  plt.figure()
  ax1=plt.subplot(211)
  ax2=plt.subplot(212)
  x=np.arange(0,nx)*1.0
  phi=0.05*np.random.random(nx)+np.sin(twopi/nx*4*x)
  #kx=2*np.pi*np.arange(0.,nx)/nx
  kx=np.fft.fftfreq(nx)*twopi
  print(kx[0:5])
  print(kx[-5:])
  phik=np.fft.fft(phi)
  ax1.plot(phi)
  ax1.grid()
  #ax1.plot(phik.real)
  #ax2.plot(phik.imag)
  print(phik[0:5])
  print(phik[-5:])
  ek=1j*kx*phik
  e=np.fft.ifft(ek)
  ax2.plot(e.real)
  ax2.plot(e.imag,'--')
  ax2.grid()
  for ix in [5,13,22]:
    # the real part should match, and imag part neglegible
    print(e[ix])
    print((phi[ix+1]-phi[ix-1])/2)
  plt.show()
else:
  # real 3D data
  bx=np.fromfile('../../data/bx_%d.gda'%it,dtype=np.float32)
  bx=np.reshape(bx,(nx,ny,nz),order='F')

  by=np.fromfile('../../data/by_%d.gda'%it,dtype=np.float32)
  by=np.reshape(by,(nx,ny,nz),order='F')

  bz=np.fromfile('../../data/bz_%d.gda'%it,dtype=np.float32)
  bz=np.reshape(bz,(nx,ny,nz),order='F')

  kx,ky,kz=np.meshgrid(np.fft.fftfreq(nx),np.fft.fftfreq(ny),np.fft.fftfreq(ny),indexing='ij')
  kx*=2*np.pi
  ky*=2*np.pi
  kz*=2*np.pi
  bxk=np.fft.fftn(bx)
  byk=np.fft.fftn(by)
  bzk=np.fft.fftn(bz)
  jxk=1j*(ky*bzk-kz*byk)
  jyk=1j*(kz*bxk-kx*bzk)
  jzk=1j*(kx*byk-ky*bxk)
  #jxk=bxk
  #jyk=byk
  #jzk=bzk
  jx=np.fft.ifftn(jxk)
  jy=np.fft.ifftn(jyk)
  jz=np.fft.ifftn(jzk)

  ix,iy,iz=5,200,13
  plt.figure(figsize=(12,4))
  ax1=plt.subplot(111)
  #check the calculation
  ax1.plot(jx[:,iy,iz],'-k')
  ax1.plot((bz[:,iy+1,iz]-bz[:,iy-1,iz])/2-(by[:,iy,iz+1]-by[:,iy,iz-1])/2,'r')
  ax1.plot((bz[:,iy+1,iz]-bz[:,iy,iz])-(by[:,iy,iz+1]-by[:,iy,iz]),'b')
  print(jx[ix,iy,iz])
  print((bz[ix,iy+1,iz]-bz[ix,iy-1,iz])/2-(by[ix,iy,iz+1]-by[ix,iy,iz-1])/2)
  print(jy[ix,iy,iz])
  print((bx[ix,iy,iz+1]-bx[ix,iy,iz-1])/2-(bz[ix+1,iy,iz]-bz[ix-1,iy,iz])/2)
  print(jz[ix,iy,iz])
  print((by[ix+1,iy,iz]-by[ix-1,iy,iz])/2-(bx[ix,iy+1,iz]-bx[ix,iy-1,iz])/2)
  plt.show()

  np.savez('j_%d.npz'%it,jx=jx,jy=jy,jz=jz)

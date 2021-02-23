''' Helmholtz Decomposition using FFT,
    To-do: parallelization '''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import os

def read_binary(filename,shape):
  fd=open(filename,'rb')
 # log.update('Reading '+filename)
  # read Fortran binary array
  data=np.fromfile(filename,dtype='f4').reshape((nx,ny,nz),order='F')
  fd.close()
  return data

def calc_irrotational_field(fxk,fyk,fzk):
  phik=1j*(kx*fxk+ky*fyk+kz*fzk)/k2
  #phik=1j*(kx*fxk+ky*fyk)/k2
  #print(abs(phik).max())
  fcxk=-1j*kx*phik
  fcyk=-1j*ky*phik
  fczk=-1j*kz*phik
  fcx=np.fft.ifftn(fcxk)
  fcy=np.fft.ifftn(fcyk)
  fcz=np.fft.ifftn(fczk)
  return fcx,fcy,fcz

def calc_solenoidal_field(fxk,fyk,fzk):
  axk=1j*(ky*fzk-kz*fyk)/k2
  ayk=1j*(kz*fxk-kx*fzk)/k2
  azk=1j*(kx*fyk-ky*fxk)/k2

  fsxk=1j*(ky*azk-kz*ayk)
  fsyk=1j*(kz*axk-kx*azk)
  fszk=1j*(kx*ayk-ky*axk)

  fsx=np.fft.ifftn(fsxk)
  fsy=np.fft.ifftn(fsyk)
  fsz=np.fft.ifftn(fszk)
  return fsx,fsy,fsz

nx=32
ny=32
nz=32
nt=1
dt=0.025
exec(open('./para.out').read())
shape=(nx,ny,nz)
pi=np.pi

kx,ky,kz=np.meshgrid(np.fft.fftfreq(nx),np.fft.fftfreq(ny),np.fft.fftfreq(nz),indexing='ij')
kx*=2*pi/xmax*nx
ky*=2*pi/ymax*ny
kz*=2*pi/zmax*nz
k2=kx**2+ky**2+kz**2
k2[0,0,0]=1 # to avoid division by zero

ix=54
iy=21
iz=189
compressible_flow=np.zeros(nt)
solenoidal_flow=np.zeros(nt)
snaps=np.zeros(nt)
for it in range(nt):  
  snap= it*tinterval
  if it==0:
    snap=1
  snaps[it]=snap
  print('t=',snap)
  vx=read_binary('../../data/vix_%d.gda'%snap,shape)
  vy=read_binary('../../data/viy_%d.gda'%snap,shape)
  vz=read_binary('../../data/viz_%d.gda'%snap,shape)
  vxk=np.fft.fftn(vx)
  vyk=np.fft.fftn(vy)
  vzk=np.fft.fftn(vz)
  fcx,fcy,fcz=calc_irrotational_field(vxk,vyk,vzk)
  sz=1
  for i in shape:
    sz*=i
  compressible_flow[it]=np.sum(fcx**2+fcy**2+fcz**2)*0.5/sz

  fsx,fsy,fsz=calc_solenoidal_field(vxk,vyk,vzk)
  solenoidal_flow[it]=np.sum(fsx**2+fsy**2+fsz**2)*0.5/sz
np.savetxt('decomp_%s.dat'%run_name,(snaps,compressible_flow,solenoidal_flow))


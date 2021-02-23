import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

tracking_width=14
cova=100
dt=0.025

tags=[98,88,1583,483,403,123,225]
tags=[260,293,235,405]
tags=[260]

ts=36000
te=48000
ts=16000
te=28000
ts=0
te=80000
ts=28000
te=40000
time=np.arange(ts,te)*dt
for i in tags:
  f=FortranFile('../../data/p_%04d.dat'%i)
  d=[]
  while True:
    try:
      tmp=f.read_reals(np.float64)
    except:
      break
    d.append(tmp)
  data=np.array(d).reshape(-1,tracking_width) 
  print(data.shape)
  x=data[ts:te,0]
  y=data[ts:te,1]
  z=data[ts:te,2]
  vx=data[ts:te,3]
  vy=data[ts:te,4]
  vz=data[ts:te,5]
  v=np.sqrt(vx**2+vy**2+vz**2)
  ex=data[ts:te,8]
  ey=data[ts:te,9]
time=np.arange(ts,te)*dt
for i in tags:
  f=FortranFile('../../data/p_%04d.dat'%i)
  d=[]
  while True:
    try:
      tmp=f.read_reals(np.float64)
    except:
      break
    d.append(tmp)
  data=np.array(d).reshape(-1,tracking_width) 
  print(data.shape)
  x=data[ts:te,0]
  y=data[ts:te,1]
  z=data[ts:te,2]
  vx=data[ts:te,3]
  vy=data[ts:te,4]
  vz=data[ts:te,5]
  v=np.sqrt(vx**2+vy**2+vz**2)
  ex=data[ts:te,8]
  ey=data[ts:te,9]
  ez=data[ts:te,10]
  bx=data[ts:te,11]
  by=data[ts:te,12]
  bz=data[ts:te,13]
  b=np.sqrt(bx**2+by**2+bz**2)
  vdotE=(vx*ex+vy*ey+vz*ez)
  vpara=(vx*bx+vy*by+vz*bz)/b
  vperp=np.sqrt(v**2-vpara**2)
  # v_perp vector
  vpx=vx-vpara*bx/b
  vpy=vy-vpara*by/b
  vpz=vz-vpara*bz/b
  # E_perp vector
  epara=(ex*bx+ey*by+ez*bz)/b
  epx=ex-epara*bx/b
  epy=ey-epara*by/b
  epz=ez-epara*bz/b
  vpEp=vpx*epx+vpy*epy+vpz*epz

  plt.figure(figsize=(12,12))
  ax1=plt.subplot(6,1,1)
  ax2=plt.subplot(6,1,2)
  ax3=plt.subplot(6,1,3)
  ax4=plt.subplot(6,1,4)
  ax5=plt.subplot(6,1,5)
  ax6=plt.subplot(6,1,6)
  ax1.plot(time,x,'-',label='x')
  ax1.plot(time,y,'-',label='y')
  ax1.plot(time,z,'-',label='z')

  ax2.plot(time,vx,'-',label='vx')
  ax2.plot(time,vy,'-',label='vy')
  ax2.plot(time,vz,'-',label='vz')

  ax3.plot(time,ex,'-',label='Ex')
  ax3.plot(time,ey,'-',label='Ey')
  ax3.plot(time,ez,'-',label='Ez')

  ax4.plot(time,bx,'-',label='Bx')
  ax4.plot(time,by,'-',label='By')
  ax4.plot(time,bz,'-',label='Bz')

  ax5.plot(time,v,'-',label=r'$|{\bf v}|$')
  ax5.plot(time,vpara,'-',label=r'$v_\parallel$')
  #ax3.plot(time,vperp,'-',label=r'$v_\perp$')

  #ax4.plot(time,vdotE,'-',label=r'${\bf v}\cdot {\bf E}$')
  #ax6.plot(time,vx*ex,'-',label=r'$ v_xE_x$')
  #ax6.plot(time,vy*ey,'-',label=r'$ v_yE_y$')
  #ax6.plot(time,vz*ez,'-',label=r'$ v_zE_z$')
  ax6.plot(time,vpEp,'-',label=r'$ {\bf v}_\perp\cdot {\bf E}_\perp$')
  ax6.plot(time,vdotE-vpEp,'-',label=r'$ {v}_\parallel{E}_\parallel$')

  ax1.legend()
  ax1.grid()
  ax1.set_xticklabels([])
  ax2.legend()
  ax2.grid()
  ax2.set_xticklabels([])
  ax3.legend()
  ax3.grid()
  ax3.set_xticklabels([])
  ax4.legend()
  ax4.grid()
  ax4.set_xticklabels([])
  ax5.legend()
  ax5.grid()
  ax5.set_xticklabels([])
  ax6.legend()
  ax6.grid()
  ax6.set_xlabel(r'$t\Omega_i$')
  plt.suptitle('Particle #%d'%i)
  plt.savefig('partacc_%04d_%d_%d.png'%(i,ts,te))

  if False:
    plt.figure(figsize=(8,10))
    fs=0
    fe=400
    ax1=plt.subplot(2,1,1)
    for f in ['Ex','Ey','Ez']:
      if f=='Ex':
        field=ex-ex.mean()
      elif f=='Ey':
        field=ey-ey.mean()
      else:
        field=ez-ez.mean()
      freq=2*np.pi/(time[-1]-time[0])*np.arange(len(field))
      power=abs(np.fft.rfft(field))**2
      ax1.loglog(freq[fs:fe],power[fs:fe],label=f)
    ax1.set_xticklabels([])
    ax1.legend()
    ax1.grid()
    ax1.set_ylim(ymin=1e-8)

    ax2=plt.subplot(2,1,2)
    for f in ['Bx','By','Bz']:
      if f=='Bx':
        field=bx-bx.mean()
      elif f=='By':
        field=by-by.mean()
      else:
        field=bz-bz.mean()
      freq=2*np.pi/(time[-1]-time[0])*np.arange(len(field))
      power=abs(np.fft.rfft(field))**2
      ax2.loglog(freq[fs:fe],power[fs:fe],label=f)
    ax2.set_xlabel(r'$\omega/\Omega_i$')
    ax2.legend()
    ax2.grid()
    ax2.set_ylim(ymin=1e-4)
    plt.suptitle('Particle #%d'%i)
    plt.savefig('partacc_%04d_f.png'%i)
plt.show()

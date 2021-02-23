import numpy as np
import matplotlib.pyplot as plt
import struct
import sys

names=['H','4He','3He','O','Fe']
frac=[0.9,0.09985,3e-5,5e-5,7e-5] #density ratio in the simulation
nppc=[64,8,8,8,8]
abnd=[12.,10.99,7.70,8.93,7.67] # soalr photosphere abundance (log nH=12)
ncores=64
#for it in [1000,2000,3000,4000]:
for it in [10000]:
  plt.figure(figsize=(8,8))
  ax=plt.subplot(211)
  ax2=plt.subplot(212)
  f=open('../../data/particle_%d.bin'%it,'rb')
  recl=40
  ntotal,=struct.unpack("i",f.read(recl)[0:4])
  print('n_total=',ntotal)
  npes=7500
  exec(open('./para.out').read())

  n=0
  data=()
  #for i in range(npes):
  for i in range(ncores):
    n_in_volume,=struct.unpack("i",f.read(recl)[0:4])
    print('i=',i,'n_in_volume=',n_in_volume)
    data+=struct.unpack("10f"*n_in_volume,f.read(recl*n_in_volume))
    n+=n_in_volume
  f.close()
  data=np.array(data).reshape(n,10)
  #qp=frac/NPPC
  for s in range(5):
    specname=names[s]
    index=np.where(abs(data[:,6]-frac[s]/nppc[s])<1e-7)
    data2=data[index]
    energy=data2[:,3]**2+data2[:,4]**2+data2[:,5]**2
    hist,bins=np.histogram(energy,bins=np.logspace(np.log10(1e-7),np.log10(2e-4),100),normed=True)
    ax.plot(bins[:-1],hist,label='%s'%specname)
    ax2.plot(bins[:-1],hist*np.power(10,abnd[s]-abnd[0]),label='%s'%specname)
  ax.grid()
  ax.legend(loc='best')
  ax.set_ylabel(r'$f(E)$')
  ax.set_xscale('log')
  ax.set_yscale('log')
  ax2.grid()
  ax2.legend(loc='best')
  ax2.set_xlabel(r'$E=v^2/c^2$')
  ax2.set_xscale('log')
  ax2.set_yscale('log')
  ax2.set_ylabel(r'$f(E)\times$ abundance')
    #plt.setp(ax.get_xticklabels(),rotation=30)
    #plt.xlim(-0.004,0.004)
  plt.suptitle(r't=%d'%(it*dt))
  plt.subplots_adjust(left=0.08,right=0.95,top=0.95,bottom=0.1,wspace=0.2,hspace=0.2)
  plt.savefig('dist2_%d.png'%(it*dt))
plt.show()
  

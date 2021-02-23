import numpy as np
import matplotlib.pyplot as plt
import struct
import sys

plt.figure(figsize=(8,6))
ax=plt.subplot(111)
for it in [10000,20000,30000,40000]:
  f=open('../../data/particle_%d.bin'%it,'rb')
  recl=40
  ntotal,=struct.unpack("i",f.read(recl)[0:4])
  print('n_total=',ntotal)
  npes=7500
  exec(open('./para.out').read())
  if len(sys.argv)>1:
    s=int(sys.argv[1])

  n=0
  data=()
  #for i in range(npes):
  for i in range(80):
    n_in_volume,=struct.unpack("i",f.read(recl)[0:4])
    print('i=',i,'n_in_volume=',n_in_volume)
    data+=struct.unpack("10f"*n_in_volume,f.read(recl*n_in_volume))
    n+=n_in_volume
  f.close()
  data=np.array(data).reshape(n,10)
  #qp=frac/NPPC
  if s==1:
    specname='H'
    index=np.where(abs(data[:,6]-0.0140625)<1e-7)
  elif s==2:
    specname='4He'
    index=np.where(abs(data[:,6]-0.01248125)<1e-7)
  elif s==3:
    specname='3He'
    index=np.where(abs(data[:,6]-3.75e-6)<1e-8)
  elif s==4:
    specname='O'
    index=np.where(abs(data[:,6]-6.25e-6)<1e-8)
  elif s==5:
    specname='Fe'
    index=np.where(abs(data[:,6]-8.75e-6)<1e-8)
  else:
    print('wrong argument:',s)
  data2=data[index]
  energy=data2[:,3]**2+data2[:,4]**2+data2[:,5]**2

  plt.hist(energy,bins=np.logspace(np.log10(1e-7),np.log10(2e-4),100),normed=True,histtype='step',label='t=%d'%(it*dt))
plt.grid()
plt.legend(loc='best')
plt.xlabel('E')
plt.ylabel('f(E)')
ax.set_xscale('log')
ax.set_yscale('log')
  #plt.setp(ax.get_xticklabels(),rotation=30)
  #plt.xlim(-0.004,0.004)
plt.suptitle(r'%s'%(specname))
plt.subplots_adjust(left=0.08,right=0.95,top=0.95,bottom=0.1,wspace=0.2,hspace=0.38)
plt.savefig('%s_dist_%d.png'%(specname,it))
plt.show()
  

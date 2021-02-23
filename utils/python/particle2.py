import numpy as np
import matplotlib.pyplot as plt
import struct

exec(open('./para.out').read())
plt.figure(figsize=(6,9))
ax1=plt.subplot(3,1,1)
ax2=plt.subplot(3,1,2)
ax3=plt.subplot(3,1,3)
axes=[ax1,ax2,ax3]
names=[r'$v_x/v_A$',r'$v_y/v_A$',r'$v_z/v_A$']
v0=6000. # c/v_A
#for it in [20000,40000,60000,80000]:
for it in [10000,20000]:
  f=open('../../data/particle_%d.bin'%it,'rb')
  recl=40
  ntotal,=struct.unpack("i",f.read(recl)[0:4])
  print('n_total=',ntotal)
  npes=7500
  n=0
  data=()
  #for i in range(npes):
  for i in range(0,npes,191):
    n_in_volume,=struct.unpack("i",f.read(recl)[0:4])
    print('i=',i,'n_in_volume=',n_in_volume)
    data+=struct.unpack("10f"*n_in_volume,f.read(recl*n_in_volume))
    n+=n_in_volume
  f.close()
  data=np.array(data).reshape(n,10)
  for i in range(3):
    ax=axes[i]
    hist,bin_edges=np.histogram(data[:,i+3]*v0,bins='auto',range=(-20,20),density=True) 
    ax.plot(bin_edges[:-1],hist,label=r'$t\Omega_i=%d$'%(it*dt))

for i in range(3):
  ax=axes[i]
  #plt.setp(ax.get_xticklabels(),rotation=30)
  #ax.set_xlim(-1e-4,1e-4)
  #ax.set_xlim(-0.5,0.5)
  ax.set_xlabel(names[i])
  ax.grid()
  ax.legend(loc='best')
plt.tight_layout()
plt.savefig('vel_dist.png')
plt.show()
  

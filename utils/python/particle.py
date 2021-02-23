import numpy as np
import matplotlib.pyplot as plt
import struct
import sys

it=2000
it=4000
it=6000
it=1
f=open('../../data/particle_%d.bin'%it,'rb')
recl=40
ntotal,=struct.unpack("i",f.read(recl)[0:4])
print('n_total=',ntotal)
npes=7500
exec(open('./para.out').read())
if len(sys.argv)>1:
  s=int(sys.argv[1])
else:
  s=1

n=0
data=()
for i in range(npes):
#for i in range(5):
  n_in_volume,=struct.unpack("i",f.read(recl)[0:4])
  print('i=',i,'n_in_volume=',n_in_volume)
  data+=struct.unpack("10f"*n_in_volume,f.read(recl*n_in_volume))
  n+=n_in_volume
f.close()
data=np.array(data).reshape(n,10)
print(data[0:10,:])
plt.figure(figsize=(8,4))
#qp=frac/NPPC
if s==1:
  specname='H'
  index=np.where(abs(data[:,6]-0.0140625)<1e-7)
elif s==2:
  specname='4He'
  index=np.where(abs(data[:,6]-6.227375e-3)<1e-7)
elif s==3:
  specname='3He'
  index=np.where(abs(data[:,6]-3.75e-7)<1e-8)
elif s==4:
  specname='O'
  index=np.where(abs(data[:,6]-6.25e-7)<1e-8)
elif s==5:
  specname='Fe'
  index=np.where(abs(data[:,6]-8.75e-7)<1e-8)
else:
  print('wrong argument:',s)
data2=data[index]
labels=['x','y','z','vx','vy','vz']
ax1=plt.subplot(121)
ax2=plt.subplot(122)
for i in range(3):
  #plt.subplot(3,2,i*2+1)
  ax1.hist(data2[:,i],50,histtype='step',label=labels[i])
  #plt.xlabel(labels[i])
  #ax=plt.subplot(3,2,i*2+2)
  ax2.hist(data2[:,i+3],50,histtype='step',label=labels[i+3])
plt.setp(ax2.get_xticklabels(),rotation=30)
ax2.set_xlim(-0.004,0.004)
ax1.grid()
ax2.grid()
ax1.legend()
ax2.legend()
plt.suptitle(r'%s, $t\Omega_i=%d$'%(specname,it*dt))
plt.subplots_adjust(left=0.08,right=0.95,top=0.90,bottom=0.1,wspace=0.2,hspace=0.38)
plt.savefig('%s_%d.png'%(specname,it))
#plt.show()
  

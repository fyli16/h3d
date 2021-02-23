''' Fourier analysis for h3d '''
import numpy as np
import matplotlib.pyplot as plt

def cal_growth(x,y,start,end,ax,color,override=None,factor=1.):
  x=x[start:end]
  y=np.log(y[start:end])
  p=np.polyfit(x,y,1)
  print("Fitting: %s"%p)
  if override != None:
    p[0]=override
  y=np.exp(p[0]*x+p[1])*factor
  ax.plot(x,y,'--',label=r'$\gamma=%.4f$'%(p[0]),color=color)

dt=0.025
nt=101
exec(open('./para.out').read())
time=np.arange(nt)*dt*1000
den=np.load('mode_history_den.npz')
bx=np.load('mode_history_bx.npz')
modes1=bx['modes']
power1=bx['power']
modes2=den['modes']
power2=den['power']
selected=[(0,0,36)]
selected=[(0,1,21)]
selected=[(0,1,2)]
f=plt.figure(figsize=(8,6))
ax1=f.add_subplot(111)
for s in selected:
  for i in range(len(modes1)):
    if s==tuple(modes1[i]):
      line,=ax1.semilogy(time,power1[:,i],label='Bx')
      line,=ax1.semilogy(time,power2[:,i],label='density')
      break
#ax1.set_xlim(0,3000)
ax1.set_ylim(1e5)
ax1.set_ylabel('Fourier power')
ax1.set_title('mode '+str(s))
ax1.grid()
ax1.legend(loc='best')
plt.savefig('mode_%d_%d_%d.png'%(s))
plt.show()


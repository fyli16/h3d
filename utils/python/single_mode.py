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
selected=[(0,0,10),(0,1,11),(0,0,17),(0,0,7)]
selected=[(0,0,10),(0,0,16),(0,0,17),(0,0,18)]
selected=[(0,0,10),(0,0,17),(0,0,7)]
f=plt.figure(figsize=(11,9))
ax1=f.add_subplot(211)
ax2=f.add_subplot(212)
for s in selected:
  for i in range(len(modes1)):
    if s==tuple(modes1[i]):
      line,=ax1.semilogy(time,power1[:,i],label=str(s))
      if s==(0,0,7):
        cal_growth(time,power1[:,i],25,50,ax1,line.get_color(),factor=3.)
      line,=ax2.semilogy(time,power2[:,i],label=str(s))
      if s==(0,0,17):
        cal_growth(time,power2[:,i],50,80,ax2,line.get_color(),factor=2.)
      break
ax1.set_title('Bx')
ax1.set_ylabel('Fourier power')
ax1.grid()
ax1.legend(loc='best')
ax2.set_title('density')
ax2.set_ylabel('Fourier power')
ax2.set_xlabel(r'$t\Omega_i$')
ax2.grid()
ax2.legend(loc='best')
plt.savefig('modes.png')
plt.show()


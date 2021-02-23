import numpy as np

import matplotlib.pyplot as plt

run_name='055'
dt=0.025
nt=20
exec(open('./para.out').read())
time=np.arange(nt)*dt*1000

f=plt.figure(figsize=(8,8))
ax=plt.subplot(111)
f2=plt.figure(figsize=(8,8))
ax2=plt.subplot(111)
ref_time=8
for s in range(1,nspec+1):
  tpara=wspec[s-1]*np.fromfile('tpar_%s_h3d-%s.dat'%(s,run_name))
  tperp=wspec[s-1]*np.fromfile('tperp_%s_h3d-%s.dat'%(s,run_name))
  line,=ax.plot(time,tperp[:nt],label=r'$T_{%s\perp}$'%s)
  ax.plot(time,tpara[:nt],ls='--',color=line.get_color(),label=r'$T_{%s\parallel}$'%s)
  ax2.plot(time[ref_time:nt],tperp[ref_time:nt]/tperp[ref_time],label=r'$T_{%s\perp}$'%s)
ax.grid()
ax.set_xlabel(r'$t\Omega_i$')
ax.legend(loc='best')
ax.set_yscale('log')
ax.set_ylabel(r'$T_{s,\perp}/T_{s,\perp(0)}$')
f.savefig('temp.png')
f.show()
ax2.grid()
ax2.set_xlabel(r'$t\Omega_i$')
ax2.legend(loc='best')
ax2.set_title(r'$T_{s,\perp}/T_{s,\perp}(t=%d)$'%time[ref_time])
ax2.set_yscale('log')
ax2.set_xlim(xmin=0)
f2.savefig('temp2.png')
f2.show()




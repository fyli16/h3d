from common_modules import *
tweak_rcParams(fs=12, home='/home1/07918/fyli')
mpl.use('TKAgg')

# import numpy as np
# import matplotlib.pyplot as plt

# run_name='h3d-106'
# dt=0.025
# nt=20
# ifshow=True
exec(open('./para.out').read())
time = np.arange(nt)*dt*tinterval
cpath = 'cdata'  # compiled data path


print('plotting E/B energy ...')
fig, ax1 = plt.subplots(1,1, figsize=[5, 3.5])
b = np.fromfile('{}/magnetic_energy_{}.dat'.format(cpath, run_name) )
e = np.fromfile('{}/electric_energy_{}.dat'.format(cpath, run_name) )
ax1.plot(time, b, label=r'$\delta B^2$')
ax1.plot(time, e, label=r'$\delta E^2$')
ax1.grid()
ax1.legend()
ax1.set_xlabel(r'$t\Omega_{ci}$')
plt.tight_layout()
plt.savefig(cpath+'/fields.png')
if ifshow: plt.show()
plt.close('all')


print('plotting temperature ...')
fig, ax1 = plt.subplots(1,1, figsize=[5, 3.5])
# plt.figure(figsize=(8,8))
# ax=plt.subplot(111)
for s in range(1, nspec+1):
  tpara = np.fromfile('{}/tpar_{}_{}.dat'.format(cpath, s, run_name) )
  tperp = np.fromfile('{}/tperp_{}_{}.dat'.format(cpath, s, run_name) )
  line, = ax1.plot(time, tperp, label=r'$T_{%s\perp}$'%s)
  ax1.plot(time, tpara, ls='--', color=line.get_color(), label=r'$T_{%s\parallel}$'%s)
ax1.grid()
ax1.set_xlabel(r'$t\Omega_{ci}$')
ax1.legend(loc='best')
#ax1.set_yscale('log')
plt.tight_layout()
plt.savefig(cpath+'/temp.png')
if ifshow: plt.show()


print('plotting den_fluc ...')
den_fluc = np.fromfile('{}/den_fluc_{}.dat'.format(cpath, run_name) )
# plt.cla()
fig, ax1 = plt.subplots(1,1, figsize=[5, 3.5])
ax1.plot(time, den_fluc, label=r'$\delta\rho^2$')
ax1.legend()
ax1.grid()
ax1.set_xlabel(r'$t\Omega_{ci}$')
plt.tight_layout()
plt.savefig(cpath+'/den_fluc.png')
if ifshow: plt.show()
plt.close('all')

print ('done!')


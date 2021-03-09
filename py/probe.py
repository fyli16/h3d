from h3d import *

def load_input(path):
    global denmin, nspec, nx, ny, nz, xmax, ymax, zmax, \
        tmax, dt, nwrtdata, y, z, \
        ndumps, timesteps, times 

    nml = f90nml.read( join(path,'input.f90') )
    denmin = nml['datum']['denmin']
    nspec = nml['datum']['nspec']
    nx = nml['datum']['nx']
    ny = nml['datum']['ny']
    nz = nml['datum']['nz']
    xmax = nml['datum']['xmax']
    ymax = nml['datum']['ymax']
    zmax = nml['datum']['zmax']
    tmax = nml['datum']['tmax']
    dt = nml['datum']['dtwci']
    nwrtdata = nml['datum']['nwrtdata']
    # box axes
    y = np.linspace(0, ymax, ny)
    z = np.linspace(0, zmax, nz)
    # timesteps
    ndumps = int(tmax/(dt*nwrtdata)+1)
    timesteps = np.zeros(ndumps)
    for i in range(ndumps):
        timesteps[i]=1 if i==0 else i*nwrtdata
    times = timesteps * dt
    return

# it, ex, ey, ez, bx, by, bz
# 0,  1,  2,  3,  4,  5,  6

path='1d-b0.5'
load_input(path)
myid=28
data = np.loadtxt(path+'/data/probes_%.04d.dat'%myid)
it = data[:,0]; t = it*dt
ex, ey, ez = data[:,1], data[:,2], data[:,3]
bx, by, bz = data[:,4], data[:,5], data[:,6]
RX = 0.5*(ex+by)
LX = 0.5*(ex-by)

w0wci=5*twopi/224
coeff = twopi/w0wci

print ('plotting ...')
fig, axes = plt.subplots(3,2, figsize=[9,5.5],)#sharex=True)
ax1 = axes[0,0]
ax1.plot(t, ex, '--', linewidth=0.5, label='ex')
ax1.plot(t, by, '--', linewidth=0.5, label='by')
ax1.legend(loc='upper right')
# ax1.set_xlabel(r'$\omega_{ci}t$')
# ax1.set_ylabel(r'$b_x$')
ax1 = axes[0,1]
f, F = FFT(t, bx)
ax1.plot(f*coeff, F)
ax1.set_xlim(0,2)

ax1 = axes[1,0]
ax1.plot(t, RX)
ax1.set_ylabel(r'$(e_x+b_y)/2$')
ax1 = axes[1,1]
f, F = FFT(t, RX)
ax1.plot(f*coeff, F)
# idx1, idx2 = int(0/dt), int(448/dt)
# f, F = FFT(t[idx1:idx2], RX[idx1:idx2])
# ax1.plot(f*coeff, F)
ax1.set_xlim(0,2)
ax1.set_ylabel('spec. amp.')

ax1 = axes[2,0]
ax1.plot(t, LX)
ax1.set_ylabel(r'$(e_x-b_y)/2$')
ax1.set_xlabel(r'$\omega_{ci}t$')
ax1 = axes[2,1]
f, F = FFT(t, LX)
ax1.plot(f*coeff, F)
# idx1, idx2 = int(1552/dt), int(2000/dt)
# f, F = FFT(t[idx1:idx2], LX[idx1:idx2])
# ax1.plot(f*coeff, F)
ax1.set_xlim(0,2)
ax1.set_xlabel(r'$\omega/\omega_0$')

plt.tight_layout()
plt.show()
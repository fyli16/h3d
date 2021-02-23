from h3d import *

def load_input(path):
    global denmin, nspec, nx, ny, nz, xmax, ymax, zmax, \
        tmax, dt, nwrtdata, y, z, \
        ndumps, timesteps, times 

    nml = f90nml.read(path + '/data/finput.f90')
    denmin = nml['datum']['denmin']
    nspec = nml['datum']['nspec']
    nx = nml['datum']['nx']
    ny = nml['datum']['ny']
    nz = nml['datum']['nz']
    xmax = nml['datum']['xmax']
    ymax = nml['datum']['ymax']
    zmax = nml['datum']['zmax']
    tmax = nml['datum']['maximum_simulation_time']
    dt = nml['datum']['dtwci']
    nwrtdata = nml['datum']['nwrtdata']
    #@ box axes
    y = np.linspace(0, ymax, ny)
    z = np.linspace(0, zmax, nz)
    #@ timesteps
    ndumps = int(tmax/(dt*nwrtdata)+1)
    timesteps = np.zeros(ndumps)
    for i in range(ndumps):
        timesteps[i]=1 if i==0 else i*nwrtdata
    times = timesteps * dt
    return


def get_snapshot_data(ts, path, field):
    fname = '%s_%d.gda' % (field, ts)
    fname = join(path, 'data', fname)
    data = np.fromfile(fname, dtype=np.float32)
    data = data.reshape(nz, ny).transpose()
    return (data)


path='.'
field = 'den'
load_input(path)


# data = np.zeros( (ndumps, len(path)) )
# mean = np.zeros( (ndumps, len(path)) )
# for i, p in enumerate(path):
#     for j, ts in enumerate(timesteps):
#         showProgressBar(ndumps-1, j)
#         # d = get_snapshot_data(p, field, ts)
#         d = get_snapshot_data(ts, p, field)
#         # data[j,i] = np.sqrt(np.sum((d-np.mean(d))**2.)/np.size(d))
#         data[j,i] = np.sqrt(np.sum((d-1)**2.)/np.size(d))
#         mean[j,i] = np.mean(d)

# fig, ax1 = plt.subplots(1,1, figsize=[6, 4])
# for i in range(len(path)):
#     ax1.plot(times, data[:,i], label='run #%d'%i)
# ax1.legend()
# ax1.set_xlabel(r'$\omega_{ci} t$')
# # ax1.set_ylabel(r'$\delta\rho^2$')
# # ax1.set_ylabel(r'$\sqrt{\langle(\rho-\langle\rho\rangle)^2\rangle}$')
# ax1.set_ylabel(r'$\sqrt{\langle(\rho-1)^2\rangle}$')
# plt.tight_layout()
# plt.show()
# plt.close('all')

# fig, ax1 = plt.subplots(1,1, figsize=[6, 4])
# for i in range(len(path)):
#     ax1.plot(times, mean[:,i], label='run #%d'%i)
# ax1.legend()
# ax1.set_xlabel(r'$\omega_{ci} t$')
# # ax1.set_ylabel(r'$\delta\rho^2$')
# # ax1.set_ylabel(r'$\sqrt{\langle(\rho-\langle\rho\rangle)^2\rangle}$')
# ax1.set_ylabel(r'$\langle\rho\rangle$')
# plt.tight_layout()
# plt.show()
# plt.close('all')

        



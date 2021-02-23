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
    # fname = a.path+'/data/'+a.field+'_'+str(ts)+'.gda'
    fname = '%s_%d.gda' % (field, ts)
    fname = join(path, 'data', fname)
    data = np.fromfile(fname, dtype=np.float32)
    # data = data.reshape(nz, ny).transpose()
    return (data)


load_input('.')
#            0      1      2    3     4     5     6       7        8
field_list=['den', 'bx', 'by', 'bz', 'ex', 'ey', 'ez', 'tpar_1', 'tperp_1']
dsize = ny*nz
data = np.zeros( (ndumps, dsize, len(field_list)) )
# tpar = np.zeros( ndumps )
# tperp = np.zeros( ndumps )
# e_ene = np.zeros( ndumps )
# b_ene = np.zeros( ndumps )
mkdir('hist')
datafile = join('hist', 'hist.npy')
if not exists(datafile):
    for i, field in enumerate(field_list):
        for j, ts in enumerate(timesteps):
            showProgressBar(ndumps-1, j)
            data[j,:,i] = get_snapshot_data(ts, '.', field)
            # data[j,i] = np.sqrt(np.sum((d-np.mean(d))**2.)/np.size(d))
            # data[j,i] = np.sqrt(np.sum((d-1)**2.)/np.size(d))
            # mean[j,i] = np.mean(d)
    # data.tofile(datafile)
    np.save(datafile, data)
else:
    # data = np.fromfile(datafile)
    data = np.load(datafile)

#               0      1      2    3     4     5     6       7        8
# field_list=['den', 'bx', 'by', 'bz', 'ex', 'ey', 'ez', 'tpar_1', 'tperp_1']
delta_rho = np.zeros(ndumps) 
bperp2 = np.zeros(ndumps)
delta_bz2 =  np.zeros(ndumps)
eperp2 = np.zeros(ndumps)
ez2= np.zeros(ndumps)
tpar = np.zeros(ndumps)
tperp = np.zeros(ndumps)
for i in range(ndumps):
    delta_rho[i]=np.sqrt(np.sum((data[i,:,0]-1.)**2.)/dsize)
    bperp2[i]=0.5*np.sum(data[i,:,1]**2+data[i,:,2]**2)/dsize
    delta_bz2[i]=0.5*np.sum((data[i,:,3]-1)**2)/dsize
    eperp2[i]=0.5*np.sum(data[i,:,4]**2+data[i,:,5]**2)/dsize
    ez2[i]=0.5*np.sum(data[i,:,6]**2)/dsize
    tpar[i]=np.mean(data[i,:,7])
    tperp[i]=np.mean(data[i,:,8])

fig, axes = plt.subplots(3,1, figsize=[5,7], sharex=True)
# # ax1.set_title(a.fieldname)
# # ax1.plot(times, delta_rho, label=r'$\delta\rho$')
# ax1.plot(times, bperp2, )
# # ax1.legend()
# ax1.set_xlabel(r'$\omega_{ci} t$')
# # ax1.set_ylabel(r'$\sqrt{\langle(\rho-1)^2\rangle}$')
ax1=axes[0]
ax1.plot(times, delta_rho, label=r'$\sqrt{(\delta\rho)^2}$'); 
ax1.legend()#loc='lower right')
ax1=axes[1]
ax1.plot(times, tpar, label=r'$\langle T_\parallel\rangle$'); 
ax1.plot(times, tperp, label=r'$\langle T_\perp\rangle$'); 
ax1.legend()#loc='lower right')
ax1=axes[2]
ax1.plot(times, bperp2, label=r'$\frac{1}{2}\langle B_\perp^2\rangle$'); 
# ax1.plot(times, delta_bz2, label=r'$(\delta B_z)^2$'); 
# ax1.legend(loc='lower left')
# ax1=axes[2]; 
ax1.plot(times, eperp2, label=r'$\frac{1}{2}\langle E_\perp^2\rangle$'); 
# ax1.plot(times, ez2, label=r'$E_z^2$'); 
ax1.legend()#loc='upper right')
ax1.set_ylim(2e-3, 6e-3)
ax1.set_xlabel(r'$\omega_{ci} t$')
plt.tight_layout()
plt.show()


    



from h3d import *

def load_input(path):
    global denmin, nspec, nx, ny, nz, xmax, ymax, zmax, \
        tmax, dt, nwrtdata, y, z, \
        ndumps, timesteps, times 

    nml = f90nml.read(path + '/input.f90')
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

# read sim data
def read_sim_data(path):
    load_input(path)
    dsize = ny*nz
    data = np.zeros( (ndumps, dsize, len(field_list)) )
    subpath = join(path, 'hist')
    mkdir(subpath)
    datafile = join(subpath, 'hist.npy')
    if not exists(datafile):
        for i, field in enumerate(field_list):
            for j, ts in enumerate(timesteps):
                showProgressBar(ndumps-1, j)
                data[j,:,i] = get_snapshot_data(ts, path, field)
        np.save(datafile, data)
    else:
        data = np.load(datafile)
    print ('sim data extracted')
    return (data)

def get_spectral_max(path, spec_rho, spec_rho_max, k_pos):
    data = read_sim_data(path)
    for m in range(ndumps):
        rho=data[m,:,0] 
        rho=rho.reshape(nz,ny).transpose()
        mean=np.mean(rho,axis=0)  # averaged over y
        f, F = FFT(z, mean-1)
        spec_rho[m,:] = F
        spec_rho_max[m] = F.max()
        k_pos[m]=f[np.argmax(F)]
    return (spec_rho, spec_rho_max, k_pos)

# path_list = ['test', '1d-resis0.001']
# label_list = ['resis=1e-6', 'resis=1e-3']

# path_list = ['1d-b0.05', '1d-b0.1', '1d-b0.2', '1d-b0.5']
# label_list = ['b0=0.05', 'b0=0.1', 'b0=0.2', 'b0=0.5']
b0_list=[0.001, 0.01, 0.05, 0.1, 0.2, 0.5]
b0_list2=[0.05, 0.1, 0.2, 0.5]

# path_list = ['run1', 'run2', 'run3']
# label_list = path_list

#            0      1      2    3     4     5     6       7        8
field_list=['den', 'bx', 'by', 'bz', 'ex', 'ey', 'ez', 'tpar_1', 'tperp_1']

fig, axes = plt.subplots(1,1, figsize=[5,3.8])
ax1=axes

for i, b0 in enumerate(b0_list2):
    path = '1d-b%s'%b0
    print(path)
    load_input(path)
    spec_rho=np.zeros( (ndumps, int(nz/2)) )
    spec_rho_max =  np.zeros( ndumps )
    k_pos =  np.zeros( ndumps )
    spec_rho, spec_rho_max, k_pos = get_spectral_max(path, 
                            spec_rho, spec_rho_max, k_pos)
    ax1.semilogy(times, spec_rho_max, 'o', markersize=3, 
                markerfacecolor='none', label=r'$b_0=%s$'%b0_list2[i])
ax1.legend()
ax1.set_xlabel(r'$\omega_{ci}t$')
ax1.set_ylabel(r'$F(\delta\rho)$')

# # exponential fit
# t0, t1 = 300., 700.
# id1, id2 = int(t0/(dt*nwrtdata)), int(t1/(dt*nwrtdata))
# pidx=np.polyfit(times[id1:id2],np.log(spec_rho_max[id1:id2]),1)
# yfit=np.exp(times[id1:id2]*pidx[0]+pidx[1])
# ax1.semilogy(times[id1:id2], yfit, 'r-', label=r'$F_{max}$')

plt.tight_layout()
plt.show()
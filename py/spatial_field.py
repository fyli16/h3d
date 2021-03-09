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

def get_snapshot_data(ts, path, field):
    # fname = a.path+'/data/'+a.field+'_'+str(ts)+'.gda'
    fname = '%s_%d.gda' % (field, ts)
    fname = join(path, 'data', fname)
    data = np.fromfile(fname, dtype=np.float32)
    # data = data.reshape(nz, ny).transpose()
    return (data)

# read sim data
def read_sim_data(path):
    print (path)
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
    return (data)
    
    # for m in range(ndumps):
    #     rho=data[m,:,0] 
    #     rho=rho.reshape(nz,ny).transpose()
    #     mean=np.mean(rho,axis=0)  # averaged over y
    #     f, F = FFT(z, mean-1)
    #     spec_rho[m,:,k] = F
    #     spec_rho_max[m,k] = F.max()
    #     k_pos[m,k]=f[np.argmax(F)]

def get_field_mean(path, fieldname, t):
    load_input(path)
    idx = np.argmin(np.abs(times-t))
    data = read_sim_data(path)
    d_ = data[idx,:,]



#            0      1      2    3     4     5     6       7        8
field_list=['den', 'bx', 'by', 'bz', 'ex', 'ey', 'ez', 'tpar_1', 'tperp_1']

plot_type = 'snaptshot'

if plot_type=='snapshot':
    path='1d-b0.1'
    load_input(path)
    t=100
    idx = np.argmin(np.abs(times-t))
    data = read_sim_data(path)
    ex_ = data[idx,:,4]
    ex_ = ex_.reshape(nz,ny).transpose()
    ex = np.mean(ex_,axis=0)
    by_ = data[idx,:,2]
    by_ = by_.reshape(nz,ny).transpose()
    by = np.mean(by_,axis=0)
    RX = 0.5*(ex+by)
    LX = 0.5*(ex-by)

    w0wci=5*twopi/224
    coeff = twopi/w0wci

    print ('plotting ...')
    fig, axes = plt.subplots(3,2, figsize=[9,5.5],)#sharex=True)
    ax1 = axes[0,0]
    ax1.plot(z, ex, '--', linewidth=0.5, label='ex')
    ax1.plot(z, by, '--', linewidth=0.5, label='by')
    ax1.legend(loc='upper right')
    # ax1.set_xlabel(r'$\omega_{ci}t$')
    # ax1.set_ylabel(r'$b_x$')
    ax1 = axes[0,1]
    f, F = FFT(z, by)
    ax1.plot(f*coeff, F)
    ax1.set_xlim(0,2)

    ax1 = axes[1,0]
    ax1.plot(z, RX)
    ax1.set_ylabel(r'$(e_x+b_y)/2$')
    ax1 = axes[1,1]
    f, F = FFT(z, RX)
    ax1.plot(f*coeff, F)
    # idx1, idx2 = int(0/dt), int(448/dt)
    # f, F = FFT(t[idx1:idx2], RX[idx1:idx2])
    # ax1.plot(f*coeff, F)
    ax1.set_xlim(0,2)
    ax1.set_ylabel('spec. amp.')

    ax1 = axes[2,0]
    ax1.plot(z, LX)
    ax1.set_ylabel(r'$(e_x-b_y)/2$')
    ax1.set_xlabel(r'$z/d_i$')
    ax1 = axes[2,1]
    f, F = FFT(z, LX)
    ax1.plot(f*coeff, F)
    # idx1, idx2 = int(1552/dt), int(2000/dt)
    # f, F = FFT(t[idx1:idx2], LX[idx1:idx2])
    # ax1.plot(f*coeff, F)
    ax1.set_xlim(0,2)
    ax1.set_xlabel(r'$k/k_0$')

    plt.tight_layout()
    plt.show()
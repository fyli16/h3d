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

def get_field_idx(fieldname):
    idx=0
    for i, field in enumerate(field_list):
        if field==fieldname:
            idx=i
            break
    return (idx)

def get_field_mean(data, fieldname, idump):
    d_ = data[idump,:,get_field_idx(fieldname)]
    d_ = d_.reshape(nz,ny).transpose()
    d = np.mean(d_, axis=0)
    return (d)

def get_spectral_max(path, spec, spec_max, k_pos):
    data = read_sim_data(path)
    for m in range(ndumps):
        ex = get_field_mean(data, 'ex', m)
        by = get_field_mean(data, 'by', m)
        # RX = 0.5*(ex+by)
        LX = 0.5*(ex-by)
        f, F = FFT(z, LX)
        spec[m,:] = F
        spec_max[m] = F.max()
        k_pos[m] = f[np.argmax(F)]
    return (spec, spec_max, k_pos)



#            0      1      2    3     4     5     6       7        8
field_list=['den', 'bx', 'by', 'bz', 'ex', 'ey', 'ez', 'tpar_1', 'tperp_1']

# plot_type = 'snapshot'
plot_type = 'growth'

path='1d-b0.1'; t=100
load_input(path)
w0wci = 5*twopi/224
coeff = twopi/w0wci

if plot_type=='snapshot':
    ex = get_field_mean(path, 'ex', t)
    by = get_field_mean(path, 'by', t)
    RX = 0.5*(ex+by)
    LX = 0.5*(ex-by)

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

elif plot_type == 'growth':
    b0_list = [0.05, 0.1, 0.2, 0.5]

    tends=np.array([[800,1300],[350,550],[200,300],[100,200]])
    growth = np.zeros(len(b0_list))
    growth[:]=np.nan

    fig, axes = plt.subplots(1,1, figsize=[5,3.8])
    ax1=axes

    for i, b0 in enumerate(b0_list):
        path = '1d-b%s'%b0
        load_input(path)
        spec = np.zeros( (ndumps, int(nz/2)) )
        spec_max =  np.zeros( ndumps )
        k_pos =  np.zeros( ndumps )
        spec, spec_max, k_pos = get_spectral_max(path, 
                                spec, spec_max, k_pos)
        line, = ax1.semilogy(times, spec_max, 'o', markersize=3, 
                markerfacecolor='none', label=r'$b_0=%s$'%b0_list[i])
        
        t0, t1 = tends[i,0], tends[i,1]
        if not (t0==0 and t1==0):
            id1, id2 = int(t0/(dt*nwrtdata)), int(t1/(dt*nwrtdata))
            pidx=np.polyfit(times[id1:id2],np.log(spec_max[id1:id2]),1)
            yfit=np.exp(times[id1:id2]*pidx[0]+pidx[1])
            growth[i] = pidx[0]
            ax1.semilogy(times[id1:id2], yfit, '-', color=line.get_color())
        else:
            growth[i] = np.nan


    ax1.legend()
    ax1.set_xlabel(r'$\omega_{ci}t$')
    ax1.set_ylabel(r'$F(L_x)$')

    np.save('growth_rates_Lx.npy', growth)

    plt.tight_layout()
    plt.show()
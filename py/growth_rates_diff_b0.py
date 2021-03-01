from h3d import *

def load_input(path):
    global denmin, nspec, nx, ny, nz, xmax, ymax, zmax, \
        tmax, dt, nwrtdata, y, z, \
        ndumps, timesteps, times 

    nml = f90nml.read(path + '/data/input.f90')
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

b0_list=[0.001, 0.01, 0.05, 0.1, 0.2, 0.5]
load_input('test')  # get ndumps, nz
spec_rho=np.zeros( (ndumps, int(nz/2), len(b0_list)) )
spec_rho_max =  np.zeros( (ndumps, len(b0_list)) )
k_pos =  np.zeros( (ndumps, len(b0_list)) )
#            0      1      2    3     4     5     6       7        8
field_list=['den', 'bx', 'by', 'bz', 'ex', 'ey', 'ez', 'tpar_1', 'tperp_1']

# read sim data
for k, b0 in enumerate(b0_list):
    path='1d-b%s'%b0
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
    
    for m in range(ndumps):
        rho=data[m,:,0] 
        rho=rho.reshape(nz,ny).transpose()
        mean=np.mean(rho,axis=0)  # averaged over y
        f, F = FFT(z, mean-1)
        spec_rho[m,:,k] = F
        spec_rho_max[m,k] = F.max()
        k_pos[m,k]=f[np.argmax(F)]

#==========================================================#
#>>   measured growth against theory
#==========================================================#
tends=np.array([[0,0],[0,0],[500,1200],[200,550],[100,300],[0,150]])
growth = np.zeros(len(b0_list))
growth[:]=np.nan
print ('plotting ...')
fig, axes = plt.subplots(1,2, figsize=[10,3.8])
ax1=axes[0]
for i in range(len(b0_list)):
    line,= ax1.semilogy(times, spec_rho_max[:,i], 'o', markersize=3, 
                markerfacecolor='none', label=r'$b_0=%s$'%b0_list[i])
    ax1.set_xlabel(r'$\omega_{ci}t$')
    ax1.set_ylabel('Spectral max.')
    t0, t1 = tends[i,0], tends[i,1]
    if not (t0==0 and t1==0):
        print (b0_list[i])
        id1, id2 = int(t0/(dt*nwrtdata)), int(t1/(dt*nwrtdata))
        pidx=np.polyfit(times[id1:id2],np.log(spec_rho_max[id1:id2,i]),1)
        yfit=np.exp(times[id1:id2]*pidx[0]+pidx[1])
        growth[i] = pidx[0]
        ax1.semilogy(times[id1:id2], yfit, '-', color=line.get_color())
ax1.legend(loc='upper right')
#@ append analytical results
ax1 = axes[1]
gmax=np.load('test/h3dtest-gmax.npy')
b0_theory=np.load('test/h3dtest-b0.npy')
from math import e
convert=5*twopi/224
ax1.plot(np.array(b0_list), growth/convert, 'ko', 
        markersize=7, markerfacecolor='none',label='sim.')
# ax2 = ax1.twinx()
ax1.plot(b0_theory[1:], gmax[1:], 'r-', label='theory')
ax1.legend()
ax1.set_xlim(0,.6)
ax1.set_xlabel(r'$b_0$')
ax1.set_ylabel(r'$\gamma_{max}/\omega_0$')
# ax1.set_ylabel(r'analytical $\gamma_{max}$')
# ax2.yaxis.label.set_color('r')
# ax2.tick_params(axis='y',colors='r')
plt.tight_layout()
plt.show()


#==========================================================#
#>>   history of quantities
#==========================================================#
#               0      1      2    3     4     5     6       7        8
# field_list=['den', 'bx', 'by', 'bz', 'ex', 'ey', 'ez', 'tpar_1', 'tperp_1']
# delta_rho = np.zeros(ndumps) 
# bperp2 = np.zeros(ndumps)
# delta_bz2 =  np.zeros(ndumps)
# eperp2 = np.zeros(ndumps)
# ez2= np.zeros(ndumps)
# tpar = np.zeros(ndumps)
# tperp = np.zeros(ndumps)
# for i in range(ndumps):
#     delta_rho[i]=np.sqrt(np.sum((data[i,:,0]-1.)**2.)/dsize)
#     bperp2[i]=0.5*np.sum(data[i,:,1]**2+data[i,:,2]**2)/dsize
#     delta_bz2[i]=0.5*np.sum((data[i,:,3]-1)**2)/dsize
#     eperp2[i]=0.5*np.sum(data[i,:,4]**2+data[i,:,5]**2)/dsize
#     ez2[i]=0.5*np.sum(data[i,:,6]**2)/dsize
#     tpar[i]=np.mean(data[i,:,7])
#     tperp[i]=np.mean(data[i,:,8])

# fig, axes = plt.subplots(3,1, figsize=[5,7], sharex=True)
# ax1=axes[0]
# ax1.plot(times, delta_rho, label=r'$\sqrt{(\delta\rho)^2}$'); 
# ax1.legend()#loc='lower right')
# ax1=axes[1]
# ax1.plot(times, tpar, label=r'$\langle T_\parallel\rangle$'); 
# ax1.plot(times, tperp, label=r'$\langle T_\perp\rangle$'); 
# ax1.legend()#loc='lower right')
# ax1=axes[2]
# ax1.plot(times, bperp2, label=r'$\frac{1}{2}\langle B_\perp^2\rangle$'); 
# ax1.plot(times, eperp2, label=r'$\frac{1}{2}\langle E_\perp^2\rangle$'); 
# ax1.legend()#loc='upper right')
# ax1.set_ylim(2e-3, 6e-3)
# ax1.set_xlabel(r'$\omega_{ci} t$')
# plt.tight_layout()
# plt.show()

# spec_rho=np.zeros( (ndumps, int(nz/2)) )
# spec_rho_max =  np.zeros(ndumps)
# k_pos =  np.zeros(ndumps)
# for i in range(ndumps):
#     showProgressBar(ndumps-1, i)
#     rho=data[i,:,0]
#     rho=rho.reshape(nz,ny).transpose()
#     mean=np.mean(rho,axis=0)
#     f, F = FFT(z, mean-1)
#     spec_rho[i,:] = F
#     spec_rho_max[i] = F.max()
#     k_pos[i]=f[np.argmax(F)]

if 0:
    fig, axes = plt.subplots(2,1, figsize=[5,7])
    # ax1.set_title('%s, t = %.2f'%(a.fieldname, ts*dt))
    ax1=axes[0]
    ax1.loglog(times, spec_rho_max, 'k', label=r'$F_{max}$')
    ax2 =ax1.twinx()
    ax2.semilogx(times, k_pos, 'r', label=r'$k_{max}$')
    ax1.set_xlabel(r'$\omega_{ci}t$')
    ax1.set_ylabel('Spectral max.')
    ax2.set_ylabel('k of spectral max.')
    ax2.yaxis.label.set_color('r')
    ax2.tick_params(axis='y',colors='r')
    # ax1.set_xscale('log')
    # ax1.set_yscale('log')
    ax1.set_xlim(100,)

    ax1=axes[1]
    im1 = ax1.imshow(spec_rho, aspect='auto', origin='lower', cmap=jbew,
            extent=[0,f.max(), 0, tmax])
    # plt.colorbar(im1, ax=ax1, pad=0.02)
    cbax = ax1.inset_axes([1.01, 0, 0.03, 1])
    cb = plt.colorbar(im1, cax=cbax, )
    ax1.set_xlabel('k')
    ax1.set_ylabel(r'$\omega_{ci}t$')

    plt.tight_layout()
    plt.show()

if 0:
    # fig, axes = plt.subplots(2,1, figsize=[5,7])
    fig, axes = plt.subplots(1,1, figsize=[5,3.8])
    # ax1.set_title('%s, t = %.2f'%(a.fieldname, ts*dt))
    ax1=axes
    ax1.semilogy(times, spec_rho_max, 'ko', markersize=3, markerfacecolor='none', label=r'$F_{max}$')
    # ax2 =ax1.twinx()
    # ax2.semilogx(times, k_pos, 'r', label=r'$k_{max}$')
    ax1.set_xlabel(r'$\omega_{ci}t$')
    ax1.set_ylabel('Spectral max.')
    t0, t1 = 300., 700.
    id1, id2 = int(t0/(dt*nwrtdata)), int(t1/(dt*nwrtdata))
    print (times.shape)
    print(id1, id2)
    p=np.polyfit(times[id1:id2],np.log(spec_rho_max[id1:id2]),1)
    yfit=np.exp(times[id1:id2]*p[0]+p[1])
    ax1.semilogy(times[id1:id2], yfit, 'r-', label=r'$F_{max}$')
    # ax2.set_ylabel('k of spectral max.')
    # ax2.yaxis.label.set_color('r')
    # ax2.tick_params(axis='y',colors='r')
    # ax1.set_xscale('log')
    # ax1.set_yscale('log')
    # ax1.set_xlim(100,)
    plt.tight_layout()
    plt.show()







    



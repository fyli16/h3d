from common_modules import *
tweak_rcParams(fs=14, home='/home1/07918/fyli')
mpl.use('TKAgg')

## input parameters
# exec(open('./para.out').read())
import f90nml
fpath = '../data/'
nml = f90nml.read(fpath + 'finput.dat')
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
tinterval = nml['datum']['nwrtdata']
nt = int(tmax/(dt*tinterval)+1)
print ('total number of dump steps: ', nt)
time = np.arange(nt)*dt*tinterval


time = np.arange(nt) * dt * tinterval
y = np.linspace(0, ymax, ny)
z = np.linspace(0, zmax, nz)
fpath = '../data/'

def plot_den_single(i=0):
    tframe = 1 if i == 0 else i * tinterval
    print (i, tframe)
    fname = fpath + 'den_' + str(tframe) + '.gda'
    rho = np.fromfile(fname, dtype=np.float32)
    # sz, = rho.shape
    # rho_avg = np.mean(rho)
    # den_fluc[ct] = np.sum((rho - rho_avg)**2) / sz
    # rho = rho.reshape(ny, nz)

    fig, axes = plt.subplots(2,1, figsize=[10,4], sharex=True)
    ax1 = axes[0]
    rho1 = rho.reshape(nz, ny).transpose()
    im1 = ax1.imshow(rho1, aspect='auto', origin='lower', cmap='jet', extent=[0, zmax, 0, ymax])
    # plt.colorbar(im1, ax=ax1, pad=0.02)
    cbax = ax1.inset_axes([1.01, 0, 0.02, 1])
    cb = plt.colorbar(im1, cax=cbax, )#orientation="vertical")
    # ax1.set_xlabel('y')
    # ax1.set_ylabel('z')

    ax1 = axes[1]
    rho_trans_avg = np.mean(rho1, axis=0)
    ax1.plot(z, rho_trans_avg, '-')

    plt.tight_layout()
    plt.show()


def get_trans_mean(field, tframe, spec='none'):
    if spec=='none':
        fname = fpath + field + '_' + str(tframe) + '.gda'
    else:
        fname = fpath + field + '_' + str(spec) + '_' + str(tframe) + '.gda'
    data = np.fromfile(fname, dtype=np.float32)
    data = data.reshape(nz, ny).transpose()
    trans_avg = np.mean(data, axis=0)
    return (trans_avg)

def ebden_eta_single(i):
    imag_path = 'ebden_eta'
    mkdir(imag_path)
    step = 1 if i == 0 else i * tinterval
    print (step)
    bx = get_trans_mean('bx', step)
    by = get_trans_mean('by', step)
    bz = get_trans_mean('bz', step)
    ex = get_trans_mean('ex', step)
    ey = get_trans_mean('ey', step)
    ez = get_trans_mean('ez', step)
    den = get_trans_mean('den', step)
    eta = get_trans_mean('eta', step)

    fig, axes = plt.subplots(4,2, figsize=[8,7], sharex=True)
    fig.suptitle('t={:.1f}'.format(step*dt), fontsize=16)
    ax1 = axes[0,0]; ax1.plot(bx, 'C0'); ax1.set_ylabel(r'$B_x$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1 = axes[1,0]; ax1.plot(by, 'C0'); ax1.set_ylabel(r'$B_y$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1 = axes[2,0]; ax1.plot(bz, 'C0'); ax1.set_ylabel(r'$B_z$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1 = axes[0,1]; ax1.plot(ex, 'C0'); ax1.set_ylabel(r'$E_x$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1 = axes[1,1]; ax1.plot(ey, 'C0'); ax1.set_ylabel(r'$E_y$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1 = axes[2,1]; ax1.plot(ez, 'C0'); ax1.set_ylabel(r'$E_z$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1 = axes[3,0]; ax1.plot(den, 'C0'); ax1.set_ylabel(r'$n$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1.set_xlabel('z'); ax1.set_xlim(0, zmax)
    ax1 = axes[3,1]; ax1.plot(eta, 'C0'); ax1.set_ylabel(r'$\eta$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1.set_xlabel('z'); ax1.set_xlim(0, zmax)
    plt.tight_layout()
    plt.savefig(imag_path+'/ebden_eta_{:08d}'.format(step))
    plt.close('all')

def ebden_eta_movie():
    processor = functools.partial(ebden_eta_single)
    dump_num = np.arange(201)
    with mp.Pool(processes=20) as worker_pool:
        frames = worker_pool.map(processor, dump_num)
        worker_pool.close()
        worker_pool.join()
    makeMovie('ebden_eta', 'ebden_eta.mp4', fps=3)
        


def pt_single(i):
    imag_path = 'pt'
    mkdir(imag_path)
    step = 1 if i == 0 else i * tinterval
    print (step)
    pxx = get_trans_mean('p-xx', step, 1)
    pxy = get_trans_mean('p-xy', step, 1)
    pxz = get_trans_mean('p-xz', step, 1)
    pyy = get_trans_mean('p-yy', step, 1)
    pyz = get_trans_mean('p-yz', step, 1)
    pzz = get_trans_mean('p-zz', step, 1)
    tpar = get_trans_mean('tpar', step, 1)
    tperp = get_trans_mean('tperp', step, 1)

    fig, axes = plt.subplots(4,2, figsize=[8,7], sharex=True)
    fig.suptitle('t={:.1f}'.format(step*dt), fontsize=16)
    ax1 = axes[0,0]; ax1.plot(pxx, 'C0'); ax1.set_ylabel(r'$p-xx$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1 = axes[1,0]; ax1.plot(pxy, 'C0'); ax1.set_ylabel(r'$p-xy$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1 = axes[2,0]; ax1.plot(pxz, 'C0'); ax1.set_ylabel(r'$p-xz$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1 = axes[0,1]; ax1.plot(pyy, 'C0'); ax1.set_ylabel(r'$p-yy$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1 = axes[1,1]; ax1.plot(pyz, 'C0'); ax1.set_ylabel(r'$p-yz$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1 = axes[2,1]; ax1.plot(pzz, 'C0'); ax1.set_ylabel(r'$p-zz$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1 = axes[3,0]; ax1.plot(tpar, 'C0'); ax1.set_ylabel(r'$T_\parallel^1$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1.set_xlabel('z'); ax1.set_xlim(0, zmax)
    ax1 = axes[3,1]; ax1.plot(tperp, 'C0'); ax1.set_ylabel(r'$T_\perp^1$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1.set_xlabel('z'); ax1.set_xlim(0, zmax)
    plt.tight_layout()
    plt.savefig(imag_path+'/pt_{:08d}'.format(step))
    plt.close('all')

def pt_movie():
    processor = functools.partial(pt_single)
    dump_num = np.arange(201)
    with mp.Pool(processes=20) as worker_pool:
        frames = worker_pool.map(processor, dump_num)
        worker_pool.close()
        worker_pool.join()
    makeMovie('pt', 'pt.mp4', fps=3)


def fv_single(i):
    imag_path = 'fv'
    mkdir(imag_path)
    step = 1 if i == 0 else i * tinterval
    print (step)
    fox = get_trans_mean('fox', step)
    foy = get_trans_mean('foy', step)
    foz = get_trans_mean('foz', step)
    vix = get_trans_mean('vix', step)
    viy = get_trans_mean('viy', step)
    viz = get_trans_mean('viz', step)
    vxs = get_trans_mean('vxs', step, 1)
    vys = get_trans_mean('vys', step, 1)
    vzs = get_trans_mean('vzs', step, 1)


    fig, axes = plt.subplots(3,3, figsize=[12,5.5], sharex=True)
    fig.suptitle('t={:.1f}'.format(step*dt), fontsize=16)
    ax1 = axes[0,0]; ax1.plot(fox, 'C0'); ax1.set_ylabel(r'$fox$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1 = axes[1,0]; ax1.plot(foy, 'C0'); ax1.set_ylabel(r'$foy$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1 = axes[2,0]; ax1.plot(foz, 'C0'); ax1.set_ylabel(r'$foz$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1.set_xlabel('z'); ax1.set_xlim(0, zmax)
    ax1 = axes[0,1]; ax1.plot(vix, 'C0'); ax1.set_ylabel(r'$vix$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1 = axes[1,1]; ax1.plot(viy, 'C0'); ax1.set_ylabel(r'$viy$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1 = axes[2,1]; ax1.plot(viz, 'C0'); ax1.set_ylabel(r'$viz$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1.set_xlabel('z'); ax1.set_xlim(0, zmax)
    ax1 = axes[0,2]; ax1.plot(vxs, 'C0'); ax1.set_ylabel(r'$vxs$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1 = axes[1,2]; ax1.plot(vys, 'C0'); ax1.set_ylabel(r'$vys$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1 = axes[2,2]; ax1.plot(vzs, 'C0'); ax1.set_ylabel(r'$vzs$'); ax1.yaxis.set_major_formatter(sci_format())
    ax1.set_xlabel('z'); ax1.set_xlim(0, zmax)

    plt.tight_layout()
    plt.savefig(imag_path+'/fv_{:08d}'.format(step))
    plt.close('all')
        
def fv_movie():
    processor = functools.partial(fv_single)
    dump_num = np.arange(201)
    with mp.Pool(processes=20) as worker_pool:
        frames = worker_pool.map(processor, dump_num)
        worker_pool.close()
        worker_pool.join()
    makeMovie('fv', 'fv.mp4', fps=3)

if __name__ == "__main__":
    t_start = systime.time()
    plot_den_single(int(sys.argv[1]))
    # ebden_eta_single(200)
    # ebden_eta_movie()
    # pt_single(200)
    # pt_movie()
    # fv_single(200)
    # fv_movie()
    print ('time used: {} s'.format(systime.time()-t_start))

from h3d import *
home='/home1/07918/fyli'

def load_input(path):
  global node_conf, nprocs, denmin, nspec, nx, ny, nz, \
      xmax, ymax, zmax, tmax, dt, dB_B0, n_diag_mesh, y, z, \
      ndumps, timesteps, times, mask_zs, mask_r, \
      n_wave_cycles, lambda0, wpiwci

  print ('loading input file of \"%s\"'%path)
  nml = f90nml.read(path + '/input.f90')
  
  node_conf = nml['input']['node_conf']
  denmin = nml['input']['denmin']
  nspec = nml['input']['nspec']
  nx = nml['input']['nx']
  ny = nml['input']['ny']
  nz = nml['input']['nz']
  xmax = nml['input']['xmax']
  ymax = nml['input']['ymax']
  zmax = nml['input']['zmax']
  tmax = nml['input']['tmax']
  dt = nml['input']['dtwci']
  mask_zs = nml['input']['mask_zs']
  mask_r = nml['input']['mask_r']
  dB_B0 = nml['input']['dB_B0']
  n_wave_cycles = nml['input']['n_wave_cycles']
#   n_wave_cycles = nml['input']['num_wave_cycles']
  wpiwci = nml['input']['wpiwci']
  n_diag_mesh = nml['input']['n_diag_mesh']

  # derived
  nprocs = int(node_conf[0]*node_conf[1])
  y = np.linspace(0, ymax, ny)
  z = np.linspace(0, zmax, nz)
  lambda0 = zmax/n_wave_cycles # wavelength
  ndumps = int(tmax/(dt*n_diag_mesh)+1)
  timesteps = np.arange(ndumps)*n_diag_mesh
  times = timesteps * dt
  
  return

def get_snapshot_data(path, field, ts):
  subfolder = field
  if '_' in field:
    subfolder = field[:field.index('_')]
  fname = '%s/%s_%d.gda' % (subfolder, field, ts)
  fpath = join(path, 'data', fname)
  if not exists(fpath): 
    data = np.zeros(int(ny*nz))
  else:
    data = np.fromfile(fpath, dtype=np.float32)
    # data = data.reshape(nz, ny).transpose()
  return (data)

def get_snapshot_data_mean(path, field, ts, axis=0):
    d = get_snapshot_data(path, field, ts)
    d = d.reshape(nz,ny).transpose()
    d_mean = np.mean(d, axis=axis)  # averaged over y
    return (d_mean)

#            0      1      2    3     4     5     6       7        8
field_list=['den', 'bx', 'by', 'bz', 'ex', 'ey', 'ez', 'tpar_1', 'tperp_1']

def read_sim_data(path):
#     load_input(path)
    dsize = ny*nz
    data = np.zeros( (ndumps, dsize, len(field_list)) )
    datafile = join(path, 'hist.npy')
    if not exists(datafile):
        for i, field in enumerate(field_list):
            for j, ts in enumerate(timesteps):
                showProgressBar(ndumps-1, j)
                data[j,:,i] = get_snapshot_data(path, field, ts)
        np.save(datafile, data)
    else:
        data = np.load(datafile)
    print ('sim data extracted')
    return (data)

def get_spectral_max(path, field, spec, spec_max, k_pos):
    data = read_sim_data(path)
    for m in range(ndumps):
        idx = field_list.index(field)
        d = data[m,:,idx] 
        d = d.reshape(nz,ny).transpose()
        d_mean = np.mean(d, axis=0)  # averaged over y
        if field=='den': d_mean -= 1.
        f, F = FFT(z, d_mean)
        spec[m,:] = F
        spec_max[m] = F.max()
        k_pos[m]=f[np.argmax(F)]
    return (spec, spec_max, k_pos)


def get_field_mean(data, fieldname, idump):
    idx = field_list.index(field)
    d_ = data[idump,:,idx]
    d_ = d_.reshape(nz,ny).transpose()
    d = np.mean(d_, axis=0)
    return (d)

def get_spectral_max_spatial_field(path, spec, spec_max, k_pos):
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
  


def spt(path0):
  # path0 = 'mask/case_z30'
  path = join(home, 'pdi', path0)
  load_input(path)

  fld = 'bx'
  # ndumps =800
  spt = np.zeros((ndumps, nz))
  for i in range(ndumps):
    showProgressBar(ndumps-1, i)
    spt[i,:] = get_snapshot_data_mean(path, fld, timesteps[i], axis=0)

  fig, axes = plt.subplots(1,1, figsize=[5,4], dpi=150)
  ax1 = axes
  img1 = ax1.imshow(spt, aspect='auto', origin='lower', cmap='bwr',
                  extent=[0, nz/lambda0, 0, ndumps*n_diag_mesh*dt/lambda0])
  img1.set_clim(-dB_B0, dB_B0)
  cbar = plt.colorbar(img1, ax=ax1, pad=0.02)
  cbar.set_label(fld)
  ax1.set_xlabel(r'$z/\lambda_0$')
  ax1.set_ylabel(r'$t/T_0$')

  plt.tight_layout()
  plt.show()



def snap(path0, t):
  path = join(home, 'pdi', path0)

  load_input(path)

  fig, axes = plt.subplots(3,1, figsize=[7, 7.5], dpi=120)

  # t = 400
  # t = 330
  fld_list = ['den', 'bx', 'ex']

  idx = np.argmin(np.abs(times-t))
  ts = timesteps[idx]
  print('t=', ts*dt)

  for i, field in enumerate(fld_list):
    ax1=axes[i]
  #   ax1.set_title(fld_list[i])
    d = get_snapshot_data_mean(path, fld_list[i], ts, axis=0)
    line,=ax1.plot(z/lambda0, d, 'C0', linewidth=0.5)
    ax1.set_xlim(0, zmax/lambda0)
  #   ax1.set_xlim(0, 2)
  #   ax1.set_ylabel(fld_list[i])

    if i==0: 
      ax1.set_ylabel(fld_list[i])
      ax1.set_title(r'$t=%.2fT_0$'%(ts*dt/lambda0))
      ax1.set_ylim(0.8,1.2)
      d = get_snapshot_data_mean(path, 'ez', ts, axis=0)
      print('mean of ez: ', np.mean(d))
      ax2 = ax1.twinx()
      line,=ax2.plot(z/lambda0, d,'C1', linewidth=0.5)
      ax2.yaxis.label.set_color(line.get_color())
      ax2.tick_params(axis='y',colors=line.get_color())
      ax2.set_ylabel('ez')
  #     ax2.set_ylim(-2e-3,2e-3)
      ax2.set_ylim(-2e-3,2e-3)

    if i==1:
      ax1.set_ylabel('bx; by')
      d = get_snapshot_data_mean(path, 'by', ts, axis=0)
      line,=ax1.plot(z/lambda0, d, 'C1', linewidth=0.5)
  #     ax1.plot([0, zmax/lambda0], [0.1,0.1], 'k--', linewidth=0.5)
      # ax1.set_xlabel(r'$z/\lambda_0$')
      ax1.set_ylim(-2*dB_B0,2*dB_B0)
      # plot fm
      d = np.zeros(nz)
      z_cell = np.arange(nz)+1
      d[z_cell<=zmax-mask_zs] = 1.0
      d[z_cell>zmax-mask_zs] = 1-(mask_r*(z_cell[z_cell>nz-mask_zs]-nz+mask_zs)/mask_zs)**2.
      d[z_cell<=mask_zs] = 1-(mask_r*(z_cell[z_cell<=mask_zs]-mask_zs)/mask_zs)**2.
      ax2 = ax1.twinx()
      line,=ax2.plot(z/lambda0, d,'C2', linewidth=0.5)
      ax2.yaxis.label.set_color(line.get_color())
      ax2.tick_params(axis='y',colors=line.get_color())
      ax2.set_ylabel(r'$f_M$')
      ax2.set_ylim(0,4/3.)

    if i==2:
      ax1.set_ylabel('ex; ey')
      d = get_snapshot_data_mean(path, 'ey', ts, axis=0)
      line,=ax1.plot(z/lambda0, d, 'C1', linewidth=0.5)
  #     ax1.plot([0, zmax/lambda0], [0.1,0.1], 'k--', linewidth=0.5)
      ax1.set_xlabel(r'$z/\lambda_0$')
      ax1.set_ylim(-2*dB_B0,2*dB_B0)
      # plot fm
      d = np.zeros(nz)
      z_cell = np.arange(nz)+1
      d[z_cell<=zmax-mask_zs] = 1.0
      d[z_cell>zmax-mask_zs] = 1-(mask_r*(z_cell[z_cell>nz-mask_zs]-nz+mask_zs)/mask_zs)**2.
      d[z_cell<=mask_zs] = 1-(mask_r*(z_cell[z_cell<=mask_zs]-mask_zs)/mask_zs)**2.
      ax2 = ax1.twinx()
      line,=ax2.plot(z/lambda0, d,'C2', linewidth=0.5)
      ax2.yaxis.label.set_color(line.get_color())
      ax2.tick_params(axis='y',colors=line.get_color())
      ax2.set_ylabel(r'$f_M$')
      ax2.set_ylim(0,4/3.)
      
  plt.tight_layout(pad=2.)
  plt.show()



# energy evolution (comparison)
def ene(path0):
  fld_list = ['it', 'efld', 'bfld', 'efluid', 'ethermal', 'eptcl']

  def get_bfld_trans():
    path = join(home, 'pdi', path0)
    load_input(path)
    Bz_ene = (1/wpiwci)**2.*nx*ny*nz*0.5
    print('Bz_ene: ', Bz_ene)
    print('nx, ny, nz = ', nx, ny, nz)

    energy = np.loadtxt(join(path, 'data', 'energy.dat'))
    print ('energy.shape = ', energy.shape)

    ene_dict = {}
    for i in range(6):
      ene_dict[fld_list[i]] = energy[:,i]
    taxis = ene_dict['it']*dt
    ene_dict['bfld'] = ene_dict['bfld'] - Bz_ene
    # taxis = energy[:,0]*dt
    # bfld_trans = energy[:,2]-Bz_ene
    # efld = energy[:,1]
    # efluid = energy[:,3]
    # ethermal = energy[:,4]
    # eptcl = energy[:,5]
    return (taxis, ene_dict)

  taxis, ene_dict = get_bfld_trans()

  fig, axes = plt.subplots(5,1,figsize=[8,8], sharex=True)
  for i in range(1,6):
    ax1 = axes[i-1]
    # ax1.set_title('mask5')
    # taxis, bfld_trans = get_bfld_trans('mask/case_waveabsorb')
    ax1.plot(taxis/lambda0, ene_dict[fld_list[i]])#label=fld_list[2]+' (trans)')
    # ax1.plot(taxis/lambda0, efld, )
    # ax1.legend()
    ax1.set_xlim(0,)
    if i<=3:
      ax1.set_ylim(0,)

    if i==5: 
      ax1.set_xlabel(r'$t/T_0$')

    ax1.set_ylabel(fld_list[i])
    if i==2:
      ax1.set_ylabel('bfld (perp)')
    
    ax1.yaxis.set_major_formatter(sci_format())

  plt.tight_layout()
  plt.show()



if __name__ == "__main__":
  import sys
  if sys.argv[2] == 'spt':
    spt(sys.argv[1])
  elif sys.argv[2] == 'snap':
    snap(sys.argv[1], t=int(sys.argv[3]))
  elif sys.argv[2] == 'ene':
    ene(sys.argv[1])
  else:
    raise Exception("arguments parsed incorrectly!")
  
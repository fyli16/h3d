"""
Analysis procedures for h3d
"""

from common_modules import *
tweak_rcParams(fs=12, home='/home1/07918/fyli')
mpl.use('TKAgg')
# mpl.use('qt5agg')

# import collections
# import math
# import os.path
# import struct

# import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib import rc
# from matplotlib.colors import LogNorm
# from matplotlib.ticker import MaxNLocator
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# from mpl_toolkits.mplot3d import Axes3D


def calc_electric_energy():
    electric_energy = np.zeros(nt)
    for ct in range(nt):
        showProgressBar(nt-1, ct)
        tframe = 1 if ct == 0 else ct * tinterval
        fname = fpath + 'ex_' + str(tframe) + '.gda'
        ex = np.fromfile(fname, dtype=np.float32)
        fname = fpath + 'ey_' + str(tframe) + '.gda'
        ey = np.fromfile(fname, dtype=np.float32)
        fname = fpath + 'ez_' + str(tframe) + '.gda'
        ez = np.fromfile(fname, dtype=np.float32)
        sz, = ex.shape
        electric_energy[ct] = np.sum(ex**2 + ey**2 + ez**2) * 0.5/sz
    oname = cpath + '/electric_energy.dat'
    electric_energy.tofile(oname)


def calc_magnetic_energy():
    magnetic_energy = np.zeros(nt)
    for ct in range(nt):
        showProgressBar(nt-1, ct)
        tframe = 1 if ct == 0 else ct * tinterval
        fname = fpath + 'bx_' + str(tframe) + '.gda'
        bx = np.fromfile(fname, dtype=np.float32)
        fname = fpath + 'by_' + str(tframe) + '.gda'
        by = np.fromfile(fname, dtype=np.float32)
        fname = fpath + 'bz_' + str(tframe) + '.gda'
        bz = np.fromfile(fname, dtype=np.float32)
        sz, = bx.shape
        magnetic_energy[ct] = np.sum(bx**2 + by**2 + (bz-1.)**2) * 0.5/sz
    oname = cpath + '/magnetic_energy.dat'
    magnetic_energy.tofile(oname)


def calc_temperature():
    tpara_avg = np.zeros(nt)
    tperp_avg = np.zeros(nt)
    for s in range(1,nspec+1):
      print('species', s)
      for ct in range(nt):
          showProgressBar(nt-1, ct)
          tframe = 1 if ct == 0 else ct * tinterval
          fname = fpath + 'tpar_' +str(s)+'_'+ str(tframe) + '.gda'
          tpara = np.fromfile(fname, dtype=np.float32)
          fname = fpath + 'tperp_' +str(s)+'_'+ str(tframe) + '.gda'
          tperp = np.fromfile(fname, dtype=np.float32)
          tpara_avg[ct] = np.mean(tpara)
          tperp_avg[ct] = np.mean(tperp)
      oname = '{}/tpar_{}.dat'.format(cpath, s)
      tpara_avg.tofile(oname)
      oname = '{}/tperp_{}.dat'.format(cpath, s)
      tperp_avg.tofile(oname)


def calc_density_fluctuation():
    den_fluc = np.zeros(nt)
    for ct in range(nt):
        showProgressBar(nt-1, ct)
        tframe = 1 if ct == 0 else ct * tinterval
        fname = fpath + 'den_' + str(tframe) + '.gda'
        rho = np.fromfile(fname, dtype=np.float32)
        sz, = rho.shape
        rho_avg = np.sum(rho) / sz
        den_fluc[ct] = np.sum((rho - rho_avg)**2) / sz
    oname = cpath + '/den_fluc.dat'
    den_fluc.tofile(oname)

def plot_history(ifshow):
    fig, axes = plt.subplots(1,3, figsize=[12, 3.5])
    ax1 = axes[0]
    b_ene = np.fromfile('{}/magnetic_energy.dat'.format(cpath) )
    e_ene = np.fromfile('{}/electric_energy.dat'.format(cpath) )
    ax1.plot(time, b_ene, label=r'$\delta B^2$')
    ax1.plot(time, e_ene, label=r'$\delta E^2$')
    ax1.legend()
    ax1.grid()
    ax1.set_xlabel(r'$t\Omega_{ci}$')

    ax1 = axes[1]
    for s in range(1, nspec+1):
        tpara = np.fromfile('{}/tpar_{}.dat'.format(cpath, s) )
        tperp = np.fromfile('{}/tperp_{}.dat'.format(cpath, s) )
        line, = ax1.plot(time, tperp, label=r'$T_{%s\perp}$'%s)
        ax1.plot(time, tpara, ls='--', color=line.get_color(), label=r'$T_{%s\parallel}$'%s)
    ax1.grid()
    ax1.set_xlabel(r'$t\Omega_{ci}$')
    ax1.legend()

    ax1 = axes[2]
    den_fluc = np.fromfile('{}/den_fluc.dat'.format(cpath) )
    ax1.plot(time, den_fluc, label=r'$\delta\rho^2$')
    ax1.legend()
    ax1.grid()
    ax1.set_xlabel(r'$t\Omega_{ci}$')

    plt.tight_layout()
    plt.savefig(cpath+'/history.png')
    if ifshow: plt.show()
    plt.close('all')

def load_input():
    import f90nml
    global denmin, nspec, nx, ny, nz, xmax, ymax, zmax, tmax, dt, tinterval
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
    return 


if __name__ == "__main__":
    fpath = '../data/'

    ## input parameters
    # exec(open('./para.out').read())
    load_input()
    # import f90nml
    # nml = f90nml.read(fpath + 'finput.dat')
    # denmin = nml['datum']['denmin']
    # nspec = nml['datum']['nspec']
    # nx = nml['datum']['nx']
    # ny = nml['datum']['ny']
    # nz = nml['datum']['nz']
    # xmax = nml['datum']['xmax']
    # ymax = nml['datum']['ymax']
    # zmax = nml['datum']['zmax']
    # tmax = nml['datum']['maximum_simulation_time']
    # dt = nml['datum']['dtwci']
    # tinterval = nml['datum']['nwrtdata']
    nt = int(tmax/(dt*tinterval)+1)
    print ('total number of dump steps: ', nt)
    time = np.arange(nt)*dt*tinterval

    # print (nx, ny, nz, time)

    ## extract and save history data
    cpath = 'cdata'  # compiled data path
    # mkdir(cpath)
    # print ('\ncalculating electric energy ...')
    # calc_electric_energy()
    # print ('\ncalculating magnetic energy ...')
    # calc_magnetic_energy()
    # print ('\ncalculating temperature ...')
    # calc_temperature()
    # print ('\ncalculating density fluctuation ...')
    # calc_density_fluctuation()

    ## plot history
    print ('\nplotting history ...')
    plot_history(ifshow=1)

    print ('\ndone!')

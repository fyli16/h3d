"""
Analysis procedures for h3d
"""
# import collections
# import math
# import os.path
# import struct

from common_modules import *
tweak_rcParams(fs=12, home='/home1/07918/fyli')
mpl.use('TKAgg')
# mpl.use('qt5agg')
# import matplotlib.pyplot as plt
# import numpy as np

# from matplotlib import rc
# from matplotlib.colors import LogNorm
# from matplotlib.ticker import MaxNLocator
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# from mpl_toolkits.mplot3d import Axes3D

# rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
# mpl.rc('text', usetex=True)
# mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]

# font = {
#     'family': 'serif',
#     #'color'  : 'darkred',
#     'color': 'black',
#     'weight': 'normal',
#     'size': 24,
# }


def calc_electric_energy(tinterval, nt, fpath, run_name):
    """Calculate electric energy

    Args:
        tinterval: time step interval
        nt: number of time frames
        fpath: data file path
    """
    electric_energy = np.zeros(nt)
    for ct in range(nt):
        showProgressBar(nt-1, ct)
        tframe = 1 if ct == 0 else ct * tinterval
        # print(tframe)
        fname = fpath + 'ex_' + str(tframe) + '.gda'
        ex = np.fromfile(fname, dtype=np.float32)
        fname = fpath + 'ey_' + str(tframe) + '.gda'
        ey = np.fromfile(fname, dtype=np.float32)
        fname = fpath + 'ez_' + str(tframe) + '.gda'
        ez = np.fromfile(fname, dtype=np.float32)
        sz, = ex.shape
        electric_energy[ct] = np.sum(ex**2 + ey**2 + ez**2) * 0.5/sz
    # plt.plot(electric_energy)
    # if ifshow: plt.show()
    oname = '{}/electric_energy_{}.dat'.format(cpath, run_name)
    electric_energy.tofile(oname)


def calc_magnetic_energy(tinterval, nt, fpath, run_name):
    """Calculate magnetic energy

    Args:
        tinterval: time step interval
        nt: number of time frames
        fpath: data file path
    """
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
    # plt.plot(magnetic_energy)
    # if ifshow: plt.show()
    oname = '{}/magnetic_energy_{}.dat'.format(cpath, run_name)
    magnetic_energy.tofile(oname)


def calc_temperature(tinterval, nt, nspec, fpath, run_name):
    """Calculate temperature

    Args:
        tinterval: time step interval
        nt: number of time frames
        fpath: data file path
        spec: number of species 
    """
    tpara_avg = np.zeros(nt)
    tperp_avg = np.zeros(nt)
    for s in range(1,nspec+1):
      print('species',s)
      for ct in range(nt):
          showProgressBar(nt-1, ct)
          tframe = 1 if ct == 0 else ct * tinterval
        #   print(tframe)
          #fname = fpath + 'den_' + str(tframe) + '.gda'
          #rho = np.fromfile(fname, dtype=np.float32)
          fname = fpath + 'tpar_' +str(s)+'_'+ str(tframe) + '.gda'
          tpara = np.fromfile(fname, dtype=np.float32)
          fname = fpath + 'tperp_' +str(s)+'_'+ str(tframe) + '.gda'
          tperp = np.fromfile(fname, dtype=np.float32)
        #   print(tperp[-1])
          tpara_avg[ct] = np.mean(tpara)
          tperp_avg[ct] = np.mean(tperp)
    #   plt.plot(tpara_avg,label='T %d para'%s)
    #   plt.plot(tperp_avg,label='T %d perp'%s)
      oname = '{}/tpar_{}_{}.dat'.format(cpath, s, run_name)
      tpara_avg.tofile(oname)
      oname = '{}/tperp_{}_{}.dat'.format(cpath, s, run_name)
      tperp_avg.tofile(oname)
    # plt.legend(loc='best')
    # if ifshow: plt.show()


def calc_density_fluctuation(tinterval, nt, fpath, run_name):
    """Calculate density fluctuation

    Args:
        tinterval: time step interval
        nt: number of time frames
        fpath: data file path
    """
    den_fluc = np.zeros(nt)
    for ct in range(nt):
        showProgressBar(nt-1, ct)
        tframe = 1 if ct == 0 else ct * tinterval
        # print(tframe)
        fname = fpath + 'den_' + str(tframe) + '.gda'
        rho = np.fromfile(fname, dtype=np.float32)
        sz, = rho.shape
        rho_avg = np.sum(rho) / sz
        den_fluc[ct] = np.sum((rho - rho_avg)**2) / sz
    # plt.plot(den_fluc)
    # if ifshow: plt.show()
    oname = '{}/den_fluc_{}.dat'.format(cpath, run_name)
    den_fluc.tofile(oname)


if __name__ == "__main__":
    # run_name = 'pdi-001'
    # tinterval = 1000
    # nt = 20
    # nspec = 1
    # ifshow=True
    exec(open('./para.out').read())
    fpath = '../../../' + run_name + '/data/'
    cpath = 'cdata'  # compiled data path
    mkdir(cpath)

    print ('\ncalculating electric energy ...')
    calc_electric_energy(tinterval, nt, fpath, run_name)
    print ('\ncalculating magnetic energy ...')
    calc_magnetic_energy(tinterval, nt, fpath, run_name)
    print ('\ncalculating temperature ...')
    calc_temperature(tinterval, nt, nspec, fpath, run_name)
    print ('\ncalculating density fluctuation ...')
    calc_density_fluctuation(tinterval, nt, fpath, run_name)
    print ('\ndone!')

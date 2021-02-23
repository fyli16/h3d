"""
Analysis procedures to calculate power spectrum
"""
import collections
import math
import multiprocessing
import os.path
import struct

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed
from matplotlib import rc
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D

import palettable
import pic_information
from contour_plots import plot_2d_contour, read_2d_fields
from energy_conversion import read_data_from_json
from shell_functions import mkdir_p

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]

font = {
    'family': 'serif',
    #'color'  : 'darkred',
    'color': 'black',
    'weight': 'normal',
    'size': 24,
}

colors = palettable.colorbrewer.qualitative.Set1_9.mpl_colors

# colors = palettable.colorbrewer.qualitative.Dark2_8.mpl_colors


def calc_power_spectrum_mag(pic_info,
                            ct,
                            run_name,
                            shock_pos,
                            base_dir='../../'):
    """calculate power spectrum using magnetic fields

    Args:
        pic_info: namedtuple for the PIC simulation information.
        ct: current time frame.
    """
    xmin, xmax = 0, pic_info.lx_di
    xmin, xmax = 0, 105
    zmin, zmax = -0.5 * pic_info.lz_di, 0.5 * pic_info.lz_di
    kwargs = {
        "current_time": ct,
        "xl": xmin,
        "xr": xmax,
        "zb": zmin,
        "zt": zmax
    }
    fname = base_dir + 'data1/vex.gda'
    x, z, vel = read_2d_fields(pic_info, fname, **kwargs)
    nx, = x.shape
    nz, = z.shape
    # data_cum = np.sum(vel, axis=0) / nz
    # data_grad = np.abs(np.gradient(data_cum))
    # xs = 5
    # max_index = np.argmax(data_grad[xs:])
    # xm = x[max_index]
    xm = x[shock_pos]

    xmin, xmax = 0, xm
    fname = base_dir + 'data1/bx.gda'
    x, z, bx = read_2d_fields(pic_info, fname, **kwargs)
    fname = base_dir + 'data1/by.gda'
    x, z, by = read_2d_fields(pic_info, fname, **kwargs)
    fname = base_dir + 'data1/bz.gda'
    x, z, bz = read_2d_fields(pic_info, fname, **kwargs)
    smime = math.sqrt(pic_info.mime)
    lx = np.max(x) - np.min(x)
    lz = np.max(z) - np.min(z)

    bx_k = np.fft.rfft2(bx)
    by_k = np.fft.rfft2(by)
    bz_k = np.fft.rfft2(bz)
    b2_k = np.absolute(bx_k)**2 + np.absolute(by_k)**2 + np.absolute(bz_k)**2
    xstep = lx / nx
    kx = np.fft.fftfreq(nx, xstep)
    idx = np.argsort(kx)
    zstep = lz / nz
    kz = np.fft.fftfreq(nz, zstep)
    idz = np.argsort(kz)
    print np.min(kx), np.max(kx), np.min(kz), np.max(kz)

    kxs, kzs = np.meshgrid(kx[:nx // 2 + 1], kz)
    ks = np.sqrt(kxs * kxs + kzs * kzs)
    # kmin, kmax = np.min(ks), np.max(ks)
    # kbins = np.linspace(kmin, kmax, nx//2+1, endpoint=True)
    kmin = 1E-2
    kmax = np.max(ks)
    kmin_log, kmax_log = math.log10(kmin), math.log10(kmax)
    kbins = 10**np.linspace(kmin_log, kmax_log, 64, endpoint=True)
    ps, kbins_edges = np.histogram(
        ks, bins=kbins, weights=b2_k * ks, density=True)
    w1, h1 = 0.8, 0.8
    xs, ys = 0.15, 0.95 - h1
    fig = plt.figure(figsize=[7, 5])
    ax1 = fig.add_axes([xs, ys, w1, h1])
    for index, k in np.ndenumerate(kbins):
        pass
        # print index, k
    psm = 25
    # pindex = -5.0/3.0
    pindex = -2.0
    power_k = kbins[psm:]**pindex
    shift = 22
    ax1.loglog(kbins_edges[:-1], ps, linewidth=2)
    ax1.loglog(
        kbins[psm:psm + shift],
        power_k[:shift] * 2 / power_k[0],
        linestyle='--',
        linewidth=2,
        color='k')
    # power_index = "{%0.2f}" % pindex
    power_index = '-2.0'
    tname = r'$\sim k^{' + power_index + '}$'
    ax1.text(
        0.45,
        0.7,
        tname,
        color='black',
        fontsize=24,
        horizontalalignment='left',
        verticalalignment='center',
        transform=ax1.transAxes)
    ax1.tick_params(labelsize=16)
    ax1.set_xlabel(r'$kd_i$', fontdict=font, fontsize=20)
    ax1.set_ylabel(r'$E_B(k)$', fontdict=font, fontsize=20)
    ax1.set_xlim([1E-2, 3E1])
    ax1.set_ylim([1E-3, 3E1])

    fig_dir = '../img/img_power_spectrum/' + run_name + '/'
    mkdir_p(fig_dir)
    fname = fig_dir + '/ps_mag_' + str(ct).zfill(3) + '.jpg'
    fig.savefig(fname, dpi=300)

    # plt.show()
    plt.close()


if __name__ == "__main__":
    # base_dir = '/net/scratch3/xiaocanli/2D-90-Mach4-sheet4-multi/'
    # run_name = '2D-90-Mach4-sheet4-multi'
    base_dir = '/net/scratch2/guofan/for_Senbei/2D-90-Mach4-sheet6-2/'
    run_name = '2D-90-Mach4-sheet6-2'
    picinfo_fname = '../data/pic_info/pic_info_' + run_name + '.json'
    pic_info = read_data_from_json(picinfo_fname)
    ct = pic_info.ntf - 2
    # ct = 200
    cts = range(10, pic_info.ntf - 1)

    xmin, xmax = 0, pic_info.lx_di
    xmin, xmax = 0, 105
    zmin, zmax = -0.5 * pic_info.lz_di, 0.5 * pic_info.lz_di
    # kwargs = {"current_time":ct, "xl":xmin, "xr":xmax, "zb":zmin, "zt":zmax}
    # fname = base_dir + 'data1/vex.gda'
    # x, z, vel = read_2d_fields(pic_info, fname, **kwargs) 
    kwargs = {
        "current_time": 0,
        "xl": xmin,
        "xr": xmax,
        "zb": zmin,
        "zt": zmax
    }
    fields_interval = pic_info.fields_interval
    tframe = str(fields_interval * ct)
    fname = base_dir + 'data/vex_' + tframe + '.gda'
    x, z, pxx = read_2d_fields(pic_info, fname, **kwargs)
    fname = '../data/shock_pos/shock_pos_' + run_name + '.txt'
    shock_loc = np.genfromtxt(fname, dtype=np.int32)
    sloc = shock_loc[ct]
    xm = x[sloc]

    def processInput(ct):
        print ct
        sloc = shock_loc[ct]
        xm = x[sloc]
        # calc_avg_bfield(pic_info, ct, run_name, xm, base_dir)
        # calc_power_spectrum_mag(pic_info, ct, run_name, sloc, base_dir)
        # calc_power_spectrum_vel(pic_info, ct, 'e', run_name, sloc, base_dir)
        # calc_power_spectrum_vel_comp(pic_info, ct, 'i', run_name, sloc,
        #                              base_dir, single_file=False)
        plot_power_spectrum_vel_comp_du(
            pic_info, ct, 'i', run_name, xm, base_dir, single_file=False)

    num_cores = multiprocessing.cpu_count()
    # num_cores = 8
    Parallel(n_jobs=num_cores)(delayed(processInput)(ct) for ct in cts)
    # calc_power_spectrum_mag(pic_info, ct, run_name, sloc, base_dir)
    # calc_power_spectrum_vel(pic_info, ct, 'i', run_name, sloc, base_dir)
    # calc_avg_bfield(pic_info, ct, run_name, xm, base_dir)
    # plot_avg_bfiled(pic_info)
    # calc_power_spectrum_vel_comp(pic_info, ct, 'i', run_name, sloc, base_dir,
    #                              single_file=False)
    # plot_power_spectrum_vel_comp_du(pic_info, ct, 'i', run_name, xm,
    #                                 base_dir, single_file=False)
    # for ct in cts:
    #     print ct
    #     plot_power_spectrum_vel_comp_du(pic_info, ct, 'i', run_name, xm, base_dir)

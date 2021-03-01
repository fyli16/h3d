#!/usr/bin/env python

from common_modules import *
tweak_rcParams(fs=14, home='/home1/07918/fyli')
mpl.use('TKAgg')

import argparse
import functools

from os.path import join
import f90nml

def argument_parse():
    parser = argparse.ArgumentParser(description='Analysis of H3D data')
    parser.add_argument('-p', dest='path', type=str, default='.', 
                        metavar="path", help="path to simulation data")
    parser.add_argument('-st', dest='showtimes', action='store_true', 
                        help="show available timesteps and exit")
    parser.add_argument('-sp', dest='showplots', action='store_true', 
                        help="show plots before saving")
    parser.add_argument('-t', dest='time', type=float, metavar="time", 
                        help="single time to process. Will use the nearest simulation dumps")
    parser.add_argument('-ts', dest='timestep', type=int, metavar="timestep", 
                        help="single timestep to process. Will use the nearest simulation dumps")
    parser.add_argument('-tl', dest='timelist',type=float, nargs=3,         
                        metavar="timelist", 
                        help="time list to process: [t0, tend, dt]")
    parser.add_argument('-m', dest='movie', action='store_true', 
                        help="movie plots using all available timesteps")
    parser.add_argument('-np', dest='nprocs', type=int, default=10, 
                        help="number of processes to use, default 10")
    parser.add_argument('-fs',  dest='figsize', type=float, nargs=2, default=[5,3.8],
                        help='figure size')
    parser.add_argument('-xr', dest='xr', type=float, nargs=2, 
                        help='x-axis range')
    parser.add_argument('-yr', dest='yr', type=float, nargs=2, 
                        help='y-axis range')
    parser.add_argument('-log', dest='logscale', default='none', 
                        help="set logscale for axes: 'x' for x axis, 'y' for y axis, and 'xy' for both axes")
    parser.add_argument('-dir', dest='dir', type=str,metavar="dir", 
                        help="directory to save plots")
    parser.add_argument('-fps', dest='fps', type=int, default=2, 
                        help="frames per second for making movies")
    parser.add_argument('-test', action='store_true', 
                        help='test arguments and return')
    
    #@ subparser handlers
    subparsers = parser.add_subparsers(dest='mode', 
                            help='Modes for different diagnostics')
    
    #@ field
    field = subparsers.add_parser('field', help="plot a given field quantity (e.g., e, b, v)")
    field.set_defaults(plot_func=field_plot)
    field.add_argument('-fn', dest='fieldname', default='bx', 
                        help="field name")
    field.add_argument("-type", dest="type", #action="store_true",
                        help="type of plots to make: snap, spec, or hist; they correspond to snapshot of a given time frame, spectral of a given frame, or history evolution of all time frames, respectively") 
    # field.add_argument("-spec", dest="spectrum", action="store_true",
    #                     help="obtain spectrum of data") 
    # field.add_argument("-env", dest="envelope", action="store_true",
    #                     help="extract envelope of data") 
    # field.add_argument('-hist', dest='history', action="store_true", 
    #                     help="plot history evolution")

    return (parser.parse_args())


def get_field(ts):
    # fname = a.path+'/data/'+a.fieldname+'_'+str(ts)+'.gda'
    fname = '%s_%d.gda' % (a.fieldname, ts)
    fname = join(a.path, 'data', fname)
    data = np.fromfile(fname, dtype=np.float32)
    data = data.reshape(nz, ny).transpose()
    return (data)

def snapshot(ts, data):
    fig, axes = plt.subplots(2,1, figsize=a.figsize, sharex=True)
    ax1 = axes[0]
    ax1.set_title('%s, t = %.2f'%(a.fieldname, ts*dt))
    im1 = ax1.imshow(data, aspect='auto', origin='lower', cmap='jet', extent=[0, zmax, 0, ymax])
    # plt.colorbar(im1, ax=ax1, pad=0.02)
    cbax = ax1.inset_axes([1.01, 0, 0.02, 1])
    cb = plt.colorbar(im1, cax=cbax, )#orientation="vertical")
    ax1.set_ylabel('y')

    ax1 = axes[1]
    mean = np.mean(data, axis=0)
    central_lineout = data[int(ny/2),:]
    ax1.plot(z, mean, '-', label='averaged lineout')
    ax1.plot(z, central_lineout, '-', label='central lineout')
    ax1.legend(loc='upper right')
    ax1.set_xlabel('z')
    # ax1.set_ylabel('z')

    plt.tight_layout()
    if a.showplots: plt.show()
    path = join(a.dir, '%s_%06d.png'%(figname, ts) )
    plt.savefig(path, pad_inches=0.01)

def spectrum(ts, data):
    mean = np.mean(data, axis=0)
    if 'den' in a.fieldname: mean-=1.
    f, F = FFT(z, mean)
    fig, ax1 = plt.subplots(1,1, figsize=a.figsize)
    ax1.set_title('%s, t = %.2f'%(a.fieldname, ts*dt))
    ax1.plot(f, F)
    if 'x' in a.logscale: ax1.set_xscale('log')
    if 'y' in a.logscale: ax1.set_yscale('log')
    ax1.set_xlabel(r'$k_z$')
    ax1.set_ylabel('Spec. Amp.')
    plt.tight_layout()
    plt.show()

def history():
    # subpath = 'hist'
    # subpath = join(outpath, subpath)
    # mkdir(subpath)
    datafile = join(a.dir, '%s.dat'% (a.fieldname) )
    if not exists(datafile):
        data = np.zeros(ndumps)
        for i, ts in enumerate(timesteps):
            showProgressBar(ndumps-1, i)
            d = get_field(ts)
            # data[i] = np.sum((d-np.mean(d))**2.)/np.size(d)
            data[i] = np.sqrt(np.sum((d-1.)**2.)/np.size(d))
        data.tofile(datafile)
    else:
        data = np.fromfile(datafile)
    
    fig, ax1 = plt.subplots(1,1, figsize=a.figsize)
    ax1.set_title(a.fieldname)
    ax1.plot(times, data)
    ax1.set_xlabel(r'$\omega_{ci} t$')
    ax1.set_ylabel(r'$\sqrt{\langle(\rho-1)^2\rangle}$')

    plt.tight_layout()
    plt.show()


def field_plot(ts):
    print ('t = %.2f' % (ts*dt) )
    data = get_field(ts)
    if a.type=='snap': snapshot(ts, data)
    elif a.type=='spec': spectrum(ts, data)
    elif a.type=='hist': history()
    return
    


# def get_timesteps():
#     global ndumps, timesteps, times 
#     ndumps = int(tmax/(dt*nwrtdata)+1)
#     timesteps = np.zeros(ndumps)
#     for i in range(ndumps):
#         timesteps[i]=1 if i==0 else i*nwrtdata
#     times = timesteps * dt
#     return 

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

# class load_input():
#     def __init__(self, path):
#         self.path = path
#         nml = f90nml.read(self.path + '/data/finput.f90')
#         self.denmin = nml['datum']['denmin']
#         self.nspec = nml['datum']['nspec']
#         self.nx = nml['datum']['nx']
#         self.ny = nml['datum']['ny']
#         self.nz = nml['datum']['nz']
#         self.xmax = nml['datum']['xmax']
#         self.ymax = nml['datum']['ymax']
#         self.zmax = nml['datum']['zmax']
#         self.tmax = nml['datum']['maximum_simulation_time']
#         self.dt = nml['datum']['dtwci']
#         self.nwrtdata = nml['datum']['nwrtdata']
#         ## box axes
#         self.y = np.linspace(0, self.ymax, self.ny)
#         self.z = np.linspace(0, self.zmax, self.nz)
#         ## timesteps
#         self.ndumps = int(self.tmax/(self.dt*self.nwrtdata)+1)
#         self.timesteps = np.zeros(self.ndumps)
#         for i in range(self.ndumps):
#             self.timesteps[i]=1 if i==0 else i*self.nwrtdata
#         self.times = self.timesteps * self.dt


def main():
    t_start = systime.time()

    global a, outpath, \
        denmin, nspec, nx, ny, nz, xmax, ymax, zmax, \
        tmax, dt, nwrtdata, y, z, \
        ndumps, timesteps, times, \
        figname
    a = argument_parse(); print ('arguments parsed')
    print ('simulation path = %s' % a.path)
    if a.test:
        print(a); return

    load_input(a.path)
    # b = load_input(a.path)
    # denmin = b.denmin; nspec=b.nspec; nx=b.nx; ny=b.ny; nz=b.nz
    # xmax=b.xmax; ymax=b.ymax; zmax=b.zmax; tmax=b.tmax 
    # dt=b.dt; nwrtdata=b.nwrtdata
    # y=b.y; z=b.z
    # ndumps=b.ndumps; timesteps=b.timesteps; times=b.times
    print ('input parameters loaded')
    # get_timesteps()

    # output path
    outpath = join(a.path, 'figs'); mkdir(outpath)
    if a.mode=='field':
        figname = '_'.join([a.fieldname,a.type])
    if a.dir: 
        a.dir=join(outpath, a.dir); mkdir(a.dir)
    else:
        a.dir=join(outpath, a.mode); mkdir(a.dir)
        if a.mode=='field':
            a.dir=join(a.dir, figname); mkdir(a.dir)
    print('Output figures directory: ', a.dir)
    

    # analyze a single time, timestep, or make movie plots using a timelist or all available steps
    if a.time:
        id_min=np.abs(timesteps-int(a.time/dt)).argmin()
        a.plot_func(timesteps[id_min])
    elif a.timestep:
        a.plot_func(a.timesteps)
    elif a.movie:
        if a.timelist: 
            timesteps_load = timesteps.copy()
            times = np.arange(a.timelist[0], a.timelist[1]+a.timelist[2]/2.0, a.timelist[2])
            timesteps=np.zeros(len(times))
            for i, j in enumerate(times):
                id_min=np.abs(timesteps_load-int(j/dt)).argmin()
                timesteps[i]=timesteps_load[id_min]
        print (times)
        processor = functools.partial(a.plot_func)
        with mp.Pool(processes=a.nprocs) as worker_pool:
            frames = worker_pool.map(processor, timesteps)
            worker_pool.close()
            worker_pool.join()
        makeMovie(a.dir, '%s.mp4'%figname, fps=a.fps)

    print ('time used: %.2fs'% (systime.time()-t_start) )
    print ('====================================================')
    return



if __name__ == '__main__':
    main()

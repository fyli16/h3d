''' Fourier analysis for h3d, for publication '''
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import os

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = \
[r"\usepackage{amsmath, bm}",
 r"\DeclareMathAlphabet{\mathsfit}{\encodingdefault}{\sfdefault}{m}{sl}",
 r"\SetMathAlphabet{\mathsfit}{bold}{\encodingdefault}{\sfdefault}{bx}{sl}",
 r"\newcommand{\tensorsym}[1]{\bm{\mathsfit{#1}}}"]

def read_binary(filename,shape):
  fd=open(filename,'rb')
 # log.update('Reading '+filename)
  # read Fortran binary array
  data=np.fromfile(filename,dtype='f4').reshape((nx,ny,nz),order='F')
  fd.close()
  return data

def fmt(x,pos):
  # formatting log scale colorbar legend
  a,b='{:.1e}'.format(x).split('e')
  b=int(b)
  return r'${}\times 10^{{{}}}$'.format(a,b)

def spectrum(field,fname,time):
  # calculate and plot 2D power spectrum of a field array
  field-=np.mean(field) #remove the average
  p=abs(np.fft.rfftn(field))**2
  #p=2*np.log(p)
  #p=abs(np.fft.rfftn(ex))**2
  levels= np.max(p)*np.power(10,np.linspace(-4,0,20))
  print(levels)
  norm=colors.BoundaryNorm(boundaries=levels,ncolors=256)
  plt.contourf(p[0:31,0:31,0],levels=levels,norm=norm,cmap='jet')
  plt.colorbar(format=ticker.FuncFormatter(fmt))
  #plt.xlim(0,20)
  #plt.ylim(0,20)
  plt.title(fname+', t=%.1f'%time)
  plt.ylabel(r'$k_xL_x/2\pi$', fontsize=14)
  # plt.grid()

def get_angle(x,y):
  if x==0:
    if y>0:
      th=pi/2
    elif y<0:
      th=3*pi/2
    else:
      print('WARNING: ill defined')
      th=None
  else:
    th=np.arctan(y/x)
    if x<0:
      th+=pi
    if th<0:
      th+=2*pi
  return th

def spectrum_1D(field):
  # calculate power spectrum vs k_para and k_perp
  # k_perp=k_z
  field-=np.mean(field) #remove the average
  p=abs(np.fft.rfftn(field))**2
  print('p.shape=',p.shape)
  p_para=np.zeros(p.shape[2])
  for k in range(p.shape[2]):
    for i in range(p.shape[0]):
      for j in range(p.shape[1]):
        p_para[k]+=p[i,j,k]
  p_para/=p.shape[0]*p.shape[1]
  k_para=2*pi/zmax*np.arange(p.shape[2])

  nkx=int(p.shape[0]/2)
  nky=int(p.shape[1]/2)
  xs=2*pi/xmax*np.arange(nkx)
  ys=2*pi/ymax*np.arange(nky)
  x,y=np.meshgrid(xs,ys)
  f=np.zeros((nkx,nky))
  for i in range(nkx):
    for j in range(nky):
      for k in range(p.shape[2]):
        f[i,j]+=p[i,j,k]
  f/=p.shape[2]
  # map x,y to kperp, theta
  kp=[]
  th=[]
  g=[]
  for i in range(len(xs)):
    for j in range(len(ys)):
      angle=get_angle(x[i,j],y[i,j])
      if angle!=None:
        kp.append(np.sqrt(x[i,j]**2+y[i,j]**2))
        th.append(get_angle(x[i,j],y[i,j]))
        g.append(f[i,j])
  dkp=np.sqrt((x[0,1]-x[0,0])**2+(y[1,0]-y[0,0])**2)
  nkp=int(max(kp)/dkp)+1
  nth_max=int((nkx*2+nky*2)*1.4)
  g2=np.zeros((nkp,nth_max))
  th2=np.zeros((nkp,nth_max))
  ind=np.zeros(nkp,dtype=np.int)
  # bin kperp
  for i in range(len(kp)):
    ik=int(kp[i]/dkp)
    #print(i,kp[i],th[i],ik)
    g2[ik,ind[ik]]=g[i]
    th2[ik,ind[ik]]=th[i]
    ind[ik]+=1
  h=np.zeros(nkp)
  # average over irregular theta's
  for ik in range(nkp):
    sort=np.argsort(th2[ik,0:ind[ik]])
    #print(sort)
    th=np.array(th2[ik,sort])
    g=np.array(g2[ik,sort])
    #print(th)
    #print(g)
    for j in range(ind[ik]-1):
        h[ik]+=(g[j]+g[j+1])*(th[j+1]-th[j])
    h[ik]+=(g[-1]+g[0])*(2*pi+th[0]-th[-1])
  h/=4*pi
  k_perp=np.arange(nkp)*dkp
  return k_para,p_para,k_perp,h

def fit_power(x,y):
  # fit a power law y=x^p
  x2=np.log(x)
  y2=np.log(y)
  p=np.polyfit(x2,y2,1)
  yfit=np.exp(x2*p[0]+p[1])
  return p,yfit


def single_mode(field,nx,ny,nz):
  field-=np.mean(field) #remove the average
  p=abs(np.fft.rfftn(field))**2
  return p[nx,ny,nz]

def mode_power(field,wavenumber):
  ''' Return the wave power of selected modes'''
  field-=np.mean(field) #remove the average
  p=abs(np.fft.rfftn(field))**2
  output=[]
  for nx,ny,nz in wavenumber:
    output.append(p[nx,ny,nz])
  return np.array(output)


def spectrum_sum(names):
    print(names)
    fullpath='../../data/%s.gda'%(names[0])
    b=read_binary(fullpath,shape)
    k_para,p_para,k_perp,p_perp=spectrum_1D(b)
    for name in names[1:]:
      fullpath='../../data/%s.gda'%(name)
      b=read_binary(fullpath,shape)
      k_para,p_para1,k_perp,p_perp1=spectrum_1D(b)
      p_para+=p_para1
      p_perp+=p_perp1
    return k_para,p_para,k_perp,p_perp

def plot_1D_spec(fname,snaps,title,fit=None,kperponly=False):
  # averaged over a few snapshots
  if fit==None:
    fit=True
  names=[]
  if fname=='bp':
    for i in range(len(snaps)):
      names.append('bx_%d'%snaps[i])
      names.append('by_%d'%snaps[i])
  elif fname=='ep':
    for i in range(len(snaps)):
      names.append('ex_%d'%snaps[i])
      names.append('ey_%d'%snaps[i])
  else:
    for i in range(len(snaps)):
      names.append('%s_%d'%(fname,snaps[i]))
  # Load or process data
  filename='%s_%d_%d.npz'%(fname,snaps[0],snaps[-1])
  if os.path.exists(filename):
    print('Loading data file: '+filename)
    data=np.load(filename)
    k_para=data['k_para']
    k_perp=data['k_perp']
    p_para=data['p_para']
    p_perp=data['p_perp']
  else:
    k_para,p_para,k_perp,p_perp=spectrum_sum(names)
    p_perp/=len(snaps)
    p_para/=len(snaps)
    np.savez(filename,k_para=k_para,k_perp=k_perp, p_para=p_para, p_perp=p_perp)
  line,=ax1.loglog(k_perp,k_perp*p_perp,'-',label=r'%s$(k_\perp)$'\
      %(title))
  if not kperponly:
    ax1.loglog(k_para,p_para,'--',lw=1.5,color=line.get_color(),label=\
        r'%s$(k_\parallel)$'%(title))
    if fit:
      # fit k_parallel to a power law
      nstart=2
      nend=30
      p,yfit=fit_power(k_para[nstart:nend],p_para[nstart:nend])
      print(p)
      ax1.plot(k_para[nstart:nend],yfit,':',color=line.get_color())
      ax1.text(k_para[int((nstart+nend)/2)],yfit[int(len(yfit)/2)]*2,r'$k_\parallel^{%.1f}$'%p[0],color=line.get_color())
  if fit:
    # fit k_perp to a power law
    nstart=2
    nend=30
    p,yfit=fit_power(k_perp[nstart:nend],p_perp[nstart:nend])
    print(p)
    ax1.plot(k_perp[nstart:nend],yfit,':',color=line.get_color())
    ax1.text(k_perp[int((nstart+nend)/2)],yfit[int(len(yfit)/2)]*2,r'$k_\perp^{%.1f}$'%p[0],color=line.get_color())
  print(max(p_perp))
  #plt.savefig('fft1d_%s.eps'%fname)
  #plt.show()

nx=32
ny=32
nz=32
nt=1
dt=0.025
exec(open('./para.out').read())
shape=(nx,ny,nz)
pi=np.pi

if True:
  fit=False
  #snaps=[20000]
  snaps=[20000,19000,21000]
#  snaps=[20000,20040,20080,20120,20160,20200]
  
  plt.figure(figsize=(9,4))
  if True:
    #plot_1D_spec('bx',snaps,r'$B_x$',fit=fit)
    #plot_1D_spec('by',snaps,r'$B_y$',fit=fit,kperponly=True)
    ax1=plt.subplot(121)
    plot_1D_spec('bx',snaps,r'$\delta B_x$',fit=fit)
    #plot_1D_spec('bp',snaps,r'$\delta B_\perp$',fit=fit)
    plot_1D_spec('by',snaps,r'$\delta B_y$',fit=fit)
    plot_1D_spec('bz',snaps,r'$\delta B_z$',fit=fit) #ax1=plt.subplot(121)
    #plot_1D_spec('den',snaps,r'$\delta n$',fit=fit) #ax1=plt.subplot(121)
    k=np.array([0.08,0.6])
    ind=-5./3
    ind=-2.8
    y=k**ind*2e3
    ax1.plot(k,y,':',color='k')
    ax1.text(k[0]*1.2,y[0]*1.2,r'$k^{%.1f}$'%ind, fontsize=14)
    ax1.text(0.15,0.9, '(a)', transform=ax1.transAxes,fontsize=16)
    # ax1.grid()
    ax1.legend(loc='best')
    ax1.set_ylabel('Power', fontsize=14)
    ax1.set_xlabel(r'$kd_i$', fontsize=14)
    ax1.set_ylim(1e2,1e8)
    ax1.tick_params(bottom=True, top=True, left=True, right=True)
    ax1.tick_params(axis='x', which='minor', direction='in')
    ax1.tick_params(axis='x', which='major', direction='in')
    ax1.tick_params(axis='y', which='minor', direction='in', right=True)
    ax1.tick_params(axis='y', which='major', direction='in')
  
  if True:
    ax1=plt.subplot(122)
    plot_1D_spec('ex',snaps,r'$\delta E_x$',fit=fit)
    #plot_1D_spec('bp',snaps,r'$\delta B_\perp$',fit=fit)
    plot_1D_spec('ey',snaps,r'$\delta E_y$',fit=fit)
    plot_1D_spec('ez',snaps,r'$\delta E_z$',fit=fit) #ax1=plt.subplot(121)
    k=np.array([0.08,0.6])
    ind=-5./3
    ind=-2.0
    y=k**ind*2e4
    ax1.plot(k,y,':',color='k')
    ax1.text(k[0]*1.2,y[0]*1.2,r'$k^{%.1f}$'%ind, fontsize=14)
    ax1.text(0.15,0.9, '(b)', transform=ax1.transAxes,fontsize=16)
    # ax1.grid()
    ax1.legend(loc='best')
    ax1.set_ylabel('Power', fontsize=14)
    ax1.set_xlabel(r'$kd_i$', fontsize=14)
    ax1.set_ylim(1e2, 1e8)
    ax1.tick_params(bottom=True, top=True, left=True, right=True)
    ax1.tick_params(axis='x', which='minor', direction='in')
    ax1.tick_params(axis='x', which='major', direction='in')
    ax1.tick_params(axis='y', which='minor', direction='in', right=True)
    ax1.tick_params(axis='y', which='major', direction='in')
  
  #plt.suptitle(r'$t\Omega_i=%d$'%(snaps[0]*dt))
  plt.tight_layout()
  # plt.savefig('fft1d_p.png')
  plt.savefig('img_pdf/fft1d_p.pdf')
  if True:
    plt.figure(figsize=(5,4))
    ax1=plt.subplot(111)
    plot_1D_spec('bz',snaps,r'$\delta B_z$',fit=fit) #ax1=plt.subplot(121)
    plot_1D_spec('den',snaps,r'$\delta n$',fit=fit) #ax1=plt.subplot(121)

    k=np.array([0.10,0.6])
    ind=-5./3
    y=k**ind*4e4
    ax1.plot(k,y,':',color='k')
    ax1.text(k[0]*2,y[0]*0.6,r'$k^{%.1f}$'%ind, fontsize=14)

    # ax1.grid()
    ax1.legend(loc='best')
    ax1.set_ylabel('Power', fontsize=14)
    ax1.set_xlabel(r'$kd_i$', fontsize=14)
    ax1.set_ylim(1e2, 1e8)
    ax1.tick_params(bottom=True, top=True, left=True, right=True)
    ax1.tick_params(axis='x', which='minor', direction='in')
    ax1.tick_params(axis='x', which='major', direction='in')
    ax1.tick_params(axis='y', which='minor', direction='in', right=True)
    ax1.tick_params(axis='y', which='major', direction='in')
    plt.tight_layout()
    # plt.savefig('fft1d_p.png')
    plt.savefig('img_pdf/fft1d_n.pdf')
  


  plt.show()


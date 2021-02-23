''' Fourier analysis for h3d '''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker

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
  plt.ylabel(r'$k_xL_x/2\pi$')
  plt.grid()

def spectrum_x(field,fname,ax):
  # calculate and plot 1D power spectrum of a field array along x
  field-=np.mean(field) #remove the average
  p=abs(np.fft.rfftn(field))**2
  print('p.shape=',p.shape)
  ps0=int(p.shape[0]/2)
  p_para=np.zeros(ps0)
  for i in range(ps0):
    for j in range(p.shape[1]):
      for k in range(p.shape[2]):
        p_para[i]+=p[i,j,k]
  p_para/=p.shape[1]*p.shape[2]
  k_para=2*pi/xmax*np.arange(ps0)

  ax.semilogy(k_para,p_para,label=fname)
  #plt.xlim(0,20)
  #plt.ylim(0,20)

def spectrum_y(field,fname,ax):
  # calculate and plot 1D power spectrum of a field array along y
  field-=np.mean(field) #remove the average
  p=abs(np.fft.rfftn(field))**2
  print('p.shape=',p.shape)
  ps0=int(p.shape[1]/2)
  p_para=np.zeros(ps0)
  for j in range(ps0):
    for i in range(p.shape[0]):
      for k in range(p.shape[2]):
        p_para[i]+=p[i,j,k]
  p_para/=p.shape[0]*p.shape[2]
  k_para=2*pi/xmax*np.arange(ps0)

  ax.semilogy(k_para,p_para,label=fname)
  #plt.xlim(0,20)
  #plt.ylim(0,20)

def spectrum_z(field,fname,ax):
  # calculate and plot 1D power spectrum of a field array along z
  field-=np.mean(field) #remove the average
  p=abs(np.fft.rfftn(field))**2
  print('p.shape=',p.shape)
  ps0=int(p.shape[2])
  p_para=np.zeros(ps0)
  for k in range(ps0):
    for i in range(p.shape[0]):
      for j in range(p.shape[1]):
        p_para[k]+=p[i,j,k]
  p_para/=p.shape[0]*p.shape[1]
  k_para=2*pi/zmax*np.arange(ps0)

  ax.semilogy(k_para,p_para,label=fname)
  #plt.xlim(0,20)
  #plt.ylim(0,20)

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

def plot_1D_spec(fname,snaps,title,fits=None,kperponly=False):
  plt.figure()
  ax1=plt.subplot(111)
  if fits==None:
    fits=np.repeat(True,len(snaps))
  for i in range(len(snaps)):
    it=snaps[i]
    print('it=%d'%it)
    fullpath='../../data/%s_%d.gda'%(fname,it)
    bx=read_binary(fullpath,shape)
    k_para,p_para,k_perp,p_perp=spectrum_1D(bx)
    #print(k_para)
    #print(p_para)
    #print(k_perp)
    #print(p_perp)
    line,=ax1.loglog(k_perp,p_perp,'-',label=r'$g(k_\perp), t\Omega_i=%d$'%(it*dt))
    if not kperponly:
      ax1.loglog(k_para,p_para,'--',lw=1.5,color=line.get_color(),label=r'$f(k_\parallel), t\Omega_i=%d$'%(it*dt))
      if fits[i]:
        # fit k_parallel to a power law
        nstart=2
        nend=30
        p,yfit=fit_power(k_para[nstart:nend],p_para[nstart:nend])
        print(p)
        ax1.plot(k_para[nstart:nend],yfit,':',color=line.get_color())
        ax1.text(k_para[int((nstart+nend)/2)],yfit[int(len(yfit)/2)]*2,r'$k_\parallel^{%.1f}$'%p[0],color=line.get_color())

      
    if fits[i]:
      # fit k_perp to a power law
      nstart=2
      nend=30
      p,yfit=fit_power(k_perp[nstart:nend],p_perp[nstart:nend])
      print(p)
      ax1.plot(k_perp[nstart:nend],yfit,':',color=line.get_color())
      ax1.text(k_perp[int((nstart+nend)/2)],yfit[int(len(yfit)/2)]*2,r'$k_\perp^{%.1f}$'%p[0],color=line.get_color())
  ax1.grid()
  ax1.legend(loc='best')
  ax1.set_ylabel('Power')
  ax1.set_title(title)
  ax1.set_xlabel(r'$kd_i$')
  ax1.set_ylim(ymin=max(min(p_perp),max(p_perp)*1e-6))
  print(max(p_perp))
  plt.savefig('fft1d_%s.png'%fname)
  plt.savefig('fft1d_%s.eps'%fname)
  plt.show()

nx=32
ny=32
nz=32
nt=1
dt=0.025
exec(open('./para.out').read())
shape=(nx,ny,nz)
pi=np.pi

if False:
  snaps=[1,24000,46000,50000,70000,100000,143000,145000,160000,200000]
  snaps=[70000]
  snaps=[143000,145000,160000,200000]
  snaps=[1,20000,40000,42000,60000]
  snaps=[6000,8000,10000]
  snaps=[1000,2000,4000,5000]
  snaps=[6000,7000,8000]
  plt.figure(figsize=(8,8))
  for it in snaps:
    print('it=%d'%it)
    plt.clf()
    plt.subplot(211)
    fullpath='../../data/'+'bz_%d.gda'%it
    bx=read_binary(fullpath,shape)
    spectrum(bx,'Bz',it*dt)

    plt.subplot(212)
    fullpath='../../data/'+'ex_%d.gda'%it
    den=read_binary(fullpath,shape)
    spectrum(den,'Ex',it*dt)

    plt.xlabel(r'$k_yL_y/2\pi$')

    plt.savefig('fft_%d.png'%it)
  plt.show()

if False:
  amp=[]
  time=[]
  modes=[]
  for i in range(2):
    for j in range(2):
      for k in range(40):
        modes.append((i,j,k))
  for i in range(30):
    for j in range(30):
      for k in range(10):
        modes.append((i,j,k))
  mode_history=np.zeros((nt,len(modes)))
  for name in ['bx','den']:
    print(name)
    for i in range(nt):
      it=1 if i==0 else i*1000
      print('it= %d'%it)
      fullpath='../../data/'+'%s_%d.gda'%(name,it)
      den=read_binary(fullpath,shape)
      mode_history[i,:]=mode_power(den,modes)
      time.append(it*dt)
    np.savez('mode_history_%s.npz'%name,modes=modes,power=mode_history)

if False:
  snaps=[1000,50000,100000]
  snaps=[50000]
  snaps=[24000, 48000]
  snaps=[10000, 20000]
  snaps=[1000, 20000,40000]
  snaps=[1000, 10000,20000]
  snaps=[1000, 10000,20000]
  snaps=[10000,15000,20000]
  snaps=[2000,4000,6000]
  snaps=[20000,40000,60000]
  fits=[False,True,False]
  fits=[False,False,False]
  
  #plot_1D_spec('bx',snaps,r'$B_x$',fits=fits)
  #plot_1D_spec('by',snaps,r'$B_y$',fits=fits,kperponly=True)
  plot_1D_spec('bx',snaps,r'$B_x$',fits=fits)
  plot_1D_spec('by',snaps,r'$B_y$',fits=fits)
  plot_1D_spec('bz',snaps,r'$B_z$',fits=fits)
  plot_1D_spec('ex',snaps,r'$E_x$',fits=fits)
  plot_1D_spec('ey',snaps,r'$E_y$',fits=fits)
  plot_1D_spec('ez',snaps,r'$E_z$',fits=fits)
  #plot_1D_spec('bz',snaps,r'$B_z$',fits=fits,kperponly=True)
  #plot_1D_spec('ex',snaps,r'$E_x$')
  #plot_1D_spec('ey',snaps,r'$E_y$')
  #plot_1D_spec('ez',snaps,r'$E_z$')

if False:
  snaps=[20000,40000,60000,80000]
  snaps=[1000,5000,10000]
  snaps=[6000,8000,10000]
  snaps=range(1000,11000,1000)
  for it in snaps:
    plt.figure(figsize=(8,8))
    ax1=plt.subplot(211)
    ax2=plt.subplot(212)
    print('it=%d'%it)
    fullpath='../../data/'+'bx_%d.gda'%it
    by=read_binary(fullpath,shape)
    spectrum_x(by,'Bx',ax1)

    fullpath='../../data/'+'by_%d.gda'%it
    by=read_binary(fullpath,shape)
    spectrum_x(by,'By',ax1)

    fullpath='../../data/'+'bz_%d.gda'%it
    bz=read_binary(fullpath,shape)
    spectrum_x(bz,'Bz',ax1)

    fullpath='../../data/'+'ex_%d.gda'%it
    ey=read_binary(fullpath,shape)
    spectrum_x(ey,'Ex',ax2)

    fullpath='../../data/'+'ey_%d.gda'%it
    ey=read_binary(fullpath,shape)
    spectrum_x(ey,'Ey',ax2)

    fullpath='../../data/'+'ez_%d.gda'%it
    ez=read_binary(fullpath,shape)
    spectrum_x(ez,'Ez',ax2)

    ax1.grid()
    ax2.grid()
    ax1.legend()
    ax2.legend()
    ax1.set_xlim(0,13)
    ax2.set_xlim(0,13)
    ax1.set_title(r'$t\Omega_i=%d$'%(it*dt))
    ax2.set_xlabel(r'$k_xd_i$')
    plt.savefig('fft_x_%d.png'%it)
    plt.show()

if False:
  snaps=[6000,8000,10000]
  snaps=range(1000,11000,1000)
  for it in snaps:
    plt.figure(figsize=(8,8))
    ax1=plt.subplot(211)
    ax2=plt.subplot(212)
    print('it=%d'%it)
    fullpath='../../data/'+'bx_%d.gda'%it
    by=read_binary(fullpath,shape)
    spectrum_x(by,'Bx',ax1)

    fullpath='../../data/'+'by_%d.gda'%it
    by=read_binary(fullpath,shape)
    spectrum_x(by,'By',ax1)

    fullpath='../../data/'+'bz_%d.gda'%it
    bz=read_binary(fullpath,shape)
    spectrum_x(bz,'Bz',ax1)

    fullpath='../../data/'+'ex_%d.gda'%it
    ey=read_binary(fullpath,shape)
    spectrum_x(ey,'Ex',ax2)

    fullpath='../../data/'+'ey_%d.gda'%it
    ey=read_binary(fullpath,shape)
    spectrum_x(ey,'Ey',ax2)

    fullpath='../../data/'+'ez_%d.gda'%it
    ez=read_binary(fullpath,shape)
    spectrum_x(ez,'Ez',ax2)

    ax1.grid()
    ax2.grid()
    ax1.legend()
    ax2.legend()
    ax1.set_xlim(0,13)
    ax2.set_xlim(0,13)
    ax1.set_title(r'$t\Omega_i=%d$'%(it*dt))
    ax2.set_xlabel(r'$k_yd_i$')
    plt.savefig('fft_y_%d.png'%it)
    plt.show()

if False:
  snaps=[6000,8000,10000]
  snaps=range(1000,11000,1000)
  for it in snaps:
    plt.figure(figsize=(8,8))
    ax1=plt.subplot(211)
    ax2=plt.subplot(212)
    print('it=%d'%it)
    fullpath='../../data/'+'bx_%d.gda'%it
    by=read_binary(fullpath,shape)
    spectrum_z(by,'Bx',ax1)

    fullpath='../../data/'+'by_%d.gda'%it
    by=read_binary(fullpath,shape)
    spectrum_z(by,'By',ax1)

    fullpath='../../data/'+'bz_%d.gda'%it
    bz=read_binary(fullpath,shape)
    spectrum_z(bz,'Bz',ax1)

    fullpath='../../data/'+'ex_%d.gda'%it
    ey=read_binary(fullpath,shape)
    spectrum_z(ey,'Ex',ax2)

    fullpath='../../data/'+'ey_%d.gda'%it
    ey=read_binary(fullpath,shape)
    spectrum_z(ey,'Ey',ax2)

    fullpath='../../data/'+'ez_%d.gda'%it
    ez=read_binary(fullpath,shape)
    spectrum_z(ez,'Ez',ax2)

    ax1.grid()
    ax2.grid()
    ax1.legend()
    ax2.legend()
    ax1.set_xlim(0,13)
    ax2.set_xlim(0,13)
    ax1.set_title(r'$t\Omega_i=%d$'%(it*dt))
    ax2.set_xlabel(r'$k_zd_i$')
    plt.savefig('fft_z_%d.png'%it)
    plt.show()

if True:
  snaps=[1000,40000]
  plt.figure(figsize=(8,8))
  ax1=plt.subplot(211)
  ax2=plt.subplot(212)
  for it in snaps:
    print('it=%d'%it)
    fullpath='../../data/'+'bx_%d.gda'%it
    by=read_binary(fullpath,shape)
    spectrum_z(by,'Bx, t=%d'%(it*dt),ax1)

    fullpath='../../data/'+'by_%d.gda'%it
    by=read_binary(fullpath,shape)
    spectrum_z(by,'By, t=%d'%(it*dt),ax1)

    fullpath='../../data/'+'bz_%d.gda'%it
    bz=read_binary(fullpath,shape)
    spectrum_z(bz,'Bz, t=%d'%(it*dt),ax1)

    fullpath='../../data/'+'ex_%d.gda'%it
    ey=read_binary(fullpath,shape)
    spectrum_z(ey,'Ex, t=%d'%(it*dt),ax2)

    fullpath='../../data/'+'ey_%d.gda'%it
    ey=read_binary(fullpath,shape)
    spectrum_z(ey,'Ey, t=%d'%(it*dt),ax2)

    fullpath='../../data/'+'ez_%d.gda'%it
    ez=read_binary(fullpath,shape)
    spectrum_z(ez,'Ez, t=%d'%(it*dt),ax2)

  ax1.grid()
  ax2.grid()
  ax1.legend()
  ax2.legend()
  #ax1.set_xlim(0,13)
  #ax2.set_xlim(0,13)
  #ax1.set_title(r'$t\Omega_i=%d$'%(it*dt))
  ax2.set_xlabel(r'$k_zd_i$')
  plt.savefig('fft_z_comp_%d.png'%it)
  plt.show()

if False:
  # generate w-kz plot
  its=range(6000,10000,25)
  field_zt=np.zeros((len(its),nz))
  plt.figure()
  for ii in range(len(its)):
    it=its[ii]
    print('it=',it)
    fullpath='../../data/'+'by_%d.gda'%it
    field=read_binary(fullpath,shape)
    field_z=np.zeros(nz)
    for i in range(nx):
      for j in range(ny):
        field_z+=field[i,j,:]
    field_z/=nx*ny
    field_zt[ii,:]=field_z
  field_zt-=np.mean(field_zt) #remove the average
  p=abs(np.fft.rfftn(field_zt))**2
  print('p.shape=',p.shape)
  w=2*np.pi/((its[-1]-its[0])*dt)*np.arange(p.shape[0])
  k=2*np.pi/zmax*np.arange(p.shape[1])
  levels= np.max(p)*np.power(10,np.linspace(-5,0,20))
  print(levels)
  norm=colors.BoundaryNorm(boundaries=levels,ncolors=256)
  plt.contourf(k,w,p,norm=norm,levels=levels)
  plt.colorbar(format=ticker.FuncFormatter(fmt))
  plt.xlabel(r'$k_\parallel d_i$')
  plt.ylabel(r'$\omega/\Omega_p$')
  plt.xlim(0,2.)
  plt.ylim(0,4.)
  plt.grid()
  plt.savefig('dispersion.png')
  plt.show()


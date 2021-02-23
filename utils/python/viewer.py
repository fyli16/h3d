''' GUI post-processing tool for H3D simulations 
    by Xiangrong "Sean" Fu, 2/6/2017
'''
import os
import numpy as np
from tkinter import *
import tkinter.scrolledtext as tkst
import matplotlib
matplotlib.use("TkAgg")
# from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk as NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from tkinter.filedialog import askopenfilename,askdirectory

class Logger(tkst.ScrolledText):
  def __init__(self, widget):
    tkst.ScrolledText.__init__(self)
    self.parent = widget
    self.config(state=DISABLED)
  def update(self, line):
    self.config(state=NORMAL)
    self.insert(END, line+'\n')
    self.config(state=DISABLED)
    print(line)

def callback():
  ''' dummy function '''
  print("called the callback!")

def process_ts():
  ''' Process time series '''
  global nt,times,den_history, b_history, e_history
  filename=fullpath+'/history.npz'
  if not os.path.exists(filename):
    log.update('Processing time series')
    b_history=np.zeros(nt)
    e_history=np.zeros(nt)
    den_history=np.zeros(nt)
    w=1./(nx*ny*nz)
    for i in range(nt):
      den_history[i]=np.sum((dict_field['den'][i,:,:,:]-1.)**2)*w
      b_history[i]=np.sum(dict_field['bx'][i,...]**2+dict_field['by'][i,...]**2+\
          (dict_field['bz'][i,...]-1.)**2)*w
      e_history[i]=np.sum(dict_field['ex'][i,...]**2+dict_field['ey'][i,...]**2+\
          dict_field['ez'][i,...]**2)*w
    np.savez(filename,times=times,den=den_history,b=b_history,e=e_history)
  else:
    log.update('Loading '+filename)
    data=np.load(filename)
    den_history=data['den']
    b_history=data['b']
    e_history=data['e']
    times=data['times']
    nt=len(times)
  show_history()

def check_time_window(event):
  if t1.get()<t0.get():
    t1.set(t0.get())

def cal_growth():
  global curve
  start=t0.get()
  end=t1.get()
  if start==end:
    log.update("WARNING: t0=t1!")
  else:
    x=times[start:end]*dt
    y=np.log(den_history[start:end])/2
    p=np.polyfit(x,y,1)
    log.update("Fitting: %s"%p)
    y=np.exp(2*(p[0]*x+p[1]))
    try:
      ax3.lines.remove(curve)
    except:
      log.update("Curve does not exist!")


    curve,=ax3.plot(x,y,'--r',label='gam=%.4f'%p[0])
      
    ax3.legend()
    canvas1.draw()

def show_history():
  global t0,t1,ax3,canvas1
  win1=Toplevel()
  win1.wm_title('Time History')
  fig = Figure(figsize=(8,8))
  ax1 = fig.add_subplot(311)
  ax2 = fig.add_subplot(312)
  ax3 = fig.add_subplot(313)

  canvas1 = FigureCanvasTkAgg(fig, win1)
  canvas1.show()
  canvas1.get_tk_widget().pack(side=TOP)

  toolbar = NavigationToolbar2TkAgg(canvas1, win1)
  toolbar.update()
  canvas1._tkcanvas.pack(side=TOP)

  toolbar=Frame(win1)
  toolbar.pack(side=TOP)
  t0=IntVar(win1)
  t1=IntVar(win1)
  sl1 = Scale(toolbar, label='t start',length=200, from_=0, to=nt-1, variable=t0, \
      orient=HORIZONTAL,command=check_time_window)
  sl1.pack(side=LEFT)
  sl2 = Scale(toolbar, label='t end',length=200, from_=0, to=nt-1, variable=t1,\
      orient=HORIZONTAL, command=check_time_window)
  sl2.pack(side=LEFT)
  bt1= Button(toolbar,text='Cal Growth',command=cal_growth)
  bt1.pack(side=LEFT)

  t=dt*np.array(times)
  ax1.semilogy(t,b_history,'-b',label='dB^2')
  ax1.legend()
  ax1.grid()
  ax2.semilogy(t,e_history,'-g',label='E^2')
  ax2.legend()
  ax2.grid()
  ax3.semilogy(t,den_history,label='rho^2')
  ax3.legend()
  ax3.grid()
  ax3.set_xlabel('t')
  canvas1.draw()

def fourier():
  ''' Fourier Analysis '''
  print("FFT!")


def quit():
  import sys;
  sys.exit()

def change_plane(event):
  p=plane.get()
  #disable_slides(p)

def disable_slides(p):
  ''' disable the two axis in the plane '''
  if p=='x-y':
    sl_x.config(state=DISABLED)
    sl_y.config(state=DISABLED)
    sl_z.config(state=NORMAL)
  elif p=='x-z':
    sl_x.config(state=DISABLED)
    sl_y.config(state=NORMAL)
    sl_z.config(state=DISABLED)
  else:
    sl_x.config(state=NORMAL)
    sl_y.config(state=DISABLED)
    sl_z.config(state=DISABLED)

def load_single_file():
  global bx,by,bz,den
  fullpath=askopenfilename(initialdir='../../data')
  if len(fullpath)>1:
    head,tail=os.path.split(fullpath)
    datadir.set(head)
    filename=tail
    #status.set('Loading '+filename+' ...')
    if filename[0:2]=='bx':
      bx=read_binary(fullpath,shape)
      field.set('bx')
    elif filename[0:2]=='by':
      by=read_binary(fullpath,shape)
      field.set('by')
    elif filename[0:2]=='bz':
      bz=read_binary(fullpath,shape)
      field.set('bz')
    elif filename[0:3]=='den':
      den=read_binary(fullpath,shape)
      field.set('den')
    #status.set('Ready')

def load_all():
  global times,nt,fullpath
  p=askdirectory(initialdir=fullpath)
  if not p=='':
    fullpath=p
    datadir.set(p)
    flist=os.listdir(fullpath)
    ts=[]
    fname=dict_field.keys()[0]
    log.update('Probing %s_*.gda'%fname)
    for f in flist:
      if f[0:len(fname)]==fname:
        t=int(f.lstrip(fname+'_').rstrip('.gda'))
        ts.append(t)
    if len(ts)==0:
      log.update('No data file found!')
    else:
      ts.sort()
      times=ts
      nt=len(times)
      log.update('number of snapshots = %d'%nt)
      sl_t.config(to=nt-1)
      for fname in dict_field:
        ff=dict_field[fname]=np.zeros((nt,)+shape)
        for j in range(nt):
          #status.set('Reading %d'%j)
          filename=fname+'_%d.gda'%times[j]
          ff[j,...]=read_binary(fullpath+'/'+filename,shape)
      #status.set('Ready')

def read_binary(filename,shape):
  fd=open(filename,'rb')
  log.update('Reading '+filename)
  # read Fortran binary array
  data=np.fromfile(filename,dtype='f4').reshape((nx,ny,nz),order='F')
  fd.close()
  return data

def plot():
  ax.clear()
  ax_r.clear()
  ax_b.clear()
  fname=field.get()
  i=ix.get()
  j=iy.get()
  k=iz.get()
  ff=dict_field[fname]
  p=plane.get()
  t=it.get()
  # contour plots
  if p=='x-y':
    ax.contourf(ff[t,:,:,k].transpose())
    ax.set_xlabel('x')
    ax.set_ylabel('y')
  elif p=='x-z':
    ax.contourf(ff[t,:,j,:].transpose())
    ax.set_xlabel('x')
    ax.set_ylabel('z')
  else: # y-z
    ax.contourf(ff[t,i,:,:].transpose())
    ax.set_xlabel('y')
    ax.set_ylabel('z')
  ax.set_title(fname+', t=%d'%(times[it.get()]*dt))
  # cuts
  c=cut.get()

  # horizontal
  if p=='y-z':
    if c=='slice':
      ax_b.plot(ff[t,i,:,k])
    elif c=='average':
      ax_b.plot(np.average(ff[t,i,:,:],axis=1))
  else: # x-y or x-z
    if c=='slice':
      ax_b.plot(ff[t,:,j,k])
  ax_b.set_xticklabels([])
  ax_b.set_xlim(ax.get_xlim())
  ax_b.grid()
  if c=='slice':
    if p=='x-y':
    #  ax_b.set_title('y=%d'%j)
      ax.axhline(j,color='w',ls='--',lw=2)
    else:
    #  ax_b.set_title('z=%d'%k)
      ax.axhline(k,color='w',ls='--',lw=2)

  # vertical
  if p=='x-y':
    ax_r.plot(ff[t,i,:,k],range(ny))
  else: # x-z or y-z
    ax_r.plot(ff[t,i,j,:],range(nz))
  ax_r.set_yticklabels([])
  ax_r.set_ylim(ax.get_ylim())
  for tick in ax_r.get_xticklabels():
    tick.set_rotation(270)
  ax_r.grid()
  if c=='slice':
    if p=='y-z':
    #  ax_r.set_title('y=%d'%k)
      ax.axvline(j,color='w',ls='--',lw=2)
    else: # x-y or x-z
    #  ax_r.set_title('x=%d'%i)
      ax.axvline(i,color='w',ls='--',lw=2)
  canvas.draw()

nx=32
ny=32
nz=32
nt=1
dt=0.025
exec(open('./para.out').read())

times=[0]

shape=(nx,ny,nz)
shape2=(nt,)+shape
empty=np.zeros(shape2)
dict_field={'bx':empty,'by':empty,'bz':empty,'ex':empty,'ey':empty,'ez':empty,'den':empty,\
    'vix':empty,'viy':empty,'viz':empty}
fullpath='../../data'
den_history=np.zeros(nt)
b_history=np.zeros(nt)
e_history=np.zeros(nt)

root = Tk() # root window
root.wm_title('H3D Viewer')

datadir=StringVar(root)
datadir.set(os.path.abspath(fullpath))
ix=IntVar(root)
ix.set(0)
iy=IntVar(root)
iy.set(0)
iz=IntVar(root)
iz.set(0)
it=IntVar(root)
it.set(0)

menu = Menu(root)
root.config(menu=menu)

filemenu = Menu(menu)
menu.add_cascade(label="File", menu=filemenu)
filemenu.add_command(label="Set Direcotry...", command=callback)
filemenu.add_separator()
filemenu.add_command(label="Exit", command=quit)

procmenu = Menu(menu)
menu.add_cascade(label="Process", menu=procmenu)
procmenu.add_command(label="Time Series...", command=process_ts)
procmenu.add_command(label="Fourier Analysis...", command=fourier)

helpmenu = Menu(menu)
menu.add_cascade(label="Help", menu=helpmenu)
helpmenu.add_command(label="About...", command=callback)

toolbar1 = Frame(root)

bt1 = Button(toolbar1, text="Load All", width=4, command=load_all)
bt1.pack(side=LEFT, padx=2, pady=2)

#bt2 = Button(toolbar1, text="Load", width=4, command=load_single_file)
#bt2.pack(side=LEFT, padx=2, pady=2)


bt4 = Button(toolbar1, text="Export", width=4, command=callback)
bt4.pack(side=LEFT, padx=2, pady=2)

field = StringVar(root)
l=dict_field.keys()
#l.sort()
field.set(list(l)[0]) # default value
#om1 = apply(OptionMenu, (toolbar1, field)+tuple(l))
om1 = OptionMenu(toolbar1, field,tuple(l))
om1.config(width=5)
om1.pack(side=LEFT,padx=2,pady=2)

plane = StringVar(root)
if nx==1:
  plane.set("y-z")
elif ny==1:
  plane.set("x-z")
else:
  plane.set("x-y") # default value
om2 = OptionMenu(toolbar1, plane, "x-y", "x-z", "y-z",command=change_plane)
om2.config(width=4)
om2.pack(side=LEFT,padx=2,pady=2)

cut = StringVar(root)
cut.set("slice") # default value
om3 = OptionMenu(toolbar1, cut,"slice","average")
om3.config(width=5)
om3.pack(side=LEFT,padx=2,pady=2)

bt3 = Button(toolbar1, text="PLOT", width=6, command=plot)
bt3.pack(side=LEFT, padx=4, pady=2,expand=TRUE)

toolbar2 = Frame(root)
sl_x = Scale(toolbar2, label='x index',from_=0, to=nx-1, variable=ix, orient=HORIZONTAL)
sl_x.pack(side=LEFT)

sl_y = Scale(toolbar2, label='y index',from_=0, to=ny-1, variable=iy, orient=HORIZONTAL)
sl_y.pack(side=LEFT)

sl_z = Scale(toolbar2, label='z index',from_=0, to=nz-1, variable=iz, orient=HORIZONTAL)
sl_z.pack(side=LEFT)

sl_t = Scale(toolbar2, label='t index',from_=0, to=nt-1, variable=it,\
    length=200, orient=HORIZONTAL)
sl_t.pack(side=LEFT)
#disable_slides(plane.get())

toolbar3 = Frame(root)
lb_dir=Label(toolbar3, text="Current Direcotry")
lb_dir.pack(side=LEFT)
tx_dir=Label(toolbar3,textvariable=datadir,bg='white')
tx_dir.pack(side=LEFT,padx=2)

toolbar1.pack(side=TOP, fill=X)
toolbar3.pack(side=TOP, fill=X)
toolbar2.pack(side=TOP, fill=X)

fig = Figure(figsize=(8,8))
gs=GridSpec(3,3,hspace=0.3)
ax = fig.add_subplot(gs[0:2,0:2])
ax_r = fig.add_subplot(gs[0:2,2])
ax_b = fig.add_subplot(gs[2,0:2])
ax_r.set_yticklabels([])
ax_b.set_xticklabels([])

canvas = FigureCanvasTkAgg(fig, root)
print (dir(canvas))
# canvas.show()
canvas.draw()
canvas.get_tk_widget().pack(side=TOP)

toolbar = NavigationToolbar2TkAgg(canvas, root)
toolbar.update()
canvas._tkcanvas.pack(side=TOP)

#status=StringVar(root)
#status.set('Ready')
#sb = Label(root, textvariable=status, bd=1, relief=SUNKEN, anchor=W)
#sb.pack(side=BOTTOM, fill=X)

log=Logger(root)
log.pack(side=TOP, fill=X)
#log.config(state=DISABLED)

root.mainloop()

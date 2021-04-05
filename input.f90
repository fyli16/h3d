&datum

! ------------------ global simulation info -----------!
tmax=1000.0, ! max sim. time, in units of 1/wci
t_begin=0.0, t_end=2000.0, 
dtwci=0.01,  ! value of dt*wci
restart=.false.,  ! whether to restart from 'restart' directory
quota=24.0,  ! walltime quota (in unit of hours)
MPI_IO_format = .true. ! use MPI IO instead of traditional binary output

! ------------------ simulation domain ----------------!
nx=1, ny=4, nz=224,  ! number of cells along each dim
xmax=1., ymax=4., zmax=224.,  ! max lengths of each dim
npx=10, npy=40, npz=2240,  ! number of particles along each dim over full length (not cell)

! number of nodes (cores) along y, z; (no decompostion along x)
! make sure npy, npz can be divided by nodey, nodez, respectively
nodey=2, nodez=16  

! boundaries of the uniform region
! setting xbb/ybb/zbb to xmax/ymax/zmax would leave only the uniform region to be simulated
xaa=0., xbb=1., nax=0, nbx=1
yaa=0., ybb=4., nay=0, nby=4
zaa=0., zbb=224., naz=0, nbz=224

uniform_loading_in_logical_grid = .false., ! used in loading particles? see 'init waves'

buffer_zone=0., ! 'epsilon=buffer_zone' in parmov.f90 ??
moat_zone=3., ! seems not implemented
profile_power=0,  ! seems not implemented

! ------------------ field solver ----------------!
n_subcycles=0
nskipx=1,  nskipy=1,  nskipz=1, ! not implemented?

iterb=5,  ! ion push can use a larger step than field advance
testorbt=.false., 
norbskip=1,

! ------------------ plasma setup ----------------!
nspec=1,  ! number of ion species, maximum 5
qspec=1., ! charge of each ion species (use array if nspec>1 and the same for rest)
wspec=1., ! mass of each ion species

frac=1.0 ! means normalized to n0 (associated with wpi)
denmin=0.05,  ! density floor. When density is smaller than this value, force it to this value to avoid divergence in calculating E field
wpiwci=400., ! ratio of ion plasma frequency to ion cyclotron frequency
btspec=0.01, ! beta of each ion species 
bete=0.01,  ! beta of electrons

! resistivity 
ieta=0,  ! other models include ieta=1,2,3,4,5; see 'etacal.f90'
resis=1.e-6,  ! ieta=0 model; constant resisitivity, i.e., eta=resis
netax=10, netay=2, netaz=5, ! used in ieta=1 model; tho netaz seems not used
etamin=1.0e-6, etamax=5.0e-5,  ! used in ieta>1 models
eta_par=0, ! parallel resisitivity? sth used in 'field.f90'

! anisotropy in velocity
anisot=1.0, ! anisotropy of velocity for each species, used in 'init waves'
gama=1.66667, ! gamma factor in EoS

ave1=100.0, ave2=50.,
phib=180.0,

! performing density/velocity smoothing
smoothing=.true., 
smooth_coef=0. ! seems not implemented

! ---------------------- init waves --------------------!
dB_B0=0.1
num_cycles=5

! ------------------ diagnostic control ----------------!
nprint=100,  ! frequency at which to print simulation information
nwrtdata=1000, ! frequency at which to write data into files
nwrtparticle=4000,  ! frequency at which to write particles within a box range

restrt_write=1,  ! whether to write restart files
nwrtrestart=20000000, ! frequency at which to write restart files

! box range within which particles will be dumped
xbox_l=0., xbox_r=1.0, 
ybox_l=0., ybox_r=1.0, 
zbox_l=0., zbox_r=2.24,

! ------------------------- others ---------------------!
! these seem not really implemented
Yee=.false., ! ?
global=.true., ! ?
harris=.false., ! ?
! sth used in 'io.f90'
fxsho=1.0,  ! seems not implemented
nxcel=4,  ! seems not implemented

/

&datum

! ------------------ global simulation info -----------!
tmax=100.0, ! max sim. time, in units of 1/wci
dtwci=0.01,  ! value of dt*wci
restart=.false.,  ! whether to restart from 'restart' directory
MPI_IO_format =.true. ! use MPI IO (one file only) instead of traditional binary output

! ------------------ simulation domain ----------------!
nx=1, ny=4, nz=224,  ! number of cells along each dim
xmax=1., ymax=4., zmax=224.,  ! max lengths of each dim
npx=10,10, npy=40,40, npz=2240,2240,  ! number of particles along each dim over full length (not cell)

! number of nodes (cores) along y, z; (no decompostion along x)
! make sure npy, npz can be divided by nodey, nodez, respectively
nodey=2, nodez=16  

! boundaries of the uniform region
! setting xbb/ybb/zbb to xmax/ymax/zmax would leave only the uniform region to be simulated
xaa=0., xbb=1., nax=0, nbx=1
yaa=0., ybb=4., nay=0, nby=4
zaa=0., zbb=224., naz=0, nbz=224

uniform_loading_in_logical_grid =.false., ! used in loading particles? see 'init waves'

! ------------------ field solver ----------------!
n_subcycles=0
nskipx=1,  nskipy=1,  nskipz=1, ! not implemented?

iterb=5,  ! ion push can use a larger step than field advance
testorbt=.false.,  ! test orbit?
norbskip=1,

! ------------------ plasma setup ----------------!
nspec=1,  ! number of ion species, maximum 5
qspec=1., ! charge of each ion species (use array if nspec>1 and the same for rest)
wspec=1., ! mass of each ion species

frac=1.0 ! means normalized to n0 (associated with wpi)
denmin=0.05,  ! when density is smaller than this value, force it to this value to avoid divergence in calculating E field
wpiwci=400., ! ratio of ion plasma frequency to ion cyclotron frequency
btspec=0.01, ! beta of each ion species 
bete=0.01,  ! beta of electrons

! resistivity 
ieta=6,  ! other models ieta=1,2,3,4,5,6; see 'etacal.f90'
resis=1.e-1,  ! ieta=0 model; constant resisitivity, i.e., eta=resis
netax=10, netay=2 ! used in ieta=1 model
etamin=1.0e-6, etamax=5.0e-5,  ! used in ieta>1 models
eta_par=0, ! parallel resisitivity? sth used in 'field.f90'
eta_zs=28, ! scale length of resistive layer along z (in unit of cell size)

! anisotropy in velocity
anisot=1.0, ! anisotropy of velocity for each species, used in 'init waves'
gamma=1.66667, ! gamma factor in EoS

ave1=100.0, ave2=50.,
phib=180.0,

! performing density/velocity smoothing
smoothing=.true., 
smooth_coef=0. ! seems not implemented

! ---------------------- init waves --------------------!
dB_B0=0.1,
num_cycles=5,

! ------------------ diagnostic control ----------------!
n_print=100,  ! frequency at which to print simulation information
n_write_data=1000, ! frequency at which to write data into files
n_write_particle=400000,  ! frequency at which to write particles within a box range
n_write_restart=2000000, ! frequency at which to write restart files

tracking_binary=.true. ! write tracking data in binary (unformatted) or formatted form
tracking_mpi=.true. ! write tracking files by mpi rank

! box range within which particles will be dumped
xbox_l=0., xbox_r=1.0, 
ybox_l=0., ybox_r=1.0, 
zbox_l=0., zbox_r=2.24,

! ------------------------- others ---------------------!
! sth used in 'io.f90'
fxsho=1.0,  ! seems not implemented
nxcel=4,  ! seems not implemented

/

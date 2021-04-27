&datum

! ------------------ global simulation info -----------!
tmax = 10.0, ! max sim. time, in units of 1/wci
dtwci = 0.01,  ! value of dt*wci
restart = .false.,  ! whether to restart from 'restart' directory
MPI_IO_format = .true. ! use MPI IO (one file only) instead of traditional binary output

! MPI nodes(ranks) configuration along y, z (no decompostion along x)
! and whether the ranks are treated periodic in both directions
node_conf(:) = 2, 16,
periods(:) = .true., .true.,

! ------------------ simulation domain ----------------!
nx = 1, ny = 4, nz = 2240,  ! number of cells along each dim
xmax = 1., ymax = 4., zmax = 2240.,  ! max lengths of each dim
npx(1:5) = 10, ! number of particles along x over full length (not one cell) for maximum 5 ion species
npy(1:5) = 40, 
npz(1:5) = 22400,  

! boundaries of the uniform region
! setting xbb/ybb/zbb to xmax/ymax/zmax would leave only the uniform region to be simulated
xaa = 0., xbb = 1., nax = 0, nbx = 1
yaa = 0., ybb = 4., nay = 0, nby = 4
zaa = 0., zbb = 2240., naz = 0, nbz = 2240

! uniform loading in logical space
! used in loading particles? see 'init waves'
uniform_load_logical = .false., 

! ------------------ field solver ----------------!
n_subcycles = 0
nskipx = 1,  nskipy = 1,  nskipz = 1, ! not implemented?

iterb = 5,  ! ion push can use a larger step than field advance
eta_par = 0, ! parallel resisitivity? options: 0, 1, 2

! ------------------ plasma setup ----------------!
nspec = 1,  ! number of ion species, maximum 5
qspec(1:5) = 1., ! charge of each ion species (use array if nspec>1 and the same for rest)
wspec(1:5) = 1., ! mass of each ion species
frac(1:5) = 1., ! density normalized to n0 (associated with wpi)

denmin = 0.05,  ! when density is smaller than this value, force it to this value to avoid divergence in calculating E field
wpiwci = 400., ! ratio of ion plasma frequency to ion cyclotron frequency
beta_spec(1:5) = 0.01, ! beta of each ion species 
beta_e = 0.01, ! beta of electrons
n_sort = 10, ! frequency at which to sort particles

! resistivity 
ieta = 6,  ! other models ieta=1,2,3,4,5,6; see 'etacal.f90'
resis = 1.e-3,  ! ieta=0 model; constant resisitivity, i.e., eta=resis
netax = 10, netay = 2 ! used in ieta=1 model
etamin = 1.0e-6, etamax = 5.0e-5,  ! used in ieta>1 models
eta_zs = 280, ! scale length of resistive layer along z (in unit of cell size)

! anisotropy in velocity
anisot(1:5) = 1.0, ! anisotropy of velocity for each species, used in 'init waves'
gamma = 1.66667, ! gamma factor in EoS

ave1 = 100.0, ave2 = 50.,
phib = 180.0,

! performing density/velocity smoothing
smoothing = .true., 
smooth_pass = 1, 

! ---------------------- init waves --------------------!
dB_B0 = 0.1,
num_wave_cycles = 32,

! ------------------ diagnostic control ----------------!
n_print = 100,  ! frequency at which to print simulation progression

n_write_mesh = 1000, ! frequency at which to write mesh data 
n_write_energy = 100, ! frequency at which to write integrated energy data
n_write_probes = 0, ! frequency at which to write field probe data
n_write_tracking = 0, ! frequency at which to write tracking particle data
n_write_particle = 0, ! frequency at which to write particles within a volume
n_write_restart = 0, ! frequency at which to write restart files

tracking_binary = .true. ! write tracking data in binary (unformatted) or formatted form
tracking_mpi = .true. ! write tracking data by mpi rank

! volume within which particles will be dumped
xbox_l = 0., xbox_r = 1.0,
ybox_l = 0., ybox_r = 1.0, 
zbox_l = 0., zbox_r = 2.24,

! ------------------------- others ---------------------!
! sth used in 'io.f90'
fxsho = 1.0,  ! seems not implemented
nxcel = 4,  ! seems not implemented

/

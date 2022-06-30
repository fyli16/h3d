&input

! ------------------ global simulation info -----------!
tmax = 1200.0, ! max sim. time, in units of 1/wci
dtwci = 0.01,  ! value of dt*wci
restart = .false.,  ! whether to restart from 'restart' directory
MPI_IO_format = .true., ! use MPI IO (one file only) instead of traditional binary output

! MPI nodes(ranks) configuration along y, z (no decompostion along x)
! and whether the ranks are treated periodic in either directions
node_conf(:) = 2, 28,
periods(:) = .true., .true.,

! simulation domain
nx = 1, ny = 4, nz = 1200,  ! total number of cells along each dim
xmax = 1., ymax = 4., zmax = 1200.,  ! max lengths of each dim

! uniform loading in logical space
! used in loading particles? see 'init waves'
uniform_load_logical = .false., 

! ------------------ field solver ----------------!
n_sub_b = 5, ! number of subcycles for advancing B field
eta_par = 0, ! ? options: 0, 1, 2; used in 'ecal'

! field masking
mask = .true., ! if to perform field masking
mask_zs = 200, ! scale length (in cell) of field masking in z
mask_r = 1., ! factor r in field masking, which controls the slope of mask function

! initial waves
dB_B0 = 0.01, ! Alfven wave amplitude
n_wave_cycles = 30.0, ! number of wave cycles that would occcupy the full box dimension in z
! upramp, flat, and downramp of the wave envelope
! these only function when mask==.true.
wave_upramp = 200,  ! wave upramp length (in cell)
wave_flat = 200,  ! wave central flat length (in cell)
wave_downramp = 200, ! wave downramp length (in cell)
! sign of cos\theta which determines wave propagation direction;
! +1: along B0; -1: opposite to B0
sign_cos = 1.0, 

! ------------------ plasma setup ----------------!
nspec = 1,  ! number of ion species, maximum 5
qspec(1:5) = 1., ! charge of each ion species 
wspec(1:5) = 1., ! mass of each ion species
frac(1:5) = 1., ! density normalized to n0 (associated with wpi)
beta_spec(1:5) = 0.01, ! beta of each ion species 
beta_elec = 0.01, ! beta of electrons

ppcx(1:5) = 4, ! number of particles per cell along x 
ppcy(1:5) = 4, ! number of particles per cell along y 
ppcz(1:5) = 4, ! number of particles per cell along z 

wpiwci = 400., ! ratio of ion plasma frequency to ion cyclotron frequency 
denmin = 0.05,  ! force density lower than this to this value
n_sort = 10, ! frequency at which to sort particles

! resistivity 
ieta = 0,  ! available models ieta=1,2,3,4,5,6
resis = 1.e-6,  ! constant resisitivity (for ieta=0)
netax = 10, netay = 2 ! (for ieta=1)
etamin = 1.0e-6, etamax = 5.0e-5,  ! (for ieta>0)
eta_zs = 200, ! scale length (in cell) of resistive layer in z (for ieta=6)

! anisotropy in velocity
anisot(1:5) = 1.0, ! anisotropy of velocity for each species
gamma = 1.66667, ! gamma factor in EoS

! density/velocity smoothing
smoothing = .true., 
smooth_pass = 1, 

! ------------------ diagnostic control ----------------!
n_print = 100,  ! frequency at which to print simulation progression

n_diag_mesh = 100, ! frequency at which to write mesh data 
n_diag_energy = 100, ! frequency at which to write integrated energy data

n_diag_probe = 100, ! frequency at which to write field probe data
probe_x = 2, ! x location (in cell) of the probes

n_diag_tracking = 0, ! frequency at which to write tracking particle data
n_diag_particle = 0, ! frequency at which to write particles within a volume

n_write_restart = 0, ! frequency at which to write restart files

tracking_binary = .true., ! write tracking data in binary (unformatted) or formatted form
tracking_mpi = .true., ! write tracking data by mpi rank

! volume within which particles will be dumped
xbox_l = 0., xbox_r = 1.0,
ybox_l = 0., ybox_r = 1.0, 
zbox_l = 0., zbox_r = 2.0,

/

&input

! ------------------ global simulation info -----------!
tmax = 2000.0, ! max sim. time, in units of 1/wci
dtwci = 0.01,  ! value of dt*wci
restart = .false.,  ! whether to restart from 'restart' directory
MPI_IO_format = .true., ! use MPI IO (one file only) instead of traditional binary output

! MPI nodes(ranks) configuration along y, z (no decompostion along x)
! and whether the ranks are treated periodic in either directions
node_conf(:) = 16, 28,
periods(:) = .true., .true.,

! simulation domain
nx = 80, ny = 80, nz = 180,  ! total number of cells along each dim
xmax = 2., ymax = 2., zmax = 90.,  ! max lengths of each dim

! uniform loading in logical space
! used in loading particles? see 'init waves'
uniform_load_logical = .false., 

! ------------------ field solver ----------------!
n_sub_b = 5, ! number of subcycles for advancing B field
eta_par = 0, ! ? options: 0, 1, 2; used in 'ecal'

! field masking
mask = .true., ! if to perform field masking
mask_zs = 60, ! scale length (in cell) of field masking in z
mask_r = 0.1, ! factor r in field masking, which controls the slope of mask function
mask_B0_fac = 1.0, ! reduction factor for B0 in the mask layers (purpose: slow down wave & better absorption) 

! load waves initially inside the box
dB_B0         = 0.0, ! Alfven wave amplitude
wave_cycles   = 9.0, ! number of wave cycles that would fill the box in z
! upramp, flat, and downramp of the wave envelope
! these only function when mask==.true.
wave_upramp   = 200,  ! wave upramp length (in cell)
wave_flat     = 200,  ! wave central flat length (in cell)
wave_downramp = 200, ! wave downramp length (in cell)
sign_cos      = 1.0, ! sign of cos\theta which determines wave propagation direction; +1: along B0; -1: opposite to B0

! ------------------ inject waves  ----------------!
inj_waves_b = .false., ! inject waves via B field
inj_waves_bv = .false.,  ! inject waves via BV field
inj_waves_e = .false.,  ! inject waves via E field

inj_waves_b_rmf = .true.,  ! inject 3d RMF antenna waves via B field; all waves share the same loop radius inj_wave_radius(1)
inj_rmf_ampl_corr = 0.31892948790361647, ! amplitude correction to let amplitude at wave center = inj_dB_B0

! inject max. 4 waves; waves that are injected at the same location should be placed next to 
!    each other in the array, so their V, B will be added up during injection
inj_dB_B0(1:4)       = 1e-2,    1e-2,    0.0,    0.0,   ! injection wave amplitude
inj_wave_cycles(1:4) = 9.0,     9.0,     30.0,   30.0,  ! number of wave cycles used to determine kz
inj_sign_cos(1:4)    = 1.0,     1.0,     -1.0,    -1.0,  ! sign of cos\theta which determines wave propagation dir.
inj_wave_pol(1:4)    = 0,         1,       0,      1, ! polarization; 0:x, 1:y (left-hand), -1:y (right-hand)
inj_wave_radius(1:4) = 10,        10,       0,       0, ! raidus of wave injection (in cell)

! injection properties
inj_z_pos(1:4)       = 70,    70,      0,      0,   ! injection z position (in cell)
inj_t_upramp(1:4)    = 50.0   50.0,  200.0,  200.0,  ! injection upramp time (in 1/wci)
inj_t_flat(1:4)      = 1e8,   1e8,  200.0,  200.0,  ! injection flat time (in 1/wci)
inj_t_downramp(1:4)  = 200.0, 200.0,  200.0,  200.0,  ! injection downramp time (in 1/wci)


! ------------------ plasma setup ----------------!
nspec          = 1,  ! number of ion species, maximum 5
qspec(1:5)     = 1., ! charge of each ion species 
wspec(1:5)     = 1., ! mass of each ion species
frac(1:5)      = 1., ! density normalized to n0 (associated with wpi)
beta_spec(1:5) = 1e-4, ! beta of each ion species 
beta_elec      = 4e-4, ! beta of electrons

ppcx(1:5) = 4, ! number of particles per cell along x 
ppcy(1:5) = 4, ! number of particles per cell along y 
ppcz(1:5) = 4, ! number of particles per cell along z 

wpiwci = 300., ! ratio of ion plasma frequency to ion cyclotron frequency 
denmin = 0.05,  ! force density lower than this to this value
n_sort = 10, ! frequency at which to sort particles

! resistivity 
ieta   = 0,      ! available models ieta=1,2,3,4,5,6
resis  = 1.0e-5,  ! for ieta=0: constant resisitivity 
netax  = 10,     ! for ieta=1   
netay  = 2 
etamin = 1.0e-6, ! for ieta>0
etamax = 5.0e-5,
eta_zs = 200,    ! for ieta=6: z scale length (in cell) of the resistive layers

! anisotropy in velocity
anisot(1:5) = 1.0, ! anisotropy of velocity for each species
gamma       = 1.66667, ! gamma factor in EoS

! density/velocity smoothing
smoothing   = .true., 
smooth_pass = 1, 

! ------------------ diagnostic control ----------------!
n_print = 100,  ! frequency at which to print simulation progression

n_diag_mesh = 1000, ! frequency at which to write mesh data 
n_diag_energy = 100, ! frequency at which to write integrated energy data

n_diag_probe = 100, ! frequency at which to write field probe data
probe_x = 40, ! x location (in cell) of the probes

n_diag_tracking = 0, ! frequency at which to write tracking particle data
n_diag_particle = 0, ! frequency at which to write particles within a volume

n_write_restart = 0, ! frequency at which to write restart files

tracking_binary = .true., ! write tracking data in binary (unformatted) or formatted form
tracking_mpi = .true., ! write tracking data by mpi rank

! volume within which particles will be dumped
xbox_l = 0., xbox_r = 1.0,
ybox_l = 0., ybox_r = 4, 
zbox_l = 0., zbox_r = 90.0,

/

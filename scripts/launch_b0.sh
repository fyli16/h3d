#!/bin/sh


for i in 0.3 0.4
do 
    echo "for b=${i} ..."
    mkdir 1d-b${i}
    cd 1d-b${i}
    git clone git@github.com:fyli16/h3d.git .
    git checkout dev

    #sed -i "s/#SBATCH -N 1/#SBATCH -N 2/" scripts/job_tacc
    #sed -i "s/#SBATCH -n 56/#SBATCH -n 112/" scripts/job_tacc 
    #sed -i "s/#SBATCH -t 01:00:00/#SBATCH -t 02:00:00/" scripts/job_tacc 

    sed -i "s/dB_B0=0.1/dB_B0=${i}/" input.f90
    #sed -i "s/num_cycles=5/num_cycles=24/" input.f90

    # sed -i "s/btspec=0.01/btspec=${j}/" input.f90
    # sed -i "s/bete=0.01/bete=${j}/" input.f90

    #sed -i "s/nz=224/nz=2240/" input.f90
    #sed -i "s/zmax=224/zmax=2240/" input.f90
    #sed -i "s/npz=2240/npz=22400/" input.f90
    #sed -i "s/zbb=224/zbb=2240/" input.f90
    #sed -i "s/nbz=224/nbz=2240/" input.f90
    # sed -i "s/npx=10/npx=12/" input.f90
    # sed -i "s/npy=40/npy=48/" input.f90
    # sed -i "s/npz=2240/npz=2688/" input.f90

    #sed -i "s/nodez=28/nodez=56/" input.f90

    #sed -i "s/resis=1.e-6/resis=0./" input.f90

    sbatch scripts/job_nersc
    cd ..
done


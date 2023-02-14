# Multi-echo Spectrum Imaging on Diffusion MRI

Wu, Y., Liu, X.M., Ahmad, S., & Yap, P. T. (2023). Multi-echo Spectrum Imaging on Diffusion MRI In *Medical Image Computing and Computer Assisted Intervention–MICCAI. (To be submission)



## ME_SMSI.m

###### Multi-echo spherical mean spectrum imaging

**Input**

- Multi-echo diffusion MRI
- Echo times

**Output**

- Diffusion relaxation (restricted, hindere, isotropic)
- Non-diffusion relaxation spectrum coefficient (restricted, hindere, isotropic)



## ME_SI.m

###### Multi-echo spectrum imaging

**Input**

- Multi-echo diffusion MRI
- Echo times
- Diffusion relaxation (restricted, hindere, isotropic), estimated by ME_SMSI
- Non-diffusion relaxation spectrum coefficient (restricted, hindere, isotropic), estimated by ME_SMSI

**Output**

- Non-diffusion relaxation multi-compartment intra-voxel architecture



## ME_MICRO.m

###### Relaxation-diffusion features estimation

**Input**

- Diffusion relaxation (restricted, hindere, isotropic), estimated by ME_SMSI
- Non-diffusion relaxation spectrum coefficient (restricted, hindere, isotropic), estimated by ME_SMSI
- Non-diffusion relaxation intra-voxel architecture struture

**Output**

- Non-diffusion relaxation microarchitecture (Huynh et al., 2020)
- Non-diffusion relaxation microstructure (White et al., 2013)
  - N0 - The relative amount of cell bodies within a voxe, such as the densities of neurons, astrocytes, and oligodendrocytes. 
  - ND - The relative amount of tube-like structures within a voxel, such as axons and dendrites. 
  - NF - The relative amount of free water outside of cell structures.


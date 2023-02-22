# Multi-echo Spectrum Imaging on Diffusion MRI

Wu, Y., Liu, X.M., Zhang, X.Y., Huynh, K.H., Ahmad, S., & Yap, P. T. (2023). Multi-echo Spectrum Imaging on Diffusion MRI. MICCAI. (Submitted)


## ME_SMSI

###### Multi-echo spherical mean spectrum imaging

**Input**

- Multi-echo diffusion MRI
- Echo times

**Output**

- Diffusion relaxation (restricted, hindered, isotropic)
- Non-diffusion relaxation spectrum coefficient (restricted, hindered, isotropic)


## ME_SI

###### Multi-echo spectrum imaging

**Input**

- Multi-echo diffusion MRI
- Echo times
- Diffusion relaxation (restricted, hindered, isotropic), estimated by ME_SMSI
- Non-diffusion relaxation spectrum coefficient (restricted, hindered, isotropic), estimated by ME_SMSI

**Output**

- Non-diffusion relaxation multi-compartment intra-voxel architecture


## SE_SMSI

###### Single-echo spherical mean spectrum imaging (Huynh et al., 2020)

**Input**

- Single-echo diffusion MRI
- Echo times

**Output**

- With-diffusion relaxation spectrum coefficient (restricted, hindered, isotropic)


## SE_SI 

###### Single-echo spectrum imaging (White et al., 2013)

**Input**

- Multi-echo diffusion MRI
- Echo times

**Output**

- With-diffusion relaxation multi-compartment intra-voxel architecture


## ME_MICRO

###### Relaxation-diffusion features estimation

**Input**

- Diffusion relaxation (restricted, hindered, isotropic), estimated by ME_SMSI
- Non-diffusion relaxation spectrum coefficient (restricted, hindered, isotropic), estimated by ME_SMSI
- Non-diffusion relaxation intra-voxel architecture struture, estimated by ME_SI

**Output**

- Non-diffusion relaxation microarchitecture (Huynh et al., 2020)
- Non-diffusion relaxation microstructure (White et al., 2013)
  - N0 - The relative amount of cell bodies within a voxe, such as the densities of neurons, astrocytes, and oligodendrocytes. 
  - ND - The relative amount of tube-like structures within a voxel, such as axons and dendrites. 
  - NF - The relative amount of free water outside of cell structures.
  - 

## ME_dMICRO

###### Directional relaxation-diffusion features estimation

clear all; clc;
addpath('/home/wuye/Softwares/REDIM');
addpath('/usr/local/freesurfer/7.2.0/matlab');
dataFolder = '/home/wuye/D_disk/MTE_dMRI/axon-relaxation';
TE = [35.5 45.5];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(x),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.nii.gz')),TE,'UniformOutput',false);
fbval   = fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(TE(1)),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.bval'));
fbvec   = fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(TE(1)),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.bvec'));
fmask   = fullfile(dataFolder,'mask.nii');
mkdir(fullfile(dataFolder,'REDIM'));
redim_pipe(fdwi, fullfile(dataFolder,'REDIM','REDIM_l'), fmask, fbvec, fbval, TE,'l');
redim_pipe(fdwi, fullfile(dataFolder,'REDIM','REDIM_wl'), fmask, fbvec, fbval, TE,'wl');
redim_pipe(fdwi, fullfile(dataFolder,'REDIM','REDIM_wn'), fmask, fbvec, fbval, TE,'wn');
redim_pipe(fdwi, fullfile(dataFolder,'REDIM','REDIM_n'), fmask, fbvec, fbval, TE,'n');


clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/axon-relaxation';
fdwi    = fullfile(dataFolder,'REDIM','REDIM_l_RelaxRegressed_dwi.nii.gz');
fbval   = fullfile(dataFolder,'REDIM','REDIM_l_bvals.txt');
fbvec   = fullfile(dataFolder,'REDIM','REDIM_l_bvecs.txt');
fmask   = fullfile(dataFolder,'mask.nii');
SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'REDIM','REDIM_l'));

F_restricted = fullfile(dataFolder,'REDIM','REDIM_l','VF_restricted.nii.gz');
F_hindered   = fullfile(dataFolder,'REDIM','REDIM_l','VF_hindered.nii.gz');
F_isotropic  = fullfile(dataFolder,'REDIM','REDIM_l','VF_free.nii.gz'); 
ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'REDIM','REDIM_l'));

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/axon-relaxation';
fdwi    = fullfile(dataFolder,'REDIM','REDIM_wl_RelaxRegressed_dwi.nii.gz');
fbval   = fullfile(dataFolder,'REDIM','REDIM_wl_bvals.txt');
fbvec   = fullfile(dataFolder,'REDIM','REDIM_wl_bvecs.txt');
fmask   = fullfile(dataFolder,'mask.nii');
SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'REDIM','REDIM_wl'));

F_restricted = fullfile(dataFolder,'REDIM','REDIM_wl','VF_restricted.nii.gz');
F_hindered   = fullfile(dataFolder,'REDIM','REDIM_wl','VF_hindered.nii.gz');
F_isotropic  = fullfile(dataFolder,'REDIM','REDIM_wl','VF_free.nii.gz'); 
ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'REDIM','REDIM_wl'));

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/axon-relaxation';
fdwi    = fullfile(dataFolder,'REDIM','REDIM_wn_RelaxRegressed_dwi.nii.gz');
fbval   = fullfile(dataFolder,'REDIM','REDIM_wn_bvals.txt');
fbvec   = fullfile(dataFolder,'REDIM','REDIM_wn_bvecs.txt');
fmask   = fullfile(dataFolder,'mask.nii');
SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'REDIM','REDIM_wn'));

F_restricted = fullfile(dataFolder,'REDIM','REDIM_wn','VF_restricted.nii.gz');
F_hindered   = fullfile(dataFolder,'REDIM','REDIM_wn','VF_hindered.nii.gz');
F_isotropic  = fullfile(dataFolder,'REDIM','REDIM_wn','VF_free.nii.gz'); 
ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'REDIM','REDIM_wn'));

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/axon-relaxation';
fdwi    = fullfile(dataFolder,'REDIM','REDIM_n_RelaxRegressed_dwi.nii.gz');
fbval   = fullfile(dataFolder,'REDIM','REDIM_n_bvals.txt');
fbvec   = fullfile(dataFolder,'REDIM','REDIM_n_bvecs.txt');
fmask   = fullfile(dataFolder,'mask.nii');
SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'REDIM','REDIM_n'));

F_restricted = fullfile(dataFolder,'REDIM','REDIM_n','VF_restricted.nii.gz');
F_hindered   = fullfile(dataFolder,'REDIM','REDIM_n','VF_hindered.nii.gz');
F_isotropic  = fullfile(dataFolder,'REDIM','REDIM_n','VF_free.nii.gz'); 
ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'REDIM','REDIM_n'));
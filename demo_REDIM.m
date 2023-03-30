% clear all; clc;
% addpath('/home/wuye/Softwares/REDIM');
% addpath('/usr/local/freesurfer/7.2.0/matlab');
% dataFolder = '/home/wuye/D_disk/MTE_dMRI/CW/proc';
% TE = [85 95 105 115 125 135];
% fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
% fbval   = fullfile(dataFolder,strcat('MTE_',num2str(TE(1)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
% fbvec   = fullfile(dataFolder,strcat('MTE_',num2str(TE(1)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
% fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
% mkdir(fullfile(dataFolder,'REDIM'));
% redim_pipe(fdwi, fullfile(dataFolder,'REDIM','REDIM_l'), fmask, fbvec, fbval, TE,'l');
% redim_pipe(fdwi, fullfile(dataFolder,'REDIM','REDIM_wl'), fmask, fbvec, fbval, TE,'wl');
% redim_pipe(fdwi, fullfile(dataFolder,'REDIM','REDIM_wn'), fmask, fbvec, fbval, TE,'wn');
% redim_pipe(fdwi, fullfile(dataFolder,'REDIM','REDIM_n'), fmask, fbvec, fbval, TE,'n');
% 
% 
% clear all; clc;
% dataFolder = '/home/wuye/D_disk/MTE_dMRI/WDZ/proc';
% TE = [75 85 95 105 115 125 135];
% addpath('/home/wuye/Softwares/REDIM');
% addpath('/usr/local/freesurfer/7.2.0/matlab');
% fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
% fbval   = fullfile(dataFolder,strcat('MTE_',num2str(TE(1)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
% fbvec   = fullfile(dataFolder,strcat('MTE_',num2str(TE(1)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
% fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
% mkdir(fullfile(dataFolder,'REDIM'));
% redim_pipe(fdwi, fullfile(dataFolder,'REDIM','REDIM_l'), fmask, fbvec, fbval, TE,'l');
% redim_pipe(fdwi, fullfile(dataFolder,'REDIM','REDIM_wl'), fmask, fbvec, fbval, TE,'wl');
% redim_pipe(fdwi, fullfile(dataFolder,'REDIM','REDIM_wn'), fmask, fbvec, fbval, TE,'wn');
% redim_pipe(fdwi, fullfile(dataFolder,'REDIM','REDIM_n'), fmask, fbvec, fbval, TE,'n');
% 
% clear all; clc;
% dataFolder = '/home/wuye/D_disk/MTE_dMRI/ZJ/proc';
% TE = [75 85 95 105 115 125 135];
% addpath('/home/wuye/Softwares/REDIM');
% addpath('/usr/local/freesurfer/7.2.0/matlab');
% fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
% fbval   = fullfile(dataFolder,strcat('MTE_',num2str(TE(1)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
% fbvec   = fullfile(dataFolder,strcat('MTE_',num2str(TE(1)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
% fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
% mkdir(fullfile(dataFolder,'REDIM'));
% redim_pipe(fdwi, fullfile(dataFolder,'REDIM','REDIM_l'), fmask, fbvec, fbval, TE,'l');
% redim_pipe(fdwi, fullfile(dataFolder,'REDIM','REDIM_wl'), fmask, fbvec, fbval, TE,'wl');
% redim_pipe(fdwi, fullfile(dataFolder,'REDIM','REDIM_wn'), fmask, fbvec, fbval, TE,'wn');
% redim_pipe(fdwi, fullfile(dataFolder,'REDIM','REDIM_n'), fmask, fbvec, fbval, TE,'n');


clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
fdwi    = fullfile(dataFolder,'REDIM','REDIM_l_RelaxRegressed_dwi.nii.gz');
fbval   = fullfile(dataFolder,'REDIM','REDIM_l_bvals.txt');
fbvec   = fullfile(dataFolder,'REDIM','REDIM_l_bvecs.txt');
fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'REDIM','REDIM_l'));

F_restricted = fullfile(dataFolder,'REDIM','REDIM_l','VF_restricted.nii.gz');
F_hindered   = fullfile(dataFolder,'REDIM','REDIM_l','VF_hindered.nii.gz');
F_isotropic  = fullfile(dataFolder,'REDIM','REDIM_l','VF_free.nii.gz'); 
ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'REDIM','REDIM_l'));

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
fdwi    = fullfile(dataFolder,'REDIM','REDIM_wl_RelaxRegressed_dwi.nii.gz');
fbval   = fullfile(dataFolder,'REDIM','REDIM_wl_bvals.txt');
fbvec   = fullfile(dataFolder,'REDIM','REDIM_wl_bvecs.txt');
fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'REDIM','REDIM_wl'));

F_restricted = fullfile(dataFolder,'REDIM','REDIM_wl','VF_restricted.nii.gz');
F_hindered   = fullfile(dataFolder,'REDIM','REDIM_wl','VF_hindered.nii.gz');
F_isotropic  = fullfile(dataFolder,'REDIM','REDIM_wl','VF_free.nii.gz'); 
ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'REDIM','REDIM_wl'));

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
fdwi    = fullfile(dataFolder,'REDIM','REDIM_wn_RelaxRegressed_dwi.nii.gz');
fbval   = fullfile(dataFolder,'REDIM','REDIM_wn_bvals.txt');
fbvec   = fullfile(dataFolder,'REDIM','REDIM_wn_bvecs.txt');
fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'REDIM','REDIM_wn'));

F_restricted = fullfile(dataFolder,'REDIM','REDIM_wn','VF_restricted.nii.gz');
F_hindered   = fullfile(dataFolder,'REDIM','REDIM_wn','VF_hindered.nii.gz');
F_isotropic  = fullfile(dataFolder,'REDIM','REDIM_wn','VF_free.nii.gz'); 
ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'REDIM','REDIM_wn'));

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
fdwi    = fullfile(dataFolder,'REDIM','REDIM_n_RelaxRegressed_dwi.nii.gz');
fbval   = fullfile(dataFolder,'REDIM','REDIM_n_bvals.txt');
fbvec   = fullfile(dataFolder,'REDIM','REDIM_n_bvecs.txt');
fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'REDIM','REDIM_n'));

F_restricted = fullfile(dataFolder,'REDIM','REDIM_n','VF_restricted.nii.gz');
F_hindered   = fullfile(dataFolder,'REDIM','REDIM_n','VF_hindered.nii.gz');
F_isotropic  = fullfile(dataFolder,'REDIM','REDIM_n','VF_free.nii.gz'); 
ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'REDIM','REDIM_n'));


clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/HTY/proc';
fdwi    = fullfile(dataFolder,'REDIM','REDIM_l_RelaxRegressed_dwi.nii.gz');
fbval   = fullfile(dataFolder,'REDIM','REDIM_l_bvals.txt');
fbvec   = fullfile(dataFolder,'REDIM','REDIM_l_bvecs.txt');
fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'REDIM','REDIM_l'));

F_restricted = fullfile(dataFolder,'REDIM','REDIM_l','VF_restricted.nii.gz');
F_hindered   = fullfile(dataFolder,'REDIM','REDIM_l','VF_hindered.nii.gz');
F_isotropic  = fullfile(dataFolder,'REDIM','REDIM_l','VF_free.nii.gz'); 
ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'REDIM','REDIM_l'));

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/HTY/proc';
fdwi    = fullfile(dataFolder,'REDIM','REDIM_wl_RelaxRegressed_dwi.nii.gz');
fbval   = fullfile(dataFolder,'REDIM','REDIM_wl_bvals.txt');
fbvec   = fullfile(dataFolder,'REDIM','REDIM_wl_bvecs.txt');
fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'REDIM','REDIM_wl'));

F_restricted = fullfile(dataFolder,'REDIM','REDIM_wl','VF_restricted.nii.gz');
F_hindered   = fullfile(dataFolder,'REDIM','REDIM_wl','VF_hindered.nii.gz');
F_isotropic  = fullfile(dataFolder,'REDIM','REDIM_wl','VF_free.nii.gz'); 
ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'REDIM','REDIM_wl'));

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/HTY/proc';
fdwi    = fullfile(dataFolder,'REDIM','REDIM_wn_RelaxRegressed_dwi.nii.gz');
fbval   = fullfile(dataFolder,'REDIM','REDIM_wn_bvals.txt');
fbvec   = fullfile(dataFolder,'REDIM','REDIM_wn_bvecs.txt');
fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'REDIM','REDIM_wn'));

F_restricted = fullfile(dataFolder,'REDIM','REDIM_wn','VF_restricted.nii.gz');
F_hindered   = fullfile(dataFolder,'REDIM','REDIM_wn','VF_hindered.nii.gz');
F_isotropic  = fullfile(dataFolder,'REDIM','REDIM_wn','VF_free.nii.gz'); 
ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'REDIM','REDIM_wn'));

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/HTY/proc';
fdwi    = fullfile(dataFolder,'REDIM','REDIM_n_RelaxRegressed_dwi.nii.gz');
fbval   = fullfile(dataFolder,'REDIM','REDIM_n_bvals.txt');
fbvec   = fullfile(dataFolder,'REDIM','REDIM_n_bvecs.txt');
fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'REDIM','REDIM_n'));

F_restricted = fullfile(dataFolder,'REDIM','REDIM_n','VF_restricted.nii.gz');
F_hindered   = fullfile(dataFolder,'REDIM','REDIM_n','VF_hindered.nii.gz');
F_isotropic  = fullfile(dataFolder,'REDIM','REDIM_n','VF_free.nii.gz'); 
ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'REDIM','REDIM_n'));
clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
TE = [75 85 95 105 115];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
fmask   = fullfile(dataFolder,'MTE_mask_ROI.nii.gz');

ME_SMSIx(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSIx_1'),spectrum='scheme/new_spectrum_1.mat');
ME_SMSIx(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSIx_2'),spectrum='scheme/new_spectrum_2.mat');
% ME_SMSIx(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSIx_3'),spectrum='scheme/new_spectrum_3.mat');
ME_SMSIx(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSIx_4'),spectrum='scheme/new_spectrum_4.mat');




clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
F_restricted = fullfile(dataFolder,'ME_SMSIx_1','VF_restricted.nii.gz');
F_hindered   = fullfile(dataFolder,'ME_SMSIx_1','VF_hindered.nii.gz');
F_isotropic  = fullfile(dataFolder,'ME_SMSIx_1','VF_free.nii.gz'); 

fmask   = fullfile(dataFolder,'MTE_mask_ROI.nii.gz');
ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'ME_SMSIx_1'),spectrum='scheme/new_spectrum_1.mat');

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
F_restricted = fullfile(dataFolder,'ME_SMSIx_2','VF_restricted.nii.gz');
F_hindered   = fullfile(dataFolder,'ME_SMSIx_2','VF_hindered.nii.gz');
F_isotropic  = fullfile(dataFolder,'ME_SMSIx_2','VF_free.nii.gz'); 

fmask   = fullfile(dataFolder,'MTE_mask_ROI.nii.gz');
ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'ME_SMSIx_2'),spectrum='scheme/new_spectrum_2.mat');

% clear all; clc;
% dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
% F_restricted = fullfile(dataFolder,'ME_SMSIx_3','VF_restricted.nii.gz');
% F_hindered   = fullfile(dataFolder,'ME_SMSIx_3','VF_hindered.nii.gz');
% F_isotropic  = fullfile(dataFolder,'ME_SMSIx_3','VF_free.nii.gz'); 
% 
% fmask   = fullfile(dataFolder,'MTE_mask_ROI.nii.gz');
% ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'ME_SMSIx_3'),spectrum='scheme/new_spectrum_3.mat');

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
F_restricted = fullfile(dataFolder,'ME_SMSIx_4','VF_restricted.nii.gz');
F_hindered   = fullfile(dataFolder,'ME_SMSIx_4','VF_hindered.nii.gz');
F_isotropic  = fullfile(dataFolder,'ME_SMSIx_4','VF_free.nii.gz'); 

fmask   = fullfile(dataFolder,'MTE_mask_ROI.nii.gz');
ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'ME_SMSIx_4'),spectrum='scheme/new_spectrum_4.mat');

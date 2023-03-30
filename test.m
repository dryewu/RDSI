clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/CW/proc/REDIM';
fdwi    = fullfile(dataFolder,'REDIM_RelaxRegressed_dwi.nii.gz');
fbval   = fullfile(dataFolder,'REDIM_bvals.txt');
fbvec   = fullfile(dataFolder,'REDIM_bvecs.txt');
fmask   = fullfile(dataFolder,'../MTE_mask.nii.gz');
SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'SE_SMSI'));

F_restricted = fullfile(dataFolder,'SE_SMSI','VF_restricted.nii.gz');
F_hindered   = fullfile(dataFolder,'SE_SMSI','VF_hindered.nii.gz');
F_isotropic  = fullfile(dataFolder,'SE_SMSI','VF_free.nii.gz'); 
ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'SE_SMSI'));
% 
% dataFolder = '/home/wuye/D_disk/MTE_dMRI/HTY/proc';
% F_restricted = fullfile(dataFolder,'ME_SMSIx_default','VF_restricted.nii.gz');
% F_hindered   = fullfile(dataFolder,'ME_SMSIx_default','VF_hindered.nii.gz');
% F_isotropic  = fullfile(dataFolder,'ME_SMSIx_default','VF_free.nii.gz'); 
% ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'ME_SMSIx_default'));
% 
% TE = [85 95 105 115 125 135];
% dataFolder = '/home/wuye/D_disk/MTE_dMRI/HTY/proc/SE_SMSI_default';
% for i = 1:length(TE)
%     F_restricted = fullfile(dataFolder,num2str(TE(i)),'VF_restricted.nii.gz');
%     F_hindered   = fullfile(dataFolder,num2str(TE(i)),'VF_hindered.nii.gz');
%     F_isotropic  = fullfile(dataFolder,num2str(TE(i)),'VF_free.nii.gz'); 
%     ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,num2str(TE(i))));
% end
%         
%     
%     
%     
% 

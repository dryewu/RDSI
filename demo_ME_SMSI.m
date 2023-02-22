% clear all; clc;
% dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
% TE = [75 85 95 105 115];
% fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
% fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
% fmask   = fullfile(dataFolder,'MTE_mask_ROI.nii.gz');
% ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI'));

% clear all; clc;
% dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
% fdwi    = fullfile(dataFolder,'DKI_AP_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz');
% fbval   = fullfile(dataFolder,'DKI_AP_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval');
% fmask   = fullfile(dataFolder,'MTE_mask_ROI.nii.gz');
% SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'SE_SMSI'));

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
TE = [75 85 95 105 115];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbvec   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
ft2r    = cellfun(@(x)fullfile(dataFolder,'ME_SMSI',strcat('T2_',num2str(x),'.nii.gz')),{'free','hindered','restricted'},'UniformOutput',false);
fmask   = fullfile(dataFolder,'MTE_mask_ROI.nii.gz');
ME_SI(fdwi,fbvec,fbval,ft2r,fmask,TE,fullfile(dataFolder,'ME_SI'));


% 
% clear all; clc;
% dataFolder = '/home/wuye/D_disk/MTE_dMRI/HTY/proc';
% TE = [85 95 105 115 125 135];
% fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
% fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
% fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
% ME_SMSI(fdwi,fbval,fmask,TE);
% 
% 
% clear all; clc;
% dataFolder = '/home/wuye/D_disk/MTE_dMRI/WDZ/proc';
% TE = [75 85 95 105 115 125 135];
% fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
% fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
% fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
% ME_SMSI(fdwi,fbval,fmask,TE);
% 
% 
% clear all; clc;
% dataFolder = '/home/wuye/D_disk/MTE_dMRI/ZJ/proc';
% TE = [75 85 95 105 115 125 135];
% fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
% fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
% fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
% ME_SMSI(fdwi,fbval,fmask,TE);
% 
% 
% clear all; clc;
% dataFolder = '/home/wuye/D_disk/MTE_dMRI/CW/proc';
% TE = [85 95 105 115 125 135];
% fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
% fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
% fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
% ME_SMSI(fdwi,fbval,fmask,TE);
% 
% 
% clear all; clc;
% dataFolder = '/home/wuye/D_disk/MTE_dMRI/LXF/proc';
% TE = [75 85 95 105 115 125 135];
% fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
% fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
% fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
% ME_SMSI(fdwi,fbval,fmask,TE);


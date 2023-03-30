clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
TE = [75 85 95 105 115];
for i = 1:length(TE)
    fdwi    = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz'));
    fbval   = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
    fbvec   = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
    fmask   = fullfile(dataFolder,'MTE_mask_ROI.nii.gz');
%     SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'SE_SMSI_ROI_default',num2str(TE(i))));
    SE_SI(fdwi,fbvec,fbval,fmask,fullfile(dataFolder,'SE_SI_ROI_default',num2str(TE(i))));
end

% clear all; clc;
% dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
% TE = [75 85 95 105 115];
% for i = 1:length(TE)
%     fdwi    = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz'));
%     fbval   = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
%     fbvec   = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
%     fmask   = fullfile(dataFolder,'MTE_mask_ROI.nii.gz');
%     SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'SE_SMSI_ROI',num2str(TE(i))),spectrum='scheme/new_spectrum.mat');
%     SE_SI(fdwi,fbvec,fbval,fmask,fullfile(dataFolder,'SE_SI_ROI',num2str(TE(i))),spectrum='scheme/new_spectrum.mat');
% end
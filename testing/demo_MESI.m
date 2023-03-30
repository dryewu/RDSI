clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
TE = [75 85 95 105 115];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbvec   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
ft2r    = cellfun(@(x)fullfile(dataFolder,'ME_SMSI',strcat('T2_',num2str(x),'.nii.gz')),{'restricted','hindered','free'},'UniformOutput',false);
fmask   = fullfile(dataFolder,'MTE_mask_ROI.nii.gz');
tic;ME_SI(fdwi,fbvec,fbval,ft2r,fmask,TE,fullfile(dataFolder,'ME_SI_ROI'));toc

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
TE = [75 85 95 105 115];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbvec   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
ft2r    = cellfun(@(x)fullfile(dataFolder,'ME_SMSI',strcat('T2_',num2str(x),'.nii.gz')),{'restricted','hindered','free'},'UniformOutput',false);
fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
tic;ME_SI(fdwi,fbvec,fbval,ft2r,fmask,TE,fullfile(dataFolder,'ME_SI'));toc

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/HTY/proc';
TE = [85 95 105 115 125 135];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbvec   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
ft2r    = cellfun(@(x)fullfile(dataFolder,'ME_SMSI',strcat('T2_',num2str(x),'.nii.gz')),{'restricted','hindered','free'},'UniformOutput',false);
fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
tic;ME_SI(fdwi,fbvec,fbval,ft2r,fmask,TE,fullfile(dataFolder,'ME_SI'));toc

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/ZJ/proc';
TE = [75 85 95 105 115 125 135];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbvec   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
ft2r    = cellfun(@(x)fullfile(dataFolder,'ME_SMSI',strcat('T2_',num2str(x),'.nii.gz')),{'restricted','hindered','free'},'UniformOutput',false);
fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
tic;ME_SI(fdwi,fbvec,fbval,ft2r,fmask,TE,fullfile(dataFolder,'ME_SI'));toc

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/CW/proc';
TE = [85 95 105 115 125 135];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbvec   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
ft2r    = cellfun(@(x)fullfile(dataFolder,'ME_SMSI',strcat('T2_',num2str(x),'.nii.gz')),{'restricted','hindered','free'},'UniformOutput',false);
fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
tic;ME_SI(fdwi,fbvec,fbval,ft2r,fmask,TE,fullfile(dataFolder,'ME_SI'));toc
% 
% 
% clear all; clc;
% dataFolder = '/home/wuye/D_disk/MTE_dMRI/LXF/proc';
% TE = [75 85 95 105 115 125 135];
% fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
% fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
% fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
% tic;ME_SI(fdwi,fbvec,fbval,ft2r,fmask,TE,fullfile(dataFolder,'ME_SI'));toc
% 

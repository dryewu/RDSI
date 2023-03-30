clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/CW/proc';
TE = [85 95 105 115 125 135];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','400'),useBshell=[0,400]);
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','800'),useBshell=[0,800]);
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','1600'),useBshell=[0,1600]);
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','3200'),useBshell=[0,3200]);

TE = [85 135];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','TE_85_135'));

TE = [85 105 125];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','TE_85_105_125'));

TE = [95 115 135];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','TE_95_115_135'));

TE = [85 95 105 115];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','TE_85_95_105_115'));

TE = [105 115 125 135];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','TE_105_115_125_135'));

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
TE = [75 85 95 105 115];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','400'),useBshell=[0,400]);
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','800'),useBshell=[0,800]);
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','1600'),useBshell=[0,1600]);
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','3200'),useBshell=[0,3200]);

TE = [75 115];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','TE_75_115'));

TE = [75 95 115];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','TE_75_95_115'));

TE = [75 85 95];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','TE_75_85_95'));

TE = [95 105 115];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','TE_95_105_115'));







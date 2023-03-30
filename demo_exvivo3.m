
clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/axon-relaxation';
TE = [35.5 45.5];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(x),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(x),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.bval')),TE,'UniformOutput',false);
fbvec   = arrayfun(@(x)fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(x),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.bvec')),TE,'UniformOutput',false);
fmask   = fullfile(dataFolder,'mask.nii.gz');
ME_SMSIx(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSIx_1_50'),spectrum='scheme/exvivo_spectrum_1.mat',x0=50);
ME_SMSIx(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSIx_2_50'),spectrum='scheme/exvivo_spectrum_2.mat',x0=50);

ME_SMSIx(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSIx_1_35'),spectrum='scheme/exvivo_spectrum_1.mat',x0=35);
ME_SMSIx(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSIx_2_35'),spectrum='scheme/exvivo_spectrum_2.mat',x0=35);


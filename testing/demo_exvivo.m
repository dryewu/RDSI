clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/axon-relaxation';
TE = [35.5 45.5];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(x),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(x),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.bval')),TE,'UniformOutput',false);
fmask   = fullfile(dataFolder,'mask.nii.gz');

ME_SMSIx(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI'),spectrum='scheme/exvivo_spectrum.mat');
SE_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'SE_SMSI'),spectrum='scheme/exvivo_spectrum.mat');


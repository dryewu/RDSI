
clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/axon-relaxation';
TE = [35.5 45.5];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(x),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(x),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.bval')),TE,'UniformOutput',false);
fbvec   = arrayfun(@(x)fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(x),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.bvec')),TE,'UniformOutput',false);
fmask   = fullfile(dataFolder,'mask.nii.gz');
ME_SMSIx(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSIx_3_70'),spectrum='scheme/exvivo_spectrum_3.mat',x0=70);
ME_SMSIx(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSIx_4_70'),spectrum='scheme/exvivo_spectrum_4.mat',x0=70);

ME_SMSIx(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSIx_3_40'),spectrum='scheme/exvivo_spectrum_3.mat',x0=40);
ME_SMSIx(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSIx_4_40'),spectrum='scheme/exvivo_spectrum_4.mat',x0=40);

ME_SMSIx(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSIx_3_15'),spectrum='scheme/exvivo_spectrum_3.mat',x0=15);
ME_SMSIx(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSIx_4_15'),spectrum='scheme/exvivo_spectrum_4.mat',x0=15);

ME_SMSIx(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSIx_3_55'),spectrum='scheme/exvivo_spectrum_3.mat',x0=55);
ME_SMSIx(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSIx_4_55'),spectrum='scheme/exvivo_spectrum_4.mat',x0=55);

ME_SMSIx(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSIx_3_35'),spectrum='scheme/exvivo_spectrum_3.mat',x0=35);
ME_SMSIx(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSIx_4_35'),spectrum='scheme/exvivo_spectrum_4.mat',x0=35);


% 
% clear all; clc;
% dataFolder = '/home/wuye/D_disk/MTE_dMRI/axon-relaxation';
% TE = [35.5 45.5];
% fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(x),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.nii.gz')),TE,'UniformOutput',false);
% fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(x),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.bval')),TE,'UniformOutput',false);
% fbvec   = arrayfun(@(x)fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(x),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.bvec')),TE,'UniformOutput',false);
% fmask   = fullfile(dataFolder,'mask.nii.gz');
% ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI_1_70'),spectrum='scheme/exvivo_spectrum_3.mat',x0=70);
% ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI_2_70'),spectrum='scheme/exvivo_spectrum_4.mat',x0=70);
% 
% ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI_1_40'),spectrum='scheme/exvivo_spectrum_3.mat',x0=40);
% ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI_2_40'),spectrum='scheme/exvivo_spectrum_4.mat',x0=40);
% 
% ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI_1_15'),spectrum='scheme/exvivo_spectrum_3.mat',x0=15);
% ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI_2_15'),spectrum='scheme/exvivo_spectrum_4.mat',x0=15);
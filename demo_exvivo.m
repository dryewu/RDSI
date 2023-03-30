% clear all; clc;
% dataFolder = '/home/wuye/D_disk/MTE_dMRI/axon-relaxation';
% addpath('/home/wuye/Softwares/REDIM');
% addpath('/usr/local/freesurfer/7.2.0/matlab');
% TE = [35.5 45.5];
% mkdir(fullfile(dataFolder,'ME_SMSI'))
% fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(x),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.nii.gz')),TE,'UniformOutput',false);
% fmask   = fullfile(dataFolder,'mask.nii');
% fbval   = fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(TE(1)),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.bval'));
% fbvec   = fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(TE(1)),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.bvec'));
% mkdir(fullfile(dataFolder,'REDIM'));
% redim_pipe(fdwi, fullfile(dataFolder,'REDIM','REDIM'), fmask, fbvec, fbval, TE);

% clear all; clc;
% dataFolder = '/home/wuye/D_disk/MTE_dMRI/axon-relaxation';
% TE = [35.5 45.5];
% mkdir(fullfile(dataFolder,'SE_SMSI'))
% for i = 1:length(TE)
%     fdwi    = fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(TE(i)),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.nii.gz'));
%     fbvec   = fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(TE(i)),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.bvec'));
%     fbval   = fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(TE(i)),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.bval'));
%     fmask   = fullfile(dataFolder,'mask.nii.gz');
%     SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'SE_SMSI',num2str(TE(i))),spectrum='scheme/exvivo_spectrum.mat');
%     SE_SI(fdwi,fbvec,fbval,fmask,fullfile(dataFolder,'SE_SI',num2str(TE(i))));
% end
% 
% clear all; clc;
% dataFolder = '/home/wuye/D_disk/MTE_dMRI/axon-relaxation';
% TE = [35.5 45.5];
% fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(x),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.nii.gz')),TE,'UniformOutput',false);
% fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(x),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.bval')),TE,'UniformOutput',false);
% fbvec   = arrayfun(@(x)fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(x),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.bvec')),TE,'UniformOutput',false);
% fmask   = fullfile(dataFolder,'mask.nii.gz');
% ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI'),spectrum='scheme/exvivo_spectrum.mat');
% ME_SMSIx(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSIx'),spectrum='scheme/exvivo_spectrum.mat');
% 
% ft2r    = fullfile(dataFolder,'ME_SMSIx','T2.nii.gz');
% ME_SIx(fdwi,fbvec,fbval,fmask,TE,fullfile(dataFolder,'ME_SIx'),spectrum='scheme/exvivo_spectrum.mat');
% 

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/axon-relaxation';
TE = [35.5 45.5];
fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(x),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(x),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.bval')),TE,'UniformOutput',false);
fbvec   = arrayfun(@(x)fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(x),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.bvec')),TE,'UniformOutput',false);
fmask   = fullfile(dataFolder,'mask.nii.gz');
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI_1_70'),spectrum='scheme/exvivo_spectrum_1.mat',x0=70);
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI_2_70'),spectrum='scheme/exvivo_spectrum_2.mat',x0=70);

ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI_1_40'),spectrum='scheme/exvivo_spectrum_1.mat',x0=40);
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI_2_40'),spectrum='scheme/exvivo_spectrum_2.mat',x0=40);

ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI_1_15'),spectrum='scheme/exvivo_spectrum_1.mat',x0=15);
ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI_2_15'),spectrum='scheme/exvivo_spectrum_2.mat',x0=15);
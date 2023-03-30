clear all; clc
addpath('/home/wuye/Softwares/REDIM')
addpath('/usr/local/freesurfer/7.2.0/matlab');

dataFolder = '/home/wuye/D_disk/MTE_dMRI/axon-relaxation';
TE = [35.5 45.5];
Input_data    = arrayfun(@(x)fullfile(dataFolder,strcat('M0683_TR_3500_TE_',num2str(x),'_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.nii.gz')),TE,'UniformOutput',false);
fn_bval = fullfile(dataFolder,'M0683_TR_3500_TE_45.5_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.bval');
fn_bvec = fullfile(dataFolder,'M0683_TR_3500_TE_45.5_d_09.6_D_17.5_dirs_96_D_G_wMax_06_axes_02.bvec');
fn_mask   = fullfile(dataFolder,'mask.nii');
Output_prefix = fullfile(dataFolder,strcat('REDIM'),'REDIM');
redim_pipe(Input_data, Output_prefix, fn_mask, fn_bvec, fn_bval, TE);


clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
TE = [75 85 95 105 115];
Input_data    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fn_bval   = fullfile(dataFolder,strcat('MTE_',num2str(TE(1)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
fn_bvec   = fullfile(dataFolder,strcat('MTE_',num2str(TE(1)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
fn_mask   = fullfile(dataFolder,'MTE_mask.nii.gz');
Output_prefix = fullfile(dataFolder,strcat('REDIM'),'REDIM');
redim_pipe(Input_data, Output_prefix, fn_mask, fn_bvec, fn_bval, TE);

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/HTY/proc';
TE = [85 95 105 115 125 135];
Input_data    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fn_bval   = fullfile(dataFolder,strcat('MTE_',num2str(TE(1)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
fn_bvec   = fullfile(dataFolder,strcat('MTE_',num2str(TE(1)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
fn_mask   = fullfile(dataFolder,'MTE_mask.nii.gz');
Output_prefix = fullfile(dataFolder,strcat('REDIM'),'REDIM');
redim_pipe(Input_data, Output_prefix, fn_mask, fn_bvec, fn_bval, TE);



clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/WDZ/proc';
TE = [75 85 95 105 115 125 135];
Input_data    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fn_bval   = fullfile(dataFolder,strcat('MTE_',num2str(TE(1)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
fn_bvec   = fullfile(dataFolder,strcat('MTE_',num2str(TE(1)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
fn_mask   = fullfile(dataFolder,'MTE_mask.nii.gz');
Output_prefix = fullfile(dataFolder,strcat('REDIM'),'REDIM');
redim_pipe(Input_data, Output_prefix, fn_mask, fn_bvec, fn_bval, TE);



clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/ZJ/proc';
TE = [75 85 95 105 115 125 135];
Input_data    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fn_bval   = fullfile(dataFolder,strcat('MTE_',num2str(TE(1)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
fn_bvec   = fullfile(dataFolder,strcat('MTE_',num2str(TE(1)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
fn_mask   = fullfile(dataFolder,'MTE_mask.nii.gz');
Output_prefix = fullfile(dataFolder,strcat('REDIM'),'REDIM');
redim_pipe(Input_data, Output_prefix, fn_mask, fn_bvec, fn_bval, TE);



clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/CW/proc';
TE = [85 95 105 115 125 135];
Input_data    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fn_bval   = fullfile(dataFolder,strcat('MTE_',num2str(TE(1)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
fn_bvec   = fullfile(dataFolder,strcat('MTE_',num2str(TE(1)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
fn_mask   = fullfile(dataFolder,'MTE_mask.nii.gz');
Output_prefix = fullfile(dataFolder,strcat('REDIM'),'REDIM');
redim_pipe(Input_data, Output_prefix, fn_mask, fn_bvec, fn_bval, TE);



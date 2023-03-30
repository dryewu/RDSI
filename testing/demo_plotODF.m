clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
F_mask   = fullfile(dataFolder,'MTE_mask.nii.gz');
F_img = fullfile(dataFolder,'ME_SI_ROI','free.nii.gz');
F_SH = arrayfun(@(x)fullfile(dataFolder,'ME_SI_ROI',strcat('FOD_hindered_',num2str(x),'.nii.gz')),[1:28],'UniformOutput',false);

ROI = [85 100; 105 120; 48 48];
plotODF(F_img,F_SH,F_mask,ROI)
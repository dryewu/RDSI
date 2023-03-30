clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
TE = [75 85 95 105 115];
mkdir(fullfile(dataFolder,'SE_SI_default'))
mkdir(fullfile(dataFolder,'SE_SI'))
for i = 1:length(TE)
    fdwi    = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz'));
    fbval   = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
    fbvec  = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
    fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
    SE_SI(fdwi,fbvec,fbval,fmask,fullfile(dataFolder,'SE_SI_default',num2str(TE(i))));
    SE_SI(fdwi,fbvec,fbval,fmask,fullfile(dataFolder,'SE_SI',num2str(TE(i))),spectrum='scheme/new_spectrum.mat');
end

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/HTY/proc';
TE = [85 95 105 115 125 135];
mkdir(fullfile(dataFolder,'SE_SI_default'))
mkdir(fullfile(dataFolder,'SE_SI'))
for i = 1:length(TE)
    fdwi    = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz'));
    fbval   = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
    fbvec  = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
    fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
    SE_SI(fdwi,fbvec,fbval,fmask,fullfile(dataFolder,'SE_SI_default',num2str(TE(i))));
    SE_SI(fdwi,fbvec,fbval,fmask,fullfile(dataFolder,'SE_SI',num2str(TE(i))),spectrum='scheme/new_spectrum.mat');
end


clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/WDZ/proc';
TE = [75 85 95 105 115 125 135];
mkdir(fullfile(dataFolder,'SE_SI_default'))
mkdir(fullfile(dataFolder,'SE_SI'))
for i = 1:length(TE)
    fdwi    = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz'));
    fbval   = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
    fbvec  = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
    fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
    SE_SI(fdwi,fbvec,fbval,fmask,fullfile(dataFolder,'SE_SI_default',num2str(TE(i))));
    SE_SI(fdwi,fbvec,fbval,fmask,fullfile(dataFolder,'SE_SI',num2str(TE(i))),spectrum='scheme/new_spectrum.mat');
end


clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/ZJ/proc';
TE = [75 85 95 105 115 125 135];
mkdir(fullfile(dataFolder,'SE_SI_default'))
mkdir(fullfile(dataFolder,'SE_SI'))
for i = 1:length(TE)
    fdwi    = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz'));
    fbval   = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
    fbvec  = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
    fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
    SE_SI(fdwi,fbvec,fbval,fmask,fullfile(dataFolder,'SE_SI_default',num2str(TE(i))));
    SE_SI(fdwi,fbvec,fbval,fmask,fullfile(dataFolder,'SE_SI',num2str(TE(i))),spectrum='scheme/new_spectrum.mat');
end


clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/CW/proc';
TE = [85 95 105 115 125 135];
mkdir(fullfile(dataFolder,'SE_SI_default'))
mkdir(fullfile(dataFolder,'SE_SI'))
for i = 1:length(TE)
    fdwi    = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz'));
    fbval   = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
    fbvec  = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
    fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
    SE_SI(fdwi,fbvec,fbval,fmask,fullfile(dataFolder,'SE_SI_default',num2str(TE(i))));
    SE_SI(fdwi,fbvec,fbval,fmask,fullfile(dataFolder,'SE_SI',num2str(TE(i))),spectrum='scheme/new_spectrum.mat');
end



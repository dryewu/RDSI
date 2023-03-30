clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/CW/proc';
TE = [85 95 105 115 125 135];
F_restricted = fullfile(dataFolder,'ME_SMSI','VF_restricted.nii.gz');
F_hindered   = fullfile(dataFolder,'ME_SMSI','VF_hindered.nii.gz');
F_isotropic  = fullfile(dataFolder,'ME_SMSI','VF_free.nii.gz'); 
T_restricted = fullfile(dataFolder,'ME_SMSI','T2_restricted.nii.gz');
T_hindered   = fullfile(dataFolder,'ME_SMSI','T2_hindered.nii.gz');
T_isotropic  = fullfile(dataFolder,'ME_SMSI','T2_free.nii.gz'); 

fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'ME_SMSI'),...
    ME=true, TE=TE, T2_restricted=T_restricted, T2_hindered=T_hindered, T2_isotropic=T_isotropic);




clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/CW/proc';
TE = [85 95 105 115 125 135];
for i = 1:length(TE)
    F_restricted = fullfile(dataFolder,'SE_SMSI',num2str(TE(i)),'VF_restricted.nii.gz');
    F_hindered   = fullfile(dataFolder,'SE_SMSI',num2str(TE(i)),'VF_hindered.nii.gz');
    F_isotropic  = fullfile(dataFolder,'SE_SMSI',num2str(TE(i)),'VF_free.nii.gz'); 
    fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
    ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'SE_SMSI',num2str(TE(i))));
end
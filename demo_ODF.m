
clear all; clc;
clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
TE = [75 85 95 105 115];
for i = 1:length(TE)
    F_restricted = fullfile(dataFolder,'SE_SMSI_ROI_default',num2str(TE(i)),'VF_restricted.nii.gz');
    F_hindered   = fullfile(dataFolder,'SE_SMSI_ROI_default',num2str(TE(i)),'VF_hindered.nii.gz');
    F_isotropic  = fullfile(dataFolder,'SE_SMSI_ROI_default',num2str(TE(i)),'VF_free.nii.gz'); 
    
    fmask   = fullfile(dataFolder,'MTE_mask_ROI.nii.gz');
    ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,fullfile(dataFolder,'SE_SMSI_ROI_default',num2str(TE(i))));

end



clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
TE = [75 85 95 105 115];
ROI = [35 125; 25 135; 48 48];
for i = 1:5
    F_SH = arrayfun(@(x)fullfile(dataFolder,'SE_SI_ROI_default',num2str(TE(i)),strcat('FOD_restricted_',num2str(x),'.nii.gz')),[1:28],'UniformOutput',false);
    F_mask   = fullfile(dataFolder,'MTE_mask_ROI.nii.gz');
    plotODFs(F_SH,F_mask,ROI)
end


clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
TE = [75 85 95 105 115];
ROI = [35 125; 25 135; 48 48];

F_SH = arrayfun(@(x)fullfile(dataFolder,'ME_SIx_ROI_default',strcat('FOD_restricted_',num2str(x),'.nii.gz')),[1:28],'UniformOutput',false);
F_mask   = fullfile(dataFolder,'MTE_mask_ROI.nii.gz');
plotODFs(F_SH,F_mask,ROI,scale=0.75)

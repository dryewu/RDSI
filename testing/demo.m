clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
TE = [75 85 95 105 115];
for i = 1:length(TE)
    fdwi    = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz'));
    fbval   = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
    fbvec   = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
    fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
    try 
        SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'SE_SMSI',num2str(TE(i)))); 
    catch 
        continue; 
    end

%     try 
%         SE_SI(fdwi,fbval,fbvec,fmask,fullfile(dataFolder,'SE_SI',num2str(TE(i))));
%     catch 
%         continue; 
%     end  
end

fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI'));
end  
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','400'),useBshell=[0,400]);
end  
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','800'),useBshell=[0,800]);
end
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','1600'),useBshell=[0,1600]);
end
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','3200'),useBshell=[0,3200]);
end
try 
    ME_SMT(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMT'));
end

%
clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/HTY/proc';
TE = [85 95 105 115 125 135];
for i = 1:length(TE)
    fdwi    = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz'));
    fbval   = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
    fbvec   = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
    fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
    try 
        SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'SE_SMSI',num2str(TE(i)))); 
    catch 
        continue; 
    end

%     try 
%         SE_SI(fdwi,fbval,fbvec,fmask,fullfile(dataFolder,'SE_SI',num2str(TE(i))));
%     catch 
%         continue; 
%     end  
end

fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI'));
end  
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','400'),useBshell=[0,400]);
end  
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','800'),useBshell=[0,800]);
end
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','1600'),useBshell=[0,1600]);
end
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','3200'),useBshell=[0,3200]);
end
try 
    ME_SMT(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMT'));
end

%
clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/WDZ/proc';
TE = [75 85 95 105 115 125 135];
for i = 1:length(TE)
    fdwi    = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz'));
    fbval   = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
    fbvec   = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
    fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
    try 
        SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'SE_SMSI',num2str(TE(i)))); 
    catch 
        continue; 
    end

%     try 
%         SE_SI(fdwi,fbval,fbvec,fmask,fullfile(dataFolder,'SE_SI',num2str(TE(i))));
%     catch 
%         continue; 
%     end  
end

fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI'));
end  
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','400'),useBshell=[0,400]);
end  
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','800'),useBshell=[0,800]);
end
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','1600'),useBshell=[0,1600]);
end
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','3200'),useBshell=[0,3200]);
end
try 
    ME_SMT(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMT'));
end

%
clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/ZJ/proc';
TE = [75 85 95 105 115 125 135];
for i = 1:length(TE)
    fdwi    = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz'));
    fbval   = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
    fbvec   = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
    fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
    try 
        SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'SE_SMSI',num2str(TE(i)))); 
    catch 
        continue; 
    end

%     try 
%         SE_SI(fdwi,fbval,fbvec,fmask,fullfile(dataFolder,'SE_SI',num2str(TE(i))));
%     catch 
%         continue; 
%     end  
end

fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI'));
end  
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','400'),useBshell=[0,400]);
end  
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','800'),useBshell=[0,800]);
end
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','1600'),useBshell=[0,1600]);
end
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','3200'),useBshell=[0,3200]);
end
try 
    ME_SMT(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMT'));
end

%
clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/CW/proc';
TE = [85 95 105 115 125 135];
for i = 1:length(TE)
    fdwi    = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz'));
    fbval   = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
    fbvec   = fullfile(dataFolder,strcat('MTE_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
    fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
    try 
        SE_SMSI(fdwi,fbval,fmask,fullfile(dataFolder,'SE_SMSI',num2str(TE(i)))); 
    catch 
        continue; 
    end

%     try 
%         SE_SI(fdwi,fbval,fbvec,fmask,fullfile(dataFolder,'SE_SI',num2str(TE(i))));
%     catch 
%         continue; 
%     end  
end

fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI'));
end  
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','400'),useBshell=[0,400]);
end  
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','800'),useBshell=[0,800]);
end
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','1600'),useBshell=[0,1600]);
end
try 
    ME_SMSI(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMSI','3200'),useBshell=[0,3200]);
end
try 
    ME_SMT(fdwi,fbval,fmask,TE,fullfile(dataFolder,'ME_SMT'));
end
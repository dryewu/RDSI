function ME_MICRO(dataFolder,TE)
    %{
    ░█▀▀█ ░█▀▀▀█ ░█▀▄▀█ ░█▀▀▀ ░█▀▀▄ ▀█▀
    ░█─── ░█──░█ ░█░█░█ ░█▀▀▀ ░█─░█ ░█─
    ░█▄▄█ ░█▄▄▄█ ░█──░█ ░█▄▄▄ ░█▄▄▀ ▄█▄

    Parameter estimation with multi-echo spectrum Imaging (ME-SMSI)


        Created by Ye Wu, PhD (dr.yewu@outlook.com)

        - Nanjing University of Science and Technology
        - University of North Carolina at Chapel Hill
        
    %}

    %% load multi-echo dMRI dataset
%     dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
%     TE = [75 85 95 105 115];
    
    fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
    fvf     = fullfile(dataFolder,'VF.nii.gz');     % estimated by ME_SMSI
    ft2r    = fullfile(dataFolder,'T2r.nii.gz');    % estimated by ME_SMSI
    fparams = fullfile(dataFolder,'params.mat');    % used in ME_SMSI

    ME_mask_info    = niftiinfo(fmask);
    ME_vf_info      = niftiinfo(fvf);
    ME_t2r_info     = niftiinfo(ft2r);
    ME_mask         = round(niftiread(ME_mask_info));
    ME_vf           = niftiread(ME_vf_info);
    ME_t2r          = niftiread(ME_t2r_info);
    params = load(fparams);

    adc_restricted      = params.adc_restricted;
    adc_hindered        = params.adc_hindered;
    adc_isotropic       = params.adc_isotropic;

    num_restricted      = params.num_restricted;  
    num_hindered        = params.num_hindered;    
    num_isotropic       = params.num_isotropic;   

    ind_restricted      = 1:num_restricted;
    ind_hindered        = 1+num_restricted:num_restricted+num_hindered;
    ind_anisotropic     = 1:num_restricted+num_hindered;
    ind_isotropic       = 1+num_restricted+num_hindered : num_restricted+num_hindered+num_isotropic;

    fodf_restricted_info= arrayfun(@(x)fullfile(dataFolder,strcat('FOD_restricted_',num2str(x),'.nii')),1:num_restricted,'UniformOutput',false);
    fodf_hindered_info  = arrayfun(@(x)fullfile(dataFolder,strcat('FOD_hindered_',num2str(x),'.nii')),1:num_hindered,'UniformOutput',false);
    
    fodf_restricted     = cellfun(@(x)niftiread(x),fodf_restricted_info,'UniformOutput',false);
    fodf_hindered       = cellfun(@(x)niftiread(x),fodf_hindered_info,'UniformOutput',false);

    lmax = 6; 
    nmax = lmax2nsh(lmax);

    %% WOT: without T2
    % Anisotropic VF
    WOT_AVF = sum(ME_vf(:,:,:,ind_anisotropic),4);

    % Intra-cellular VF 
    WOT_ICVF = sum(ME_vf(:,:,:,ind_restricted),4) ./ WOT_AVF;

    % Extra-cellular VF 
    WOT_ECVF = sum(ME_vf(:,:,:,ind_hindered),4) ./ WOT_AVF;

    % Isotropic VF
    WOT_IVF = sum(ME_vf(:,:,:,ind_isotropic),4);

    % Microscopic AD
    WOT_uAD = sum(ME_vf.*[adc_restricted(:,1); adc_hindered(:,1); adc_isotropic(:,1)]',4) ./ sum(ME_vf,4);

    % Microscopic RD
    WOT_uRD = sum(ME_vf.*[adc_restricted(:,2); adc_hindered(:,2); adc_isotropic(:,2)]',4) ./ sum(ME_vf,4);

    % Microscopic MD
    WOT_uMD = (WOT_uAD + WOT_uRD * 2) ./ 3;

    % Microscopic FA
    WOT_uFA = (WOT_uAD - WOT_uRD) ./ sqrt(WOT_uAD.^2 + 2*WOT_uRD.^2);

    % Intra-cellular AD
    WOT_uICAD = sum(ME_vf(:,:,:,ind_restricted).*adc_restricted(:,1)',4) ./ sum(ME_vf(:,:,:,ind_restricted),4);
    
    % Intra-cellular RD
    WOT_uICRD = sum(ME_vf(:,:,:,ind_restricted).*adc_restricted(:,2)',4) ./ sum(ME_vf(:,:,:,ind_restricted),4);
    
    % Extra-cellular AD
    WOT_uECAD = sum(ME_vf(:,:,:,ind_hindered).*adc_hindered(:,1)',4) ./ sum(ME_vf(:,:,:,ind_hindered),4);
    
    % Extra-cellular RD
    WOT_uECRD = sum(ME_vf(:,:,:,ind_hindered).*adc_hindered(:,2)',4) ./ sum(ME_vf(:,:,:,ind_hindered),4);
       
    % Microscopic sphericity
    WOT_uCS= WOT_uRD./WOT_uAD;

    % Microscopic linearity
    WOT_uCL= (WOT_uAD-WOT_uMD)./WOT_uMD;

    % Microscopic plane
    WOT_uCP = (WOT_uMD-WOT_uRD)./WOT_uAD;

    % diffusion-relaxation exchange
    WT_R2EX = abs(1./(ME_t2r(:,:,:,1)) - 1./(ME_t2r(:,:,:,2)));

    % restricted isotropic diffusion 
    % WOT_RNO = cellfun(@(x)x(:,:,:,1)fodf_restricted

    % restricted directional diffusion
    % WOT_RND = 
    
    % hindered isotropic diffusion 
    % WOT_HNO = cellfun(@(x)x(:,:,:,1)fodf_restricted

    % hindered directional diffusion
    % WOT_HND = 


    info = ME_vf_info;
    info.ImageSize = size(ME_vf,[1,2,3]);

    niftiwrite(WOT_AVF.*ME_mask,fullfile(dataFolder,'WOT_AVF.nii'),info,'Compressed', true);
    niftiwrite(WOT_ICVF.*ME_mask,fullfile(dataFolder,'WOT_ICVF.nii'),info,'Compressed', true);
    niftiwrite(WOT_ECVF.*ME_mask,fullfile(dataFolder,'WOT_ECVF.nii'),info,'Compressed', true);
    niftiwrite(WOT_IVF.*ME_mask,fullfile(dataFolder,'WOT_IVF.nii'),info,'Compressed', true);
    niftiwrite(WOT_uAD.*ME_mask,fullfile(dataFolder,'WOT_uAD.nii'),info,'Compressed', true);
    niftiwrite(WOT_uMD.*ME_mask,fullfile(dataFolder,'WOT_uMD.nii'),info,'Compressed', true);
    niftiwrite(WOT_uFA.*ME_mask,fullfile(dataFolder,'WOT_uFA.nii'),info,'Compressed', true);
    niftiwrite(WOT_uICAD.*ME_mask,fullfile(dataFolder,'WOT_uICAD.nii'),info,'Compressed', true);
    niftiwrite(WOT_uECAD.*ME_mask,fullfile(dataFolder,'WOT_uECAD.nii'),info,'Compressed', true);
    niftiwrite(WOT_uICRD.*ME_mask,fullfile(dataFolder,'WOT_uICRD.nii'),info,'Compressed', true);
    niftiwrite(WOT_uECRD.*ME_mask,fullfile(dataFolder,'WOT_uECRD.nii'),info,'Compressed', true);
    niftiwrite(WOT_uCS.*ME_mask,fullfile(dataFolder,'WOT_uCS.nii'),info,'Compressed', true);
    niftiwrite(WOT_uCL.*ME_mask,fullfile(dataFolder,'WOT_uCL.nii'),info,'Compressed', true);
    niftiwrite(WOT_uCP.*ME_mask,fullfile(dataFolder,'WOT_uCP.nii'),info,'Compressed', true);
    niftiwrite(WT_R2EX.*ME_mask,fullfile(dataFolder,'WT_R2EX.nii'),info,'Compressed', true);
    niftiwrite(WOT_RNO.*ME_mask,fullfile(dataFolder,'WOT_RNO.nii'),info,'Compressed', true);
    niftiwrite(WOT_RND.*ME_mask,fullfile(dataFolder,'WOT_RND.nii'),info,'Compressed', true);
    niftiwrite(WOT_HNO.*ME_mask,fullfile(dataFolder,'WOT_HNO.nii'),info,'Compressed', true);
    niftiwrite(WOT_HND.*ME_mask,fullfile(dataFolder,'WOT_HND.nii'),info,'Compressed', true);

end


function nsh = lmax2nsh(lmax)
    nsh = (lmax+1) * (lmax+2) / 2;
end

























    
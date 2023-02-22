function ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,outpath,options)
    %{
    ░█▀▀█ ░█▀▀▀█ ░█▀▄▀█ ░█▀▀▀ ░█▀▀▄ ▀█▀
    ░█─── ░█──░█ ░█░█░█ ░█▀▀▀ ░█─░█ ░█─
    ░█▄▄█ ░█▄▄▄█ ░█──░█ ░█▄▄▄ ░█▄▄▀ ▄█▄

    Microstructure parameter estimation


        Created by Ye Wu, PhD (dr.yewu@outlook.com)

        - Nanjing University of Science and Technology
        - University of North Carolina at Chapel Hill
        
    %}

    arguments
        F_restricted            string   {mustBeNonzeroLengthText} = false  % [x,y,z,num_restricted]
        F_hindered              string   {mustBeNonzeroLengthText} = false  % [x,y,z,num_indered]
        F_isotropic             string   {mustBeNonzeroLengthText} = false  % [x,y,z,num_isotropic]
        fmask                   string   {mustBeFile}              = false           
        outpath                 string   {mustBeNonzeroLengthText} = false

        options.TE              (1,:)    {mustBeNumeric}           = false
        options.T2_restricted   string   {mustBeNonzeroLengthText} = false  % [x,y,z,1]
        options.T2_hindered     string   {mustBeNonzeroLengthText} = false  % [x,y,z,1]
        options.T2_isotropic    string   {mustBeNonzeroLengthText} = false  % [x,y,z,1]
        options.spectrum        string   {mustBeFolder} = 'scheme/default_spectrum.mat'
    end

    ME_F_restricted_info   = niftiinfo(F_restricted);
    ME_F_hindered_info     = niftiinfo(F_hindered);
    ME_F_isotropic_info    = niftiinfo(F_isotropic);
    ME_mask_info           = niftiinfo(fmask);

    ME_F_restricted        = niftiread(ME_F_restricted_info);
    ME_F_hindered          = niftiread(ME_F_hindered_info);
    ME_F_isotropic         = niftiread(ME_F_isotropic_info);
    ME_mask                = round(niftiread(ME_mask_info));


    if options.T2_restricted
        ME_T2_restricted_info   = niftiinfo(T2_restricted);
        ME_T2_restricted        = niftiread(ME_T2_restricted_info);
    end

    if options.T2_hindered
        ME_T2_hindered_info     = niftiinfo(T2_hindered);
        ME_T2_hindered          = niftiread(ME_T2_hindered_info);
    end

    if options.T2_isotropic
        ME_T2_isotropic_info    = niftiinfo(T2_isotropic);
        ME_T2_isotropic         = niftiread(ME_T2_isotropic_info);
    end

    default_spectrum = load(options.spectrum);

    adc_restricted   = default_spectrum.adc_restricted;
    adc_hindered     = default_spectrum.adc_hindered;
    adc_isotropic    = default_spectrum.adc_isotropic;

    num_restricted  = size(adc_restricted,1);  
    num_hindered    = size(adc_hindered,1);    
    num_isotropic   = size(adc_isotropic,1);  


    %% WOT: without T2
    % Anisotropic VF
    WOT_AVF = sum(ME_F_restricted+ME_F_hindered,4);

    % Intra-cellular VF 
    WOT_ICVF = sum(ME_F_restricted,4) ./ WOT_AVF;

    % Extra-cellular VF 
    WOT_ECVF = sum(ME_F_hindered,4) ./ WOT_AVF;

    % Isotropic VF
    WOT_IVF = sum(ME_F_isotropic,4);

    % Microscopic AD
    ME_vf   = cat(4,ME_F_restricted,ME_F_hindered,ME_F_isotropic);
    WOT_uAD = sum(ME_vf.*[adc_restricted(:,1); adc_hindered(:,1); adc_isotropic(:,1)]',4) ./ sum(ME_vf,4);

    % Microscopic RD
    WOT_uRD = sum(ME_vf.*[adc_restricted(:,2); adc_hindered(:,2); adc_isotropic(:,2)]',4) ./ sum(ME_vf,4);

    % Microscopic MD
    WOT_uMD = (WOT_uAD + WOT_uRD * 2) ./ 3;

    % Microscopic FA
    WOT_uFA = (WOT_uAD - WOT_uRD) ./ sqrt(WOT_uAD.^2 + 2*WOT_uRD.^2);

    % Intra-cellular AD
    WOT_uICAD = sum(ME_F_restricted.*adc_restricted(:,1)',4) ./ sum(ME_F_restricted,4);
    
    % Intra-cellular RD
    WOT_uICRD = sum(ME_F_restricted.*adc_restricted(:,2)',4) ./ sum(ME_F_restricted,4);
    
    % Extra-cellular AD
    WOT_uECAD = sum(ME_F_hindered.*adc_hindered(:,1)',4) ./ sum(ME_F_hindered,4);
    
    % Extra-cellular RD
    WOT_uECRD = sum(ME_F_hindered.*adc_hindered(:,2)',4) ./ sum(ME_F_hindered,4);
       
    % Microscopic sphericity
    WOT_uCS= WOT_uRD./WOT_uAD;

    % Microscopic linearity
    WOT_uCL= (WOT_uAD-WOT_uMD)./WOT_uMD;

    % Microscopic plane
    WOT_uCP = (WOT_uMD-WOT_uRD)./WOT_uAD;


    if ~exist(outpath,'folder')
        mkdir('outpath')
    end
    info = ME_F_restricted_info;
    info.ImageSize = size(ME_F_restricted_info,[1,2,3]);

    niftiwrite(WOT_AVF.*ME_mask,fullfile(outpath,'WOT_AVF.nii'),info,'Compressed', true);
    niftiwrite(WOT_ICVF.*ME_mask,fullfile(outpath,'WOT_ICVF.nii'),info,'Compressed', true);
    niftiwrite(WOT_ECVF.*ME_mask,fullfile(outpath,'WOT_ECVF.nii'),info,'Compressed', true);
    niftiwrite(WOT_IVF.*ME_mask,fullfile(outpath,'WOT_IVF.nii'),info,'Compressed', true);
    niftiwrite(WOT_uAD.*ME_mask,fullfile(outpath,'WOT_uAD.nii'),info,'Compressed', true);
    niftiwrite(WOT_uMD.*ME_mask,fullfile(outpath,'WOT_uMD.nii'),info,'Compressed', true);
    niftiwrite(WOT_uFA.*ME_mask,fullfile(outpath,'WOT_uFA.nii'),info,'Compressed', true);
    niftiwrite(WOT_uICAD.*ME_mask,fullfile(outpath,'WOT_uICAD.nii'),info,'Compressed', true);
    niftiwrite(WOT_uECAD.*ME_mask,fullfile(outpath,'WOT_uECAD.nii'),info,'Compressed', true);
    niftiwrite(WOT_uICRD.*ME_mask,fullfile(outpath,'WOT_uICRD.nii'),info,'Compressed', true);
    niftiwrite(WOT_uECRD.*ME_mask,fullfile(outpath,'WOT_uECRD.nii'),info,'Compressed', true);
    niftiwrite(WOT_uCS.*ME_mask,fullfile(outpath,'WOT_uCS.nii'),info,'Compressed', true);
    niftiwrite(WOT_uCL.*ME_mask,fullfile(outpath,'WOT_uCL.nii'),info,'Compressed', true);
    niftiwrite(WOT_uCP.*ME_mask,fullfile(outpath,'WOT_uCP.nii'),info,'Compressed', true);
    niftiwrite(WT_R2EX.*ME_mask,fullfile(outpath,'WT_R2EX.nii'),info,'Compressed', true);
    niftiwrite(WOT_RNO.*ME_mask,fullfile(outpath,'WOT_RNO.nii'),info,'Compressed', true);
    niftiwrite(WOT_RND.*ME_mask,fullfile(outpath,'WOT_RND.nii'),info,'Compressed', true);
    niftiwrite(WOT_HNO.*ME_mask,fullfile(outpath,'WOT_HNO.nii'),info,'Compressed', true);
    niftiwrite(WOT_HND.*ME_mask,fullfile(outpath,'WOT_HND.nii'),info,'Compressed', true);

    if options.T2_restricted && options.T2_hindered && options.T2_hindered 

        % relaxation rate on restricted
        WT_R2r = 1./(ME_T2_restricted);

        % relaxation rate on hindered
        WT_R2h = 1./(ME_T2_hindered);

        % relaxation rate on isotropic
        WT_R2f = 1./(ME_T2_hindered);

        % relaxation rate on exchange between restricted and hindered
        WT_E2rh = abs(1./(ME_T2_restricted)-(1./(ME_T2_hindered));

        % relaxation rate on exchange between restricted and isotropic
        WT_E2rf = abs(1./(ME_T2_restricted)-(1./(ME_T2_isotropic));

        % relaxation rate on exchange between hindered and isotropic
        WT_E2hf = abs(1./(ME_T2_hindered)-(1./(ME_T2_isotropic));

        niftiwrite(WT_R2r.*ME_mask,fullfile(outpath,'WT_R2r.nii'),info,'Compressed', true);
        niftiwrite(WT_R2h.*ME_mask,fullfile(outpath,'WT_R2h.nii'),info,'Compressed', true);
        niftiwrite(WT_R2f.*ME_mask,fullfile(outpath,'WT_R2f.nii'),info,'Compressed', true);
        niftiwrite(WT_E2rh.*ME_mask,fullfile(outpath,'WT_E2rh.nii'),info,'Compressed', true);
        niftiwrite(WT_E2rf.*ME_mask,fullfile(outpath,'WT_E2rf.nii'),info,'Compressed', true);
        niftiwrite(WT_E2hf.*ME_mask,fullfile(outpath,'WT_E2hf.nii'),info,'Compressed', true);


        ME_F_restricted = ME_F_restricted .* ME_T2_restricted;
        ME_F_hindered   = ME_F_hindered .* ME_T2_hindered;
        ME_F_isotropic  = ME_F_isotropic .* ME_T2_isotropic;
    
        %% WT: with T2
        % Anisotropic VF
        WT_AVF = sum(ME_F_restricted+ME_F_hindered,4);
    
        % Intra-cellular VF 
        WT_ICVF = sum(ME_F_restricted,4) ./ WT_AVF;
    
        % Extra-cellular VF 
        WT_ECVF = sum(ME_F_hindered,4) ./ WT_AVF;
    
        % Isotropic VF
        WT_IVF = sum(ME_F_isotropic,4);
    
        % Microscopic AD
        WT_uAD = sum(ME_vf.*[adc_restricted(:,1); adc_hindered(:,1); adc_isotropic(:,1)]',4) ./ sum(ME_vf,4);
    
        % Microscopic RD
        WT_uRD = sum(ME_vf.*[adc_restricted(:,2); adc_hindered(:,2); adc_isotropic(:,2)]',4) ./ sum(ME_vf,4);
    
        % Microscopic MD
        WT_uMD = (WT_uAD + WT_uRD * 2) ./ 3;
    
        % Microscopic FA
        WT_uFA = (WT_uAD - WT_uRD) ./ sqrt(WT_uAD.^2 + 2*WT_uRD.^2);
    
        % Intra-cellular AD
        WT_uICAD = sum(ME_F_restricted.*adc_restricted(:,1)',4) ./ sum(ME_F_restricted,4);
        
        % Intra-cellular RD
        WT_uICRD = sum(ME_F_restricted.*adc_restricted(:,2)',4) ./ sum(ME_F_restricted,4);
        
        % Extra-cellular AD
        WT_uECAD = sum(ME_F_hindered.*adc_hindered(:,1)',4) ./ sum(ME_F_hindered,4);
        
        % Extra-cellular RD
        WT_uECRD = sum(ME_F_hindered.*adc_hindered(:,2)',4) ./ sum(ME_F_hindered,4);
           
        % Microscopic sphericity
        WT_uCS= WT_uRD./WT_uAD;
    
        % Microscopic linearity
        WT_uCL= (WT_uAD-WT_uMD)./WT_uMD;
    
        % Microscopic plane
        WT_uCP = (WT_uMD-WT_uRD)./WT_uAD;
    
        niftiwrite(WT_AVF.*ME_mask,fullfile(outpath,'WT_AVF.nii'),info,'Compressed', true);
        niftiwrite(WT_ICVF.*ME_mask,fullfile(outpath,'WT_ICVF.nii'),info,'Compressed', true);
        niftiwrite(WT_ECVF.*ME_mask,fullfile(outpath,'WT_ECVF.nii'),info,'Compressed', true);
        niftiwrite(WT_IVF.*ME_mask,fullfile(outpath,'WT_IVF.nii'),info,'Compressed', true);
        niftiwrite(WT_uAD.*ME_mask,fullfile(outpath,'WT_uAD.nii'),info,'Compressed', true);
        niftiwrite(WT_uMD.*ME_mask,fullfile(outpath,'WT_uMD.nii'),info,'Compressed', true);
        niftiwrite(WT_uFA.*ME_mask,fullfile(outpath,'WT_uFA.nii'),info,'Compressed', true);
        niftiwrite(WT_uICAD.*ME_mask,fullfile(outpath,'WT_uICAD.nii'),info,'Compressed', true);
        niftiwrite(WT_uECAD.*ME_mask,fullfile(outpath,'WT_uECAD.nii'),info,'Compressed', true);
        niftiwrite(WT_uICRD.*ME_mask,fullfile(outpath,'WT_uICRD.nii'),info,'Compressed', true);
        niftiwrite(WT_uECRD.*ME_mask,fullfile(outpath,'WT_uECRD.nii'),info,'Compressed', true);
        niftiwrite(WT_uCS.*ME_mask,fullfile(outpath,'WT_uCS.nii'),info,'Compressed', true);
        niftiwrite(WT_uCL.*ME_mask,fullfile(outpath,'WT_uCL.nii'),info,'Compressed', true);
        niftiwrite(WT_uCP.*ME_mask,fullfile(outpath,'WT_uCP.nii'),info,'Compressed', true);
        niftiwrite(WT_R2EX.*ME_mask,fullfile(outpath,'WT_R2EX.nii'),info,'Compressed', true);
        niftiwrite(WT_RNO.*ME_mask,fullfile(outpath,'WT_RNO.nii'),info,'Compressed', true);
        niftiwrite(WT_RND.*ME_mask,fullfile(outpath,'WT_RND.nii'),info,'Compressed', true);
        niftiwrite(WT_HNO.*ME_mask,fullfile(outpath,'WT_HNO.nii'),info,'Compressed', true);
        niftiwrite(WT_HND.*ME_mask,fullfile(outpath,'WT_HND.nii'),info,'Compressed', true);

    end
end


















    
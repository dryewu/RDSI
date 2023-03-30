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
        
        options.ME              (1,1)    {mustBeNumericOrLogical}  = false
        options.TE              (1,:)    {mustBeVector}            = false
        options.T2_restricted   string   {mustBeNonzeroLengthText} = false  % [x,y,z,1]
        options.T2_hindered     string   {mustBeNonzeroLengthText} = false  % [x,y,z,1]
        options.T2_isotropic    string   {mustBeNonzeroLengthText} = false  % [x,y,z,1]
        options.spectrum        string   {mustBeFile} = 'scheme/default_spectrum.mat'
    end

    ME_F_restricted_info   = niftiinfo(F_restricted);
    ME_F_hindered_info     = niftiinfo(F_hindered);
    ME_F_isotropic_info    = niftiinfo(F_isotropic);
    ME_mask_info           = niftiinfo(fmask);

    ME_F_restricted        = niftiread(ME_F_restricted_info);
    ME_F_hindered          = niftiread(ME_F_hindered_info);
    ME_F_isotropic         = niftiread(ME_F_isotropic_info);
    ME_mask                = single(round(niftiread(ME_mask_info)));


    if options.ME
        ME_T2_restricted_info   = niftiinfo(options.T2_restricted);
        ME_T2_restricted        = niftiread(ME_T2_restricted_info);

        ME_T2_hindered_info     = niftiinfo(options.T2_hindered);
        ME_T2_hindered          = niftiread(ME_T2_hindered_info);

        ME_T2_isotropic_info    = niftiinfo(options.T2_isotropic);
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
    WOT_AVF = sum(ME_F_restricted,4)+sum(ME_F_hindered,4);

    % Intra-cellular VF 
    WOT_ICVF = sum(ME_F_restricted,4) ./ WOT_AVF;

    % Extra-cellular VF 
    WOT_ECVF = sum(ME_F_hindered,4) ./ WOT_AVF;

    % Isotropic VF
    WOT_IVF = sum(ME_F_isotropic,4);

    % Microscopic AD
    ME_vf   = cat(4,ME_F_restricted,ME_F_hindered,ME_F_isotropic);
    WOT_uAD = sum(ME_vf.*reshape([adc_restricted(:,1); adc_hindered(:,1); adc_isotropic(:,1)],1,1,1,num_restricted+num_hindered+num_isotropic),4) ./ sum(ME_vf,4);

    % Microscopic RD
    WOT_uRD = sum(ME_vf.*reshape([adc_restricted(:,2); adc_hindered(:,2); adc_isotropic(:,1)],1,1,1,num_restricted+num_hindered+num_isotropic),4) ./ sum(ME_vf,4);

    % Microscopic MD
    WOT_uMD = (WOT_uAD + WOT_uRD * 2) ./ 3;

    % Microscopic FA
    WOT_uFA = (WOT_uAD - WOT_uRD) ./ sqrt(WOT_uAD.^2 + 2*WOT_uRD.^2);

    % Intra-cellular AD
    WOT_uICAD = sum(ME_F_restricted.*reshape(adc_restricted(:,1),1,1,1,num_restricted),4) ./ sum(ME_F_restricted,4);
    
    % Intra-cellular RD
    WOT_uICRD = sum(ME_F_restricted.*reshape(adc_restricted(:,2),1,1,1,num_restricted),4) ./ sum(ME_F_restricted,4);
    
    % Extra-cellular AD
    WOT_uECAD = sum(ME_F_hindered.*reshape(adc_hindered(:,1),1,1,1,num_hindered),4) ./ sum(ME_F_hindered,4);
    
    % Extra-cellular RD
    WOT_uECRD = sum(ME_F_hindered.*reshape(adc_hindered(:,2),1,1,1,num_hindered),4) ./ sum(ME_F_hindered,4);
       
    % Microscopic sphericity
    WOT_uCS= WOT_uRD./WOT_uAD;

    % Microscopic linearity
    WOT_uCL= (WOT_uAD-WOT_uMD)./WOT_uMD;

    % Microscopic plane
    WOT_uCP = (WOT_uMD-WOT_uRD)./WOT_uAD;
    
    % mean effective neurite radius 
    WOT_NRr = mean((ME_F_restricted.*reshape(adc_restricted(:,1).*adc_restricted(:,2),1,1,1,num_restricted)).^0.25,4);
    WOT_NRh = mean((ME_F_hindered.*reshape(adc_hindered(:,1).*adc_hindered(:,2),1,1,1,num_hindered)).^0.25,4);

    % std effective neurite radius 
    WOT_SNRr = std((ME_F_restricted.*reshape(adc_restricted(:,1).*adc_restricted(:,2),1,1,1,num_restricted)).^0.25,[],4);
    WOT_SNRh = std((ME_F_hindered.*reshape(adc_hindered(:,1).*adc_hindered(:,2),1,1,1,num_hindered)).^0.25,[],4);

    % relative effective neurite radius
    WOT_RNRr = WOT_SNRr./WOT_NRr;
    WOT_RNRh = WOT_SNRh./WOT_NRh;


    if ~exist(outpath,'dir')
        mkdir('outpath')
    end
    info = ME_F_restricted_info;
    info.ImageSize = size(ME_F_restricted,[1,2,3]);
    info.PixelDimensions = info.PixelDimensions(1:3);
    niftiwrite(single(WOT_AVF.*ME_mask),fullfile(outpath,'WOT_AVF.nii'),info,'Compressed', true);
    niftiwrite(single(WOT_ICVF.*ME_mask),fullfile(outpath,'WOT_ICVF.nii'),info,'Compressed', true);
    niftiwrite(single(WOT_ECVF.*ME_mask),fullfile(outpath,'WOT_ECVF.nii'),info,'Compressed', true);
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
    niftiwrite(WOT_NRr.*ME_mask,fullfile(outpath,'WOT_NRr.nii'),info,'Compressed', true);
    niftiwrite(WOT_NRh.*ME_mask,fullfile(outpath,'WOT_NRh.nii'),info,'Compressed', true);
    niftiwrite(WOT_SNRr.*ME_mask,fullfile(outpath,'WOT_SNRr.nii'),info,'Compressed', true);
    niftiwrite(WOT_SNRh.*ME_mask,fullfile(outpath,'WOT_SNRh.nii'),info,'Compressed', true);
    niftiwrite(WOT_RNRr.*ME_mask,fullfile(outpath,'WOT_RNRr.nii'),info,'Compressed', true);
    niftiwrite(WOT_RNRh.*ME_mask,fullfile(outpath,'WOT_RNRh.nii'),info,'Compressed', true);


%     niftiwrite(WT_R2EX.*ME_mask,fullfile(outpath,'WT_R2EX.nii'),info,'Compressed', true);
%     niftiwrite(WOT_RNO.*ME_mask,fullfile(outpath,'WOT_RNO.nii'),info,'Compressed', true);
%     niftiwrite(WOT_RND.*ME_mask,fullfile(outpath,'WOT_RND.nii'),info,'Compressed', true);
%     niftiwrite(WOT_HNO.*ME_mask,fullfile(outpath,'WOT_HNO.nii'),info,'Compressed', true);
%     niftiwrite(WOT_HND.*ME_mask,fullfile(outpath,'WOT_HND.nii'),info,'Compressed', true);
    
    if options.ME
        % relaxation rate on restricted
        WT_R2r = 1./(ME_T2_restricted);

        % relaxation rate on hindered
        WT_R2h = 1./(ME_T2_hindered);

        % relaxation rate on isotropic
        WT_R2f = 1./(ME_T2_isotropic);

        % relaxation rate on exchange between restricted and hindered
        WT_E2rh = abs(1./(ME_T2_restricted)-(1./(ME_T2_hindered)));

        % relaxation rate on exchange between restricted and isotropic
        WT_E2rf = abs(1./(ME_T2_restricted)-(1./(ME_T2_isotropic)));

        % relaxation rate on exchange between hindered and isotropic
        WT_E2hf = abs(1./(ME_T2_hindered)-(1./(ME_T2_isotropic)));

        niftiwrite(single(WT_R2r.*ME_mask),fullfile(outpath,'WT_R2r.nii'),info,'Compressed', true);
        niftiwrite(single(WT_R2h.*ME_mask),fullfile(outpath,'WT_R2h.nii'),info,'Compressed', true);
        niftiwrite(single(WT_R2f.*ME_mask),fullfile(outpath,'WT_R2f.nii'),info,'Compressed', true);
        niftiwrite(single(WT_E2rh.*ME_mask),fullfile(outpath,'WT_E2rh.nii'),info,'Compressed', true);
        niftiwrite(single(WT_E2rf.*ME_mask),fullfile(outpath,'WT_E2rf.nii'),info,'Compressed', true);
        niftiwrite(single(WT_E2hf.*ME_mask),fullfile(outpath,'WT_E2hf.nii'),info,'Compressed', true);
    end
end

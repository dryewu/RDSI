function ME_dMICRO(F_restricted,F_hindered,F_isotropic,fmask,outpath,options)
    %{
    ░█▀▀█ ░█▀▀▀█ ░█▀▄▀█ ░█▀▀▀ ░█▀▀▄ ▀█▀
    ░█─── ░█──░█ ░█░█░█ ░█▀▀▀ ░█─░█ ░█─
    ░█▄▄█ ░█▄▄▄█ ░█──░█ ░█▄▄▄ ░█▄▄▀ ▄█▄

    Directional microstructure parameter estimation


        Created by Ye Wu, PhD (dr.yewu@outlook.com)

        - Nanjing University of Science and Technology, China
        - University of North Carolina at Chapel Hill, USA
        
    %}

    arguments
        F_restricted            (1,:)    {mustBeNonzeroLengthText} = false  % [x,y,z,num_restricted*nsh]
        F_hindered              (1,:)    {mustBeNonzeroLengthText} = false  % [x,y,z,num_indered*nsh]
        F_isotropic             string   {mustBeNonzeroLengthText} = false  % [x,y,z,num_isotropic*nsh]
        fmask                   string   {mustBeFile}              = false           
        outpath                 string   {mustBeNonzeroLengthText} = false

        options.TE              (1,:)    {mustBeNumeric}           = false
        options.T2_restricted   string   {mustBeNonzeroLengthText} = false  % [x,y,z,1]
        options.T2_hindered     string   {mustBeNonzeroLengthText} = false  % [x,y,z,1]
        options.T2_isotropic    string   {mustBeNonzeroLengthText} = false  % [x,y,z,1]
        options.spectrum        string   {mustBeFolder} = 'scheme/default_spectrum.mat'
    end

    ME_F_restricted_info   = cellfun(@(x)niftiinfo(x),F_restricted,'UniformOutput',false);
    ME_F_hindered_info     = cellfun(@(x)niftiinfo(x),F_hindered,'UniformOutput',false);
    ME_F_isotropic_info    = cellfun(@(x)niftiinfo(x),F_isotropic,'UniformOutput',false);
    ME_mask_info           = niftiinfo(fmask);

    ME_F_restricted        = cellfun(@(x)niftiread(x),ME_F_restricted_info,'UniformOutput',false);
    ME_F_hindered          = cellfun(@(x)niftiread(x),ME_F_hindered_info,'UniformOutput',false);
    ME_F_isotropic         = cellfun(@(x)niftiread(x),ME_F_isotropic_info,'UniformOutput',false);
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
    WOT_AVF = sum(cat(5,ME_F_restricted{:}),5) + sum(cat(5,ME_F_hindered{:}),5);

    % Intra-cellular VF 
    WOT_ICVF = sum(cat(5,ME_F_restricted{:}),5) ./ WOT_AVF;

    % Extra-cellular VF 
    WOT_ECVF = sum(cat(5,ME_F_hindered{:}),5) ./ WOT_AVF;

    % Isotropic VF
    WOT_IVF = sum(ME_F_isotropic,4);

    % Microscopic AD

    temp1_restricted = cellfun(@(x,y) x.*y,ME_F_restricted,adc_restricted(:,1), 'UniformOutput',false);
    temp1_hindered = cellfun(@(x,y) x.*y,ME_F_hindered,adc_hindered(:,1), 'UniformOutput',false);
    temp1_isotropic = cellfun(@(x,y) x.*y,ME_F_isotropic,adc_isotropic(:,1), 'UniformOutput',false);

    temp2_restricted = cellfun(@(x,y) x.*y,ME_F_restricted,adc_restricted(:,2), 'UniformOutput',false);
    temp2_hindered = cellfun(@(x,y) x.*y,ME_F_hindered,adc_hindered(:,2), 'UniformOutput',false);
    temp2_isotropic = cellfun(@(x,y) x.*y,ME_F_isotropic,adc_isotropic(:,2), 'UniformOutput',false);

    WOT_uAD = (sum(cat(5,temp1_restricted{:}),5) + sum(cat(5,temp1_hindered{:}),5) + sum(cat(5,temp1_isotropic{:}),5)) ./ (WOT_AVF + WOT_IVF);

    % Microscopic RD
    WOT_uRD = (sum(cat(5,temp2_restricted{:}),5) + sum(cat(5,temp2_hindered{:}),5) + sum(cat(5,temp2_isotropic{:}),5))./ (WOT_AVF + WOT_IVF);

    % Microscopic MD
    WOT_uMD = (WOT_uAD + WOT_uRD * 2) ./ 3;

    % Microscopic FA
    WOT_uFA = (WOT_uAD - WOT_uRD) ./ sqrt(WOT_uAD.^2 + 2*WOT_uRD.^2);

    % Intra-cellular AD
    WOT_uICAD = sum(cat(5,temp1_restricted{:}),5) ./ sum(cat(5,ME_F_restricted{:}),5);
    
    % Intra-cellular RD
    WOT_uICRD = sum(cat(5,temp2_restricted{:}),5) ./ sum(cat(5,ME_F_restricted{:}),5);
    
    % Extra-cellular AD
    WOT_uECAD = sum(cat(5,temp1_hindered{:}),5) ./ sum(cat(5,ME_F_hindered{:}),5);
    
    % Extra-cellular RD
    WOT_uECRD = sum(cat(5,temp2_hindered{:}),5) ./ sum(cat(5,ME_F_hindered{:}),5);
       
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

    niftiwrite(WOT_AVF.*ME_mask,fullfile(outpath,'dWOT_AVF.nii'),info,'Compressed', true);
    niftiwrite(WOT_ICVF.*ME_mask,fullfile(outpath,'dWOT_ICVF.nii'),info,'Compressed', true);
    niftiwrite(WOT_ECVF.*ME_mask,fullfile(outpath,'dWOT_ECVF.nii'),info,'Compressed', true);
    niftiwrite(WOT_IVF.*ME_mask,fullfile(outpath,'dWOT_IVF.nii'),info,'Compressed', true);
    niftiwrite(WOT_uAD.*ME_mask,fullfile(outpath,'dWOT_uAD.nii'),info,'Compressed', true);
    niftiwrite(WOT_uMD.*ME_mask,fullfile(outpath,'dWOT_uMD.nii'),info,'Compressed', true);
    niftiwrite(WOT_uFA.*ME_mask,fullfile(outpath,'dWOT_uFA.nii'),info,'Compressed', true);
    niftiwrite(WOT_uICAD.*ME_mask,fullfile(outpath,'dWOT_uICAD.nii'),info,'Compressed', true);
    niftiwrite(WOT_uECAD.*ME_mask,fullfile(outpath,'dWOT_uECAD.nii'),info,'Compressed', true);
    niftiwrite(WOT_uICRD.*ME_mask,fullfile(outpath,'dWOT_uICRD.nii'),info,'Compressed', true);
    niftiwrite(WOT_uECRD.*ME_mask,fullfile(outpath,'dWOT_uECRD.nii'),info,'Compressed', true);
    niftiwrite(WOT_uCS.*ME_mask,fullfile(outpath,'dWOT_uCS.nii'),info,'Compressed', true);
    niftiwrite(WOT_uCL.*ME_mask,fullfile(outpath,'dWOT_uCL.nii'),info,'Compressed', true);
    niftiwrite(WOT_uCP.*ME_mask,fullfile(outpath,'dWOT_uCP.nii'),info,'Compressed', true);
    niftiwrite(WT_R2EX.*ME_mask,fullfile(outpath,'dWT_R2EX.nii'),info,'Compressed', true);
    niftiwrite(WOT_RNO.*ME_mask,fullfile(outpath,'dWOT_RNO.nii'),info,'Compressed', true);
    niftiwrite(WOT_RND.*ME_mask,fullfile(outpath,'dWOT_RND.nii'),info,'Compressed', true);
    niftiwrite(WOT_HNO.*ME_mask,fullfile(outpath,'dWOT_HNO.nii'),info,'Compressed', true);
    niftiwrite(WOT_HND.*ME_mask,fullfile(outpath,'dWOT_HND.nii'),info,'Compressed', true);

end


function nsh = lmax2nsh(lmax)
    nsh = (lmax+1) * (lmax+2) / 2;
end

















    
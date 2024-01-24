function ME_MICRO(F_restricted,F_hindered,F_isotropic,fmask,outpath,options)
    %{
    ░█▀▀█ ░█▀▀▀█ ░█▀▄▀█ ░█▀▀▀ ░█▀▀▄ ▀█▀
    ░█─── ░█──░█ ░█░█░█ ░█▀▀▀ ░█─░█ ░█─
    ░█▄▄█ ░█▄▄▄█ ░█──░█ ░█▄▄▄ ░█▄▄▀ ▄█▄

    Microstructure parameter estimation


        Created by Ye Wu, PhD (dr.yewu@outlook.com)

        - Nanjing University of Science and Technology, China
        - University of North Carolina at Chapel Hill, USA
        
    %}

    arguments
        F_restricted            string   {mustBeNonzeroLengthText} = false  % [x,y,z,num_restricted]
        F_hindered              string   {mustBeNonzeroLengthText} = false  % [x,y,z,num_indered]
        F_isotropic             string   {mustBeNonzeroLengthText} = false  % [x,y,z,num_isotropic]
        fmask                   string   {mustBeFile}              = false           
        outpath                 string   {mustBeNonzeroLengthText} = false

        options.TE              (1,:)    {mustBeVector}            = false
        options.F_T2            string   {mustBeNonzeroLengthText} = false  % [x,y,z,nb*3]

        options.spectrum        string   {mustBeFile} = which('default_spectrum.mat')
        options.index           (1,:)    {mustBeNonzeroLengthText} = {'WOT_AVF','WOT_ICVF','WOT_ECVF','WOT_IVF','WOT_uAD','WOT_uRD',...
                                                                        'WOT_uMD','WOT_uFA','WOT_uICAD','WOT_uICRD','WOT_uECAD','WOT_uECRD',...
                                                                        'WOT_uCS','WOT_uCL','WOT_uCP','WOT_NRr','WOT_NRh','WOT_SNRr','WOT_SNRh',...
                                                                        'WOT_RNRr','WOT_RNRh','WT_R2r','WT_R2h','WT_R2f','WT_E2rh','WT_E2rf','WT_E2hf',...
                                                                        'WT_R2rh','WT_R2hr','WT_R2rf','WT_R2fr','WT_R2hf','WT_R2fh'};
    end

        
    ME_F_restricted_info   = niftiinfo(F_restricted);
    ME_F_hindered_info     = niftiinfo(F_hindered);
    ME_F_isotropic_info    = niftiinfo(F_isotropic);
    ME_mask_info           = niftiinfo(fmask);

    ME_F_restricted        = niftiread(ME_F_restricted_info);
    ME_F_hindered          = niftiread(ME_F_hindered_info);
    ME_F_isotropic         = niftiread(ME_F_isotropic_info);
    ME_mask                = single(round(niftiread(ME_mask_info)));


    if exist(options.F_T2,'file')
        ME_T2_info   = niftiinfo(options.F_T2);
        ME_T2        = niftiread(ME_T2_info);
        [~,~,~,nm] = size(ME_T2);

        ME_T2_restricted        = ME_T2(:,:,:,1:nm/3);
        ME_T2_hindered          = ME_T2(:,:,:,1+nm/3 : 2*nm/3);
        ME_T2_isotropic         = ME_T2(:,:,:,1+2*nm/3:end);

    end

    default_spectrum = load(options.spectrum);

    adc_restricted   = default_spectrum.adc_restricted;
    adc_hindered     = default_spectrum.adc_hindered;
    adc_isotropic    = default_spectrum.adc_isotropic;

    num_restricted  = size(adc_restricted,1);  
    num_hindered    = size(adc_hindered,1);    
    num_isotropic   = size(adc_isotropic,1);  

    if ~exist(outpath,'dir')
        mkdir('outpath')
    end

    %% WOT: without T2
    for i = 1:length(options.index)
        switch options.index{i}
            case 'WOT_AVF'  % Anisotropic VF 
                ME_vf   = cat(4,ME_F_restricted,ME_F_hindered,ME_F_isotropic);
                data = (sum(ME_F_restricted,4)+sum(ME_F_hindered,4))./ sum(ME_vf,4);

            case 'WOT_ICVF' % Intra-cellular VF 
                WOT_AVF = sum(ME_F_restricted,4)+sum(ME_F_hindered,4);
                data = sum(ME_F_restricted,4) ./ WOT_AVF;

            case 'WOT_ECVF' % Extra-cellular VF 
                WOT_AVF = sum(ME_F_restricted,4)+sum(ME_F_hindered,4);
                data = sum(ME_F_hindered,4) ./ WOT_AVF;

            case 'WOT_IVF' % Isotropic VF
                ME_vf   = cat(4,ME_F_restricted,ME_F_hindered,ME_F_isotropic);
                data = sum(ME_F_isotropic,4)./ sum(ME_vf,4);

            case 'WOT_uAD' % Microscopic AD
                ME_vf   = cat(4,ME_F_restricted,ME_F_hindered,ME_F_isotropic);
                data = sum(ME_vf.*reshape([adc_restricted(:,1); adc_hindered(:,1); adc_isotropic(:,1)],1,1,1,num_restricted+num_hindered+num_isotropic),4) ./ sum(ME_vf,4);

            case 'WOT_uRD' % Microscopic RD
                ME_vf   = cat(4,ME_F_restricted,ME_F_hindered,ME_F_isotropic);
                data = sum(ME_vf.*reshape([adc_restricted(:,2); adc_hindered(:,2); adc_isotropic(:,1)],1,1,1,num_restricted+num_hindered+num_isotropic),4) ./ sum(ME_vf,4);

            case 'WOT_uMD' % Microscopic MD
                ME_vf   = cat(4,ME_F_restricted,ME_F_hindered,ME_F_isotropic);
                WOT_uAD = sum(ME_vf.*reshape([adc_restricted(:,1); adc_hindered(:,1); adc_isotropic(:,1)],1,1,1,num_restricted+num_hindered+num_isotropic),4) ./ sum(ME_vf,4);
                WOT_uRD = sum(ME_vf.*reshape([adc_restricted(:,2); adc_hindered(:,2); adc_isotropic(:,1)],1,1,1,num_restricted+num_hindered+num_isotropic),4) ./ sum(ME_vf,4);
                data = (WOT_uAD + WOT_uRD * 2) ./ 3;

            case 'WOT_uFA' % Microscopic FA
                ME_vf   = cat(4,ME_F_restricted,ME_F_hindered,ME_F_isotropic);
                WOT_uAD = sum(ME_vf.*reshape([adc_restricted(:,1); adc_hindered(:,1); adc_isotropic(:,1)],1,1,1,num_restricted+num_hindered+num_isotropic),4) ./ sum(ME_vf,4);
                WOT_uRD = sum(ME_vf.*reshape([adc_restricted(:,2); adc_hindered(:,2); adc_isotropic(:,1)],1,1,1,num_restricted+num_hindered+num_isotropic),4) ./ sum(ME_vf,4);
                data = (WOT_uAD - WOT_uRD) ./ sqrt(WOT_uAD.^2 + 2*WOT_uRD.^2);

            case 'WOT_uICAD' % Intra-cellular AD
                data = sum(ME_F_restricted.*reshape(adc_restricted(:,1),1,1,1,num_restricted),4) ./ sum(ME_F_restricted,4);
                
            case 'WOT_uICRD' % Intra-cellular RD
                data = sum(ME_F_restricted.*reshape(adc_restricted(:,2),1,1,1,num_restricted),4) ./ sum(ME_F_restricted,4);

            case 'WOT_uECAD' % Extra-cellular AD
                data = sum(ME_F_hindered.*reshape(adc_hindered(:,1),1,1,1,num_hindered),4) ./ sum(ME_F_hindered,4);

            case 'WOT_uECRD' % Extra-cellular RD
                data = sum(ME_F_hindered.*reshape(adc_hindered(:,2),1,1,1,num_hindered),4) ./ sum(ME_F_hindered,4);

            case 'WOT_uCS' % Microscopic sphericity
                ME_vf   = cat(4,ME_F_restricted,ME_F_hindered,ME_F_isotropic);
                WOT_uAD = sum(ME_vf.*reshape([adc_restricted(:,1); adc_hindered(:,1); adc_isotropic(:,1)],1,1,1,num_restricted+num_hindered+num_isotropic),4) ./ sum(ME_vf,4);
                WOT_uRD = sum(ME_vf.*reshape([adc_restricted(:,2); adc_hindered(:,2); adc_isotropic(:,1)],1,1,1,num_restricted+num_hindered+num_isotropic),4) ./ sum(ME_vf,4);
                data= WOT_uRD./WOT_uAD;

            case 'WOT_uCL' % Microscopic linearity
                ME_vf   = cat(4,ME_F_restricted,ME_F_hindered,ME_F_isotropic);
                WOT_uAD = sum(ME_vf.*reshape([adc_restricted(:,1); adc_hindered(:,1); adc_isotropic(:,1)],1,1,1,num_restricted+num_hindered+num_isotropic),4) ./ sum(ME_vf,4);
                WOT_uRD = sum(ME_vf.*reshape([adc_restricted(:,2); adc_hindered(:,2); adc_isotropic(:,1)],1,1,1,num_restricted+num_hindered+num_isotropic),4) ./ sum(ME_vf,4);
                WOT_uMD = (WOT_uAD + WOT_uRD * 2) ./ 3;
                data= (WOT_uAD-WOT_uMD)./WOT_uMD;

            case 'WOT_uCP' % Microscopic plane
                ME_vf   = cat(4,ME_F_restricted,ME_F_hindered,ME_F_isotropic);
                WOT_uAD = sum(ME_vf.*reshape([adc_restricted(:,1); adc_hindered(:,1); adc_isotropic(:,1)],1,1,1,num_restricted+num_hindered+num_isotropic),4) ./ sum(ME_vf,4);
                WOT_uRD = sum(ME_vf.*reshape([adc_restricted(:,2); adc_hindered(:,2); adc_isotropic(:,1)],1,1,1,num_restricted+num_hindered+num_isotropic),4) ./ sum(ME_vf,4);
                WOT_uMD = (WOT_uAD + WOT_uRD * 2) ./ 3;
                data = (WOT_uMD-WOT_uRD)./WOT_uAD;
                
            case 'WOT_NRr' % mean effective Intra-neurite radius 
                data = mean((ME_F_restricted.*reshape(adc_restricted(:,1).*adc_restricted(:,2),1,1,1,num_restricted)).^0.25,4);

            case 'WOT_NRh' % mean effective Extra-neurite radius 
                data = mean((ME_F_hindered.*reshape(adc_hindered(:,1).*adc_hindered(:,2),1,1,1,num_hindered)).^0.25,4);

            case 'WOT_SNRr' % std effective Intra-neurite radius 
                data = std((ME_F_restricted.*reshape(adc_restricted(:,1).*adc_restricted(:,2),1,1,1,num_restricted)).^0.25,[],4);

            case 'WOT_SNRh' % std effective Extra-neurite radius 
                data = std((ME_F_hindered.*reshape(adc_hindered(:,1).*adc_hindered(:,2),1,1,1,num_hindered)).^0.25,[],4);

            case 'WOT_RNRr' % relative effective Intra-neurite radius
                WOT_SNRr = std((ME_F_restricted.*reshape(adc_restricted(:,1).*adc_restricted(:,2),1,1,1,num_restricted)).^0.25,[],4);
                WOT_NRr = mean((ME_F_restricted.*reshape(adc_restricted(:,1).*adc_restricted(:,2),1,1,1,num_restricted)).^0.25,4);
                data = WOT_SNRr./WOT_NRr;

            case 'WOT_RNRh' % relative effective Extra-neurite radius
                WOT_NRh = mean((ME_F_hindered.*reshape(adc_hindered(:,1).*adc_hindered(:,2),1,1,1,num_hindered)).^0.25,4);
                WOT_SNRh = std((ME_F_hindered.*reshape(adc_hindered(:,1).*adc_hindered(:,2),1,1,1,num_hindered)).^0.25,[],4);
                data = WOT_SNRh./WOT_NRh;

            case 'WT_R2r' % relaxation rate on restricted
                data = 1./(ME_T2_restricted);
            
            case 'WT_R2h' % relaxation rate on hindered
                data = 1./(ME_T2_hindered);
            
            case 'WT_R2f' % relaxation rate on isotropic
                data = 1./(ME_T2_isotropic);
            
            case 'WT_E2rh' % relaxation rate on exchange between restricted and hindered
                data = abs(1./(ME_T2_restricted-ME_T2_hindered));
                % data = abs(1./(ME_T2_restricted)-(1./(ME_T2_hindered)));
            
            case 'WT_E2rf' % relaxation rate on exchange between restricted and isotropic
                data = abs(1./(ME_T2_restricted-ME_T2_isotropic));
                % data = abs(1./(ME_T2_restricted)-(1./(ME_T2_isotropic)));
            
            case 'WT_E2hf' % relaxation rate on exchange between hindered and isotropic
                data = abs(1./(ME_T2_hindered-ME_T2_isotropic));
                % data = abs(1./(ME_T2_hindered)-(1./(ME_T2_isotropic)));

            case 'WT_R2rh' % relational rate on exchange between restricted and hindered
                data = ME_T2_restricted./ME_T2_hindered;
            
            case 'WT_R2hr' % relational rate on exchange between restricted and hindered
                data = ME_T2_hindered./ME_T2_restricted;

            case 'WT_R2rf' % relational rate on exchange between restricted and isotropic
                data = ME_T2_restricted./ME_T2_isotropic;

            case 'WT_R2fr' % relational rate on exchange between restricted and isotropic
                data = ME_T2_isotropic./ME_T2_restricted;

            case 'WT_R2hf' % relational rate on exchange between hindered and isotropic
                data = ME_T2_hindered./ME_T2_isotropic;

            case 'WT_R2fh' % relational rate on exchange between hindered and isotropic
                data = ME_T2_isotropic./ME_T2_hindered;
        end
        info = ME_F_restricted_info;
        info.ImageSize = size(data);
        info.PixelDimensions = info.PixelDimensions(1:length(size(data)));
        niftiwrite(single(max(0,data).*ME_mask),fullfile(outpath,strcat(options.index{i},'.nii')),info,'Compressed', true);
    end
end

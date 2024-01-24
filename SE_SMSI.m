function SE_SMSI(fdwi,fbval,fmask,outpath,options)
    %{
    ░█▀▀█ ░█▀▀▀█ ░█▀▄▀█ ░█▀▀▀ ░█▀▀▄ ▀█▀
    ░█─── ░█──░█ ░█░█░█ ░█▀▀▀ ░█─░█ ░█─
    ░█▄▄█ ░█▄▄▄█ ░█──░█ ░█▄▄▄ ░█▄▄▀ ▄█▄

    Single-echo spherical mean spectrum imaging (SE-SMSI)


        Created by Ye Wu, PhD (dr.yewu@outlook.com)

        - Nanjing University of Science and Technology, China
        - University of North Carolina at Chapel Hill, USA
        
    %}
    
    arguments
        fdwi                    string   {mustBeFile}
        fbval                   string   {mustBeFile}
        fmask                   string   {mustBeFile}
        outpath                 string   {mustBeNonzeroLengthText}

        options.normalizeToS0   (1,1)    {mustBeNumericOrLogical} = true
        options.useBshell       (1,:)    {mustBeNumericOrLogical} = false
        options.spectrum        string   {mustBeFile} = 'scheme/default_spectrum.mat'
    end

    addpath('third/osqp');
    addpath('scheme');

    assert(exist(fdwi,'file'),'Input mask %s does not exist', fdwi);
    assert(exist(fbval,'file'),'Input mask %s does not exist', fbval);
    assert(exist(fmask,'file'),'Input mask %s does not exist', fmask);
    assert(exist(options.spectrum,'file'),'Input spectrum %s does not exist', options.spectrum);

    %% load single-echo dMRI dataset
    SE_dwi_info     = niftiinfo(fdwi);
    SE_dwi          = niftiread(SE_dwi_info);
    SE_bval         = round(importdata(fbval)'/100)*100;
    SE_mask_info    = niftiinfo(fmask);
    SE_mask         = round(niftiread(SE_mask_info));
    clear fdwi fmask SE_mask_info;

    if options.useBshell
        ind         = ismember(fbval,options.useBshell);
        SE_dwi      = SE_dwi(:,:,:,ind);
        SE_bval     = SE_bval(ind);
        clear ind
    end

    SE_bshell       = unique(SE_bval);
    SE_dwi_mean     = single(zeros([size(SE_dwi,[1,2,3]),length(SE_bshell)]));
    for j = 1:length(SE_bshell)
        ind = SE_bval == SE_bshell(j);
        SE_dwi_mean(:,:,:,j) = mean(SE_dwi(:,:,:,ind),4);
    end

    %% Normalization S/S0
    if options.normalizeToS0
        SE_dwi_mean = SE_dwi_mean(:,:,:,2:end)./(SE_dwi_mean(:,:,:,1)+eps);
        SE_bshell   = SE_bshell(2:end);
    end

    %% kernel
    default_spectrum = load(options.spectrum);

    adc_restricted   = default_spectrum.adc_restricted;
    adc_hindered     = default_spectrum.adc_hindered;
    adc_isotropic    = default_spectrum.adc_isotropic;

    num_restricted  = size(adc_restricted,1);  
    num_hindered    = size(adc_hindered,1);    
    num_isotropic   = size(adc_isotropic,1);   
    
    kernel_restricted = cell(1,num_restricted);
    kernel_hindered   = cell(1,num_hindered);
    kernel_isotropic  = cell(1,num_isotropic);

    for j = 1:num_restricted
        kernel_restricted{1,j} = single(zeros(length(SE_bshell),1));
        for k = 1:length(SE_bshell)
            R = sqrt(SE_bshell(k)*(adc_restricted(j,1)-adc_restricted(j,2)));
            kernel_restricted{1,j}(k,1) = (sqrt(pi)*exp(-SE_bshell(k)*adc_restricted(j,2))*erf(R))/(2*R);
        end
    end

    for j = 1:num_hindered
        kernel_hindered{1,j} = single(zeros(length(SE_bshell),1));
        for k = 1:length(SE_bshell)
            R = sqrt(SE_bshell(k)*(adc_hindered(j,1)-adc_hindered(j,2)));
            kernel_hindered{1,j}(k,1) = (sqrt(pi)*exp(-SE_bshell(k)*adc_hindered(j,2))*erf(R))/(2*R);
        end
    end

    for j = 1:num_isotropic
        kernel_isotropic{1,j} = single(zeros(length(SE_bshell),1));
        for k = 1:length(SE_bshell)
            kernel_isotropic{1,j}(k,1) = exp(-SE_bshell(k)*adc_isotropic(j)); 
        end
    end
    
    blk_kernel = cell2mat([kernel_restricted kernel_hindered kernel_isotropic]);


    %% Vectorization & Masked & arrayed
    SE_mask_ind = find(SE_mask>0.5);
    SE_dwi_mean_array = reshape(SE_dwi_mean,size(SE_dwi_mean,1)*size(SE_dwi_mean,2)*size(SE_dwi_mean,3),size(SE_dwi_mean,4));
    SE_dwi_mean_array = SE_dwi_mean_array(SE_mask_ind,:)';
    SE_dwi_mean_array_ind = all(SE_dwi_mean_array);
    SE_dwi_mean_array = SE_dwi_mean_array(:,SE_dwi_mean_array_ind);
    clear SE_dwi_mean

    %% sum(total vf) = 1 if normalizeToS0 is true
    alpha_coef = zeros(num_restricted+num_hindered+num_isotropic,size(SE_dwi_mean_array,2));
    alpha_nmse = zeros(1,size(SE_dwi_mean_array,2));

    H   = double(blk_kernel'*blk_kernel);
    K   = double(-blk_kernel'*SE_dwi_mean_array);
    A1  = [diag(ones(size(blk_kernel,2),1)); ones(1,size(blk_kernel,2))];
    A2  = [zeros(size(blk_kernel,2),1);1];
    A3  = [ones(size(blk_kernel,2),1);1];

    %% optimization
    parfor i = 1:size(SE_dwi_mean_array,2)
        
        f = K(:,i);
        S = SE_dwi_mean_array(:,i);
        prob = osqp;
        prob.setup(H,f,A1,A2,A3,'alpha',0.1,'verbose',0);
        res = prob.solve();
        alpha_coef(:,i) = res.x;
        alpha_nmse(1,i) = norm(blk_kernel * res.x - S) ./ norm(S);
    end
    
    %% save results
    if ~exist(outpath,'dir')
        mkdir(outpath)
    end

    % save VF restricted
    temp2 = single(zeros(num_restricted,length(SE_dwi_mean_array_ind)));
    temp2(:,SE_dwi_mean_array_ind) = alpha_coef(1:num_restricted,:);
    temp = single(zeros(num_restricted,size(SE_mask,1)*size(SE_mask,2)*size(SE_mask,3)));
    temp(:,SE_mask_ind) = temp2;
    vf = reshape(temp',size(SE_mask,1),size(SE_mask,2),size(SE_mask,3),num_restricted);

    info_vf= SE_dwi_info;
    info_vf.Datatype = 'single';
    info_vf.ImageSize = size(vf);
    info_vf.PixelDimensions = info_vf.PixelDimensions(1:length(size(vf)));
    niftiwrite(vf,fullfile(outpath,'VF_restricted.nii'),info_vf,'Compressed', true);

    % save VF hindered
    temp2 = single(zeros(num_hindered,length(SE_dwi_mean_array_ind)));
    temp2(:,SE_dwi_mean_array_ind) = alpha_coef(1+num_restricted:num_restricted+num_hindered,:);
    temp = single(zeros(num_hindered,size(SE_mask,1)*size(SE_mask,2)*size(SE_mask,3)));
    temp(:,SE_mask_ind) = temp2;
    vf = reshape(temp',size(SE_mask,1),size(SE_mask,2),size(SE_mask,3),num_hindered);

    info_vf= SE_dwi_info;
    info_vf.Datatype = 'single';
    info_vf.ImageSize = size(vf);
    info_vf.PixelDimensions = info_vf.PixelDimensions(1:length(size(vf)));
    info_vf.PixelDimensions = info_vf.PixelDimensions(1:length(size(vf)));
    niftiwrite(vf,fullfile(outpath,'VF_hindered.nii'),info_vf,'Compressed', true);

    % save VF free water
    temp2 = single(zeros(num_isotropic,length(SE_dwi_mean_array_ind)));
    temp2(:,SE_dwi_mean_array_ind) = alpha_coef(1+num_restricted+num_hindered:num_restricted+num_hindered+num_isotropic,:);
    temp = single(zeros(num_isotropic,size(SE_mask,1)*size(SE_mask,2)*size(SE_mask,3)));
    temp(:,SE_mask_ind) = temp2;
    vf = reshape(temp',size(SE_mask,1),size(SE_mask,2),size(SE_mask,3),num_isotropic);

    info_vf= SE_dwi_info;
    info_vf.Datatype = 'single';
    info_vf.ImageSize = size(vf);
    info_vf.PixelDimensions = info_vf.PixelDimensions(1:length(size(vf)));
    niftiwrite(vf,fullfile(outpath,'VF_free.nii'),info_vf,'Compressed', true);

    % save NMSE
    temp2 = single(zeros(1,length(SE_dwi_mean_array_ind)));
    temp2(:,SE_dwi_mean_array_ind) = alpha_nmse;
    temp = single(zeros(1,size(SE_mask,1)*size(SE_mask,2)*size(SE_mask,3)));
    temp(:,SE_mask_ind) = temp2;
    vf = reshape(temp',size(SE_mask,1),size(SE_mask,2),size(SE_mask,3));

    info_vf= SE_dwi_info;
    info_vf.Datatype = 'single';
    info_vf.ImageSize = size(vf);
    info_vf.PixelDimensions = info_vf.PixelDimensions(1:length(size(vf)));
    niftiwrite(vf,fullfile(outpath,'NMSE.nii'),info_vf,'Compressed', true);
end
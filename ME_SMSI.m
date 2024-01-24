function ME_SMSI(fdwi,fbval,fmask,TE,outpath,options)
    %{
    ░█▀▀█ ░█▀▀▀█ ░█▀▄▀█ ░█▀▀▀ ░█▀▀▄ ▀█▀
    ░█─── ░█──░█ ░█░█░█ ░█▀▀▀ ░█─░█ ░█─
    ░█▄▄█ ░█▄▄▄█ ░█──░█ ░█▄▄▄ ░█▄▄▀ ▄█▄

    Multi-echo spherical mean spectrum imaging (ME-SMSI)


        Created by Ye Wu, PhD (dr.yewu@outlook.com)

        - Nanjing University of Science and Technology, China
        - University of North Carolina at Chapel Hill, USA
        
    %}
    
    arguments
        fdwi                    (1,:)    {mustBeNonzeroLengthText}
        fbval                   (1,:)    {mustBeNonzeroLengthText}
        fmask                   string   {mustBeFile}
        TE                      (1,:)    {mustBeNumeric}
        outpath                 string   {mustBeNonzeroLengthText}

        options.normalizeToS0   (1,1)    {mustBeNumericOrLogical} = true
        options.useBshell       (1,:)    {mustBeNumericOrLogical} = false
        options.lambda          (1,1)    {mustBeNumeric} = 0.015
        options.x0              (1,1)    {mustBeNumeric} = 100
        options.spectrum        string   {mustBeFile} = 'scheme/default_spectrum.mat'
    end

    addpath('third/osqp');
    addpath('scheme');

    cellfun(@(x)assert(exist(x,'file'),'Input DWI %s does not exist', x),fdwi,'UniformOutput',false);
    cellfun(@(x)assert(exist(x,'file'),'Input Bval %s does not exist', x),fbval,'UniformOutput',false);
    assert(exist(fmask,'file'),'Input mask %s does not exist', fmask);
    assert(exist(options.spectrum,'file'),'Input spectrum %s does not exist', options.spectrum);

    %% load multi-echo dMRI dataset
    ME_dwi_info     = cellfun(@(x)niftiinfo(x),fdwi,'UniformOutput',false);
    ME_dwi          = cellfun(@(x)niftiread(x),ME_dwi_info,'UniformOutput',false);
    ME_bval         = cellfun(@(x)round(importdata(x)'/100)*100,fbval,'UniformOutput',false); 
    ME_mask_info    = niftiinfo(fmask);
    ME_mask         = round(niftiread(ME_mask_info));
    clear ind fdwi fmask;

    if options.useBshell
        ind         = cellfun(@(x)find(ismember(x,options.useBshell)),fbval,'UniformOutput',false); 
        ME_dwi      = cellfun(@(x,y)x(:,:,:,y),ME_dwi,ind,'UniformOutput',false);
        ME_bval     = cellfun(@(x,y)x(y),ME_bval,ind,'UniformOutput',false); 
        clear ind
    end

    ME_bshell       = cellfun(@(x)unique(x),ME_bval,'UniformOutput',false); 
    ME_dwi_mean     = cell(size(ME_dwi));
    for i = 1:length(TE)
        ME_dwi_mean{i} = single(zeros([size(ME_dwi{i},[1,2,3]),length(ME_bshell{i})]));
        for j = 1:length(ME_bshell{i})
            ind = ME_bval{i} == ME_bshell{i}(j);
            ME_dwi_mean{i}(:,:,:,j) = mean(ME_dwi{i}(:,:,:,ind),4);
        end
    end

    %% Normalization S/S0
    if options.normalizeToS0
        ME_dwi_mean = cellfun(@(x)x(:,:,:,2:end)./(x(:,:,:,1)+eps),ME_dwi_mean,'UniformOutput',false);
        ME_bshell   = cellfun(@(x)x(2:end),ME_bshell,'UniformOutput',false);
    end

    %% kernel
    default_spectrum = load(options.spectrum);

    adc_restricted   = default_spectrum.adc_restricted;
    adc_hindered     = default_spectrum.adc_hindered;
    adc_isotropic    = default_spectrum.adc_isotropic;

    num_restricted  = size(adc_restricted,1);  
    num_hindered    = size(adc_hindered,1);    
    num_isotropic   = size(adc_isotropic,1);   
    
    kernel_restricted = cell(length(TE),num_restricted);
    kernel_hindered   = cell(length(TE),num_hindered);
    kernel_isotropic  = cell(length(TE),num_isotropic);


    for i = 1:length(TE)
        for j = 1:num_restricted
            kernel_restricted{i,j} = single(zeros(length(ME_bshell{i}),1));
            for k = 1:length(ME_bshell{i})
                R = sqrt(ME_bshell{i}(k)*(adc_restricted(j,1)-adc_restricted(j,2)));
                kernel_restricted{i,j}(k,1) = (sqrt(pi)*exp(-ME_bshell{i}(k)*adc_restricted(j,2))*erf(R))/(2*R);
            end
        end

        for j = 1:num_hindered
            kernel_hindered{i,j} = single(zeros(length(ME_bshell{i}),1));
            for k = 1:length(ME_bshell{i})
                R = sqrt(ME_bshell{i}(k)*(adc_hindered(j,1)-adc_hindered(j,2)));
                kernel_hindered{i,j}(k,1) = (sqrt(pi)*exp(-ME_bshell{i}(k)*adc_hindered(j,2))*erf(R))/(2*R);
            end
        end

        for j = 1:num_isotropic
            kernel_isotropic{i,j} = single(zeros(length(ME_bshell{i}),1));
            for k = 1:length(ME_bshell{i})
                kernel_isotropic{i,j}(k,1) = exp(-ME_bshell{i}(k)*adc_isotropic(j)); 
            end
        end
    end

    blk_kernel_restricted = [];
    for i = 1:num_restricted
        blk_kernel_restricted = [blk_kernel_restricted blkdiag(kernel_restricted{:,i})];
    end
    
    blk_kernel_hindered = [];
    for i = 1:num_hindered
        blk_kernel_hindered = [blk_kernel_hindered blkdiag(kernel_hindered{:,i})];
    end

    blk_kernel_isotropic = [];
    for i = 1:num_isotropic
        blk_kernel_isotropic = [blk_kernel_isotropic blkdiag(kernel_isotropic{:,i})];
    end
    
    blk_kernel = [blk_kernel_restricted blk_kernel_hindered blk_kernel_isotropic];


    %% Vectorization & Masked & arrayed
    ME_mask_ind = find(ME_mask>0.5);
    ME_dwi_mean_array = cellfun(@(x)reshape(x,size(x,1)*size(x,2)*size(x,3),size(x,4)),ME_dwi_mean,'UniformOutput',false);
    ME_dwi_mean_array = cellfun(@(x)x(ME_mask_ind,:)',ME_dwi_mean_array,'UniformOutput',false);
    ME_dwi_mean_array = cell2mat(ME_dwi_mean_array');
    ME_dwi_mean_array_ind = all(ME_dwi_mean_array);
    ME_dwi_mean_array = ME_dwi_mean_array(:,ME_dwi_mean_array_ind);
    clear ME_dwi_mean

    %% sum(total vf) = 1 if normalizeToS0 is true
    if options.normalizeToS0
        Aeq     = [zeros(1,3) ones(1,num_restricted+num_hindered+num_isotropic)];
        beq     = 1;
    else
        Aeq     = [];
        beq     = [];
    end
    lb = zeros(num_restricted+num_hindered+num_isotropic+3,1);
    option = optimoptions(@fmincon,'Display','off');
    alpha_coef = zeros(num_restricted+num_hindered+num_isotropic+3,size(ME_dwi_mean_array,2));
    alpha_init = [options.x0*ones(3,1);rand(num_restricted+num_hindered+num_isotropic,1)];

    I = [blk_kernel;options.lambda*diag(ones(1,size(blk_kernel,2)))];
    H   = double(I'*I);
    K   = double(-I'*[ME_dwi_mean_array;zeros(size(blk_kernel,2),size(ME_dwi_mean_array,2))]);
    A1  = diag(ones(size(blk_kernel,2),1));
    A2  = zeros(size(blk_kernel,2),1);
    A3  = ones(size(blk_kernel,2),1);

    A1_new  = diag(ones(num_restricted+num_hindered+num_isotropic,1));
    A2_new  = zeros(num_restricted+num_hindered+num_isotropic,1);
    A3_new  = ones(num_restricted+num_hindered+num_isotropic,1);

    clear blk_kernel;

    %% optimization
    parfor i = 1:size(ME_dwi_mean_array,2)
        
        f = K(:,i);
        prob = osqp;
        prob.setup(H,f,A1,A2,A3,'alpha',0.1,'verbose',0);
        res = prob.solve();
        fun = @(x)echoFunc(x,res.x,TE,num_restricted,num_hindered,num_isotropic);
        beta = max(0,fmincon(fun,alpha_init,[],[],Aeq,beq,lb,[],[],option));

        beta_restricted  = exp(-TE'./beta(1));
        beta_hindered    = exp(-TE'./beta(2));
        beta_isotropic   = exp(-TE'./beta(3));

        kernel_new = [reshape(repmat(beta_restricted',length(ME_bshell{1}),1),[],1) .* cell2mat(kernel_restricted) ...
                        reshape(repmat(beta_hindered',length(ME_bshell{1}),1),[],1) .* cell2mat(kernel_hindered) ...
                        reshape(repmat(beta_isotropic',length(ME_bshell{1}),1),[],1) .* cell2mat(kernel_isotropic)];

        f_new = ME_dwi_mean_array(:,i);
        I_new = [kernel_new;options.lambda*diag(ones(1,size(kernel_new,2)))];
        H_new   = double(I_new'*I_new);
        K_new   = double(-I_new'*[f_new;zeros(size(kernel_new,2),1)]);

        prob = osqp;
        prob.setup(H_new,K_new,A1_new,A2_new,A3_new,'alpha',0.1,'verbose',0);
        res = prob.solve();

        alpha_coef(:,i) = [beta(1:3,1);res.x];

    end
    

    %% save results
    if ~exist(outpath,'dir')
        mkdir(outpath)
    end

    % save T2 restricted
    temp2 = single(zeros(1,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(1,:);
    temp = single(zeros(1,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3));

    info_t2r = ME_dwi_info{1};
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    info_t2r.PixelDimensions = info_t2r.PixelDimensions(1:length(size(T2r)));
    niftiwrite(single(T2r),fullfile(outpath,'T2_restricted.nii'),info_t2r,'Compressed', true);

    % save T2 hindered
    temp2 = single(zeros(1,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(2,:);
    temp = single(zeros(1,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3));

    info_t2r = ME_dwi_info{1};
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    info_t2r.PixelDimensions = info_t2r.PixelDimensions(1:length(size(T2r)));
    niftiwrite(single(T2r),fullfile(outpath,'T2_hindered.nii'),info_t2r,'Compressed', true);

    % save T2 free water
    temp2 = single(zeros(1,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(3,:);
    temp = single(zeros(1,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3));

    info_t2r = ME_dwi_info{1};
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    info_t2r.PixelDimensions = info_t2r.PixelDimensions(1:length(size(T2r)));
    niftiwrite(single(T2r),fullfile(outpath,'T2_free.nii'),info_t2r,'Compressed', true);

    % save VF restricted
    temp2 = single(zeros(num_restricted,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(4:3+num_restricted,:);
    temp = single(zeros(num_restricted,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    vf = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3),num_restricted);

    info_vf= ME_dwi_info{1};
    info_vf.Datatype = 'single';
    info_vf.ImageSize = size(vf);
    info_vf.Transform = ME_dwi_info{1}.Transform;
    info_vf.PixelDimensions = info_vf.PixelDimensions(1:length(size(vf)));
    niftiwrite(single(vf),fullfile(outpath,'VF_restricted.nii'),info_vf,'Compressed', true);

    % save VF hindered
    temp2 = single(zeros(num_hindered,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(4+num_restricted:3+num_restricted+num_hindered,:);
    temp = single(zeros(num_hindered,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    vf = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3),num_hindered);

    info_vf= ME_dwi_info{1};
    info_vf.Datatype = 'single';
    info_vf.ImageSize = size(vf);
    info_vf.Transform = ME_dwi_info{1}.Transform;
    info_vf.PixelDimensions = info_vf.PixelDimensions(1:length(size(vf)));
    niftiwrite(single(vf),fullfile(outpath,'VF_hindered.nii'),info_vf,'Compressed', true);

    % save VF free water
    temp2 = single(zeros(num_isotropic,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(4+num_restricted+num_hindered:end,:);
    temp = single(zeros(num_isotropic,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    vf = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3),num_isotropic);

    info_vf= ME_dwi_info{1};
    info_vf.Datatype = 'single';
    info_vf.ImageSize = size(vf);
    info_vf.Transform = ME_dwi_info{1}.Transform;
    info_vf.PixelDimensions = info_vf.PixelDimensions(1:length(size(vf)));
    niftiwrite(single(vf),fullfile(outpath,'VF_free.nii'),info_vf,'Compressed', true);
end

function FF = echoFunc(x,alpha,TE,num_restricted,num_hindered,num_isotropic)
    % x: [ T_2^r, T_2^h, T_2^f, w_1^r, w_2^r...w_1^h, w_2^h...,w_1^f, w_2^f...]
    num = 1;
    for i = 1:num_restricted
        for j = 1:length(TE)
            F(num) = exp(-TE(j)/x(1)) * x(3+i) - alpha(num);
            num = num + 1;
        end
    end

    for i = 1:num_hindered
        for j = 1:length(TE)
            F(num) = exp(-TE(j)/x(2)) * x(3+num_restricted+i) - alpha(num);
            num = num + 1;
        end
    end

    for i = 1:num_isotropic
        for j = 1:length(TE)
            F(num) = exp(-TE(j)/x(3)) * x(3+num_restricted+num_hindered+i) - alpha(num);
            num = num + 1;
        end
    end

    FF = norm(F)^2;
end
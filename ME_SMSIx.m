function ME_SMSIx(fdwi,fbval,fmask,TE,outpath,options)
    %{
    ░█▀▀█ ░█▀▀▀█ ░█▀▄▀█ ░█▀▀▀ ░█▀▀▄ ▀█▀
    ░█─── ░█──░█ ░█░█░█ ░█▀▀▀ ░█─░█ ░█─
    ░█▄▄█ ░█▄▄▄█ ░█──░█ ░█▄▄▄ ░█▄▄▀ ▄█▄

    Multi-echo spherical mean spectrum imaging with B-dependent (ME-SMSIx)


        Created by Ye Wu, PhD (dr.yewu@outlook.com)

        - Nanjing University of Science and Technology, China, China
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
        options.lambda          (1,1)    {mustBeNumeric} = 0.01
        options.x0              (1,1)    {mustBeNumeric} = 75
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
    clear ind fdwi fbval fmask;

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
    
    kernel_restricted = zeros(length(TE),length(ME_bshell{1}),num_restricted);
    kernel_hindered   = zeros(length(TE),length(ME_bshell{1}),num_hindered);
    kernel_isotropic  = zeros(length(TE),length(ME_bshell{1}),num_isotropic);


    for i = 1:length(TE)
        for j = 1:num_restricted
            for k = 1:length(ME_bshell{i})
                R = sqrt(ME_bshell{i}(k)*(adc_restricted(j,1)-adc_restricted(j,2)));
                kernel_restricted(i,k,j) = (sqrt(pi)*exp(-ME_bshell{i}(k)*adc_restricted(j,2))*erf(R))/(2*R);
            end
        end

        for j = 1:num_hindered
            for k = 1:length(ME_bshell{i})
                R = sqrt(ME_bshell{i}(k)*(adc_hindered(j,1)-adc_hindered(j,2)));
                kernel_hindered(i,k,j) = (sqrt(pi)*exp(-ME_bshell{i}(k)*adc_hindered(j,2))*erf(R))/(2*R);
            end
        end

        for j = 1:num_isotropic
            for k = 1:length(ME_bshell{i})
                kernel_isotropic(i,k,j) = exp(-ME_bshell{i}(k)*adc_isotropic(j)); 
            end
        end
    end
   
    blk_kernel_restricted = [];
    for i = 1:num_restricted
        blk_kernel_restricted = [blk_kernel_restricted diag(reshape(kernel_restricted(:,:,i)',[],1))];
    end
    
    blk_kernel_hindered = [];
    for i = 1:num_hindered
        blk_kernel_hindered = [blk_kernel_hindered diag(reshape(kernel_hindered(:,:,i)',[],1))];
    end

    blk_kernel_isotropic = [];
    for i = 1:num_isotropic
        blk_kernel_isotropic = [blk_kernel_isotropic diag(reshape(kernel_isotropic(:,:,i)',[],1))];
    end
    
    blk_kernel = [blk_kernel_restricted blk_kernel_hindered blk_kernel_isotropic];
    nrelax = 3*size(blk_kernel,1)/length(TE);

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
        Aeq     = [zeros(1,nrelax) ones(1,num_restricted+num_hindered+num_isotropic)];
        beq     = 1;
    else
        Aeq     = [];
        beq     = [];
    end
    lb = zeros(num_restricted+num_hindered+num_isotropic+nrelax,1);
    option = optimoptions(@fmincon,'Display','off');
    alpha_coef = zeros(num_restricted+num_hindered+num_isotropic+nrelax,size(ME_dwi_mean_array,2));
    alpha_init = [options.x0*ones(nrelax,1);rand(num_restricted+num_hindered+num_isotropic,1)];

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

    kernel_restricted_new = permute(kernel_restricted,[2 1 3]);
    kernel_hindered_new = permute(kernel_hindered,[2 1 3]);
    kernel_isotropic_new = permute(kernel_isotropic,[2 1 3]);

    %% optimization
    parfor i = 1:size(ME_dwi_mean_array,2)
        
        f = K(:,i);
        prob = osqp;
        prob.setup(H,f,A1,A2,A3,'alpha',0.1,'verbose',0);
        res = prob.solve();
        fun = @(x)echoFunc(x,res.x,TE,nrelax/3,num_restricted,num_hindered,num_isotropic);
        beta = max(0,fmincon(fun,alpha_init,[],[],Aeq,beq,lb,[],[],option));

        beta_restricted = exp(-TE./beta(1:nrelax/3,1)).*kernel_restricted_new;
        beta_hindered = exp(-TE./beta(1+nrelax/3:2*nrelax/3,1)).*kernel_hindered_new;
        beta_isotropic = exp(-TE./beta(1+2*nrelax/3:nrelax,1)).*kernel_isotropic_new;

        kernel_new = reshape(cat(3,beta_restricted,beta_hindered,beta_isotropic),[],num_restricted+num_hindered+num_isotropic);
        f_new = ME_dwi_mean_array(:,i);
        I_new = [kernel_new;options.lambda*diag(ones(1,size(kernel_new,2)))];
        H_new   = double(I_new'*I_new);
        K_new   = double(-I_new'*[f_new;zeros(size(kernel_new,2),1)]);

        prob = osqp;
        prob.setup(H_new,K_new,A1_new,A2_new,A3_new,'alpha',0.1,'verbose',0);
        res = prob.solve();

        alpha_coef(:,i) = max(0,[beta(1:nrelax,1);res.x]);
    end
    
    %% save results
    if ~exist(outpath,'dir')
        mkdir(outpath)
    end

    % save T2 restricted
    temp2 = single(zeros(nrelax,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(1:nrelax,:);
    temp = single(zeros(nrelax,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3),nrelax);

    info_t2r = ME_dwi_info{1};
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    info_t2r.PixelDimensions = info_t2r.PixelDimensions(1:length(size(T2r)));
    niftiwrite(single(T2r),fullfile(outpath,'T2.nii'),info_t2r,'Compressed', true);

    % save VF restricted
    temp2 = single(zeros(num_restricted,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(nrelax+1:nrelax+num_restricted,:);
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
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(nrelax+1+num_restricted:nrelax+num_restricted+num_hindered,:);
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
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(nrelax+1+num_restricted+num_hindered:end,:);
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

function FF = echoFunc(x,alpha,TE,num_shell,num_restricted,num_hindered,num_isotropic)
    % x: [ T_2^r(b1), T_2^r(b2), T_2^h(b1), T_2^h(b2), T_2^f(b1), T_2^f(b2), w_1^r, w_2^r...w_1^h, w_2^h...,w_1^f, w_2^f...]

    num = 1;
    for i = 1:num_restricted
        for j = 1:length(TE)
            for k = 1:num_shell
                F(num) = exp(-TE(j)/x(k)) * x(3*num_shell+i) - alpha(num);
                num = num + 1;
            end
        end
    end

    for i = 1:num_hindered
        for j = 1:length(TE)
            for k = 1:num_shell
                F(num) = exp(-TE(j)/x(k+num_shell)) * x(3*num_shell+num_restricted+i) - alpha(num);
                num = num + 1;
            end
        end
    end

    for i = 1:num_isotropic
        for j = 1:length(TE)
            for k = 1:num_shell
                F(num) = exp(-TE(j)/x(k+num_shell*2)) * x(3*num_shell+num_restricted+num_hindered+i) - alpha(num);
                num = num + 1;
            end
        end
    end

    FF = norm(F)^2;
end
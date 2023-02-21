function ME_SMSI(dataFolder,TE, options)
% function ME_SMSI()
    %{
    ░█▀▀█ ░█▀▀▀█ ░█▀▄▀█ ░█▀▀▀ ░█▀▀▄ ▀█▀
    ░█─── ░█──░█ ░█░█░█ ░█▀▀▀ ░█─░█ ░█─
    ░█▄▄█ ░█▄▄▄█ ░█──░█ ░█▄▄▄ ░█▄▄▀ ▄█▄

    Multi-echo spherical mean spectrum Imaging (ME-SMSI)


        Created by Ye Wu, PhD (dr.yewu@outlook.com)

        - Nanjing University of Science and Technology
        - University of North Carolina at Chapel Hill
        
    %}
    
    arguments
        dataFolder              string   {mustBeFolder}
        TE                      (1,:)    {mustBeNumeric}
        
        options.normalizeToS0   (1,1)    {mustBeNumericOrLogical} = true
        options.useBshell       (1,:)    {mustBeNumericOrLogical} = false
        options.adc_restricted  (2,:)    {mustBeNumericOrLogical} = false
        options.adc_hindered    (2,:)    {mustBeNumericOrLogical} = false
        options.adc_isotropic   (1,:)    {mustBeNumericOrLogical} = false
        options.core            (1,1)    {mustBeInteger,mustBeNonnegative} = 0
    end

    addpath('third/osqp');
    addpath('scheme');
    %% load multi-echo dMRI dataset
    fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
    fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
    fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');

    ME_dwi_info     = cellfun(@(x)niftiinfo(x),fdwi,'UniformOutput',false);
    ME_dwi          = cellfun(@(x)niftiread(x),ME_dwi_info,'UniformOutput',false);
    ME_bval         = cellfun(@(x)round(importdata(x)'/100)*100,fbval,'UniformOutput',false); 
    ME_mask_info    = niftiinfo(fmask);
    ME_mask         = round(niftiread(ME_mask_info));
    clear ind fdwi fbval fmask ME_mask_info;

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
    default_spectrum = load('default_spectrum.mat');
    if ~options.adc_restricted
        adc_restricted = default_spectrum.adc_restricted;
    end

    if ~options.adc_hindered
        adc_hindered = default_spectrum.adc_hindered;
    end

    if ~options.adc_isotropic
        adc_isotropic = default_spectrum.adc_isotropic;
    end

    num_restricted  = size(adc_restricted,1);  kernel_restricted = cell(length(TE),num_restricted);
    num_hindered    = size(adc_hindered,1);    kernel_hindered   = cell(length(TE),num_hindered);
    num_isotropic   = size(adc_isotropic,1);   kernel_isotropic  = cell(length(TE),num_isotropic);
    
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


    %% optimization sum(sum of tissue-based vf) = 1;
%     lb      = [10;10;10;zeros(num_restricted+num_hindered+num_isotropic,1)];
%     ub      = [inf;inf;inf;ones(num_restricted+num_hindered+num_isotropic,1)];
%     Aeq     = [zeros(3,3) blkdiag(ones(1,num_restricted),ones(1,num_hindered),ones(1,num_isotropic))];
%     beq     = [1;1;1];

    %% sum(total vf) = 1;
    lb      = [10;10;10;zeros(num_restricted+num_hindered+num_isotropic,1)];
    ub      = [inf;inf;inf;ones(num_restricted+num_hindered+num_isotropic,1)];
    Aeq     = [zeros(1,3) ones(1,num_restricted+num_hindered+num_isotropic)];
    beq     = 1;
    options = optimoptions(@fmincon,'Display','off');
    alpha_coef = zeros(num_restricted+num_hindered+num_isotropic+3,size(ME_dwi_mean_array,2));
    alpha_init = [100;100;100;rand(num_restricted+num_hindered+num_isotropic,1)];

    H = double(blk_kernel'*blk_kernel);
    A1 = diag(ones(size(blk_kernel,2),1));
    A2 = zeros(size(blk_kernel,2),1);
    A3 = ones(size(blk_kernel,2),1);
    parfor i = 1:size(ME_dwi_mean_array,2)
        
        f = -double(blk_kernel'*ME_dwi_mean_array(:,i));
        prob = osqp;
        prob.setup(H,f,A1,A2,A3,'alpha',0.1,'verbose',0);
        res = prob.solve();
        beta = res.x;

        fun = @(x)echoFunc(x,beta,TE,num_restricted,num_hindered,num_isotropic);
%         alpha_coef(:,i) = fmincon(fun,alpha_init,[],[],Aeq,beq,lb,ub,[],options);
        alpha_coef(:,i) = fmincon(fun,alpha_init,[],[],Aeq,beq,[],[],[],options);
    end

    T2r = alpha_coef(1:3,:);
    vf = alpha_coef(4:end,:);

    temp2 = single(zeros(3,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = T2r;
    temp = single(zeros(3,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3),3);

    info_t2r = ME_dwi_info{1};
    info_t2r.ImageSize = size(T2r);
    niftiwrite(T2r,fullfile(dataFolder,'T2r_norm.nii'),info_t2r,'Compressed', true);

    temp2 = single(zeros(num_restricted+num_hindered+num_isotropic,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = vf;
    temp = single(zeros(num_restricted+num_hindered+num_isotropic,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    vf = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3),num_restricted+num_hindered+num_isotropic);

    info_vf= ME_dwi_info{1};
    info_vf.ImageSize = size(vf);
    niftiwrite(vf,fullfile(dataFolder,'VF_norm.nii'),info_vf,'Compressed', true);

    % save options
    save(fullfile(dataFolder,'params_norm.mat'),'num_restricted','adc_restricted',...
                                           'num_hindered','adc_hindered',...
                                           'num_isotropic','adc_isotropic',...
                                            '-v7.3');
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

    FF = norm(F);
end





















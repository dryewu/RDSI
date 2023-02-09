function demo_ME_SMSI()
    %{
    ░█▀▀█ ░█▀▀▀█ ░█▀▄▀█ ░█▀▀▀ ░█▀▀▄ ▀█▀
    ░█─── ░█──░█ ░█░█░█ ░█▀▀▀ ░█─░█ ░█─
    ░█▄▄█ ░█▄▄▄█ ░█──░█ ░█▄▄▄ ░█▄▄▀ ▄█▄

    Parameter estimation with multi-echo spherical mean spectrum Imaging (ME-SMSI)


        Created by Ye Wu, PhD (dr.yewu@outlook.com)

        - Nanjing University of Science and Technology
        - University of North Carolina at Chapel Hill
        
    %}
    
    run('third/spams/start_spams.m');
    
    %% load multi-echo dMRI dataset
    dataFolder = '/home/wuye/D_disk/MTE_dMRI/WDZ';
    TE = [75 85 95 105 115 125 135];
    
    fdwi    = arrayfun(@(x)fullfile(dataFolder,'proc',strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased.nii.gz')),TE,'UniformOutput',false);
    fbvec   = arrayfun(@(x)fullfile(dataFolder,'proc',strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased.bvec')),TE,'UniformOutput',false);
    fbval   = arrayfun(@(x)fullfile(dataFolder,'proc',strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased.bval')),TE,'UniformOutput',false);
    fmask   = fullfile(dataFolder,'proc','MTE_mask.nii.gz');

    ME_dwi_info     = cellfun(@(x)niftiinfo(x),fdwi,'UniformOutput',false);
    ME_dwi          = cellfun(@(x)niftiread(x),ME_dwi_info,'UniformOutput',false);
    ME_bval         = cellfun(@(x)round(importdata(x)'/100)*100,fbval,'UniformOutput',false); 
    ME_bvec         = cellfun(@(x)importdata(x)',fbvec,'UniformOutput',false);
    ME_mask_info    = niftiinfo(fmask);
    ME_mask         = round(niftiread(ME_mask_info));
    
    ME_bshell       = cellfun(@(x)unique(x),ME_bval,'UniformOutput',false); 
    ME_dwi_mean     = cell(size(ME_dwi));
    for i = 1:length(TE)
        ME_dwi_mean{i} = single(zeros([size(ME_dwi{i},[1,2,3]),length(ME_bshell{i})]));
        for j = 1:length(ME_bshell{i})
            ind = find(ME_bval{i} == ME_bshell{i}(j));
            ME_dwi_mean{i}(:,:,:,j) = mean(ME_dwi{i}(:,:,:,ind),4);
        end
    end

    %% Normalization S/S0
    ME_dwi_mean = cellfun(@(x)x(:,:,:,2:end)./(eps+x(:,:,:,1)),ME_dwi_mean,'UniformOutput',false);
    ME_bshell   = cellfun(@(x)x(2:end),ME_bshell,'UniformOutput',false);

    %% kernel
    adc_restricted = [];
    adc_hindered = [];
    adc_isotropic = (0 : 0.2e-3 : 3e-3)';

    for adc_long_diff = 0.5e-3 : 0.2e-3 : 1.5e-3
        for adc_trans_diff = 0e-3 : 0.2e-3 : 0.9e-3
            if adc_long_diff > adc_trans_diff * (pi/2)
                if adc_long_diff / adc_trans_diff >= (pi/2)^2
                    adc_restricted = [adc_restricted; adc_long_diff adc_trans_diff];
                else
                    adc_hindered = [adc_hindered; adc_long_diff adc_trans_diff];
                end
            end
        end
    end

    num_restricted  = size(adc_restricted,1);  kernel_restricted = cell(length(TE),num_restricted);
    num_hindered    = size(adc_hindered,1);    kernel_hindered   = cell(length(TE),num_hindered);
    num_isotropic   = size(adc_isotropic,1);   kernel_isotropic   = cell(length(TE),num_isotropic);
    
    for i = 1:length(TE)
        for j = 1:num_restricted
            kernel_restricted{i,j} = single(zeros(length(ME_bshell{i}),1));
            for k = 1:length(ME_bshell{i})
                R = sqrt((ME_bshell{i}(k)+eps)*(adc_restricted(j,1)-adc_restricted(j,2)));
                kernel_restricted{i,j}(k,1) = (sqrt(pi)*exp(-(ME_bshell{i}(k)+eps)*adc_restricted(j,2))*erf(R))/(2*R);
            end
        end

        for j = 1:num_hindered
            kernel_hindered{i,j} = single(zeros(length(ME_bshell{i}),1));
            for k = 1:length(ME_bshell{i})
                R = sqrt((ME_bshell{i}(k)+eps)*(adc_hindered(j,1)-adc_hindered(j,2)));
                kernel_hindered{i,j}(k,1) = (sqrt(pi)*exp(-(ME_bshell{i}(k)+eps)*adc_hindered(j,2))*erf(R))/(2*R);
            end
        end

        for j = 1:num_isotropic
            kernel_isotropic{i,j} = single(zeros(length(ME_bshell{i}),1));
            for k = 1:length(ME_bshell{i})
                kernel_isotropic{i,j}(k,1) = exp(-(ME_bshell{i}(k)+eps)*adc_isotropic(j)); 
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


    %% optimization
    lb      = [10;10;10;zeros(num_restricted+num_hindered+num_isotropic,1)];
    ub      = [inf;inf;inf;ones(num_restricted+num_hindered+num_isotropic,1)];
    Aeq     = [zeros(3,3) blkdiag(ones(1,num_restricted),ones(1,num_hindered),ones(1,num_isotropic))];
    beq     = [1;1;1];

    options = optimoptions(@fmincon,'Display','off');
    alpha_coef = zeros(num_restricted+num_hindered+num_isotropic+3,size(ME_dwi_mean_array,2));
    alpha_init = [100;100;100;rand(num_restricted+num_hindered+num_isotropic,1)];

    parfor i = 1:size(ME_dwi_mean_array,2)
        beta = lsqnonneg(double(blk_kernel),double(ME_dwi_mean_array(:,i)));
        fun = @(x)echoFunc(x,beta,TE,num_restricted,num_hindered,num_isotropic);
        alpha_coef(:,i) = fmincon(fun,alpha_init,[],[],Aeq,beq,lb,ub,[],options);
    end

    temp2 = single(zeros(num_restricted+num_hindered+num_isotropic+3,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef;
    temp = single(zeros(num_restricted+num_hindered+num_isotropic+3,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;

    temp = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3),num_restricted+num_hindered+num_isotropic+3);
    niftiwrite(temp,'~/temp.nii','Compressed',true);
end

function FF = echoFunc(x,alpha,TE,num_restricted,num_hindered,num_isotropic)
    % x: [ T_2^h, T_2^r, T_2^f, w_1^h, w_2^h...w_1^r, w_2^r...,w_1^f, w_2^f...]
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




















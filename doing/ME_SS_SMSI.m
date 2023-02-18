function ME_SS_SMSI(dataFolder,TE,bshell)
    %{
    ░█▀▀█ ░█▀▀▀█ ░█▀▄▀█ ░█▀▀▀ ░█▀▀▄ ▀█▀
    ░█─── ░█──░█ ░█░█░█ ░█▀▀▀ ░█─░█ ░█─
    ░█▄▄█ ░█▄▄▄█ ░█──░█ ░█▄▄▄ ░█▄▄▀ ▄█▄

    Parameter estimation with multi-echo spherical mean spectrum Imaging
    (ME-SMSI) on single-shell dMRI


        Created by Ye Wu, PhD (dr.yewu@outlook.com)

        - Nanjing University of Science and Technology
        - University of North Carolina at Chapel Hill
        
    %}
    
    addpath('third/osqp');
    %% load multi-echo dMRI dataset
    fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
    fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
    fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');

    ME_dwi_info     = cellfun(@(x)niftiinfo(x),fdwi,'UniformOutput',false);
    ME_dwi          = cellfun(@(x)niftiread(x),ME_dwi_info,'UniformOutput',false);
    ME_bval         = cellfun(@(x)round(importdata(x)'/100)*100,fbval,'UniformOutput',false); 
    ME_mask_info    = niftiinfo(fmask);
    ME_mask         = round(niftiread(ME_mask_info));

    ME_dwi_mean     = single(zeros([size(ME_dwi{1},[1,2,3],length(TE))]));
    for i = 1:length(TE)
        ind = ME_bval{i} == bshell;
        ind_b0 = ME_bval{i} == 0;
        ME_dwi_mean(:,:,:,i) = mean(ME_dwi{i}(:,:,:,ind),4) ./ (eps+mean(ME_dwi{i}(:,:,:,ind_b0),4));
    end

    %% kernel
    adc_restricted = [];
    adc_hindered = [];
    adc_isotropic = (0 : 0.1e-3 : 3e-3)';

    for adc_long_diff = 0.5e-3 : 0.1e-3 : 1.5e-3
        for adc_trans_diff = 0e-3 : 0.1e-3 : 0.9e-3
            if adc_long_diff > adc_trans_diff * (pi/2)
                if adc_long_diff / adc_trans_diff >= (pi/2)^2
                    adc_restricted = [adc_restricted; adc_long_diff adc_trans_diff];
                else
                    adc_hindered = [adc_hindered; adc_long_diff adc_trans_diff];
                end
            end
        end
    end

    num_restricted  = size(adc_restricted,1);  kernel_restricted = zeros(length(TE),num_restricted);
    num_hindered    = size(adc_hindered,1);    kernel_hindered   = zeros(length(TE),num_hindered);
    num_isotropic   = size(adc_isotropic,1);   kernel_isotropic  = zeros(length(TE),num_isotropic);
    
    for i = 1:num_restricted
        R = sqrt(bshell*(adc_restricted(i,1)-adc_restricted(i,2)));
        kernel_restricted(:,i) = (sqrt(pi)*exp(-bshell*adc_restricted(i,2))*erf(R))/(2*R);
    end
    for i = 1:num_restricted
        R = sqrt(bshell*(adc_hindered(i,1)-adc_hindered(i,2)));
        kernel_hindered(:,i) = (sqrt(pi)*exp(-bshell*adc_hindered(i,2))*erf(R))/(2*R);
    end

    for i = 1:num_isotropic
        kernel_isotropic(:,i) = exp(-bshell*adc_isotropic(i)); 
    end

    blk_kernel_restricted = [];
    for i = 1:num_restricted
        blk_kernel_restricted = [blk_kernel_restricted diag(kernel_restricted(:,i))];
    end
    
    blk_kernel_hindered = [];
    for i = 1:num_hindered
        blk_kernel_hindered = [blk_kernel_hindered diag(kernel_hindered(:,i))];
    end

    blk_kernel_isotropic = [];
    for i = 1:num_isotropic
        blk_kernel_isotropic = [blk_kernel_isotropic diag(kernel_isotropic(:,i))];
    end
    
    blk_kernel = [blk_kernel_restricted blk_kernel_hindered blk_kernel_isotropic];


    %% Vectorization & Masked & arrayed
    ME_mask_ind = find(ME_mask>0.5);
    ME_dwi_mean_array = reshape(ME_dwi_mean,size(ME_dwi_mean,1)*size(ME_dwi_mean,2)*size(ME_dwi_mean,3),length(TE))';
    ME_dwi_mean_array_ind = all(ME_dwi_mean_array);
    ME_dwi_mean_array = ME_dwi_mean_array(:,ME_dwi_mean_array_ind);


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
    niftiwrite(T2r,fullfile(dataFolder,strcat('T2r_with_b',num2str(bshell),'.nii')),info_t2r,'Compressed', true);

    temp2 = single(zeros(num_restricted+num_hindered+num_isotropic,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = vf;
    temp = single(zeros(num_restricted+num_hindered+num_isotropic,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    vf = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3),num_restricted+num_hindered+num_isotropic);

    info_vf= ME_dwi_info{1};
    info_vf.ImageSize = size(vf);
    niftiwrite(vf,fullfile(dataFolder,strcat('VF_with_b',num2str(bshell),'.nii')),info_vf,'Compressed', true);

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





















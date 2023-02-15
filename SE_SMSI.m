function SE_SMSI(dataFolder,TE)
    %{
    ░█▀▀█ ░█▀▀▀█ ░█▀▄▀█ ░█▀▀▀ ░█▀▀▄ ▀█▀
    ░█─── ░█──░█ ░█░█░█ ░█▀▀▀ ░█─░█ ░█─
    ░█▄▄█ ░█▄▄▄█ ░█──░█ ░█▄▄▄ ░█▄▄▀ ▄█▄

    Parameter estimation with single-echo spherical mean spectrum Imaging (SE-SMSI)


        Created by Ye Wu, PhD (dr.yewu@outlook.com)

        - Nanjing University of Science and Technology
        - University of North Carolina at Chapel Hill
        
    %}
    
    addpath('third/osqp');
    %% load multi-echo dMRI dataset
    fdwi    = fullfile(dataFolder,strcat('MTE_',num2str(TE),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz'));
    fbval   = fullfile(dataFolder,strcat('MTE_',num2str(TE),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
    fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');

    SE_dwi_info     = niftiinfo(fdwi);
    SE_dwi          = niftiread(SE_dwi_info);
    SE_bval         = round(importdata(fbval)'/100)*100;
    SE_mask_info    = niftiinfo(fmask);
    SE_mask         = round(niftiread(SE_mask_info));
    
    SE_bshell       = unique(SE_bval);
    SE_dwi_mean     = single(zeros([size(SE_dwi,[1,2,3]),length(SE_bshell)]));
    for j = 1:length(SE_bshell)
        ind = find(SE_bval == SE_bshell(j));
        SE_dwi_mean(:,:,:,j) = mean(SE_dwi(:,:,:,ind),4);
    end


    %% Normalization S/S0
    SE_dwi_mean = SE_dwi_mean(:,:,:,2:end)./(SE_dwi_mean(:,:,:,1)+eps);
    SE_bshell   = SE_bshell(2:end);

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

    num_restricted  = size(adc_restricted,1);  kernel_restricted = cell(1,num_restricted);
    num_hindered    = size(adc_hindered,1);    kernel_hindered   = cell(1,num_hindered);
    num_isotropic   = size(adc_isotropic,1);   kernel_isotropic  = cell(1,num_isotropic);
    
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


    %% optimization 
    vf = zeros(num_restricted+num_hindered+num_isotropic,size(SE_dwi_mean_array,2));

    H = double(blk_kernel'*blk_kernel);
    A1 = [diag(ones(size(blk_kernel,2),1));ones(1,num_restricted+num_hindered+num_isotropic)];
    A2 = [zeros(size(blk_kernel,2),1);1];
    A3 = [ones(size(blk_kernel,2),1);1];
    parfor i = 1:size(SE_dwi_mean_array,2)
        f = -double(blk_kernel'*SE_dwi_mean_array(:,i));
        prob = osqp;
        prob.setup(H,f,A1,A2,A3,'alpha',0.1,'verbose',0);
        res = prob.solve();
        vf(:,i) = res.x;
    end

    temp2 = single(zeros(num_restricted+num_hindered+num_isotropic,length(SE_dwi_mean_array_ind)));
    temp2(:,SE_dwi_mean_array_ind) = vf;
    temp = single(zeros(num_restricted+num_hindered+num_isotropic,size(SE_mask,1)*size(SE_mask,2)*size(SE_mask,3)));
    temp(:,SE_mask_ind) = temp2;
    vf = reshape(temp',size(SE_mask,1),size(SE_mask,2),size(SE_mask,3),num_restricted+num_hindered+num_isotropic);

    info_vf= SE_dwi_info{1};
    info_vf.ImageSize = size(vf);
    niftiwrite(vf,fullfile(dataFolder,'VF_SE_',num2str(TE),'_norm.nii'),info_vf,'Compressed', true);

    % save options
    save(fullfile(dataFolder,'params_SE_',num2str(TE),'_norm.mat'),'num_restricted','adc_restricted',...
                                           'num_hindered','adc_hindered',...
                                           'num_isotropic','adc_isotropic',...
                                            '-v7.3');
end


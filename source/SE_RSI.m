function SE_RSI(dataFolder,TE)
    %{
    ░█▀▀█ ░█▀▀▀█ ░█▀▄▀█ ░█▀▀▀ ░█▀▀▄ ▀█▀
    ░█─── ░█──░█ ░█░█░█ ░█▀▀▀ ░█─░█ ░█─
    ░█▄▄█ ░█▄▄▄█ ░█──░█ ░█▄▄▄ ░█▄▄▀ ▄█▄

    Parameter estimation with single-echo restricted spectrum Imaging (SE-RSI)


        Created by Ye Wu, PhD (dr.yewu@outlook.com)

        - Nanjing University of Science and Technology
        - University of North Carolina at Chapel Hill
        
    %}
    
    addpath('third/osqp');
    addpath('third/csd');
    %% load multi-echo dMRI dataset
    fdwi    = fullfile(dataFolder,strcat('MTE_',num2str(TE),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz'));
    fbval   = fullfile(dataFolder,strcat('MTE_',num2str(TE),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval'));
    fbvec   = fullfile(dataFolder,strcat('MTE_',num2str(TE),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec'));
    fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');

    SE_dwi_info     = niftiinfo(fdwi);
    SE_dwi          = niftiread(SE_dwi_info);
    SE_bval         = round(importdata(fbval)'/100)*100;
    SE_bvec         = importdata(fbvec)';
    SE_mask_info    = niftiinfo(fmask);
    SE_mask         = round(niftiread(SE_mask_info));
    
    SE_bshell       = unique(SE_bval);

    %% Normalization S/S0
    ind_S0 = SE_bval == SE_bshell(1);
    SE_dwi = SE_dwi(:,:,:,~ind_S0) ./ (mean(SE_dwi(:,:,:,ind_S0),4)+eps);
    SE_bval = SE_bval(~ind_S0,:);
    SE_bvec = SE_bvec(~ind_S0,:);

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

    lmax = 6; 
    nmax = lmax2nsh(lmax);
    scheme = gen_scheme('sphere_362_vertices.txt',lmax);

    bval    = SE_bval;
    bvec    = SE_bvec;
    bshell  = unique(bval);
    nvol    = length(bval);
    
    for j = 1:num_restricted
        kernel_restricted{1,j} = zeros(nvol,nmax);
        for k = 1:length(bshell)
            order = floor(nsh2lmax(sum(bval==bshell(k))));
            DW_scheme = gen_scheme(bvec(bval==bshell(k),:),order);
            
            R_amp = response(adc_restricted(j,1),adc_restricted(j,2),bshell(k),scheme);
            R_SH = amp2SH(R_amp, scheme);
            R_RH = SH2RH(R_SH);

            m = [];
            for l = 0:2:order
                m = [ m R_RH(l/2+1)*ones(1,2*l+1) ];
            end
            
            fconv = DW_scheme.sh .* m(ones(size(DW_scheme.sh,1),1),:);
            fconv(:,end+1:nmax) = 0;
            kernel_restricted{1,j}(bval==bshell(k),:) = fconv;
            clear DW_scheme R_amp R_SH R_RH fconv m;
        end
    end
    
    for j = 1:num_hindered
        kernel_hindered{1,j} = zeros(nvol,nmax);
        for k = 1:length(bshell)
            order = floor(nsh2lmax(sum(bval==bshell(k))));
            DW_scheme = gen_scheme(bvec(bval==bshell(k),:),order);
            
            R_amp = response(adc_hindered(j,1),adc_hindered(j,2),bshell(k),scheme);
            R_SH = amp2SH(R_amp, scheme);
            R_RH = SH2RH(R_SH);

            m = [];
            for l = 0:2:order
                m = [ m R_RH(l/2+1)*ones(1,2*l+1) ];
            end
            
            fconv = DW_scheme.sh .* m(ones(size(DW_scheme.sh,1),1),:);
            fconv(:,end+1:nmax) = 0;
            kernel_hindered{1,j}(bval==bshell(k),:) = fconv;
            clear DW_scheme R_amp R_SH R_RH fconv m;
        end
    end

    for j = 1:num_isotropic
        kernel_isotropic{1,j} = single(zeros(length(bshell),1));
        for k = 1:length(bshell)
            kernel_isotropic{1,j}(k,1) = exp(-bshell(k)*adc_isotropic(j)); 
        end
    end

    blk_kernel = cell2mat([kernel_restricted kernel_hindered kernel_isotropic]);

    %% Vectorization & Masked & arrayed
    SE_mask_ind = find(SE_mask>0.5);
    SE_dwi_array = reshape(SE_dwi_mean,size(SE_dwi,1)*size(SE_dwi,2)*size(SE_dwi,3),size(SE_dwi,4));
    SE_dwi_array = SE_dwi_array(SE_mask_ind,:)';

    H = double(blk_kernel'*blk_kernel);
    f = double(blk_kernel'*SE_dwi_array);
    alpha_coef = H\f;

    temp = single(zeros(nmax,size(SE_mask,1)*size(SE_mask,2)*size(SE_mask,3)));
    info_fod = SE_dwi_info;
    info_fod.ImageSize(4) = nmax;
    for i = 1:num_restricted+num_hindered
        temp(:,SE_mask_ind) = alpha_coef((i-1)*nmax+1:i*nmax,:);
        fod = reshape(temp',size(SE_mask,1),size(SE_mask,2),size(SE_mask,3),nmax);
        if i <= num_restricted
            niftiwrite(fod,fullfile(dataFolder,strcat('FOD_SE_',num2str(TE),'_restricted_',num2str(i),'.nii')),info_fod,'Compressed', true);
        else
            niftiwrite(fod,fullfile(dataFolder,strcat('FOD_SE_',num2str(TE),'_hindered_',num2str(i-num_restricted),'.nii')),info_fod,'Compressed', true);
        end
    end

    temp = single(zeros(num_isotropic,size(SE_mask,1)*size(SE_mask,2)*size(SE_mask,3)));
    info_fod = SE_dwi_info;
    info_fod.ImageSize(4) = num_isotropic;
    temp(:,SE_mask_ind) = alpha_coef(nmax*(num_restricted+num_hindered)+1:end,:);
    fod = reshape(temp',size(SE_mask,1),size(SE_mask,2),size(SE_mask,3),num_isotropic);
    niftiwrite(fod,fullfile(dataFolder,strcat('FOD_SE_',num2str(TE),'_isotropic_',num2str(i-num_restricted),'.nii')),info_fod,'Compressed', true);

end
    
function nsh = lmax2nsh(lmax)
    nsh = (lmax+1) * (lmax+2) / 2;
end

function lmax = nsh2lmax(nsh)
    lmax = 2*(floor((sqrt(1+8*nsh)-3)/4));
end
 
function S = response(longitudinal,transverse,b,scheme)

    D = [ transverse 0 0; 0 transverse 0; 0 0 longitudinal ];
    C = s2c([ scheme.el scheme.az 1+0*scheme.az ]);
    X = C(:,1);
    Y = C(:,2);
    Z = C(:,3);

    S = exp(-b*[X.^2 Y.^2 Z.^2 2.*X.*Y 2.*X.*Z 2.*Y.*Z] * ...
            [ D(1,1) D(2,2) D(3,3) D(1,2) D(1,3) D(2,3) ]');
end


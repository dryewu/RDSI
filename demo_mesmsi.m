function demo_mesmsi()
    %{
    ░█▀▀█ ░█▀▀▀█ ░█▀▄▀█ ░█▀▀▀ ░█▀▀▄ ▀█▀
    ░█─── ░█──░█ ░█░█░█ ░█▀▀▀ ░█─░█ ░█─
    ░█▄▄█ ░█▄▄▄█ ░█──░█ ░█▄▄▄ ░█▄▄▀ ▄█▄

    Multi-echo Spectrum Imaging (MESI)


        Created by Ye Wu, PhD (dr.yewu@outlook.com)

        - Nanjing University of Science and Technology
        - University of North Carolina at Chapel Hill
        
    %}

    
    %% load multi-echo dMRI dataset
    dataFolder = '/Users/wuye/dryewu/Research/Project/MESI/Code/Data';
    TE = [75 85 95 105 115];
    
    ME_dwi = cell(1,length(TE));
    ME_dwi_mean = cell(1,length(TE));
    ME_info = cell(1,length(TE));
    ME_bval = cell(1,length(TE));
    ME_bvec = cell(1,length(TE));

    for i = 1:length(TE)
        filename = fullfile(dataFolder,strcat('MTE_PA_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased.nii.gz'));
        info = niftiinfo(filename);
        ME_dwi{1,i} = niftiread(info);
        ME_info{1,i} = info;

        filename = fullfile(dataFolder,strcat('MTE_PA_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased.bvec'));
        ME_bvec{1,i} = importdata(filename)';

        filename = fullfile(dataFolder,strcat('MTE_PA_',num2str(TE(i)),'_align_center_denoise_unring_preproc_epi_unbiased.bval'));
        ME_bval{1,i} = importdata(filename)';

        bvals = ME_bval{1,i}; 
        bvecs = ME_bvec{1,i}; 
        dwi = ME_dwi{1,i};

        bvals = round(bvals/100)*100;
        bshell = unique(bvals);

        for j = 1:length(bshell)
            ME_dwi{j,i} = dwi(:,:,:,bvals == bshell(j));
            ME_bval{j,i} = bshell(j);
            ME_bvec{j,i} = bvecs(bvals == bshell(j),:);
            ME_dwi_mean{j,i} = mean(ME_dwi{j,i},4);
        end
    end

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

    num_restricted  = size(adc_restricted,1);
    num_hindered    = size(adc_hindered,1);
    num_isotropic   = size(adc_isotropic,1);

    bshell = unique(cell2mat(ME_bval)) + eps;
    num_bshell = length(bshell);

    kernel_restricted   = cell(1,num_restricted);
    kernel_hindered     = cell(1,num_hindered);
    kernel_isotropic    = cell(1,num_isotropic);

    for i = 1:num_restricted
        kernel_restricted{1,i} = zeros(num_bshell,1);
        for j = 1:num_bshell
            R = sqrt(bshell(j)*(adc_restricted(i,1)-adc_restricted(i,2)));
            kernel_restricted{1,i}(j,1) = (sqrt(pi)*exp(-bshell(j)*adc_restricted(i,2))*erf(R))/(2*R);
        end
    end

    for i = 1:num_hindered
        kernel_hindered{1,i} = zeros(num_bshell,1);
        for j = 1:num_bshell
            R = sqrt(bshell(j)*(adc_hindered(i,1)-adc_hindered(i,2)));
            kernel_hindered{1,i}(j,1) = (sqrt(pi)*exp(-bshell(j)*adc_hindered(i,2))*erf(R))/(2*R);
        end
    end

    for i = 1:num_isotropic
        kernel_isotropic{1,i} = zeros(num_bshell,1);
        for j = 1:num_bshell
            kernel_isotropic{1,i}(j,1) = exp(-bshell(j)*adc_isotropic(i)); 
        end
    end

    kernel = [cell2mat(kernel_restricted) cell2mat(kernel_hindered) cell2mat(kernel_isotropic)];

    
    % masked dwi
    filename = fullfile(dataFolder,'MTE_mask.nii.gz');
    info = niftiinfo(filename);
    mask = round(niftiread(info));
    ind_mask = find(mask>0.5);
    ME_dwi_mean_array = single(zeros([size(ME_dwi_mean{1}) numel(ME_dwi_mean)]));
    num = 1;
    for i = 1:size(ME_dwi_mean,2)
        for j = 1:size(ME_dwi_mean,1)
            ME_dwi_mean_array(:,:,:,num) = ME_dwi_mean{j,i};
        end
    end
    [nx,ny,nz,ns] = size(ME_dwi_mean_array);
    ME_dwi_mean_array = reshape(ME_dwi_mean_array,nx*ny*nz,ns);
    ME_dwi_mean_array = ME_dwi_mean_array(ind_mask,:)';


end

    %% construct spherical mean kernel (not rely on TE)

% 
%     
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%     
% 
% 
% 
%                 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% params = [];
% params.lambda = 0.1;                % regularization constant
% params.SH_order = 4;                % spherical harmonic order -- must be even
% params.norm_flag = 0;               % normalize data to b=0 image
% params.nonlin_flag = 0;             % use nonlinear optimization with initial parameters from linear fit
% params.b0_thresh = 10;              % threshold used for considering a b-value to be 0
% params.scalefacts_flag = 0;         % calculate scaling factors from b=0 images and apply them to all subsequent frames (for multiple acquisitions)
% 
% params.ADC_long = [];               % longitudinal ADC
% params.ADC_trans = [];              % transverse ADC
% num = 0;
% for ADC_long = 0.5e-3 : 0.2e-3 : 1.5e-3
%     for ADC_trans = 0e-3 : 0.2e-3 : 0.9e-3
%         if ADC_trans < ADC_long
%             params.ADC_long = [params.ADC_long ADC_long];
%             params.ADC_trans = [params.ADC_trans ADC_trans];
%             num = num + 1;
%         end
%     end
% end
% params.ADC_iso = 0:0.2e-3:3e-3;     % minimum isotropic ADC
% params.num_ADC_aniso = num;         % number of transverse ADC size scales
% params.num_ADC_iso = length(params.ADC_iso); % number of isotropic ADC size scales
% 
% %% load data
% dwi_info = niftiinfo(dwi_filename);
% dwi = niftiread(dwi_info);
% 
% mask_info = niftiinfo(mask_filename);
% mask = niftiread(mask_info);
% 
% bvals = importdata(bval_filename);
% bvecs = importdata(bvec_filename);
% 
% if size(bvecs,2) ~=3; bvecs = bvecs'; end
% if size(bvals,2) ~=1; bvals = bvals'; end
% 
% bvals = round(bvals/params.b0_thresh)*params.b0_thresh;
% 
% %% data header
% params.numvol = length(dwi);
% params.voxsiz = mask_info.PixelDimensions;
% params.voxdim = size(mask);
% params.i_b0 = find(bvals <= params.b0_thres);
% params.i_dwi = find(bvals > params.b0_thresh);
% 
% %% create forward matrix for fitting tensor; construct tensor forward matrix
% [params.B,params.Binv] = create_tensor_matrix(bvals,bvecs);
% 
% %% construct RSI multi-FOD forward matrix
% [params.A,params.Ainv,params.icoverts,params.beta2ico] = create_rsi_matrix(MFfit);
% 
% % linear RSI fit
% MFfit = fit_rsi(vol,MFfit,parms);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [B,Binv] = create_tensor_matrix(bvals,bvecs)
%   ndirs = length(bvals);
%   B = zeros(ndirs,7);
%   B(:,7) = 1; % S0
%   for i = 1:ndirs
%     outerprod = -bvals(i).*bvecs(i,:)'*bvecs(i,:);
%     B(i,1:3) = diag(outerprod)';
%     B(i,4) = 2*outerprod(1,2);
%     B(i,5) = 2*outerprod(1,3);
%     B(i,6) = 2*outerprod(2,3);
%   end
%   Binv = pinv(B);
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [A,Ainv,V,F] = create_rsi_matrix(bvals,bvecs,params)
%   V = importdata('sphere_362_vertices.txt');
%   Q = bvecs.*repmat(sqrt(bvals),1,3);
%   F = rsi_SH_matrix(V,params.SH_order);
%   F0 = rsi_SH_matrix(V,0);
%   A = [];
% 
%   % series of anisotropic FODs for varying size scales
%   for i = 1:params.num_ADC_aniso
%     for j = 1:params.num_ADC_aniso
%         R = rsi_FOD_matrix(Q,V,params.ADC_long(i),params.ADC_trans(j));
%         A = [A R*F];
%     end
%   end
% 
%   % series of isotropic FODs for varying size scales
%   for t=1:params.num_ADC_iso
%     R = rsi_FOD_matrix(Q,V,params.ADC_iso(t),params.ADC_iso(t));
%     A = [A R*F0];
%   end
% 
%   % compute regularized inverse
%   AtA = A'*A;
%   Ainv = (AtA+params.lambda*mean(diag(AtA))*eye(size(AtA)))\A';
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function MFfit = fit_rsi(vol,parms)
%   % initialize multi-FOD parameter volumes
%   MFfit.volMF = zeros(parms.nx,parms.ny,parms.nz,MFfit.nb);
%   % loop over slices
%   for z = 1:parms.nz
%     y = reshape(vol(:,:,z,:),[parms.nx*parms.ny,parms.nf])';
%     if parms.norm_flag % normalize data to b=0
%       y0 = repmat(reshape(MFfit.volb0(:,:,z),[parms.nx*parms.ny,1])',[parms.nf 1]);
%       y = y./y0;
%     end;
%     betas = MFfit.Ainv*y;
%     betas = reshape(betas',[parms.nx parms.ny MFfit.nb]);
%     betas(isnan(betas)) = 0;
%     MFfit.volMF(:,:,z,:) = betas;
%   end
%   % dilate brain mask
%   MFfit.volmask_dilated = dilate_brain_mask(MFfit.volmask,parms);
% return;
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% bList = unique(bVals);
% bNums = length(bList);
% gNums = length(bVals);
% b0Index = bVals == 0;
% 
% bIndex = cell(1,bNums);
% gIndex = cell(1,bNums);
% for idx = 1:bNums
%     bIndex{idx} = find(bVals==bList(idx));
%     gIndex{idx} = g_data(bIndex{idx},:);
% end























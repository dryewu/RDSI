function ME_SS_SI(dataFolder,TE,bshell)
    %{
    ░█▀▀█ ░█▀▀▀█ ░█▀▄▀█ ░█▀▀▀ ░█▀▀▄ ▀█▀
    ░█─── ░█──░█ ░█░█░█ ░█▀▀▀ ░█─░█ ░█─
    ░█▄▄█ ░█▄▄▄█ ░█──░█ ░█▄▄▄ ░█▄▄▀ ▄█▄

    Parameter estimation with multi-echo spectrum Imaging (ME-SI)
    (ME-SMSI) on single-shell dMRI


        Created by Ye Wu, PhD (dr.yewu@outlook.com)

        - Nanjing University of Science and Technology
        - University of North Carolina at Chapel Hill
        
    %}
    
    addpath('third/osqp');
    addpath('third/csd');
    %% load multi-echo dMRI dataset
    fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
    fbvec   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec')),TE,'UniformOutput',false);
    fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
    fmask   = fullfile(dataFolder,'MTE_mask.nii.gz');
    fvf     = fullfile(dataFolder,strcat('VF_with_b',num2str(bshell),'.nii'));     % estimated by ME_SS_SMSI
    ft2r    = fullfile(dataFolder,strcat('T2r_with_b',num2str(bshell),'.nii'));    % estimated by ME_SS_SMSI
    fparams = fullfile(dataFolder,strcat('param_with_b',num2str(bshell),'.mat'));    % used in ME_SS_SMSI

    ME_dwi_info     = cellfun(@(x)niftiinfo(x),fdwi,'UniformOutput',false);
    ME_dwi          = cellfun(@(x)niftiread(x),ME_dwi_info,'UniformOutput',false);
    ME_bval         = cellfun(@(x)round(importdata(x)'/100)*100,fbval,'UniformOutput',false); 
    ME_bvec         = cellfun(@(x)importdata(x)',fbvec,'UniformOutput',false);
    ME_mask_info    = niftiinfo(fmask);
    ME_vf_info      = niftiinfo(fvf);
    ME_t2r_info     = niftiinfo(ft2r);
    ME_mask         = round(niftiread(ME_mask_info));
    ME_vf           = niftiread(ME_vf_info);
    ME_t2r          = niftiread(ME_t2r_info);
    params = load(fparams);
    
    %% Normalization S/S0
    ME_dwi_norm     = cell(size(ME_dwi));
    ME_bval_norm    = cell(size(ME_bval));
    ME_bvec_norm    = cell(size(ME_bvec));

    for i = 1:length(TE)
        ind = ME_bval{i} == bshell;
        ind_b0 = ME_bval{i} == 0;
        ME_dwi_norm{i} = ME_dwi{i}(:,:,:,ind) ./ (mean(ME_dwi{i}(:,:,:,ind_b0),4)+eps);

        ME_bval_norm{i} = ME_bval{i}(ind,:);
        ME_bvec_norm{i} = ME_bvec{i}(ind,:);
    end

    ME_dwi = ME_dwi_norm;    clear ME_dwi_norm;
    ME_bval = ME_bval_norm;  clear ME_bval_norm;
    ME_bvec = ME_bvec_norm;  clear ME_bvec_norm;
    
    %% kernel
    adc_restricted      = params.adc_restricted;
    adc_hindered        = params.adc_hindered;
    adc_isotropic       = params.adc_isotropic;

    num_restricted      = params.num_restricted;  
    num_hindered        = params.num_hindered;    
    num_isotropic       = params.num_isotropic;   
    
    kernel_restricted   = cell(length(TE),num_restricted);
    kernel_hindered     = cell(length(TE),num_hindered);
    kernel_isotropic    = cell(length(TE),num_isotropic);

    lmax = 6; 
    nmax = lmax2nsh(lmax);
    scheme = gen_scheme('sphere_362_vertices.txt',lmax);

    for i = 1:length(TE)
        bval    = ME_bval{i};
        bvec    = ME_bvec{i};
        bshell  = unique(bval);
        
        for j = 1:num_restricted 
            order = floor(nsh2lmax(sum(bval==bshell)));
            DW_scheme = gen_scheme(bvec(bval==bshell,:),order);
            
            R_amp = response(adc_restricted(j,1),adc_restricted(j,2),bshell,scheme);
            R_SH = amp2SH(R_amp, scheme);
            R_RH = SH2RH(R_SH);

            m = [];
            for l = 0:2:order
                m = [ m R_RH(l/2+1)*ones(1,2*l+1) ];
            end
            
            fconv = DW_scheme.sh .* m(ones(size(DW_scheme.sh,1),1),:);
            fconv(:,end+1:nmax) = 0;
            kernel_restricted{i,j} = fconv;
            clear DW_scheme R_amp R_SH R_RH fconv m;

        end
        
        for j = 1:num_hindered
            order = floor(nsh2lmax(sum(bval==bshell)));
            DW_scheme = gen_scheme(bvec(bval==bshell,:),order);
            
            R_amp = response(adc_hindered(j,1),adc_hindered(j,2),bshell,scheme);
            R_SH = amp2SH(R_amp, scheme);
            R_RH = SH2RH(R_SH);

            m = [];
            for l = 0:2:order
                m = [ m R_RH(l/2+1)*ones(1,2*l+1) ];
            end
            
            fconv = DW_scheme.sh .* m(ones(size(DW_scheme.sh,1),1),:);
            fconv(:,end+1:nmax) = 0;
            kernel_hindered{i,j} = fconv;
            clear DW_scheme R_amp R_SH R_RH fconv m;
        end

        for j = 1:num_isotropic
            kernel_isotropic{i,j} = exp(-bshell*adc_isotropic(j)); 
        end
    end  

    %% Vectorization & Masked & arrayed
    ME_mask_ind = find(ME_mask>0.5);
    ME_dwi_array = cellfun(@(x)reshape(x,size(x,1)*size(x,2)*size(x,3),size(x,4)),ME_dwi,'UniformOutput',false);
    ME_dwi_array = cellfun(@(x)x(ME_mask_ind,:)',ME_dwi_array,'UniformOutput',false);
    ME_dwi_array = cell2mat(ME_dwi_array');

    ME_t2r_array = reshape(ME_t2r,size(ME_t2r,1)*size(ME_t2r,2)*size(ME_t2r,3),3);
    ME_t2r_array = ME_t2r_array(ME_mask_ind,:)';

    ME_vf_array = reshape(ME_vf,size(ME_vf,1)*size(ME_vf,2)*size(ME_vf,3),num_restricted + num_hindered + num_isotropic);
    ME_vf_array = ME_vf_array(ME_mask_ind,:)';

    %% remove isotropic component
    ME_vf_array = ME_vf_array./sum(ME_vf_array);
    ME_vf_array(isnan(ME_vf_array)) = 0;
    ME_vf_array_isotropic = ME_vf_array(num_restricted+num_hindered+1:end,:);
    ME_vf_array_restricted = ME_vf_array(1:num_restricted,:);
    ME_vf_array_hindered = ME_vf_array(1+num_restricted:num_restricted+num_hindered,:);
    ME_t2r_array_restricted = ME_t2r_array(1,:);
    ME_t2r_array_hindered = ME_t2r_array(2,:);
    ME_t2r_array_isotropic = ME_t2r_array(3,:);
    ME_dwi_array_anisotropic = ME_dwi_array;

    num_start = 1;
    for i = 1:length(TE)
        num_end = num_start + length(ME_bval{i}) - 1;
        dwi = ME_dwi_array_anisotropic(num_start:num_end,:);
        bval    = ME_bval{i};
        bshell  = unique(bval);

        for j = 1:num_isotropic
            for k = 1:length(bshell)
                dwi(bval==bshell,:) = dwi(bval==bshell,:) - exp(-TE(i)./ME_t2r_array_isotropic).*ME_vf_array_isotropic(j,:).*kernel_isotropic{i,j}(k,1);
            end
        end
        ME_dwi_array_anisotropic(num_start:num_end,:) = dwi;
        num_start = num_end + 1;
    end

    %% subject to
    nv = size(scheme.vert,1);
    ampbasis = repmat(scheme.sh,1,num_restricted + num_hindered);
    ampbasis = mat2cell(ampbasis,nv,repmat(nmax,1,num_restricted + num_hindered));
    ampbasis = blkdiag(ampbasis{:});

    B = repmat(ones(1,nv),1,num_restricted + num_hindered);
    B = mat2cell(B,1,repmat(nv,1,num_restricted + num_hindered));
    B = blkdiag(B{:});

    A1 = [ampbasis; B*ampbasis];
    A2 = [zeros(size(ampbasis,1),1); ones(num_restricted + num_hindered,1)];
    A3 = [inf(size(ampbasis,1),1); ones(num_restricted + num_hindered,1)];
    alpha_coef = zeros(num_restricted*nmax+num_hindered*nmax,size(ME_vf_array,2));

    ME_t2r_restricted  = exp(-TE'./ME_t2r_array_restricted);
    ME_t2r_hindered    = exp(-TE'./ME_t2r_array_hindered);

    %% optimization
    parfor i = 1:size(ME_vf_array,2)
        kernel = cell2mat([ cellfun(@(x,y) x.*y, kernel_restricted,num2cell(ME_t2r_restricted(:,i).*ME_vf_array_restricted(:,i)'), 'UniformOutput',false) ...
                            cellfun(@(x,y) x.*y, kernel_hindered,num2cell(ME_t2r_hindered(:,i).* ME_vf_array_hindered(:,i)'), 'UniformOutput',false)]);
        
        dwi = ME_dwi_array_anisotropic(:,i);
        ind = dwi > 0 & ~isnan(dwi) & ~isinf(dwi);

        if sum(ind) < 25
            continue;
        end

        try
            H = double(kernel(ind,:)'*kernel(ind,:));
            f = -double(kernel(ind,:)'*dwi(ind,1));
    
            prob = osqp;
            prob.setup(H,f,A1,A2,A3,'alpha',0.1,'verbose',0);
            res = prob.solve();
            alpha_coef(:,i) = res.x;
        catch
            continue;
        end
    end

    temp = single(zeros(nmax,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    info_fod = ME_dwi_info{1};
    info_fod.ImageSize(4) = nmax;
    for i = 1:num_restricted+num_hindered
        temp(:,ME_mask_ind) = alpha_coef((i-1)*nmax+1:i*nmax,:);
        fod = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3),nmax);
        if i <= num_restricted
            niftiwrite(fod,fullfile(dataFolder,strcat('FOD_restricted_',num2str(i),'.nii')),info_fod,'Compressed', true);
        else
            niftiwrite(fod,fullfile(dataFolder,strcat('FOD_hindered_',num2str(i-num_restricted),'.nii')),info_fod,'Compressed', true);
        end
    end
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


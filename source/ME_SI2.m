% function ME_SI2(dataFolder,TE)
function ME_SI2()
    %{
    ░█▀▀█ ░█▀▀▀█ ░█▀▄▀█ ░█▀▀▀ ░█▀▀▄ ▀█▀
    ░█─── ░█──░█ ░█░█░█ ░█▀▀▀ ░█─░█ ░█─
    ░█▄▄█ ░█▄▄▄█ ░█──░█ ░█▄▄▄ ░█▄▄▀ ▄█▄

    Parameter estimation with multi-echo spectrum Imaging (ME-SI)


        Created by Ye Wu, PhD (dr.yewu@outlook.com)

        - Nanjing University of Science and Technology
        - University of North Carolina at Chapel Hill
        
    %}
    
    addpath('third/osqp');
    addpath('third/csd');
    %% load multi-echo dMRI dataset
    dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
    TE = [75 85 95 105 115];
    fdwi    = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz')),TE,'UniformOutput',false);
    fbvec   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec')),TE,'UniformOutput',false);
    fbval   = arrayfun(@(x)fullfile(dataFolder,strcat('MTE_',num2str(x),'_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval')),TE,'UniformOutput',false);
    fmask   = fullfile(dataFolder,'MTE_mask_ROI.nii.gz');
    fvf     = fullfile(dataFolder,'VF_norm.nii.gz');     % estimated by ME_SMSI
    ft2r    = fullfile(dataFolder,'T2r_norm.nii.gz');    % estimated by ME_SMSI
    fparams = fullfile(dataFolder,'params_norm.mat');    % used in ME_SMSI

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
    
    ME_bshell       = cellfun(@(x)unique(x),ME_bval,'UniformOutput',false); 
    clear fdwi fbvec fbval fmask fvf ft2r fparams;
    clear ME_vf_info ME_t2r_info ME_mask_info

    %% Normalization S/S0
    ME_dwi_norm     = cell(size(ME_dwi));
    ME_bval_norm    = cell(size(ME_bval));
    ME_bvec_norm    = cell(size(ME_bvec));

    for i = 1:length(TE)
        ind_S0 = ME_bval{i} == ME_bshell{i}(1);
        ME_dwi_norm{i} = ME_dwi{i}(:,:,:,~ind_S0) ./ (mean(ME_dwi{i}(:,:,:,ind_S0),4)+eps);

        ME_bval_norm{i} = ME_bval{i}(~ind_S0,:);
        ME_bvec_norm{i} = ME_bvec{i}(~ind_S0,:);
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
    scheme = gen_scheme('sphere_362_vertices_yz.txt',lmax);

    for i = 1:length(TE)
        bval    = ME_bval{i};
        bvec    = ME_bvec{i};
        bshell  = unique(bval);
        nvol    = length(bval);
        
        for j = 1:num_restricted
            kernel_restricted{i,j} = zeros(nvol,nmax);
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
                kernel_restricted{i,j}(bval==bshell(k),:) = fconv;
                clear DW_scheme R_amp R_SH R_RH fconv m;
            end
        end
        
        for j = 1:num_hindered
            kernel_hindered{i,j} = zeros(nvol,nmax);
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
                kernel_hindered{i,j}(bval==bshell(k),:) = fconv;
                clear DW_scheme R_amp R_SH R_RH fconv m;
            end
        end

        for j = 1:num_isotropic
            kernel_isotropic{i,j} = zeros(nvol,1);
            for k = 1:length(bshell)
                kernel_isotropic{i,j}(bval==bshell(k),1) = exp(-bshell(k)*adc_isotropic(j)); 
            end
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

    clear params ME_dwi ME_t2r ME_vf

    %% remove isotropic component
    ME_vf_array = ME_vf_array./sum(ME_vf_array);
    ME_vf_array(isnan(ME_vf_array)) = 0;
    ME_vf_array_isotropic = ME_vf_array(num_restricted+num_hindered+1:end,:);
    ME_vf_array_restricted = ME_vf_array(1:num_restricted,:);
    ME_vf_array_hindered = ME_vf_array(1+num_restricted:num_restricted+num_hindered,:);
    ME_t2r_array_restricted = ME_t2r_array(1,:);
    ME_t2r_array_hindered = ME_t2r_array(2,:);
    ME_t2r_array_isotropic = ME_t2r_array(3,:);

    %% subject to
    nv = size(scheme.vert,1);
    ampbasis = repmat(scheme.sh,1,num_restricted + num_hindered);
    ampbasis = mat2cell(ampbasis,nv,repmat(nmax,1,num_restricted + num_hindered));
    ampbasis = blkdiag(ampbasis{:});
    A1 = blkdiag(ampbasis,diag(ones(1,num_isotropic)));

    A2 = zeros(size(A1,1),1);
    A3 = ones(size(A1,1),1);
    alpha_coef = zeros(num_restricted*nmax+num_hindered*nmax+num_isotropic,size(ME_vf_array,2));

    ME_t2r_restricted  = exp(-TE'./ME_t2r_array_restricted);
    ME_t2r_hindered    = exp(-TE'./ME_t2r_array_hindered);
    ME_t2r_isotropic    = exp(-TE'./ME_t2r_array_isotropic);
    clear ampbasis ME_t2r_array_restricted ME_t2r_array_hindered B;

    %% optimization
    parfor i = 1:size(ME_vf_array_restricted,2)
        kernel = cell2mat([ cellfun(@(x,y) x.*y, kernel_restricted,num2cell(ME_t2r_restricted(:,i).*ones(1,num_restricted)), 'UniformOutput',false) ...
                            cellfun(@(x,y) x.*y, kernel_hindered,num2cell(ME_t2r_hindered(:,i).* ones(1,num_hindered)), 'UniformOutput',false) ...
                            cellfun(@(x,y) x.*y, kernel_isotropic,num2cell(ME_t2r_isotropic(:,i).* ones(1,num_isotropic)), 'UniformOutput',false)]);

        dwi = ME_dwi_array(:,i);
        ind = ~isnan(dwi) & ~isinf(dwi);
% 
%         if sum(ind) < 45
%             continue;
%         end

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

    temp = single(zeros(num_isotropic,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    info_fod = ME_dwi_info{1};
    info_fod.ImageSize(4) = num_isotropic;
    temp(:,ME_mask_ind) = alpha_coef(end-num_isotropic+1:end,:);
    fod = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3),num_isotropic);

    niftiwrite(fod,fullfile(dataFolder,'FOD_isotropic.nii'),info_fod,'Compressed', true);

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


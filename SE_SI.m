function SE_SI(fdwi,fbvec,fbval,fmask,outpath,options)
    %{
    ░█▀▀█ ░█▀▀▀█ ░█▀▄▀█ ░█▀▀▀ ░█▀▀▄ ▀█▀
    ░█─── ░█──░█ ░█░█░█ ░█▀▀▀ ░█─░█ ░█─
    ░█▄▄█ ░█▄▄▄█ ░█──░█ ░█▄▄▄ ░█▄▄▀ ▄█▄

    Single-echo spectrum imaging (SE-SI)


        Created by Ye Wu, PhD (dr.yewu@outlook.com)

        - Nanjing University of Science and Technology, China
        - University of North Carolina at Chapel Hill, USA
        
    %}

    arguments
        fdwi                    (1,:)    {mustBeFile}
        fbvec                   (1,:)    {mustBeFile}
        fbval                   (1,:)    {mustBeFile}
        fmask                   string   {mustBeFile}
        outpath                 string   {mustBeNonzeroLengthText}

        options.normalizeToS0   (1,1)    {mustBeNumericOrLogical} = true
        options.useBshell       (1,:)    {mustBeNumericOrLogical} = false
        options.spectrum        string   {mustBeFile} = 'scheme/default_spectrum.mat'
    end

    addpath('third/osqp');
    addpath('third/csd');

    assert(exist(fdwi,'file'),'Input mask %s does not exist', fdwi);
    assert(exist(fbvec,'file'),'Input mask %s does not exist', fbvec);
    assert(exist(fbval,'file'),'Input mask %s does not exist', fbval);
    assert(exist(fmask,'file'),'Input mask %s does not exist', fmask);
    assert(exist(options.spectrum,'file'),'Input spectrum %s does not exist', options.spectrum);

    %% load single-echo dMRI dataset
    SE_dwi_info     = niftiinfo(fdwi);
    SE_dwi          = niftiread(SE_dwi_info);
    SE_bvec         = importdata(fbvec)';
    SE_bval         = round(importdata(fbval)'/100)*100;
    SE_mask_info    = niftiinfo(fmask);
    SE_mask         = round(niftiread(SE_mask_info));
    clear fdwi fbvec fbval fmask SE_mask_info;

    if options.useBshell
        ind         = ismember(fbval,options.useBshell);
        SE_dwi      = SE_dwi(:,:,:,ind);
        SE_bvec     = SE_bvec(ind,:);
        SE_bval     = SE_bval(ind);
        clear ind
    end

    %% Normalization S/S0
    if options.normalizeToS0
        ind             = ismember(SE_bval,0); 
        SE_dwi_norm     = SE_dwi(:,:,:,~ind)./(mean(SE_dwi(:,:,:,~ind),4)+eps);
        SE_bval_norm    = SE_bval(~ind,:);
        SE_bvec_norm    = SE_bvec(~ind,:);
        clear ind
    
        SE_dwi = SE_dwi_norm;    clear SE_dwi_norm;
        SE_bval = SE_bval_norm;  clear SE_bval_norm;
        SE_bvec = SE_bvec_norm;  clear SE_bvec_norm;
    end
    
    %% kernel
    default_spectrum = load(options.spectrum);

    adc_restricted   = default_spectrum.adc_restricted;
    adc_hindered     = default_spectrum.adc_hindered;
    adc_isotropic    = default_spectrum.adc_isotropic;

    num_restricted  = size(adc_restricted,1);  
    num_hindered    = size(adc_hindered,1);    
    num_isotropic   = size(adc_isotropic,1);   
    
    kernel_restricted = cell(1,num_restricted);
    kernel_hindered   = cell(1,num_hindered);
    kernel_isotropic  = cell(1,num_isotropic);

    lmax = 6; 
    nmax = lmax2nsh(lmax);
    scheme = gen_scheme('scheme/sphere_362_vertices.txt',lmax);

    bshell  = unique(SE_bval);
    nvol    = length(SE_bval);
    
    for j = 1:num_restricted
        kernel_restricted{1,j} = zeros(nvol,nmax);
        for k = 1:length(bshell)
            order = min(floor(nsh2lmax(sum(SE_bval==bshell(k)))),lmax);
            DW_scheme = gen_scheme(SE_bvec(SE_bval==bshell(k),:),order);
            
            R_amp = response(adc_restricted(j,1),adc_restricted(j,2),bshell(k),scheme);
            R_SH = amp2SH(R_amp, scheme);
            R_RH = SH2RH(R_SH);

            m = [];
            for l = 0:2:order
                m = [ m R_RH(l/2+1)*ones(1,2*l+1) ];
            end
            
            fconv = DW_scheme.sh .* m(ones(size(DW_scheme.sh,1),1),:);
            fconv(:,end+1:nmax) = 0;
            kernel_restricted{1,j}(SE_bval==bshell(k),:) = fconv;
            clear DW_scheme R_amp R_SH R_RH fconv m;
        end
    end
    
    for j = 1:num_hindered
        kernel_hindered{1,j} = zeros(nvol,nmax);
        for k = 1:length(bshell)
            order = min(floor(nsh2lmax(sum(SE_bval==bshell(k)))),lmax);
            DW_scheme = gen_scheme(SE_bvec(SE_bval==bshell(k),:),order);
            
            R_amp = response(adc_hindered(j,1),adc_hindered(j,2),bshell(k),scheme);
            R_SH = amp2SH(R_amp, scheme);
            R_RH = SH2RH(R_SH);

            m = [];
            for l = 0:2:order
                m = [ m R_RH(l/2+1)*ones(1,2*l+1) ];
            end
            
            fconv = DW_scheme.sh .* m(ones(size(DW_scheme.sh,1),1),:);
            fconv(:,end+1:nmax) = 0;
            kernel_hindered{1,j}(SE_bval==bshell(k),:) = fconv;
            clear DW_scheme R_amp R_SH R_RH fconv m;
        end
    end

    for j = 1:num_isotropic
        kernel_isotropic{1,j} = zeros(nvol,1);
        for k = 1:length(bshell)
            kernel_isotropic{1,j}(SE_bval==bshell(k),1) = exp(-bshell(k)*adc_isotropic(j)); 
        end
    end
    
    blk_kernel = cell2mat([kernel_restricted kernel_hindered kernel_isotropic]);

    
    %% Vectorization & Masked & arrayed
    SE_mask_ind = find(SE_mask>0.5);
    SE_dwi_array = reshape(SE_dwi,size(SE_dwi,1)*size(SE_dwi,2)*size(SE_dwi,3),size(SE_dwi,4));
    SE_dwi_array = SE_dwi_array(SE_mask_ind,:)';
    SE_dwi_array_ind = all(SE_dwi_array);
    SE_dwi_array = SE_dwi_array(:,SE_dwi_array_ind);
    clear SE_dwi

    %% subject to
    nv = size(scheme.vert,1);
    ampbasis = repmat(scheme.sh,1,num_restricted + num_hindered);
    ampbasis = mat2cell(ampbasis,nv,repmat(nmax,1,num_restricted + num_hindered));
    ampbasis = blkdiag(ampbasis{:});

    H  = double(blk_kernel'*blk_kernel);
    K  = double(-blk_kernel'*SE_dwi_array);
    A1 = blkdiag(ampbasis,diag(ones(1,num_isotropic)));
    A2 = zeros(size(A1,1),1);
    A3 = ones(size(A1,1),1);

    alpha_coef = zeros(num_restricted*nmax+num_hindered*nmax+num_isotropic,size(SE_dwi_array,2));
    clear ampbasis;

    %% optimization
    parfor i = 1:size(SE_dwi_array,2)

            f = K(:,i);
            prob = osqp;
            prob.setup(H,f,A1,A2,A3,'alpha',0.1,'verbose',0);
            res = prob.solve();
            alpha_coef(:,i) = res.x;
    end

    %% save results
    if ~exist(outpath,'dir')
        mkdir(outpath)
    end

    % save FOD restricted
    temp = single(zeros(nmax,size(SE_mask,1)*size(SE_mask,2)*size(SE_mask,3)));
    info_fod = SE_dwi_info;
    info_fod.Datatype = 'single';
    info_fod.ImageSize(4) = nmax;
    for i = 1:num_restricted+num_hindered
        temp(:,SE_mask_ind) = alpha_coef((i-1)*nmax+1:i*nmax,:);
        fod = reshape(temp',size(SE_mask,1),size(SE_mask,2),size(SE_mask,3),nmax);
        if i <= num_restricted
            niftiwrite(single(fod),fullfile(outpath,strcat('FOD_restricted_',num2str(i),'.nii')),info_fod,'Compressed', true);
        else
            niftiwrite(single(fod),fullfile(outpath,strcat('FOD_hindered_',num2str(i-num_restricted),'.nii')),info_fod,'Compressed', true);
        end
    end

    temp = single(zeros(num_isotropic,size(SE_mask,1)*size(SE_mask,2)*size(SE_mask,3)));
    info_fod.Datatype = 'single';
    info_fod.ImageSize(4) = num_isotropic;
    temp(:,SE_mask_ind) = alpha_coef(1+(num_restricted+num_hindered)*nmax:end,:)/sqrt(4*pi);
    fod = reshape(temp',size(SE_mask,1),size(SE_mask,2),size(SE_mask,3),num_isotropic);

    niftiwrite(fod,fullfile(outpath,'FOD_free.nii'),info_fod,'Compressed', true);
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


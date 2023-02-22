function ME_SI(fdwi,fbvec,fbval,ft2r,fmask,TE,outpath,options)
    %{
    ░█▀▀█ ░█▀▀▀█ ░█▀▄▀█ ░█▀▀▀ ░█▀▀▄ ▀█▀
    ░█─── ░█──░█ ░█░█░█ ░█▀▀▀ ░█─░█ ░█─
    ░█▄▄█ ░█▄▄▄█ ░█──░█ ░█▄▄▄ ░█▄▄▀ ▄█▄

    Parameter estimation with multi-echo spectrum Imaging (ME-SI)


        Created by Ye Wu, PhD (dr.yewu@outlook.com)

        - Nanjing University of Science and Technology
        - University of North Carolina at Chapel Hill
        
    %}

    arguments
        fdwi                    (1,:)    {mustBeNonzeroLengthText}
        fbvec                   (1,:)    {mustBeNonzeroLengthText}
        fbval                   (1,:)    {mustBeNonzeroLengthText}
        ft2r                    (1,:)    {mustBeNonzeroLengthText} % [restricted, hindered, free]
        fmask                   string   {mustBeFile}
        TE                      (1,:)    {mustBeNumeric}
        outpath                 string   {mustBeNonzeroLengthText}

        options.normalizeToS0   (1,1)    {mustBeNumericOrLogical} = true
        options.useBshell       (1,:)    {mustBeNumericOrLogical} = false
        options.spectrum        string   {mustBeFile} = 'scheme/default_spectrum.mat'
    end

    addpath('third/osqp');
    addpath('third/csd');

    cellfun(@(x)assert(exist(x,'file'),'Input DWI %s does not exist', x),fdwi,'UniformOutput',false);
    cellfun(@(x)assert(exist(x,'file'),'Input Bvec %s does not exist', x),fbvec,'UniformOutput',false);
    cellfun(@(x)assert(exist(x,'file'),'Input Bval %s does not exist', x),fbval,'UniformOutput',false);
    cellfun(@(x)assert(exist(x,'file'),'Input T2 relax %s does not exist', x),ft2r,'UniformOutput',false);
    assert(exist(fmask,'file'),'Input mask %s does not exist', fmask);
    assert(exist(options.spectrum,'file'),'Input spectrum %s does not exist', options.spectrum);

    %% load multi-echo dMRI dataset
    ME_dwi_info     = cellfun(@(x)niftiinfo(x),fdwi,'UniformOutput',false);
    ME_dwi          = cellfun(@(x)niftiread(x),ME_dwi_info,'UniformOutput',false);
    ME_bval         = cellfun(@(x)round(importdata(x)'/100)*100,fbval,'UniformOutput',false); 
    ME_bvec         = cellfun(@(x)importdata(x)',fbvec,'UniformOutput',false);
    ME_mask_info    = niftiinfo(fmask);
    ME_mask         = round(niftiread(ME_mask_info));
    ME_t2r_info     = cellfun(@(x)niftiinfo(x),ft2r,'UniformOutput',false);
    ME_t2r          = cellfun(@(x)niftiread(x),ME_t2r_info,'UniformOutput',false);
    clear fdwi fbvec fbval fmask ft2r;
    clear ME_t2r_info ME_mask_info

    if options.useBshell
        ind         = cellfun(@(x)ismember(x,options.useBshell),fbval,'UniformOutput',false); 
        ME_dwi      = cellfun(@(x,y)x(:,:,:,y),ME_dwi,ind,'UniformOutput',false);
        ME_bvec     = cellfun(@(x,y)x(y,:),ME_bvec,ind,'UniformOutput',false); 
        ME_bval     = cellfun(@(x,y)x(y),ME_bval,ind,'UniformOutput',false); 
        clear ind
    end

    %% Normalization S/S0
    if options.normalizeToS0
        ind             = cellfun(@(x)ismember(x,0),ME_bval,'UniformOutput',false); 
        ME_dwi_norm     = cellfun(@(x,y)x(:,:,:,~y)./(mean(x(:,:,:,~y),4)+eps),ME_dwi,ind,'UniformOutput',false);
        ME_bval_norm    = cellfun(@(x,y)x(~y,:),ME_bval,ind,'UniformOutput',false);
        ME_bvec_norm    = cellfun(@(x,y)x(~y,:),ME_bvec,ind,'UniformOutput',false);
        clear ind
  
        ME_dwi = ME_dwi_norm;    clear ME_dwi_norm;
        ME_bval = ME_bval_norm;  clear ME_bval_norm;
        ME_bvec = ME_bvec_norm;  clear ME_bvec_norm;
    end
    
    %% kernel
    default_spectrum = load(options.spectrum);

    adc_restricted   = default_spectrum.adc_restricted;
    adc_hindered     = default_spectrum.adc_hindered;
    adc_isotropic    = default_spectrum.adc_isotropic;

    num_restricted  = size(adc_restricted,1);  
    num_hindered    = size(adc_hindered,1);    
    num_isotropic   = size(adc_isotropic,1);   
    
    kernel_restricted = cell(length(TE),num_restricted);
    kernel_hindered   = cell(length(TE),num_hindered);
    kernel_isotropic  = cell(length(TE),num_isotropic);

    lmax = 6; 
    nmax = lmax2nsh(lmax);
    scheme = gen_scheme('scheme/sphere_362_vertices.txt',lmax);

    for i = 1:length(TE)
        bval    = ME_bval{i};
        bvec    = ME_bvec{i};
        bshell  = unique(bval);
        nvol    = length(bval);
        
        for j = 1:num_restricted
            kernel_restricted{i,j} = zeros(nvol,nmax);
            for k = 1:length(bshell)
                order = min(floor(nsh2lmax(sum(bval==bshell(k)))),lmax);
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
                order = min(floor(nsh2lmax(sum(bval==bshell(k)))),lmax);
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
    ME_t2r_array = cellfun(@(x)reshape(x,size(x,1)*size(x,2)*size(x,3),1),ME_t2r,'UniformOutput',false);
    ME_t2r_array = cellfun(@(x)x(ME_mask_ind,:)',ME_t2r_array,'UniformOutput',false);
    clear ME_dwi ME_t2r

    %% subject to
    nv = size(scheme.vert,1);
    ampbasis = repmat(scheme.sh,1,num_restricted + num_hindered);
    ampbasis = mat2cell(ampbasis,nv,repmat(nmax,1,num_restricted + num_hindered));
    ampbasis = blkdiag(ampbasis{:});
    A1 = blkdiag(ampbasis,diag(ones(1,num_isotropic)));
    A2 = zeros(size(A1,1),1);
    A3 = ones(size(A1,1),1);

    alpha_coef = zeros(num_restricted*nmax+num_hindered*nmax+num_isotropic,size(ME_dwi_array,2));

    ME_t2r_restricted  = exp(-TE'./ME_t2r_array{1});
    ME_t2r_hindered    = exp(-TE'./ME_t2r_array{2});
    ME_t2r_isotropic   = exp(-TE'./ME_t2r_array{3});

    clear ampbasis ME_t2r_array;

    %% optimization
    parfor i = 1:size(ME_dwi_array,2)
        kernel = cell2mat([ cellfun(@(x,y) x.*y, kernel_restricted,num2cell(ME_t2r_restricted(:,i).*ones(1,num_restricted)), 'UniformOutput',false) ...
                            cellfun(@(x,y) x.*y, kernel_hindered,num2cell(ME_t2r_hindered(:,i).* ones(1,num_hindered)), 'UniformOutput',false) ...
                            cellfun(@(x,y) x.*y, kernel_isotropic,num2cell(ME_t2r_isotropic(:,i).* ones(1,num_isotropic)), 'UniformOutput',false)]);

        dwi = ME_dwi_array(:,i);
        ind = ~isnan(dwi) & ~isinf(dwi);

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

    %% save results
    if ~exist(outpath,'dir')
        mkdir(outpath)
    end

    % save FOD restricted
    temp = single(zeros(nmax,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    info_fod = ME_dwi_info{1};
    info_fod.Datatype = 'single';
    info_fod.ImageSize(4) = nmax;
    for i = 1:num_restricted+num_hindered
        temp(:,ME_mask_ind) = alpha_coef((i-1)*nmax+1:i*nmax,:);
        fod = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3),nmax);
        if i <= num_restricted
            niftiwrite(single(fod),fullfile(outpath,strcat('FOD_restricted_',num2str(i),'.nii')),info_fod,'Compressed', true);
        else
            niftiwrite(single(fod),fullfile(outpath,strcat('FOD_hindered_',num2str(i-num_restricted),'.nii')),info_fod,'Compressed', true);
        end
    end

    temp = single(zeros(num_isotropic,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    info_fod.Datatype = 'single';
    info_fod.ImageSize(4) = num_isotropic;
    temp(:,ME_mask_ind) = alpha_coef(end-num_isotropic+1:end,:)/sqrt(4*pi);
    fod = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3),num_isotropic);

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


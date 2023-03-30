function ME_eSMSI(fdwi,fbval,fmask,TE,outpath,options)
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
        fdwi                    (1,:)    {mustBeNonzeroLengthText}
        fbval                   (1,:)    {mustBeNonzeroLengthText}
        fmask                   string   {mustBeFile}
        TE                      (1,:)    {mustBeNumeric}
        outpath                 string   {mustBeNonzeroLengthText}

        options.normalizeToS0   (1,1)    {mustBeNumericOrLogical} = true
        options.useBshell       (1,:)    {mustBeNumericOrLogical} = false
        options.spectrum        string   {mustBeFile} = 'scheme/exchange_spectrum.mat'
    end

    addpath('third/osqp');
    addpath('scheme');

    cellfun(@(x)assert(exist(x,'file'),'Input DWI %s does not exist', x),fdwi,'UniformOutput',false);
    cellfun(@(x)assert(exist(x,'file'),'Input Bval %s does not exist', x),fbval,'UniformOutput',false);
    assert(exist(fmask,'file'),'Input mask %s does not exist', fmask);
    assert(exist(options.spectrum,'file'),'Input spectrum %s does not exist', options.spectrum);

    %% load multi-echo dMRI dataset
    ME_dwi_info     = cellfun(@(x)niftiinfo(x),fdwi,'UniformOutput',false);
    ME_dwi          = cellfun(@(x)niftiread(x),ME_dwi_info,'UniformOutput',false);
    ME_bval         = cellfun(@(x)round(importdata(x)'/100)*100,fbval,'UniformOutput',false); 
    ME_mask_info    = niftiinfo(fmask);
    ME_mask         = round(niftiread(ME_mask_info));
    clear ind fdwi fbval fmask;

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
    default_spectrum = load(options.spectrum);

    adc_restricted   = default_spectrum.adc_restricted;
    adc_hindered     = default_spectrum.adc_hindered;
    adc_isotropic    = default_spectrum.adc_isotropic;
    adc_semi_restricted    = default_spectrum.adc_semi_restricted;
    adc_free    = default_spectrum.adc_free;


    num_restricted   = default_spectrum.num_restricted;
    num_hindered     = default_spectrum.num_hindered;
    num_isotropic    = default_spectrum.num_isotropic;
    num_semi_restricted    = default_spectrum.num_semi_restricted;
    num_free    = default_spectrum.num_free;  
    
    kernel_restricted = cell(length(TE),num_restricted);
    kernel_hindered   = cell(length(TE),num_hindered);
    kernel_isotropic  = cell(length(TE),num_isotropic);
    kernel_semi_restricted   = cell(length(TE),num_semi_restricted);
    kernel_free       = cell(length(TE),num_free);

    for i = 1:length(TE)
        for j = 1:num_restricted
            kernel_restricted{i,j} = single(zeros(length(ME_bshell{i}),1));
            for k = 1:length(ME_bshell{i})
                R = sqrt(ME_bshell{i}(k)*(adc_restricted(j,1)-adc_restricted(j,2)));
                kernel_restricted{i,j}(k,1) = (sqrt(pi)*exp(-ME_bshell{i}(k)*adc_restricted(j,2))*erf(R))/(2*R);
            end
        end

        for j = 1:num_semi_restricted
            kernel_semi_restricted{i,j} = single(zeros(length(ME_bshell{i}),1));
            for k = 1:length(ME_bshell{i})
                R = sqrt(ME_bshell{i}(k)*(adc_semi_restricted(j,1)-adc_semi_restricted(j,2)));
                kernel_semi_restricted{i,j}(k,1) = (sqrt(pi)*exp(-ME_bshell{i}(k)*adc_semi_restricted(j,2))*erf(R))/(2*R);
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

        for j = 1:num_free
            kernel_free{i,j} = single(zeros(length(ME_bshell{i}),1));
            for k = 1:length(ME_bshell{i})
                kernel_free{i,j}(k,1) = exp(-ME_bshell{i}(k)*adc_free(j)); 
            end
        end
    end

    blk_kernel_restricted = [];
    for i = 1:num_restricted
        blk_kernel_restricted = [blk_kernel_restricted blkdiag(kernel_restricted{:,i})];
    end

    blk_kernel_semi_restricted = [];
    for i = 1:num_semi_restricted
        blk_kernel_semi_restricted = [blk_kernel_semi_restricted blkdiag(kernel_semi_restricted{:,i})];
    end
    

    blk_kernel_hindered = [];
    for i = 1:num_hindered
        blk_kernel_hindered = [blk_kernel_hindered blkdiag(kernel_hindered{:,i})];
    end

    blk_kernel_isotropic = [];
    for i = 1:num_isotropic
        blk_kernel_isotropic = [blk_kernel_isotropic blkdiag(kernel_isotropic{:,i})];
    end

    blk_kernel_free = [];
    for i = 1:num_free
        blk_kernel_free = [blk_kernel_free blkdiag(kernel_free{:,i})];
    end
    
    blk_kernel = [blk_kernel_restricted blk_kernel_semi_restricted blk_kernel_hindered blk_kernel_isotropic blk_kernel_free];


    %% Vectorization & Masked & arrayed
    ME_mask_ind = find(ME_mask>0.5);
    ME_dwi_mean_array = cellfun(@(x)reshape(x,size(x,1)*size(x,2)*size(x,3),size(x,4)),ME_dwi_mean,'UniformOutput',false);
    ME_dwi_mean_array = cellfun(@(x)x(ME_mask_ind,:)',ME_dwi_mean_array,'UniformOutput',false);
    ME_dwi_mean_array = cell2mat(ME_dwi_mean_array');
    ME_dwi_mean_array_ind = all(ME_dwi_mean_array);
    ME_dwi_mean_array = ME_dwi_mean_array(:,ME_dwi_mean_array_ind);
    clear ME_dwi_mean

    %% sum(total vf) = 1 if normalizeToS0 is true
    if options.normalizeToS0
        Aeq     = [zeros(1,5) ones(1,num_restricted+num_semi_restricted+num_hindered+num_isotropic+num_free)];
        beq     = 1;
    else
        Aeq     = [];
        beq     = [];
    end

%     lb = [zeros(5,1);zeros(num_restricted+num_semi_restricted+num_hindered+num_isotropic+num_free,1)];
    %ub = [inf(5,1);ones(num_restricted+num_semi_restricted+num_hindered+num_isotropic+num_free,1)];

    options = optimoptions(@fmincon,'Display','off');
    alpha_coef = zeros(num_restricted+num_semi_restricted+num_hindered+num_isotropic+num_free+5,size(ME_dwi_mean_array,2));
    alpha_init = [100;100;100;100;100;rand(num_restricted+num_semi_restricted+num_hindered+num_isotropic+num_free,1)];

    H   = double(blk_kernel'*blk_kernel);
    K   = double(-blk_kernel'*ME_dwi_mean_array);
    A1  = diag(ones(size(blk_kernel,2),1));
    A2  = zeros(size(blk_kernel,2),1);
    A3  = ones(size(blk_kernel,2),1);
    clear blk_kernel;

    %% optimization
    parfor i = 1:size(ME_dwi_mean_array,2)
        
        f = K(:,i);
        prob = osqp;
        prob.setup(H,f,A1,A2,A3,'alpha',0.1,'verbose',0);
        res = prob.solve();
        fun = @(x)echoFunc(x,min(max(0,res.x),1),TE,num_restricted,num_semi_restricted,num_hindered,num_isotropic,num_free);
        alpha_coef(:,i) = max(0,fmincon(fun,alpha_init,[],[],Aeq,beq,[],[],[],options));
    end
    
    %% save results
    if ~exist(outpath,'dir')
        mkdir(outpath)
    end

    % save T2 restricted
    temp2 = single(zeros(1,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(1,:);
    temp = single(zeros(1,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3));

    info_t2r = ME_mask_info;
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    niftiwrite(single(T2r),fullfile(outpath,'T2_restricted.nii'),info_t2r,'Compressed', true);

    % save T2 semi_restricted
    temp2 = single(zeros(1,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(2,:);
    temp = single(zeros(1,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3));

    info_t2r = ME_mask_info;
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    niftiwrite(single(T2r),fullfile(outpath,'T2_semi_restricted.nii'),info_t2r,'Compressed', true);

    % save T2 hindered
    temp2 = single(zeros(1,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(3,:);
    temp = single(zeros(1,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3));

    info_t2r = ME_mask_info;
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    niftiwrite(single(T2r),fullfile(outpath,'T2_hindered.nii'),info_t2r,'Compressed', true);

    % save T2 isotropic
    temp2 = single(zeros(1,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(4,:);
    temp = single(zeros(1,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3));

    info_t2r = ME_mask_info;
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    niftiwrite(single(T2r),fullfile(outpath,'T2_isotropic.nii'),info_t2r,'Compressed', true);

    % save T2 free water
    temp2 = single(zeros(1,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(5,:);
    temp = single(zeros(1,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3));

    info_t2r = ME_mask_info;
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    niftiwrite(single(T2r),fullfile(outpath,'T2_free.nii'),info_t2r,'Compressed', true);

    % save VF restricted
    temp2 = single(zeros(num_restricted,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(6:5+num_restricted,:);
    temp = single(zeros(num_restricted,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    vf = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3),num_restricted);

    info_vf= ME_dwi_info{1};
    info_vf.Datatype = 'single';
    info_vf.ImageSize = size(vf);
    niftiwrite(single(vf),fullfile(outpath,'VF_restricted.nii'),info_vf,'Compressed', true);

    % save VF semi_restricted
    temp2 = single(zeros(num_semi_restricted,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(6+num_restricted:5+num_restricted+num_semi_restricted,:);
    temp = single(zeros(num_semi_restricted,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    vf = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3),num_semi_restricted);

    info_vf= ME_dwi_info{1};
    info_vf.Datatype = 'single';
    info_vf.ImageSize = size(vf);
    niftiwrite(single(vf),fullfile(outpath,'VF_semi_restricted.nii'),info_vf,'Compressed', true);

    % save VF hindered
    temp2 = single(zeros(num_hindered,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(6+num_restricted+num_semi_restricted:5+num_restricted+num_semi_restricted+num_hindered,:);
    temp = single(zeros(num_hindered,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    vf = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3),num_hindered);

    info_vf= ME_dwi_info{1};
    info_vf.Datatype = 'single';
    info_vf.ImageSize = size(vf);
    niftiwrite(single(vf),fullfile(outpath,'VF_hindered.nii'),info_vf,'Compressed', true);

    % save VF isotropic
    temp2 = single(zeros(num_isotropic,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(6+num_restricted+num_semi_restricted+num_hindered:5+num_restricted+num_semi_restricted+num_hindered+num_isotropic,:);

    temp = single(zeros(num_isotropic,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    vf = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3),num_isotropic);

    info_vf= ME_dwi_info{1};
    info_vf.Datatype = 'single';
    info_vf.ImageSize = size(vf);
    niftiwrite(single(vf),fullfile(outpath,'VF_isotropic.nii'),info_vf,'Compressed', true);

    % save VF free water
    temp2 = single(zeros(num_free,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(6+num_restricted+num_semi_restricted+num_hindered+num_isotropic:end,:);
    temp = single(zeros(num_free,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    vf = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3),num_free);

    info_vf= ME_dwi_info{1};
    info_vf.Datatype = 'single';
    info_vf.ImageSize = size(vf);
    niftiwrite(single(vf),fullfile(outpath,'VF_free.nii'),info_vf,'Compressed', true);


end

function FF = echoFunc(x,alpha,TE,num_restricted,num_semi_restricted,num_hindered,num_isotropic,num_free)
    % x: [ T_2^r, T_2^h, T_2^f, w_1^r, w_2^r...w_1^h, w_2^h...,w_1^f, w_2^f...]
%     x = max(0,x);
    num = 1;
    for i = 1:num_restricted
        for j = 1:length(TE)
            F(num) = exp(-TE(j)/x(1)) * x(5+i) - alpha(num);
            num = num + 1;
        end
    end

    for i = 1:num_semi_restricted
        for j = 1:length(TE)
            F(num) = exp(-TE(j)/x(2)) * x(5+num_restricted+i) - alpha(num);
            num = num + 1;
        end
    end

    for i = 1:num_hindered
        for j = 1:length(TE)
            F(num) = exp(-TE(j)/x(3)) * x(5+num_restricted+num_semi_restricted+i) - alpha(num);
            num = num + 1;
        end
    end

    for i = 1:num_isotropic
        for j = 1:length(TE)
            F(num) = exp(-TE(j)/x(4)) * x(5+num_restricted+num_semi_restricted+num_hindered+i) - alpha(num);
            num = num + 1;
        end
    end

    for i = 1:num_free
        for j = 1:length(TE)
            F(num) = exp(-TE(j)/x(5)) * x(5+num_restricted+num_semi_restricted+num_hindered+num_isotropic+i) - alpha(num);
            num = num + 1;
        end
    end

    FF = norm(F);
end
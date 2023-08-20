function ME_SMT(fdwi,fbval,fmask,TE,outpath,options)
    %{
    ░█▀▀█ ░█▀▀▀█ ░█▀▄▀█ ░█▀▀▀ ░█▀▀▄ ▀█▀
    ░█─── ░█──░█ ░█░█░█ ░█▀▀▀ ░█─░█ ░█─
    ░█▄▄█ ░█▄▄▄█ ░█──░█ ░█▄▄▄ ░█▄▄▀ ▄█▄

    Multi-echo spherical mean technology (ME-SMSI)


        Created by Ye Wu, PhD (dr.yewu@outlook.com)

        - Nanjing University of Science and Technology, China
        - University of North Carolina at Chapel Hill, USA
        
    %}
    
    arguments
        fdwi                    (1,:)    {mustBeNonzeroLengthText}
        fbval                   (1,:)    {mustBeNonzeroLengthText}
        fmask                   string   {mustBeFile}
        TE                      (1,:)    {mustBeNumeric}
        outpath                 string   {mustBeNonzeroLengthText}

        options.normalizeToS0   (1,1)    {mustBeNumericOrLogical} = true
        options.useBshell       (1,:)    {mustBeNumericOrLogical} = false
        options.lambda          (1,1)    {mustBeNumeric} = 0.015
    end

    addpath('third/osqp');
    addpath('scheme');

    cellfun(@(x)assert(exist(x,'file'),'Input DWI %s does not exist', x),fdwi,'UniformOutput',false);
    cellfun(@(x)assert(exist(x,'file'),'Input Bval %s does not exist', x),fbval,'UniformOutput',false);
    assert(exist(fmask,'file'),'Input mask %s does not exist', fmask);

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
    
    TE_array = [];
    B_array = [];

    for i = 1:length(TE)
        for k = 1:length(ME_bshell{i})
            TE_array = [TE_array;TE(i)];
            B_array = [B_array;ME_bshell{i}(k)];
        end
    end

    %% Vectorization & Masked & arrayed
    ME_mask_ind = find(ME_mask>0.5);
    ME_dwi_mean_array = cellfun(@(x)reshape(x,size(x,1)*size(x,2)*size(x,3),size(x,4)),ME_dwi_mean,'UniformOutput',false);
    ME_dwi_mean_array = cellfun(@(x)x(ME_mask_ind,:)',ME_dwi_mean_array,'UniformOutput',false);
    ME_dwi_mean_array = cell2mat(ME_dwi_mean_array');
    ME_dwi_mean_array_ind = all(ME_dwi_mean_array);
    ME_dwi_mean_array = ME_dwi_mean_array(:,ME_dwi_mean_array_ind);
    clear ME_dwi_mean

    %% 
    A = [0 0 0 0 0 0 -(2/pi)^2 1 0 0 0; ...
         0 0 0 0 0 0 0 0 -1/1.1 1 0; ...
         0 0 0 0 0 0 0 0 1 -(pi/2)^2 0];
    b = [0;0;0];
    Aeq = [0 0 0 1 1 1 0 0 0 0 0];
    beq = 1;
    lb = [10;10;10;0.0001*ones(8,1)];
    ub = [1500;1500;1500;1;1;1;0.001;0.001;0.001;0.001;0.001];
    option = optimoptions(@fmincon,'Display','off');
    alpha_coef = zeros(11,size(ME_dwi_mean_array,2));
    alpha_nmse = zeros(1,size(ME_dwi_mean_array,2));
    alpha_init = [100;100;100;0.3;0.3;0.4;0.0017;0.0002;0.0015;0.001;0.003];

    %% optimization
    parfor i = 1:size(ME_dwi_mean_array,2)
        S = double(ME_dwi_mean_array(:,i));
        fun = @(x)echoMCFunc(x,S,TE_array,B_array);
        beta = fmincon(fun,alpha_init,A,b,[],[],lb,ub,[],option);
        alpha_coef(:,i) = beta;
        alpha_nmse(1,i) = echoMCFunc(beta,S,TE_array,B_array) / norm(S);
    end
    
    %% save results
    if ~exist(outpath,'dir')
        mkdir(outpath)
    end
 
    % save T2 r
    temp2 = single(zeros(1,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(1,:);
    temp = single(zeros(1,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3));

    info_t2r = ME_dwi_info{1};
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    info_t2r.PixelDimensions = info_t2r.PixelDimensions(1:length(size(T2r)));
    niftiwrite(single(T2r),fullfile(outpath,'T2_restricted.nii'),info_t2r,'Compressed', true);

    % save T2 h
    temp2 = single(zeros(1,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(2,:);
    temp = single(zeros(1,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3));

    info_t2r = ME_dwi_info{1};
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    info_t2r.PixelDimensions = info_t2r.PixelDimensions(1:length(size(T2r)));
    niftiwrite(single(T2r),fullfile(outpath,'T2_hindered.nii'),info_t2r,'Compressed', true);

    % save T2 f
    temp2 = single(zeros(1,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(3,:);
    temp = single(zeros(1,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3));

    info_t2r = ME_dwi_info{1};
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    info_t2r.PixelDimensions = info_t2r.PixelDimensions(1:length(size(T2r)));
    niftiwrite(single(T2r),fullfile(outpath,'T2_free.nii'),info_t2r,'Compressed', true);

    % save vf r
    temp2 = single(zeros(1,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(4,:);
    temp = single(zeros(1,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3));

    info_t2r = ME_dwi_info{1};
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    info_t2r.PixelDimensions = info_t2r.PixelDimensions(1:length(size(T2r)));
    niftiwrite(single(T2r),fullfile(outpath,'VF_restricted.nii'),info_t2r,'Compressed', true);

    % save vf h
    temp2 = single(zeros(1,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(5,:);
    temp = single(zeros(1,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3));

    info_t2r = ME_dwi_info{1};
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    info_t2r.PixelDimensions = info_t2r.PixelDimensions(1:length(size(T2r)));
    niftiwrite(single(T2r),fullfile(outpath,'VF_hindered.nii'),info_t2r,'Compressed', true);

    % save vf f
    temp2 = single(zeros(1,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(6,:);
    temp = single(zeros(1,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3));

    info_t2r = ME_dwi_info{1};
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    info_t2r.PixelDimensions = info_t2r.PixelDimensions(1:length(size(T2r)));
    niftiwrite(single(T2r),fullfile(outpath,'VF_free.nii'),info_t2r,'Compressed', true);

    temp2 = single(zeros(1,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(7,:);
    temp = single(zeros(1,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3));

    info_t2r = ME_dwi_info{1};
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    info_t2r.PixelDimensions = info_t2r.PixelDimensions(1:length(size(T2r)));
    niftiwrite(single(T2r),fullfile(outpath,'lambda_1_restricted.nii'),info_t2r,'Compressed', true);

    temp2 = single(zeros(1,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(8,:);
    temp = single(zeros(1,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3));

    info_t2r = ME_dwi_info{1};
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    info_t2r.PixelDimensions = info_t2r.PixelDimensions(1:length(size(T2r)));
    niftiwrite(single(T2r),fullfile(outpath,'lambda_2_restricted.nii'),info_t2r,'Compressed', true);

    temp2 = single(zeros(1,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(9,:);
    temp = single(zeros(1,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3));

    info_t2r = ME_dwi_info{1};
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    info_t2r.PixelDimensions = info_t2r.PixelDimensions(1:length(size(T2r)));
    niftiwrite(single(T2r),fullfile(outpath,'lambda_1_hindered.nii'),info_t2r,'Compressed', true);

    temp2 = single(zeros(1,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(10,:);
    temp = single(zeros(1,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3));

    info_t2r = ME_dwi_info{1};
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    info_t2r.PixelDimensions = info_t2r.PixelDimensions(1:length(size(T2r)));
    niftiwrite(single(T2r),fullfile(outpath,'lambda_2_hindered.nii'),info_t2r,'Compressed', true);

    temp2 = single(zeros(1,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_coef(11,:);
    temp = single(zeros(1,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3));

    info_t2r = ME_dwi_info{1};
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    info_t2r.PixelDimensions = info_t2r.PixelDimensions(1:length(size(T2r)));
    niftiwrite(single(T2r),fullfile(outpath,'lambda_free.nii'),info_t2r,'Compressed', true);

    temp2 = single(zeros(1,length(ME_dwi_mean_array_ind)));
    temp2(:,ME_dwi_mean_array_ind) = alpha_nmse;
    temp = single(zeros(1,size(ME_mask,1)*size(ME_mask,2)*size(ME_mask,3)));
    temp(:,ME_mask_ind) = temp2;
    T2r = reshape(temp',size(ME_mask,1),size(ME_mask,2),size(ME_mask,3));

    info_t2r = ME_dwi_info{1};
    info_t2r.Datatype = 'single';
    info_t2r.ImageSize = size(T2r);
    info_t2r.Transform = ME_dwi_info{1}.Transform;
    info_t2r.PixelDimensions = info_t2r.PixelDimensions(1:length(size(T2r)));
    niftiwrite(single(T2r),fullfile(outpath,'NMSE.nii'),info_t2r,'Compressed', true);

end


function FF = echoMCFunc(x,S,TE,b)
    % x = [T2r, T2f, T2f, fr, fh, ff, r_para, r_perp, h_para, h_perp, f];

    Rr = sqrt(b*max(x(7)-x(8),0));
    Rh = sqrt(b*max(x(9)-x(10),0));

    F = x(4)*exp(-TE/x(1)).*(sqrt(pi)*exp(-b*x(8)).*erf(Rr))./(2*Rr) + ...
        x(5)*exp(-TE/x(2)).*(sqrt(pi)*exp(-b*x(10)).*erf(Rh))./(2*Rh) + ...
        x(6)*exp(-TE/x(3)).*exp(-b*x(11)) - S;

    FF = norm(F)/norm(S);
end

function demo_rsi(dwi_filename,bval_filename,bvec_filename,mask_filename)
% Restriction Spectrum Imaging (RSI)
%
% =========================================================================
% Author: Ye Wu Ph.D
%
% E-mail: dr.yewu@outlook.com
%
% Created: 2023-01-12, using Matlab 2022b
% =========================================================================

params = [];
params.lambda = 0.1;                % regularization constant
params.SH_order = 4;                % spherical harmonic order -- must be even
params.norm_flag = 0;               % normalize data to b=0 image
params.nonlin_flag = 0;             % use nonlinear optimization with initial parameters from linear fit
params.b0_thresh = 10;              % threshold used for considering a b-value to be 0
params.scalefacts_flag = 0;         % calculate scaling factors from b=0 images and apply them to all subsequent frames (for multiple acquisitions)

params.ADC_long = [];               % longitudinal ADC
params.ADC_trans = [];              % transverse ADC
num = 0;
for ADC_long = 0.5e-3 : 0.2e-3 : 1.5e-3
    for ADC_trans = 0e-3 : 0.2e-3 : 0.9e-3
        if ADC_trans < ADC_long
            params.ADC_long = [params.ADC_long ADC_long];
            params.ADC_trans = [params.ADC_trans ADC_trans];
            num = num + 1;
        end
    end
end
params.ADC_iso = 0:0.2e-3:3e-3;     % minimum isotropic ADC
params.num_ADC_aniso = num;         % number of transverse ADC size scales
params.num_ADC_iso = length(params.ADC_iso); % number of isotropic ADC size scales

%% load data
dwi_info = niftiinfo(dwi_filename);
dwi = niftiread(dwi_info);

mask_info = niftiinfo(mask_filename);
mask = niftiread(mask_info);

bvals = importdata(bval_filename);
bvecs = importdata(bvec_filename);

if size(bvecs,2) ~=3; bvecs = bvecs'; end
if size(bvals,2) ~=1; bvals = bvals'; end

bvals = round(bvals/params.b0_thresh)*params.b0_thresh;

%% data header
params.numvol = length(dwi);
params.voxsiz = mask_info.PixelDimensions;
params.voxdim = size(mask);
params.i_b0 = find(bvals <= params.b0_thres);
params.i_dwi = find(bvals > params.b0_thresh);

%% create forward matrix for fitting tensor; construct tensor forward matrix
[params.B,params.Binv] = create_tensor_matrix(bvals,bvecs);

%% construct RSI multi-FOD forward matrix
[params.A,params.Ainv,params.icoverts,params.beta2ico] = create_rsi_matrix(MFfit);

% linear RSI fit
MFfit = fit_rsi(vol,MFfit,parms);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B,Binv] = create_tensor_matrix(bvals,bvecs)
  ndirs = length(bvals);
  B = zeros(ndirs,7);
  B(:,7) = 1; % S0
  for i = 1:ndirs
    outerprod = -bvals(i).*bvecs(i,:)'*bvecs(i,:);
    B(i,1:3) = diag(outerprod)';
    B(i,4) = 2*outerprod(1,2);
    B(i,5) = 2*outerprod(1,3);
    B(i,6) = 2*outerprod(2,3);
  end
  Binv = pinv(B);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,Ainv,V,F] = create_rsi_matrix(bvals,bvecs,params)
  V = importdata('sphere_362_vertices.txt');
  Q = bvecs.*repmat(sqrt(bvals),1,3);
  F = rsi_SH_matrix(V,params.SH_order);
  F0 = rsi_SH_matrix(V,0);
  A = [];

  % series of anisotropic FODs for varying size scales
  for i = 1:params.num_ADC_aniso
    for j = 1:params.num_ADC_aniso
        R = rsi_FOD_matrix(Q,V,params.ADC_long(i),params.ADC_trans(j));
        A = [A R*F];
    end
  end

  % series of isotropic FODs for varying size scales
  for t=1:params.num_ADC_iso
    R = rsi_FOD_matrix(Q,V,params.ADC_iso(t),params.ADC_iso(t));
    A = [A R*F0];
  end

  % compute regularized inverse
  AtA = A'*A;
  Ainv = (AtA+params.lambda*mean(diag(AtA))*eye(size(AtA)))\A';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MFfit = fit_rsi(vol,parms)
  % initialize multi-FOD parameter volumes
  MFfit.volMF = zeros(parms.nx,parms.ny,parms.nz,MFfit.nb);
  % loop over slices
  for z = 1:parms.nz
    y = reshape(vol(:,:,z,:),[parms.nx*parms.ny,parms.nf])';
    if parms.norm_flag % normalize data to b=0
      y0 = repmat(reshape(MFfit.volb0(:,:,z),[parms.nx*parms.ny,1])',[parms.nf 1]);
      y = y./y0;
    end;
    betas = MFfit.Ainv*y;
    betas = reshape(betas',[parms.nx parms.ny MFfit.nb]);
    betas(isnan(betas)) = 0;
    MFfit.volMF(:,:,z,:) = betas;
  end
  % dilate brain mask
  MFfit.volmask_dilated = dilate_brain_mask(MFfit.volmask,parms);
return;



















bList = unique(bVals);
bNums = length(bList);
gNums = length(bVals);
b0Index = bVals == 0;

bIndex = cell(1,bNums);
gIndex = cell(1,bNums);
for idx = 1:bNums
    bIndex{idx} = find(bVals==bList(idx));
    gIndex{idx} = g_data(bIndex{idx},:);
end























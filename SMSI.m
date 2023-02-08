function SMSI()
Path = 'D:\Work\SMSI';
addpath('D:\Work\MNAlab\Functions');
% =========
% load data
% =========
dwi = single(niftiread(fullfile(Path,'synthetic_degeneracy.nii.gz')));
bvecs = importdata(fullfile(Path,'bvecs'));
bvals = importdata(fullfile(Path,'bvals'));
% mask = int8(niftiread(fullfile(Path,'mask.nii.gz')));
mask = int8(ones(4,1,1));
fscheme = 'sphere_362_vertices.txt';

% =========
% fix bvals and bvecs
% =========
fixparam = 100;
bvals = round(bvals/fixparam)*fixparam;

if size(bvals,2) > 1
    bvals = bvals';
end

if size(bvecs,2) > 3
    bvecs = bvecs';
end

% =========
% gradient table (remove b0)
% =========
indminb = bvals == min(bvals);
bvals(indminb) = [];
bvecs(indminb,:) = [];
list_shell = unique(bvals);
n_shell = length(list_shell);


% =========
% construct observed signal (normalized to s0)
% =========
[nx,ny,nz] = deal(size(dwi,1),size(dwi,2),size(dwi,3));
full_dwi = dwi(:,:,:,~indminb)./mean(dwi(:,:,:,indminb),4);

mean_dwi = single(zeros(nx,ny,nz,n_shell));
index_shell = cell(1,n_shell);
for i = 1:n_shell
    index_shell{1,i} = find(bvals==list_shell(i));
    mean_dwi(:,:,:,i) = mean(full_dwi(:,:,:,index_shell{1,i}),4);
end


% =========
% masking observed signal
% =========
ind = find(mask>0.5);
nmask = length(ind);
nvox = nx*ny*nz;

mean_dwi_mask = reshape(mean_dwi,nvox,n_shell)'; mean_dwi_mask = mean_dwi_mask(:,ind);
full_dwi_mask = reshape(full_dwi,nvox,[])'; full_dwi_mask = full_dwi_mask(:,ind);


% =========
% setting diffusivity
% =========
isotropy = [];
for i = 0.3e-3:0.1e-3:3e-3
    isotropy = [isotropy i];
end

longitudinal = [];
transverse = [];
for i = 1.7e-3
    for j = 0.1:0.1:0.9
            longitudinal = [longitudinal i];
            transverse = [transverse i*j];
    end
end

% for i = 0.2e-3:0.1e-3:2.5e-3
%     for j = i:0.1e-3:2.5e-3
%         if j > i*2
%             longitudinal = [longitudinal j];
%             transverse = [transverse i];
%         end
%     end
% end


% =========
% calculate kernel (in spherical harmonics space) for full signal fitting
% =========
n_ani = length(longitudinal);
n_iso = length(isotropy);
assert(isequal(length(longitudinal),length(transverse)),'Dimension between longitudinal and longitudinal must be consistent')

lmax = 6;
shorder = zeros(1,n_shell);
for i = 1:n_shell
   shorder(1,i) = floor(nsh2lmax(length(index_shell{i})));
end

lmax = min(max(shorder),lmax);
nmax = lmax2nsh(lmax);
shorder = min(shorder,lmax);

% estimate anisotropy component
scheme = gen_scheme(fscheme,lmax);
fRF = cell(1,n_ani+n_iso);

for i = 1:n_ani
   fRF{1,i} = zeros(size(bvecs,1),nmax);

   for j = 1:n_shell
       DW_scheme = gen_scheme(bvecs(index_shell{j},1:3),shorder(1,j));

       R_amp = response(longitudinal(i),transverse(i),list_shell(j),scheme);
       R_SH = amp2SH(R_amp, scheme);
       R_RH = SH2RH(R_SH);

       m = [];
       for l = 0:2:shorder(1,j)
         m = [ m R_RH(l/2+1)*ones(1,2*l+1) ];
       end

       fconv = DW_scheme.sh .* m(ones(size(DW_scheme.sh,1),1),:);
       fconv(:,end+1:nmax) = 0;

       fRF{1,i}(index_shell{j},:) = fconv;
       
       clear DW_scheme R_amp R_SH R_RH fconv m;
   end
end
           
           
% estimate isotropy component
for i = 1:n_iso
   fRF{1,n_ani+i} = zeros(size(bvecs,1),1);

   for j = 1:n_shell
       DW_scheme = gen_scheme(bvecs(index_shell{j},1:3),shorder(1,j));

       R_amp = response(isotropy(i),isotropy(i),list_shell(j),scheme);
       R_SH = amp2SH(R_amp, scheme);
       R_RH = SH2RH(R_SH(1));

       fconv = DW_scheme.sh(1) .* R_RH(ones(size(DW_scheme.sh,1),1),:);

       fRF{1,n_ani+i}(index_shell{j},:) = fconv;

       clear DW_scheme R_amp R_SH R_RH fconv;

   end
end


% =========
% calculate kernel (in spherical harmonics space) for mean signal fitting
% =========
% estimate anisotropy component            
for i = 1:n_ani
    mRF{1,i} = zeros(n_shell,1);
    for j = 1:n_shell
        R = sqrt(list_shell(j)*(longitudinal(i)-transverse(i)));
        mRF{1,i}(j,1) = (sqrt(pi)*exp(-list_shell(j)*transverse(i))*erf(R))/(2*R);
    end
end

% estimate anisotropy component
for i = 1:n_iso
    mRF{1,n_ani+i} = zeros(n_shell,1);
    for j = 1:n_shell
        mRF{1,n_ani+i}(j,1) = exp(-list_shell(j)*isotropy(i)); 
    end
end


% =========
% fitting SMSI (sparsity free)
% =========
mean_vfbasis = ones(1,n_iso+n_ani);

mean_fconv = cell2mat(mRF);

% training regularization parameter
mean_lambda = 0.01;
maxiter = 100;
error_variance = 1;

for iter = 1:maxiter
    Adj = (mean_fconv'*mean_fconv + mean_lambda * eye(n_iso+n_ani)) \ mean_fconv';
    k = abs(trace(-Adj * mean_fconv));

    if k <= 20
        mean_lambda = nmask * log(error_variance^2) + k*log(nmask);
        break;
    end
    mean_lambda = nmask * log(error_variance^2) + k*log(nmask);
end

Adj = ([mean_vfbasis;mean_fconv]'*[mean_vfbasis;mean_fconv] + mean_lambda * eye(n_iso+n_ani)) \ [mean_vfbasis;mean_fconv]';
mean_weight = max(0,Adj*[ones(1,nmask);mean_dwi_mask]);
%saving(mean_weight,ind,nx,ny,nz,'~/mean_weight.nii')
merr = norm(mean_fconv*mean_weight-mean_dwi_mask)^2;

% =========
% fitting RSI (sparsity free)
% =========
single_nsh = lmax2nsh(lmax);
nsh = n_ani * single_nsh + n_iso;
full_vfbasis = zeros(1,nsh);
full_vfbasis([0:n_ani-1]*single_nsh+1) = 1/scheme.sh(1);
full_vfbasis(n_ani*single_nsh+1:end) = 1/scheme.sh(1);

full_fconv = cell2mat(fRF);

% training regularization parameter
full_lambda = 0.001;
maxiter = 100;
error_variance = 1;

for iter = 1:maxiter
    Adj = (full_fconv'*full_fconv + full_lambda * eye(nsh)) \ full_fconv';
    k = abs(trace(-Adj * full_fconv));

    if k <= 20
        full_lambda = nmask * log(error_variance^2) + k*log(nmask);
        break;
    end
    full_lambda = nmask * log(error_variance^2) + k*log(nmask);
end
% Adj = ([full_vfbasis;full_fconv]'*[full_vfbasis;full_fconv] + full_lambda * eye(nsh)) \ [full_vfbasis;full_fconv]';
% sh = Adj*[ones(1,nmask);full_dwi_mask];
Adj = (full_fconv'*full_fconv + full_lambda * eye(nsh)) \ full_fconv';
sh = Adj*full_dwi_mask;
full_weight = max(0,sh([[0:n_ani-1]*single_nsh+1 n_ani*single_nsh+1:end],:) * scheme.sh(1));
ferr = norm(full_fconv*sh-full_dwi_mask)^2;

nv = length(scheme.vert);
A_aniso = repmat(scheme.sh,1,n_ani);
A_aniso = mat2cell(A_aniso,nv,repmat(size(scheme.sh,2),1,n_ani));
A_aniso = blkdiag(A_aniso{:});
ASH = blkdiag(A_aniso,diag(scheme.sh(1,1)*ones(1,n_iso)));
FOD = ASH * sh;
clear sh;
%saving(full_weight,ind,nx,ny,nz,'~/full_weight.nii')


% =========
% fitting SMSI (reweighted)
% =========
param.mode=1;
param.lambda=min(merr,ferr);
param.pos=true;
param.numThreads=-1;
double_weight = double(bsxfun(@power,bsxfun(@times,mean_weight,full_weight),0.5));
double_weight = double_weight./sum(double_weight);

alpha = mexLassoWeighted(double([ones(1,nmask);mean_dwi_mask]),...
                        double([mean_vfbasis;mean_fconv]),...
                        1./(double_weight+0.001),...
                        param);
weights = full(alpha);
saving(weights,ind,nx,ny,nz,'~/weights.nii')


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

function saving(masked_data,ind,nx,ny,nz,filename)
    data = single(zeros(size(masked_data,1),nx*ny*nz));
    data(:,ind) = masked_data;
    data = reshape(data',nx,ny,nz,size(masked_data,1));
    niftiwrite(data,filename);
end



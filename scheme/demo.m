clear all; clc;

scheme_hcp = importdata('HCP.txt');
bvecs = importdata('bvecs')';
bvals = importdata('bvals')';

bvals = round(bvals/10)*10;
bvecs = bvecs ./ vecnorm(bvecs,2,2);
TE = [75 85 95 105 115 125 135];

scheme = [];
for i = 1:length(TE)
    scheme = [scheme; [bvecs bvals]];
end

G = zeros(size(scheme,4),1);
G(scheme(:,4)==0) = 0;
G(scheme(:,4)==400) = 0.0275;
G(scheme(:,4)==800) = 0.0561;
G(scheme(:,4)==1600) = 0.0793;
G(scheme(:,4)==3200) = 0.0971;
scheme(:,5) = G;

scheme(:,6) = 0.0431;
scheme(:,7) = 0.0106;

for i = 1:length(TE)
    scheme((i-1)*length(bvals)+1:i*length(bvals),8) = TE(i);
end
bvals = scheme(:,4);

writematrix(scheme,'MTE.scheme','FileType','text','Delimiter','tab');
writematrix(bvals,'MTE.bvals','FileType','text','Delimiter','tab');

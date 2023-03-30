clear all; clc;
Path = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc/ME_SMSIx';
F_T2 = fullfile(Path,'T2.nii.gz');
F_Free = fullfile(Path,'VF_free.nii.gz');
F_Hindered = fullfile(Path,'VF_hindered.nii.gz');
F_Restricted = fullfile(Path,'VF_restricted.nii.gz');
info_T2 = niftiinfo(F_T2);
data_T2 = niftiread(info_T2);
B = [400 800 1600 3200];

List_compont = {'VF_restricted.nii.gz','VF_hindered.nii.gz','VF_free.nii.gz'};
for k = 1:2
    
    info_Restricted = niftiinfo(fullfile(Path,List_compont{k}));
    data_Restricted = niftiread(info_Restricted);
    
    List_tissue = {'tissue_SWM_int_rs.nii.gz','tissue_WM_int_rs.nii.gz','tissue_GM_int_rs.nii.gz','tissue_SGM_int_rs.nii.gz'};
    scheme = load('scheme/new_spectrum.mat');
    
    for i = 1:length(List_tissue)
        F_Mask = fullfile(Path,'..',List_tissue{i});
        
        info_Mask = niftiinfo(F_Mask);
        mask = round(niftiread(info_Mask));
    
        ind = mask > 0.5;
        num = 1;
        for j = (k-1)*4+1:k*4
            T2 = 1./squeeze(data_T2(:,:,:,j));
            T2 = round(T2(ind),4);
            switch k 
                case 1
                    D1 = sum(data_Restricted .* reshape(scheme.adc_restricted(:,1),1,1,1,scheme.num_restricted),4);
                    D2 = sum(data_Restricted .* reshape(scheme.adc_restricted(:,2),1,1,1,scheme.num_restricted),4);
                case 2
                    D1 = sum(data_Restricted .* reshape(scheme.adc_hindered(:,1),1,1,1,scheme.num_hindered),4);
                    D2 = sum(data_Restricted .* reshape(scheme.adc_hindered(:,2),1,1,1,scheme.num_hindered),4);
            end

            D1 = round(D1(ind),4);
            D2 = round(D2(ind),4);
    
            ind2 = ~isnan(T2) & ~isnan(D1) & ~isinf(T2) & ~isinf(D1) & ~isnan(D2) & ~isinf(D2);
    
            D1 = D1(ind2);
            D2 = D2(ind2);
            T2 = T2(ind2);
            
            ind3 = abs(zscore(T2))<1 & abs(zscore(D1))<1 & abs(zscore(D2))<1;
            T2 = T2(ind3);
            D1 = D1(ind3);
            D2 = D2(ind3);
            
            figure
            kscontour([T2 D1], 'Color', 'blue'); hold on;
            kscontour([T2 D2], 'Color', 'orange');
            fig=gcf;
            fig.Position(3:4)=[300,250];
            box on;
            set(gca,'linewidth',1)

            print(gcf,fullfile('T2_d',strcat(List_compont{k}(1:end-7),'_',List_tissue{i}(1:end-7),'_',num2str(B(num)),'.png')),'-dpng');
            num = num + 1;
            close;
        end
    end
end
            


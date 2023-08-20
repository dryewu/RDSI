function plotODF(F_img,F_SH,F_mask,ROI,options)
    %{
    ░█▀▀█ ░█▀▀▀█ ░█▀▄▀█ ░█▀▀▀ ░█▀▀▄ ▀█▀
    ░█─── ░█──░█ ░█░█░█ ░█▀▀▀ ░█─░█ ░█─
    ░█▄▄█ ░█▄▄▄█ ░█──░█ ░█▄▄▄ ░█▄▄▀ ▄█▄

    Visualization of orientation distribution function


        Created by Ye Wu, PhD (dr.yewu@outlook.com)

        - Nanjing University of Science and Technology, China
        - University of North Carolina at Chapel Hill, USA
        
    %}

    arguments
        F_img                   string   {mustBeFile}
        F_SH                    string   {mustBeNonzeroLengthText}  % [x,y,z,num*nsh]
        F_mask                  string   {mustBeFile}                         
        ROI                     (3,2)    {mustBeInteger,mustBeNonnegative}% [x_star x_end; y_star y_end, z_star,z_end]

        options.scale           (1,1)    {mustBeNonnegative}  = 1
    end

    addpath('third/csd')
    img_info = niftiinfo(F_img);
    img      = niftiread(img_info);

    mask_info = niftiinfo(F_mask);
    mask      = round(niftiread(mask_info))>0.5;

    SH_info   = cellfun(@(x)niftiinfo(x),F_SH,'UniformOutput',false);
    SH        = cellfun(@(x)niftiread(x),SH_info,'UniformOutput',false);

    x_star = ROI(1,1); x_end = ROI(1,2); 
    y_star = ROI(2,1); y_end = ROI(2,2); 
    z_star = ROI(3,1); z_end = ROI(3,2); 

    img_roi = img(x_star:x_end,y_star:y_end,z_star:z_end,1);
    SH = sum(cat(5,SH{:}),5) * options.scale;
    order = round((sqrt(8*size(SH,4)+1)-3)/2);
    scheme = gen_scheme('scheme/sphere_362_vertices.txt',order);

    figure;imagesc(img_roi);colormap(gray); hold on;

    for x = x_star:x_end
        for y = y_star:y_end
            for z = z_star:z_end
                if mask(x,y,z)
                    fod = max(0,SH2amp(squeeze(SH(x,y,z,:)),scheme));
                    trimesh(scheme.mesh, scheme.vert(:,1).*fod+x-x_star+1, scheme.vert(:,2).*fod+y-y_star+1, scheme.vert(:,3).*fod+z-z_star+1, 'edgecolor','none','FaceColor','interp','FaceVertexCData',abs(scheme.vert),'specularstrength',0.11);
                    hold on;
                end
            end
        end
    end

    axis equal; axis off; grid off; view(0,90);light('position',[0,0,1],'Style','infinite');
    lighting gouraud;light('Color',[ 0.8 0.8 0.8]);shading interp;set(gcf,'color','k');
    
end
    
    























    
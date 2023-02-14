cli=/home/wuye/Softwares/pnlNipype/scripts
Path=/home/wuye/D_disk/MTE_dMRI
for i in XYQ HTY WDZ ZJ CW LXF  
do
    for j in 75 85 95 105 115 125 135
    do
    dwifslpreproc ${Path}/${i}/proc2/MTE_${j}_align_center_denoise_unring.mif ${Path}/${i}/proc2/MTE_${j}_align_center_denoise_unring_preproc.nii.gz -export_grad_fsl ${Path}/${i}/proc2/MTE_${j}_align_center_denoise_unring_preproc.bvecs ${Path}/${i}/proc2/MTE_${j}_align_center_denoise_unring_preproc.bvals -pe_dir PA -rpe_none -eddy_options " --slm=linear" -json_import ${Path}/${i}/proc2/MTE_${j}.json
    done 
done 

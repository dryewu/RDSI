cli=/home/wuye/Softwares/pnlNipype/scripts
Path=/home/wuye/D_disk/MTE_dMRI
for i in CW  HTY  LXF  WDZ  XYQ  ZJ
do

    for j in 75 85 95 105 115 125 135
    do
		dwi2response dhollander ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.nii.gz ${Path}/${i}/proc/MTE_${j}_dhollander_rf_wm.txt ${Path}/${i}/proc/MTE_${j}_dhollander_rf_gm.txt ${Path}/${i}/proc/MTE_${j}_dhollander_rf_csf.txt -mask ${Path}/${i}/proc/MTE_mask.nii.gz -fslgrad ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.bvec ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.bval -force
		dwi2fod msmt_csd ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.nii.gz ${Path}/${i}/proc/MTE_${j}_dhollander_rf_wm.txt ${Path}/${i}/proc/MTE_${j}_dhollander_fod_wm.nii.gz ${Path}/${i}/proc/MTE_${j}_dhollander_rf_gm.txt ${Path}/${i}/proc/MTE_${j}_dhollander_fod_gm.nii.gz ${Path}/${i}/proc/MTE_${j}_dhollander_rf_csf.txt ${Path}/${i}/proc/MTE_${j}_dhollander_fod_csf.nii.gz -mask ${Path}/${i}/proc/MTE_mask.nii.gz -fslgrad ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.bvec ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.bval -force
		dwi2tensor ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.nii.gz ${Path}/${i}/proc/MTE_${j}_tensor.mif -mask ${Path}/${i}/proc/MTE_mask.nii.gz -fslgrad ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.bvec ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.bval -force
		tensor2metric ${Path}/${i}/proc/MTE_${j}_tensor.mif -fa ${Path}/${i}/proc/MTE_${j}_fa.nii.gz -vector ${Path}/${i}/proc/MTE_${j}_dec.mif -force
	done
done 

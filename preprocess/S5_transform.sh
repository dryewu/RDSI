cli=/home/wuye/Softwares/pnlNipype/scripts
Path=/home/wuye/D_disk/MTE_dMRI
for i in XYQ HTY WDZ ZJ CW LXF  
do
    rm ${Path}/${i}/proc/*.mat
    rm ${Path}/${i}/proc/*.txt
    rm ${Path}/${i}/proc/*in_T1*
    rm ${Path}/${i}/proc/t2_align_center_unring_unbiased_brain.nii.gz
	mrcalc ${Path}/${i}/proc/t2_align_center_unring_unbiased.nii.gz ${Path}/${i}/proc/brain_t2_mask.nii.gz --multiply ${Path}/${i}/proc/t2_align_center_unring_unbiased_brain.nii.gz
    for j in 75 85 95 105 115 125 135
    do
		flirt -in ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2_b0_brain.nii.gz -ref ${Path}/${i}/proc/t2_align_center_unring_unbiased_brain.nii.gz -omat ${Path}/${i}/proc/MTE_${j}_DWI2T1_flirt.mat -dof 6 -bins 256 -cost corratio -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -interp trilinear
		flirt -in ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2_b0_brain.nii.gz -ref ${Path}/${i}/proc/t2_align_center_unring_unbiased_brain.nii.gz -omat ${Path}/${i}/proc/MTE_${j}_DWI2T1_flirt_reg.mat -applyxfm -init ${Path}/${i}/proc/MTE_${j}_DWI2T1_flirt.mat -dof 6 -bins 256 -cost corratio -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -interp trilinear
		transformconvert ${Path}/${i}/proc/MTE_${j}_DWI2T1_flirt_reg.mat ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2_b0_brain.nii.gz ${Path}/${i}/proc/t2_align_center_unring_unbiased_brain.nii.gz flirt_import ${Path}/${i}/proc/MTE_${j}_DWI2T1_mrtrix.txt

		mrtransform ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2_b0_brain.nii.gz ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2_b0_brain_in_T1.nii.gz -linear ${Path}/${i}/proc/MTE_${j}_DWI2T1_mrtrix.txt
		mrtransform ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2.nii.gz ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.nii.gz -linear ${Path}/${i}/proc/MTE_${j}_DWI2T1_mrtrix.txt -fslgrad ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2.bvec  ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2.bval -export_grad_fsl ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bvec ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2_in_T1.bval
		mrtransform ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_mask.nii.gz ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_mask_in_T1.nii.gz -linear ${Path}/${i}/proc/MTE_${j}_DWI2T1_mrtrix.txt -interp nearest
    done 
    mrmath ${Path}/${i}/proc/MTE_*_align_center_denoise_unring_preproc_epi_mask_in_T1.nii.gz min ${Path}/${i}/proc/MTE_mask.nii.gz -force
done 

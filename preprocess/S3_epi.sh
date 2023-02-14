cli=/home/wuye/Softwares/pnlNipype/scripts
Path=/home/wuye/D_disk/MTE_dMRI
for i in XYQ HTY WDZ ZJ CW LXF  
do
    for j in 75 85 95 105 115 125 135
    do
    	dwiextract ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.nii.gz -fslgrad ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.bvec ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.bval - -bzero | mrmath - mean ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_b0.nii.gz -axis 3
    	nifti_bet_mask --bvals ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.bval -i ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.nii.gz -o ${Path}/${i}/proc/MTE_${j}
	mrcalc ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_b0.nii.gz ${Path}/${i}/proc/MTE_${j}_mask.nii.gz --multiply ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_b0_brain.nii.gz
	if [ -f ${Path}/${i}/proc/t2_like_align_center_unring_unbiased.nii.gz ]; then
        pnl_epi --bse ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_b0_brain.nii.gz --bvals ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.bval --bvecs ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.bvec --dwi ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.nii.gz --dwimask ${Path}/${i}/proc/MTE_${j}_mask.nii.gz -n 2 -o ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi --t2 ${Path}/${i}/proc/t2_like_align_center_unring_unbiased.nii.gz --t2mask ${Path}/${i}/proc/brain_t1_mask.nii.gz
    else
        pnl_epi --bse ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_b0_brain.nii.gz --bvals ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.bval --bvecs ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.bvec --dwi ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.nii.gz --dwimask ${Path}/${i}/proc/MTE_${j}_mask.nii.gz -n 2 -o ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi --t2 ${Path}/${i}/proc/t2_align_center_unring_unbiased.nii.gz --t2mask ${Path}/${i}/proc/brain_t2_mask.nii.gz
    fi

	dwibiascorrect ants ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi.nii.gz ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2.nii.gz -fslgrad ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi.bvec ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi.bval 
	dwigradcheck ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2.nii.gz -fslgrad ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi.bvec ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi.bval -export_grad_fsl ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2.bvec ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2.bval -mask ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_mask.nii.gz


    dwiextract ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2.nii.gz -fslgrad ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2.bvec ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2.bval - -bzero | mrmath - mean ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2_b0.nii.gz -axis 3

	mrcalc ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2_b0.nii.gz ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_mask.nii.gz --multiply ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased2_b0_brain.nii.gz
    done 
done 

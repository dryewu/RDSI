cli=/home/wuye/Softwares/pnlNipype/scripts
Path=/home/wuye/D_disk/MTE_dMRI
for i in XYQ HTY WDZ ZJ CW LXF  
do
    ${cli}/align.py --axisAlign --center -i ${Path}/${i}/proc/T1w.nii.gz -o ${Path}/${i}/proc/t1_align_center
    ${cli}/align.py --axisAlign --center -i ${Path}/${i}/proc/T2w.nii.gz -o ${Path}/${i}/proc/t2_align_center

    mrdegibbs ${Path}/${i}/proc/t1_align_center.nii.gz ${Path}/${i}/proc/t1_align_center_unring.nii.gz -axes 0,1
    mrdegibbs ${Path}/${i}/proc/t2_align_center.nii.gz ${Path}/${i}/proc/t2_align_center_unring.nii.gz -axes 0,1

    N4BiasFieldCorrection -d 3 -i ${Path}/${i}/proc/t1_align_center_unring.nii.gz -o ${Path}/${i}/proc/t1_align_center_unring_unbiased.nii.gz
    N4BiasFieldCorrection -d 3 -i ${Path}/${i}/proc/t2_align_center_unring.nii.gz -o ${Path}/${i}/proc/t2_align_center_unring_unbiased.nii.gz

    nifti_atlas -t ${Path}/${i}/proc/t1_align_center_unring_unbiased.nii.gz -o ${Path}/${i}/proc/brain_t1 -n 6 --train ~/Softwares/pnl/trainingDataT1AHCC-8141805/trainingDataT1Masks-hdr.csv
    nifti_atlas -t ${Path}/${i}/proc/t2_align_center_unring_unbiased.nii.gz -o ${Path}/${i}/proc/brain_t2 -n 6 --train ~/Softwares/pnl/trainingDataT2Masks-12a14d9/trainingDataT2Masks-hdr.csv

    for j in 75 85 95 105 115 125 135
    do
    
		${cli}/align.py --axisAlign --center -i ${Path}/${i}/proc/MTE_${j}.nii.gz --bvals ${Path}/${i}/proc/MTE_${j}.bval --bvecs ${Path}/${i}/proc/MTE_${j}.bvec -o ${Path}/${i}/proc/MTE_${j}_align_center
		dwidenoise ${Path}/${i}/proc/MTE_${j}_align_center.nii.gz ${Path}/${i}/proc/MTE_${j}_align_center_denoise.nii.gz
		mrdegibbs ${Path}/${i}/proc/MTE_${j}_align_center_denoise.nii.gz ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring.nii.gz -axes 0,1
		mrconvert ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring.nii.gz ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring.mif -fslgrad ${Path}/${i}/proc/MTE_${j}.bvec ${Path}/${i}/proc/MTE_${j}.bval -json_import ${Path}/${i}/proc/MTE_${j}.json
		dwifslpreproc ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring.mif ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc.nii.gz -export_grad_fsl ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc.bvecs ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc.bvals -pe_dir PA -rpe_none -eddy_options " --slm=linear" -json_import ${Path}/${i}/proc/MTE_${j}.json
		dwibiascorrect ants ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc.nii.gz ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.nii.gz -fslgrad ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc.bvecs ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc.bvals
		dwigradcheck ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.nii.gz -fslgrad ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc.bvecs ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc.bvals -export_grad_fsl ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.bvec ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.bval
    	nifti_bet_mask --bvals ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.bval -i ${Path}/${i}/proc/MTE_${j}_align_center_denoise_unring_preproc_epi_unbiased.nii.gz -o ${Path}/${i}/proc/MTE_${j}_2

    done 
 
    mrmath ${Path}/${i}/proc/MTE_*_mask.nii.gz min ${Path}/${i}/proc/MTE_mask.nii.gz -force
done 

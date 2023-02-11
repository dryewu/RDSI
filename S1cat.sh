cli=/home/wuye/Softwares/pnlNipype/scripts
Path=/home/wuye/D_disk/MTE_dMRI
for i in XYQ HTY WDZ ZJ CW LXF  
do
    #${cli}/align.py --axisAlign --center -i ${Path}/${i}/proc2/T1w.nii.gz -o ${Path}/${i}/proc2/t1_align_center
    #${cli}/align.py --axisAlign --center -i ${Path}/${i}/proc2/T2w.nii.gz -o ${Path}/${i}/proc2/t2_align_center

    #mrdegibbs ${Path}/${i}/proc2/t1_align_center.nii.gz ${Path}/${i}/proc2/t1_align_center_unring.nii.gz -axes 0,1
    #mrdegibbs ${Path}/${i}/proc2/t2_align_center.nii.gz ${Path}/${i}/proc2/t2_align_center_unring.nii.gz -axes 0,1

    #N4BiasFieldCorrection -d 3 -i ${Path}/${i}/proc2/t1_align_center_unring.nii.gz -o ${Path}/${i}/proc2/t1_align_center_unring_unbiased.nii.gz
    #N4BiasFieldCorrection -d 3 -i ${Path}/${i}/proc2/t2_align_center_unring.nii.gz -o ${Path}/${i}/proc2/t2_align_center_unring_unbiased.nii.gz

    #nifti_atlas -t ${Path}/${i}/proc2/t1_align_center_unring_unbiased.nii.gz -o ${Path}/${i}/proc2/brain_t1 -n 6 --train ~/Softwares/pnl/trainingDataT1AHCC-8141805/trainingDataT1Masks-hdr.csv
    #nifti_atlas -t ${Path}/${i}/proc2/t2_align_center_unring_unbiased.nii.gz -o ${Path}/${i}/proc2/brain_t2 -n 6 --train ~/Softwares/pnl/trainingDataT2Masks-12a14d9/trainingDataT2Masks-hdr.csv

    #for j in 75 85 95 105 115 125 135
    #do
    #    ${cli}/align.py --axisAlign --center -i ${Path}/${i}/proc2/MTE_${j}.nii.gz --bvals ${Path}/${i}/proc2/MTE_${j}.bval --bvecs ${Path}/${i}/proc2/MTE_${j}.bvec -o ${Path}/${i}/proc2/MTE_${j}_align_center
    #    dwidenoise ${Path}/${i}/proc2/MTE_${j}_align_center.nii.gz ${Path}/${i}/proc2/MTE_${j}_align_center_denoise.nii.gz
    #    mrdegibbs ${Path}/${i}/proc2/MTE_${j}_align_center_denoise.nii.gz ${Path}/${i}/proc2/MTE_${j}_align_center_denoise_unring.nii.gz -axes 0,1
    #    mrconvert ${Path}/${i}/proc2/MTE_${j}_align_center_denoise_unring.nii.gz ${Path}/${i}/proc2/MTE_${j}_align_center_denoise_unring.mif -fslgrad ${Path}/${i}/proc2/MTE_${j}.bvec ${Path}/${i}/proc2/MTE_${j}.bval -json_import ${Path}/${i}/proc2/MTE_${j}.json
    #    dwiextract ${Path}/${i}/proc2/MTE_${j}_align_center_denoise_unring.mif - -bzero | mrmath - mean ${Path}/${i}/proc2/MTE_${j}_align_center_denoise_unring_b0.mif -axis 3
    #done 

    #list1=''
    #list2=''
    #for j in 75 85 95 105 115 125 135
    #do
    #    if [ -f ${Path}/${i}/proc2/MTE_${j}_align_center_denoise_unring.mif ]; then
    #        list1=`echo ${list1} ${Path}/${i}/proc2/MTE_${j}_align_center_denoise_unring.mif`
    #        list2=`echo ${list2} ${Path}/${i}/proc2/MTE_${j}_align_center_denoise_unring_b0.mif`
    #    fi
    #done

    #mrcat ${list1} ${Path}/${i}/proc2/all_dwis.mif -axis 3
    #mrcat ${list2} ${Path}/${i}/proc2/all_b0s.mif -axis 3

	dwifslpreproc ${Path}/${i}/proc2/all_dwis.mif ${Path}/${i}/proc2/all_dwis_preproc.mif -pe_dir PA -rpe_none -eddy_options " --slm=linear" #-se_epi ${Path}/${i}/proc2/all_b0s.mif

done 

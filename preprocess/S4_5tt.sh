cli=/home/wuye/Softwares/pnlNipype/scripts
Path=/home/wuye/D_disk/MTE_dMRI
for i in XYQ HTY WDZ ZJ CW LXF  
do
    for j in 75 85 95 105 115 125 135
    do
	5ttgen fsl ${Path}/${i}/proc/t1_align_center_unring_unbiased.nii.gz ${Path}/${i}/proc/5TT_t1.mif -nocrop -sgm_amyg_hipp -mask ${Path}/${i}/proc/brain_t1_mask.nii.gz
	5ttgen fsl ${Path}/${i}/proc/t2_align_center_unring_unbiased.nii.gz ${Path}/${i}/proc/5TT_t2.mif -nocrop -sgm_amyg_hipp -mask ${Path}/${i}/proc/brain_t2_mask.nii.gz

    done 
done 

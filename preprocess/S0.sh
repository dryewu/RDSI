cli=/home/wuye/Softwares/mricron/Resources/dcm2niix 
Path=/home/wuye/D_disk/MTE_dMRI/XYQ
rm -r ${Path}/proc2
mkdir -p ${Path}/proc2
for i in 75 85 95 105 115 125 135
do
   ${cli} -f MTE_${i} -p y -z y -o ${Path}/proc2 ${Path}/*-PA_${i}*sDKI_opt_40
done 
${cli} -f T1w -p y -z y -o ${Path}/proc2 ${Path}/201-cs_T1W_3D_TFE_32channel
${cli} -f T2w -p y -z y -o ${Path}/proc2 ${Path}/301-3D_Brain_View_FLAIR_SHC

##

cli=/home/wuye/Softwares/mricron/Resources/dcm2niix
Path=/home/wuye/D_disk/MTE_dMRI/HTY
rm -r ${Path}/proc2
mkdir -p ${Path}/proc2
for i in 75 85 95 105 115 125 135
do
   ${cli} -f MTE_${i} -p y -z y -o ${Path}/proc2 ${Path}/*-PA_${i}*sDKI_opt_40*
done 
${cli} -f T1w -p y -z y -o ${Path}/proc2 ${Path}/301-cs_T1W_3D_TFE_32_channel
${cli} -f T2w -p y -z y -o ${Path}/proc2 ${Path}/201-cs3D_T2_flair_Brain_View

##

cli=/home/wuye/Softwares/mricron/Resources/dcm2niix
Path=/home/wuye/D_disk/MTE_dMRI/WDZ
rm -r ${Path}/proc2
mkdir -p ${Path}/proc2
for i in 75 85 95 105 115 125 135
do

   ${cli} -f MTE_${i} -p y -z y -o ${Path}/proc2 ${Path}/*-PA_${i}*sDKI_opt_40*

done  
${cli} -f T1w -p y -z y -o ${Path}/proc2 ${Path}/301-cs_T1W_3D_TFE_32_channel
${cli} -f T2w -p y -z y -o ${Path}/proc2 ${Path}/201-cs3D_T2_flair_Brain_View

##

cli=/home/wuye/Softwares/mricron/Resources/dcm2niix
Path=/home/wuye/D_disk/MTE_dMRI/ZJ
rm -r ${Path}/proc2
mkdir -p ${Path}/proc2
for i in 75 85 95 105 115 125 135
do

   ${cli} -f MTE_${i} -p y -z y -o ${Path}/proc2 ${Path}/*-PA_${i}*sDKI_opt_40*

done 
${cli} -f T1w -p y -z y -o ${Path}/proc2 ${Path}/301-cs_T1W_3D_TFE_32_channel
${cli} -f T2w -p y -z y -o ${Path}/proc2 ${Path}/201-cs3D_T2_flair_Brain_View

##

cli=/home/wuye/Softwares/mricron/Resources/dcm2niix
Path=/home/wuye/D_disk/MTE_dMRI/CW
rm -r ${Path}/proc2
mkdir -p ${Path}/proc2
for i in 75 85 95 105 115 125 135
do

   ${cli} -f MTE_${i} -p y -z y -o ${Path}/proc2 ${Path}/*-PA_${i}*sDKI_opt_40*

done 
${cli} -f T1w -p y -z y -o ${Path}/proc2 ${Path}/1301-cs_T1W_3D_TFE_new 
${cli} -f T2w -p y -z y -o ${Path}/proc2 ${Path}/201-T2W_DRIVE

##

cli=/home/wuye/Softwares/mricron/Resources/dcm2niix
Path=/home/wuye/D_disk/MTE_dMRI/LXF
rm -r ${Path}/proc2
mkdir -p ${Path}/proc2
for i in 75 85 95 105 115 125 135
do

   ${cli} -f MTE_${i} -p y -z y -o ${Path}/proc2 ${Path}/*-PA_${i}*sDKI_opt_40*

done 
${cli} -f T2w -p y -z y -o ${Path}/proc2 ${Path}/901-3D_Brain_View_FLAIR_SHC





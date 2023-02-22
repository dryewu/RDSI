clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/XYQ/proc';
for TE = [75 85 95 105 115]
    SE_RSI(dataFolder,TE);
end

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/HTY/proc';
for TE = [85 95 105 115 125 135]
    SE_RSI(dataFolder,TE);
end
  
clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/WDZ/proc';
for TE = [75 85 95 105 115 125 135]
    SE_RSI(dataFolder,TE);
end

clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/ZJ/proc';
for TE = [75 85 95 105 115 125 135]
    SE_RSI(dataFolder,TE);
end
clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/CW/proc';
for TE = [85 95 105 115 125 135]
    SE_RSI(dataFolder,TE);
end
clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/LXF/proc';
for TE = [75 85 95 105 115 125 135]
    SE_RSI(dataFolder,TE);
end
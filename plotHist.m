clear all; clc;
dataFolder = '/home/wuye/D_disk/MTE_dMRI/WDZ/proc';
fT2    = fullfile(dataFolder,'ME_SMSIx_default','T2.nii.gz');
info = niftiinfo(fT2);
data = niftiread(info);

data = reshape(data,size(data,1)*size(data,2)*size(data,3),size(data,4));
ind = all(data>20 & data<350,2);
data = data(ind,:);

for i = 1:12
    [N(i,:),edges(i,:)] = histcounts(data(:,i),100);
end

b = ones(12,100);
b([1 5 9],:) = 400;
b([1 5 9]+1,:) = 800;
b([1 5 9]+2,:) = 1600;
b([1 5 9]+3,:) = 3200;

edges = edges(:,2:end);

figure;
for i = 1:4
    plot3(edges(i,:),b(i,:),N(i,:),'LineWidth',3); hold on;
end
grid on
xticks([0 100 200 300]);
yticks([400 800 1600 3200]);
xlim([0 300]);
view(40,40);
fig=gcf;
fig.Position(3:4)=[500,400];
print(gcf,fullfile('histogram','WDZ','restricted.png'),'-dpng'); close;

figure;
for i = 5:8
    plot3(edges(i,:),b(i,:),N(i,:),'LineWidth',3); hold on;
end
grid on
xticks([0 50 100 150 200 250]);
yticks([400 800 1600 3200]);
xlim([0 250]);
view(40,40);
fig=gcf;
fig.Position(3:4)=[500,400];
print(gcf,fullfile('histogram','WDZ','hindered.png'),'-dpng'); close;


figure;
for i = 9:12
    plot3(edges(i,:),b(i,:),N(i,:),'LineWidth',3); hold on;
end
grid on
xticks([0 100 200]);
yticks([400 800 1600 3200]);
xlim([0 200]);
view(40,40);
fig=gcf;
fig.Position(3:4)=[500,400];
print(gcf,fullfile('histogram','WDZ','free.png'),'-dpng'); close;

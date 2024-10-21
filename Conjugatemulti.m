close all;
clear;
clc;

csi_trace = read_bf_file('sample_data/data_phase.dat');
size0=length(csi_trace);

% antenna1 = [];
% antenna2 = [];
% antenna3 = [];

timestamp(size0,1) = 0;

for j=1:size0
    timestamp(j,1) = csi_trace{j}.timestamp_low;
end

%%%CSI data on each receiving antenna%%%
antenna1(size0,30) = 0;
antenna2(size0,30) = 0;
antenna3(size0,30) = 0;
k=1;

 for i=1:size0
     
     csi = squeeze(get_scaled_csi(csi_trace{i}));
     
%      antenna1= [antenna1; csi(1,:)];
%      antenna2= [antenna2; csi(2,:)];
%      antenna3= [antenna3; csi(3,:)];

    antenna1(k,:)=csi(1,:);
    antenna2(k,:)=csi(2,:);
    antenna3(k,:)=csi(3,:);

    k=k+1;
     
 end

%%%%%Conjugate multiplication (denoising)%%%%
output1=antenna1.*conj(antenna2);
output2=antenna2.*conj(antenna1);
output3=antenna1.*conj(antenna3);
output4=antenna3.*conj(antenna1);
output5=antenna2.*conj(antenna3);
output6=antenna3.*conj(antenna2);

figure;
plot(smooth(angle(output1(:,1)),10,'lowess'));
figure;
plot(smooth(angle(output2(:,1)),10,'lowess'));
figure;
plot(smooth(angle(output3(:,1)),10,'lowess'));
figure;
plot(smooth(angle(output4(:,1)),10,'lowess'));
figure;
plot(smooth(angle(output5(:,1)),10,'lowess'));
figure;
plot(smooth(angle(output6(:,1)),10,'lowess'));


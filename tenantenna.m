% 计算子载波的频率间隔
sub_freq_delta = (40 * 10^6) / 30; % 40 MHz/30
fc = 5.32e9; % center frequency 5.32GHz
M = 3;    % number of rx antennas
fs = 40e6; % channel bandwidth
c = 3e8;  % speed of light
% dTx = 2.6e-2; 
SubCarrInd = [-58,-54,-50,-46,-42,-38,-34,-30,-26,-22,-18,-14,-10,-6,-2,2,6,10,14,18,22,26,30,34,38,42,46,50,54,58]; % WiFi subcarrier indices at which CSI is available
N = length(SubCarrInd); % number of subcarriers
% subCarrSize = 128;  % total number fo
fgap = 312.5e3; % frequency gap in Hz between successive subcarriers in WiFi
lambda = c/fc;  % wavelength
d = lambda/2;  % distance between adjacent antennas in the linear antenna array
T = 1; % number of transmitter antennas

% 读取CSI数据文件
csi_trace = read_bf_file('data1_1126_rx6.dat'); 
size0=length(csi_trace); % CSI中的数据包的数量
antenna1_card1(size0,30) = 0;
antenna2_card1(size0,30) = 0;
antenna3_card1(size0,30) = 0;
k=1;

 for i=1:size0
     
    csi = squeeze(get_scaled_csi(csi_trace{i})); % csi_shape:[3,30]
     
    antenna1_card1(k,:)=csi(1,:);
    antenna2_card1(k,:)=csi(2,:);
    antenna3_card1(k,:)=csi(3,:);

    k=k+1;
     
 end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
csi_trace2 = read_bf_file('data1_1126_rx5.dat');
size2=length(csi_trace2);
antenna1_card2(size2,30) = 0;
antenna2_card2(size2,30) = 0;
antenna3_card2(size2,30) = 0;
k=1;

 for i=1:size2
     
     csi = squeeze(get_scaled_csi(csi_trace2{i}));
     
    antenna1_card2(k,:)=csi(1,:);
    antenna2_card2(k,:)=csi(2,:);
    antenna3_card2(k,:)=csi(3,:);

    k=k+1;
     
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  csi_trace3 = read_bf_file('data1_1126_rx3.dat');
size3=length(csi_trace3);
antenna1_card3(size3,30) = 0;
antenna2_card3(size3,30) = 0;
antenna3_card3(size3,30) = 0;
k=1;

 for i=1:size3
     
     csi = squeeze(get_scaled_csi(csi_trace3{i}));
     
    antenna1_card3(k,:)=csi(1,:);
    antenna2_card3(k,:)=csi(2,:);
    antenna3_card3(k,:)=csi(3,:);

    k=k+1;
     
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  csi_trace4 = read_bf_file('data1_1126_rx7.dat');
size4=length(csi_trace4);
antenna1_card4(size4,30) = 0;
antenna2_card4(size4,30) = 0;
antenna3_card4(size4,30) = 0;
k=1;

 for i=1:size4
     
    csi = squeeze(get_scaled_csi(csi_trace4{i}));
     
    antenna1_card4(k,:)=csi(1,:);    
    antenna2_card4(k,:)=csi(2,:);
    antenna3_card4(k,:)=csi(3,:);

    k=k+1;
     
 end


 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSI=zeros(600,300,33);
for i=151:750
a = round(i*33.3);
tencsi=zeros(300,33);
for pack= a:(a+32)
%%
% 接收器1的三个天线的第pack行csi数据组合成一个矩阵
sample_csi_trace = [antenna1_card1(pack,:)'; antenna2_card1(pack,:)'; antenna3_card1(pack,:)'];
csi_plot = reshape(sample_csi_trace, N, M); % shape：30*3

[PhsSlope, PhsCons] = removePhsSlope(csi_plot,M,SubCarrInd,N); %去除相位斜率
ToMult = exp(1i* (-PhsSlope*repmat(SubCarrInd(:),1,M) - PhsCons*ones(N,M) ));  % 30*3
csi_plot = csi_plot.*ToMult;
relChannel_noSlope = reshape(csi_plot, N, M, T); % shape:30*3*1
sample_csi_trace_sanitized = relChannel_noSlope(:); %即将3D矩阵N*M*T变成一个一维向量。
% figure;
% plot(smooth(unwrap(angle(sample_csi_trace_sanitized(:,1)))));
%% 
% 接收器2的三个天线的第pack行csi数据组合成一个矩阵
sample_csi_trace2 = [antenna1_card2(pack,:)'; antenna2_card2(pack,:)'; antenna3_card2(pack,:)'];
csi_plot2 = reshape(sample_csi_trace2, N, M);
[PhsSlope2, PhsCons2] = removePhsSlope(csi_plot2,M,SubCarrInd,N); % ToF sanitization code (Algorithm 1 in SpotFi paper)
ToMult2 = exp(1i* (-PhsSlope2*repmat(SubCarrInd(:),1,M) - PhsCons2*ones(N,M) ));
csi_plot2 = csi_plot2.*ToMult2;
relChannel_noSlope2 = reshape(csi_plot2, N, M, T);
sample_csi_trace_sanitized2 = relChannel_noSlope2(:);
%%
% 接收器3的三个天线的第pack行csi数据组合成一个矩阵
sample_csi_trace3 = [antenna1_card3(pack,:)'; antenna2_card3(pack,:)'; antenna3_card3(pack,:)'];
csi_plot3 = reshape(sample_csi_trace3, N, M);
[PhsSlope3, PhsCons3] = removePhsSlope(csi_plot3,M,SubCarrInd,N);
ToMult3 = exp(1i* (-PhsSlope3*repmat(SubCarrInd(:),1,M) - PhsCons3*ones(N,M) ));
csi_plot3 = csi_plot3.*ToMult3;
relChannel_noSlope3 = reshape(csi_plot3, N, M, T);
sample_csi_trace_sanitized3 = relChannel_noSlope3(:);
%%
% 接收器4的三个天线的第pack行csi数据组合成一个矩阵
sample_csi_trace4 = [antenna1_card4(pack,:)'; antenna2_card4(pack,:)'; antenna3_card4(pack,:)'];
csi_plot4 = reshape(sample_csi_trace4, N, M);
[PhsSlope4, PhsCons4] = removePhsSlope(csi_plot4,M,SubCarrInd,N);
ToMult4 = exp(1i* (-PhsSlope4*repmat(SubCarrInd(:),1,M) - PhsCons3*ones(N,M) ));
csi_plot4 = csi_plot4.*ToMult4;
relChannel_noSlope4 = reshape(csi_plot4, N, M, T);
sample_csi_trace_sanitized4 = relChannel_noSlope4(:);
%%
    
  %%%%%%%%%%%%%%%%  
   csi_matrix1=sample_csi_trace_sanitized(1:30,:);
   csi_matrix2=sample_csi_trace_sanitized2(1:30,:);
  
   [newcsi_matrix, phaseshiftA, phaseshiftB] = removeTimeintroPhashift(csi_matrix1,csi_matrix2);
  
   csi_matrix2=newcsi_matrix;



    
    csi_matrix3=sample_csi_trace_sanitized(31:60,:);
    csi_matrix4=sample_csi_trace_sanitized2(31:60,:);
    csi_matrix5=sample_csi_trace_sanitized(61:90,:);
    csi_matrix6=sample_csi_trace_sanitized2(61:90,:);
    
    R4 = abs(csi_matrix4);
    phase_matrix4 = unwrap(angle(csi_matrix4));
    R6 = abs(csi_matrix6);
    phase_matrix6 = unwrap(angle(csi_matrix6));
for n=1:30
   phase_matrix4(n) = phase_matrix4(n) + SubCarrInd(n)*phaseshiftA+phaseshiftB;
   phase_matrix6(n) = phase_matrix6(n) + SubCarrInd(n)*phaseshiftA+phaseshiftB; 
end
    csi_matrix4 = R4 .* exp(1i * phase_matrix4);
    csi_matrix6 = R6 .* exp(1i * phase_matrix6);

    combinedcsi=[csi_matrix5 ;csi_matrix3 ;csi_matrix1 ;csi_matrix4 ;csi_matrix6];
%%
%     fivecsi(:,(pack-a+1)) = combinedcsi;
  %%%%%%%%%%%%%%%%  

  csi_matrix1=sample_csi_trace_sanitized3(1:30,:);
  csi_matrix2=sample_csi_trace_sanitized4(1:30,:);
  
  [newcsi_matrix, phaseshiftA, phaseshiftB] = removeTimeintroPhashift(csi_matrix1,csi_matrix2);
  
  csi_matrix2=newcsi_matrix;



    
    csi_matrix3=sample_csi_trace_sanitized3(31:60,:);
    csi_matrix4=sample_csi_trace_sanitized4(31:60,:);
    csi_matrix5=sample_csi_trace_sanitized3(61:90,:);
    csi_matrix6=sample_csi_trace_sanitized4(61:90,:);
    
    R4 = abs(csi_matrix4);
    phase_matrix4 = unwrap(angle(csi_matrix4));
    R6 = abs(csi_matrix6);
    phase_matrix6 = unwrap(angle(csi_matrix6));
for n=1:30
   phase_matrix4(n) = phase_matrix4(n) + SubCarrInd(n)*phaseshiftA+phaseshiftB;
   phase_matrix6(n) = phase_matrix6(n) + SubCarrInd(n)*phaseshiftA+phaseshiftB; 
end
    csi_matrix4 = R4 .* exp(1i * phase_matrix4);
    csi_matrix6 = R6 .* exp(1i * phase_matrix6);

    combinedcsi2=[csi_matrix5 ;csi_matrix3 ;csi_matrix1 ;csi_matrix4 ;csi_matrix6];
%%
%     fivecsi2(:,(pack-a+1)) = combinedcsi2;
    %%%%%%%%%%%%%%%% 
  csi_matrix1=combinedcsi(121:150,:);
  csi_matrix2=combinedcsi2(121:150,:);
  
  [newcsi_matrix, phaseshiftA, phaseshiftB] = removeTimeintroPhashift(csi_matrix1,csi_matrix2);
  
  csi_matrix2=newcsi_matrix;



    
    csi_matrix9=combinedcsi(1:30,:);
    csi_matrix10=combinedcsi2(1:30,:);
    csi_matrix7=combinedcsi(31:60,:);
    csi_matrix8=combinedcsi2(31:60,:);
    csi_matrix5=combinedcsi(61:90,:);
    csi_matrix6=combinedcsi2(61:90,:);
    csi_matrix3=combinedcsi(91:120,:);
    csi_matrix4=combinedcsi2(91:120,:);
    
    R4 = abs(csi_matrix4);
    phase_matrix4 = unwrap(angle(csi_matrix4));
    R6 = abs(csi_matrix6);
    phase_matrix6 = unwrap(angle(csi_matrix6));
    R8 = abs(csi_matrix8);
    phase_matrix8 = unwrap(angle(csi_matrix8));
    R10 = abs(csi_matrix10);
    phase_matrix10 = unwrap(angle(csi_matrix10));
for n=1:30
   phase_matrix4(n) = phase_matrix4(n) + SubCarrInd(n)*phaseshiftA+phaseshiftB;
   phase_matrix6(n) = phase_matrix6(n) + SubCarrInd(n)*phaseshiftA+phaseshiftB;
   phase_matrix8(n) = phase_matrix8(n) + SubCarrInd(n)*phaseshiftA+phaseshiftB;
   phase_matrix10(n) = phase_matrix10(n) + SubCarrInd(n)*phaseshiftA+phaseshiftB; 
end
    csi_matrix4 = R4 .* exp(1i * phase_matrix4);
    csi_matrix6 = R6 .* exp(1i * phase_matrix6);
    csi_matrix8 = R4 .* exp(1i * phase_matrix8);
    csi_matrix10 = R6 .* exp(1i * phase_matrix10);

    combinedCSI=[csi_matrix9; csi_matrix7; csi_matrix5; csi_matrix3; csi_matrix1; csi_matrix10; csi_matrix8; csi_matrix6; csi_matrix4; csi_matrix2];  
    tencsi(:,(pack-a+1)) = combinedCSI;
end
index=i-150;
CSI(index,:,:)=tencsi;
end
    



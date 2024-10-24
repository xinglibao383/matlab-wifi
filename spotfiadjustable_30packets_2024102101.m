clear
close all
clc
SignalEndIdx = 25;
sub_freq_delta = (20 * 10^6) / 30; 
frequency = 5.64e9;
M = 3;
fs = 40e6;
c = 3e8;
antenna_distance = 2.6e-2;
SubCarrInd = [-58,-54,-50,-46,-42,-38,-34,-30,-26,-22,-18,-14,-10,-6,-2,2,6,10,14,18,22,26,30,34,38,42,46,50,54,58];
N = length(SubCarrInd);
% fgap = 312.5e3;
lambda = c/frequency;
T = 1;
csi_trace = read_bf_file(fullfile('my_data', 'lb_01_01_01.dat'));  % 1122 * 1
size0=length(csi_trace);
antenna1_card1(size0,30) = 0;   % 1122 * 30
antenna2_card1(size0,30) = 0;   % 1122 * 30
antenna3_card1(size0,30) = 0;   % 1122 * 30

k=1;
for i=1:size0
    csi = squeeze(get_scaled_csi(csi_trace{i}));
    antenna1_card1(k,:)=csi(1,:);
    antenna2_card1(k,:)=csi(2,:);
    antenna3_card1(k,:)=csi(3,:);
    k=k+1;
end
CSI=zeros(50,90,20);
for i=1:20:1000
    tempcsi=zeros(90,20);
    for pack= i:(i+19)
        sample_csi_trace = [antenna1_card1(pack,:)'; antenna2_card1(pack,:)'; antenna3_card1(pack,:)'];
        csi_plot = reshape(sample_csi_trace, N, M);
        [PhsSlope, PhsCons] = removePhsSlope(csi_plot,M,SubCarrInd,N);
        %[PhsSlope, PhsCons] = removePhsSlope_20241021(csi_plot,M,SubCarrInd,N);
        ToMult = exp(1i* (-PhsSlope*repmat(SubCarrInd(:),1,M) - PhsCons*ones(N,M) ));
        csi_plot = csi_plot.*ToMult;
        relChannel_noSlope = reshape(csi_plot, N, M, T);
        sanitized_csi0 = relChannel_noSlope(:);
        tempcsi(:,(pack-i+1)) = sanitized_csi0;
    end
    CSI(round(i/20)+1,:,:)=tempcsi;
end
ind = 9;
x = squeeze(CSI(ind,:,:));
L = size(x);
global DEBUG_PATHS
global DEBUG_PATHS_LIGHT
global OUTPUT_AOAS
global OUTPUT_TOFS
global OUTPUT_AOA_MUSIC_PEAK_GRAPH
global OUTPUT_TOF_MUSIC_PEAK_GRAPH
global OUTPUT_AOA_TOF_MUSIC_PEAK_GRAPH
global OUTPUT_SELECTIVE_AOA_TOF_MUSIC_PEAK_GRAPH
global OUTPUT_BINARY_AOA_TOF_MUSIC_PEAK_GRAPH
global OUTPUT_SUPPRESSED
global OUTPUT_FIGURES_SUPPRESSED
data_name = '-';
R = x * x'/L(2);
[Utmp,D] = eig(R);
D = abs(D);
[Dtmp,I] = sort(diag(D), 'descend');
D = diag(Dtmp);
U = Utmp(:,I);
Qn = U(:,SignalEndIdx+1:end);
Qs = U(:,1:SignalEndIdx);
eigenvectors = Qn;
theta = -90:1:90;
tau = 0:(0.125 * 10^-9):(60 * 10^-9);
Pmusic = zeros(length(theta), length(tau));
for ii = 1:length(theta)
    for jj = 1:length(tau)
        steering_vector = compute_steering_vector90(theta(ii), tau(jj), frequency, sub_freq_delta, antenna_distance);
        PP = steering_vector' * (eigenvectors * eigenvectors') * steering_vector;
        Pmusic(ii, jj) = abs(1 /  PP);
    end
end

% 保存数据到.mat文件
save('aoa-tof.mat', 'tau', 'theta', 'Pmusic');

figure('Name', 'AoA & ToF MUSIC Peaks', 'NumberTitle', 'off')
surf(tau, theta, Pmusic)
xlabel('Time of Flight')
ylabel('Angle of Arrival in degrees')
zlabel('Spectrum Peaks')
title('AoA and ToF Estimation from Modified MUSIC Algorithm')
shading interp

function steering_vector = compute_steering_vector11(theta, tau, freq, sub_freq_delta, ant_dist)
    steering_vector = zeros(40, 1);
    k = 1;
    base_element = 1;
    for ii = 1:2
        for jj = 1:20
            steering_vector(k, 1) = base_element * omega_tof_phase11(tau, sub_freq_delta)^(jj - 1);
            k = k + 1;
        end
        base_element = base_element * phi_aoa_phase11(theta, freq, ant_dist);
    end
end

function steering_vector = compute_steering_vector90(theta, tau, freq, sub_freq_delta, ant_dist)
    steering_vector = zeros(90, 1);
    k = 1;
    base_element = 1;
    for ii = 1:3
        for jj = 1:30
            steering_vector(k, 1) = base_element * omega_tof_phase11(tau, sub_freq_delta)^(jj - 1);
            k = k + 1;
        end
        base_element = base_element * phi_aoa_phase11(theta, freq, ant_dist); 
    end
end

function time_phase = omega_tof_phase11(tau, sub_freq_delta)
    time_phase = exp(-1i * 2 * pi * sub_freq_delta * tau);
end

function angle_phase = phi_aoa_phase11(theta, frequency, d)
    c = 3.0 * 10^8;
    theta = theta / 180 * pi;
    angle_phase = exp(-1i * 2 * pi * d * sin(theta) * (frequency / c));
end

function smoothed_csi = smooth_csi11(csi)
    nn = 11;
    nnn = (30-nn+1)*2;
    smoothed_csi = zeros(nnn, 2*nn);
    for n=1:nn
        smoothed_csi(:,n)=reshape(csi(1:2,n:n+(30-nn)).',nnn,1);
        smoothed_csi(:,n+nn)=reshape(csi(2:3,n:n+(30-nn)).',nnn,1);
    end
end

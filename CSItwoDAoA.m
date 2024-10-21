clear
close all
clc

SignalEndIdx = 55; %%%number of signals, adjustable
antenna_distance = 2.6e-2;
    % frequency = 5 * 10^9;
    % frequency = 5.785 * 10^9;
    frequency = 5.32 * 10^9;
    sub_freq_delta = (40 * 10^6) / 30;

csi_trace = read_bf_file('sample_data/data_2DAoA.dat');  
if exist('OUTPUT_SUPPRESSED') == 0                           
        globals_init()
    end
    
    
        data_name = ' - ';
    num_packets = 2;
    aoa_packet_data = cell(num_packets, 1); tof_packet_data = cell(num_packets, 1);
    packet_one_phase_matrix1 = 0; packet_one_phase_matrix2 = 0; packet_one_phase_matrix3 = 0;% Do computations for packet one so the packet loop can be parallelized % Get CSI for current packet
    csi_entry = csi_trace{1410};  %packet number    adjustable         
    csi = get_scaled_csi(csi_entry);
    % Only consider measurements for transmitting on one antenna
    csi1 = csi(1, :, :);
    csi1 = squeeze(csi1);
    
    csi2 = csi(2, :, :);
    csi2 = squeeze(csi2);
    
    csi3 = csi(3, :, :);
    csi3 = squeeze(csi3);

    % Sanitize ToFs with Algorithm 1
    packet_one_phase_matrix1 = unwrap(angle(csi1), pi, 2);
    sanitized_csi1 = spotfi_algorithm_1(csi1, sub_freq_delta);
    
    packet_one_phase_matrix2 = unwrap(angle(csi2), pi, 2);
    sanitized_csi2 = spotfi_algorithm_1(csi2, sub_freq_delta);
    
    packet_one_phase_matrix3 = unwrap(angle(csi3), pi, 2);
    sanitized_csi3 = spotfi_algorithm_1(csi3, sub_freq_delta);
    
 % for abc = 1:6  
    abc = 5;
    % Acquire smoothed CSI matrix
    smoothed_sanitized_csi = smooth_csi_AoD(sanitized_csi1, sanitized_csi2, sanitized_csi3, abc);
    % Run SpotFi's AoA-ToF MUSIC algorithm on the smoothed and sanitized CSI matrix
%     [aoa_packet_data{1}, tof_packet_data{1}] = aoa_tof_music(...
%             smoothed_sanitized_csi, antenna_distance, frequency, sub_freq_delta, data_name);
x = smoothed_sanitized_csi;

  
    %% Run MUSIC algorithm with SpotFi method including ToF and AoA
% x                -- the signal matrix
% antenna_distance -- the distance between the antennas in the linear array
% frequency        -- the frequency of the signal being localized
% sub_freq_delta   -- the difference between subcarrier frequencies
% data_name        -- the name of the data file being operated on, used for labeling figures
% Return:
% estimated_aoas   -- the angle of arrivals that gave peaks from running MUSIC, as a vector
% estimated_tofs   -- the time of flights that gave peaks on the estimated_aoas from running music.
%                         This is a matrix with dimensions [length(estimated_aoas, ), length(tau)].
%                         The columns are zero padded at the ends to handle different peak counts 
%                           across different AoAs.
%                         I.E. if there are three AoAs then there will be three rows in 
%                           estimated_tofs
% function [estimated_aoas, estimated_tofs] = aoa_tof_music(x, ...
%         antenna_distance, frequency, sub_freq_delta, data_name)
 
    
    % Data covarivance matrix
    R = x * x'; 
    % Find the eigenvalues and eigenvectors of the covariance matrix
  %  Rxx = X1*X1'/n; 
% [EV,D] = eig(R); 
% [EVA,I] = sort(diag(D).'); 
% EV = fliplr(EV(:,I)); 
% Un = EV(:,SignalEndIdx+1:end); 
    
[Utmp,D] = eig(R);
D = abs(D);
[Dtmp,I] = sort(diag(D), 'descend');
D = diag(Dtmp);
U = Utmp(:,I);
Qn = U(:,SignalEndIdx+1:end);
Qs = U(:,1:SignalEndIdx);    
eigenvectors = Qn;       
    
    % Peak search
    % Angle in degrees (converts to radians in phase calculations)
    

% %     
    theta = 1:3:180; %elevation adjustable
    phi = 1:3:180; %azimuth adjustable
    tau =0:(3 * 10^-9):(150 * 10^-9);%ToF adjustable
    Pmusic = zeros(length(theta),length(phi),length(tau));
    % Angle of Arrival Loop (AoA)

    for ii = 1:length(theta)
        for kk = 1:length(phi)
        % Time of Flight Loop (ToF)
            for jj = 1:length(tau)
               % for hh = 1:length(alpha)
                steering_vector = compute_steering_vector_AoD(theta(ii), phi(kk), tau(jj), ...
                        frequency, sub_freq_delta, antenna_distance);
                PP = steering_vector' * eigenvectors * eigenvectors' * steering_vector;
                Pmusic(ii,kk,jj) = abs(1 / PP);
               % end
            end
        end
    end

%     PPmusic = zeros(length(theta),length(phi));
     PPstatic = sum(Pmusic,3);
%    

 %PPmusic = Pmusic(:,:, 30);
   % if OUTPUT_AOA_TOF_MUSIC_PEAK_GRAPH && ~OUTPUT_SUPPRESSED && ~OUTPUT_FIGURES_SUPPRESSED
        % Theta (AoA) & Tau (ToF) 3D Plot
        figure('Name', 'Azimuth & Elevation MUSIC Peaks', 'NumberTitle', 'off')
        mesh(theta,phi,PPstatic)
        xlabel('elevation')
        ylabel('azimuth')
        zlabel('Spectrum Peaks')
        title('Azimuth and Elevation Estimation from Modified MUSIC Algorithm')
        grid on
        view([0,90])
        
   % end

   

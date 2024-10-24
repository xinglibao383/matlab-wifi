clear
close all
clc

SignalEndIdx = 25; %%%number of signals

sub_freq_delta = (40 * 10^6) / 30;
frequency = 5.32e9; % center frequency
M = 3;    % number of rx antennas
fs = 40e6; % channel bandwidth
c = 3e8;  % speed of light
antenna_distance = 2.6e-2;  % distance between adjacent antennas in the linear antenna array
% dTx = 2.6e-2; 
SubCarrInd = [-58,-54,-50,-46,-42,-38,-34,-30,-26,-22,-18,-14,-10,-6,-2,2,6,10,14,18,22,26,30,34,38,42,46,50,54,58]; % WiFi subcarrier indices at which CSI is available
N = length(SubCarrInd); % number of subcarriers
% subCarrSize = 128;  % total number fo
fgap = 312.5e3; % frequency gap in Hz between successive subcarriers in WiFi
lambda = c/frequency;  % wavelength
T = 1; % number of transmitter antennas
% csi_trace = read_bf_file('data1_1126_rx3.dat'); %data_1DAoA.dat
% csi_trace = read_bf_file(fullfile('my_data', 'lb_01_01_01.dat')); %data_1DAoA.dat lb_01_01_01.dat

sample_csi_traceTmp = load('./my_data/S11_01_10'); %'sample_csi_trace'
% sample_csi_trace = sample_csi_traceTmp.sample_csi_trace;
sample_csi_trace0 = sample_csi_traceTmp.csi_out;
sample_csi_trace0 =sample_csi_trace0(1,:,:,:);
sample_csi_trace0 = reshape(sample_csi_trace0,3,30,20);
%sample_csi_trace = reshape(sample_csi_trace,90,1);
% s =size(sample_csi_trace);

size0=size(sample_csi_trace0,3);
antenna1_card1(size0,30) = 0;
antenna2_card1(size0,30) = 0;
antenna3_card1(size0,30) = 0;
k=1;

 for i=1:size0
     
     csi = sample_csi_trace0(:,:,i);
     
    antenna1_card1(k,:)=csi(1,:);
    antenna2_card1(k,:)=csi(2,:);
    antenna3_card1(k,:)=csi(3,:);

    k=k+1;
     
 end
 % pack=4500; %packet number
 %% 选择第几号包
 pack=15; %packet number

sample_csi_trace = [antenna1_card1(pack,:)'; antenna2_card1(pack,:)'; antenna3_card1(pack,:)'];
csi_plot = reshape(sample_csi_trace, N, M);
[PhsSlope, PhsCons] = removePhsSlope(csi_plot,M,SubCarrInd,N);
ToMult = exp(1i* (-PhsSlope*repmat(SubCarrInd(:),1,M) - PhsCons*ones(N,M) ));
csi_plot = csi_plot.*ToMult;
relChannel_noSlope = reshape(csi_plot, N, M, T);
sanitized_csi0 = relChannel_noSlope(:);
sanitized_csi = reshape(sanitized_csi0,30,3).';

    % Acquire smoothed CSI matrix
    smoothed_sanitized_csi = smooth_csi11(sanitized_csi);
    % Run SpotFi's AoA-ToF MUSIC algorithm on the smoothed and sanitized CSI matrix
%     [aoa_packet_data{1}, tof_packet_data{1}] = aoa_tof_music(...
%             smoothed_sanitized_csi, antenna_distance, frequency, sub_freq_delta, data_name);
x = smoothed_sanitized_csi;
L = size(x);
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
    % If OUTPUT_SUPPRESSED does not exist then initialize all the globals.
    %if exist('OUTPUT_SUPPRESSED') == 0
     %   globals_init()
    %end
    %% DEBUG AND OUTPUT VARIABLES-----------------------------------------------------------------%%
    % Debug Variables
    global DEBUG_PATHS
    global DEBUG_PATHS_LIGHT
    
    % Output Variables
    global OUTPUT_AOAS
    global OUTPUT_TOFS
    global OUTPUT_AOA_MUSIC_PEAK_GRAPH
    global OUTPUT_TOF_MUSIC_PEAK_GRAPH
    global OUTPUT_AOA_TOF_MUSIC_PEAK_GRAPH
    global OUTPUT_SELECTIVE_AOA_TOF_MUSIC_PEAK_GRAPH
    global OUTPUT_BINARY_AOA_TOF_MUSIC_PEAK_GRAPH
    global OUTPUT_SUPPRESSED
    global OUTPUT_FIGURES_SUPPRESSED
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % if nargin == 4
        data_name = '-';
    %end
    
    % Data covarivance matrix
    R = x * x'/L(2); 
    % Find the eigenvalues and eigenvectors of the covariance matrix
  %  
    
[Utmp,D] = eig(R);
D = abs(D);
[Dtmp,I] = sort(diag(D), 'descend');
D = diag(Dtmp);
U = Utmp(:,I);
Qn = U(:,SignalEndIdx+1:end);
Qs = U(:,1:SignalEndIdx);    
eigenvectors = Qn;       
    
    %% TODO: Tuning theta too??

% %     
%     theta = 0:3:360; 
%     phi = 0:3:360;
    
%     theta = 0:1:360; 
%     phi = 0:1:360;
% 
%  theta = 0:6:360; 
%     phi = 0:6:360;
    % time in milliseconds
%     %% TODO: Tuning tau....
%     %tau = 0:(1.0 * 10^-9):(50 * 10^-9);
%     tau = 0:(1 * 10^-9):(10 * 10^-9);
%     Pmusic = zeros(length(theta),length(phi), length(tau));
%     % Angle of Arrival Loop (AoA)
%     for ii = 1:length(theta)
%         for kk = 1:length(phi)
%         % Time of Flight Loop (ToF)
%             for jj = 1:length(tau)
%                 steering_vector = compute_steering_vector(theta(ii), phi(kk), tau(jj), ...
%                         frequency, sub_freq_delta, antenna_distance);
%                 PP = steering_vector' * (eigenvectors * eigenvectors') * steering_vector;
%                 Pmusic(ii,kk, jj) = abs(1 /  PP);
%             end
%         end
%     end
%     
%     PPmusic = zeros(length(theta),length(phi));
%     PPmusic = sum(Pmusic,3);
   

%  %PPmusic = Pmusic(:,:, 30);
%    % if OUTPUT_AOA_TOF_MUSIC_PEAK_GRAPH && ~OUTPUT_SUPPRESSED && ~OUTPUT_FIGURES_SUPPRESSED
%         % Theta (AoA) & Tau (ToF) 3D Plot
%         figure('Name', 'Azimuth & Elevation MUSIC Peaks', 'NumberTitle', 'off')
%         mesh(theta,phi, PPmusic)
%         xlabel('elevation')
%         ylabel('azimuth')
%         zlabel('Spectrum Peaks')
%         title('Azimuth and Elevation Estimation from Modified MUSIC Algorithm')
%         grid on
%    % end

%    newtheta = 0:3:180; 
%     newphi = 0:3:180;
%    ppp = PPmusic(61:121,61:121);
%    figure('Name', 'Azimuth & Elevation MUSIC Peaks', 'NumberTitle', 'off')
%         mesh(newphi,newtheta, ppp)
%         xlabel('azimuth')
%         ylabel('elevation')
%         zlabel('Spectrum Peaks')
%         title('Azimuth and Elevation Estimation from Modified MUSIC Algorithm')
%         grid on
   
%    
   
   
   
   %% TODO: Tuning theta too??
    theta = -90:1:90; 
    % time in milliseconds
    %% TODO: Tuning tau....
    tau = -(5 * 10^-8):(0.1 * 10^-8):(5 * 10^-8);
    %tau = 0:(100.0 * 10^-9):(3000 * 10^-9);
    Pmusic = zeros(length(theta), length(tau));
    % Angle of Arrival Loop (AoA)
    for ii = 1:length(theta)
        % Time of Flight Loop (ToF)
        for jj = 1:length(tau)
            steering_vector = compute_steering_vector11(theta(ii), tau(jj), ...
                    frequency, sub_freq_delta, antenna_distance);
            PP = steering_vector' * (eigenvectors * eigenvectors') * steering_vector;
            Pmusic(ii, jj) = abs(1 /  PP);
        end
    end
% 
%     % Convert to decibels
%     % ToF loop
%     for jj = 1:size(Pmusic, 2)
%         % AoA loop
%         for ii = 1:size(Pmusic, 1)
%             Pmusic(ii, jj) = 10 * log10(Pmusic(ii, jj));
%         end
%     end
% 
%     
        % Theta (AoA) & Tau (ToF) 3D Plot
        figure('Name', 'AoA & ToF MUSIC Peaks', 'NumberTitle', 'off')
        surf(tau, theta, Pmusic)
        xlabel('Time of Flight')
        ylabel('Angle of Arrival in degrees')
        zlabel('Spectrum Peaks')
        title('AoA and ToF Estimation from Modified MUSIC Algorithm')
        grid on
        shading interp
        view([0,90])

   
   
   
   
   
   
%     if (DEBUG_PATHS || OUTPUT_AOA_MUSIC_PEAK_GRAPH) ...
%             && ~OUTPUT_SUPPRESSED && ~OUTPUT_FIGURES_SUPPRESSED
%         % Theta (AoA)
%         figure_name_string = sprintf('%s: Number of Paths: %d', data_name, num_computed_paths);
%         figure('Name', figure_name_string, 'NumberTitle', 'off')
%         plot(theta, Pmusic(:, 1), '-k')
%         xlabel('Angle, \theta')
%         ylabel('Spectrum function P(\theta, \tau)  / dB')
%         title('AoA Estimation as a function of theta')
%         grid on
%     end
% 
%     binary_peaks_pmusic = imregionalmax(Pmusic);
%     
%     % Get AoAs that have peaks
%     % fprintf('Future estimated aoas\n')
%     aoa_indices = linspace(1, size(binary_peaks_pmusic, 1), size(binary_peaks_pmusic, 1));
%     aoa_peaks_binary_vector = any(binary_peaks_pmusic, 2);
%     estimated_aoas = theta(aoa_peaks_binary_vector);
%     
%     if OUTPUT_AOAS && ~OUTPUT_SUPPRESSED
%         fprintf('Estimated AoAs\n')
%         estimated_aoas
%     end
% 
%     aoa_peak_indices = aoa_indices(aoa_peaks_binary_vector);
%     
%     % Get ToFs that have peaks
%     time_peak_indices = zeros(length(aoa_peak_indices), length(tau));
%     % AoA loop (only looping over peaks in AoA found above)
%     for ii = 1:length(aoa_peak_indices)
%         aoa_index = aoa_peak_indices(ii);
%         binary_tof_peaks_vector = binary_peaks_pmusic(aoa_index, :);
%         matching_tofs = tau(binary_tof_peaks_vector);
%         
%         % Pad ToF rows with -1s to have non-jagged matrix
%         negative_ones_for_padding = -1 * ones(1, length(tau) - length(matching_tofs));
%         time_peak_indices(ii, :) = horzcat(matching_tofs, negative_ones_for_padding);
%     end
% 
%     
%     if OUTPUT_BINARY_AOA_TOF_MUSIC_PEAK_GRAPH && ~OUTPUT_SUPPRESSED && ~OUTPUT_FIGURES_SUPPRESSED
%         figure('Name', 'BINARY Peaks over AoA & ToF MUSIC Spectrum', 'NumberTitle', 'off')
%         mesh(tau, theta, double(binary_peaks_pmusic))
%         xlabel('Time of Flight')
%         ylabel('Angle of Arrival in degrees')
%         zlabel('Spectrum Peaks')
%         title('AoA and ToF Estimation from Modified MUSIC Algorithm')
%         grid on
%     end
% 
%     if OUTPUT_SELECTIVE_AOA_TOF_MUSIC_PEAK_GRAPH && ~OUTPUT_SUPPRESSED && ~OUTPUT_FIGURES_SUPPRESSED
%         % Theta (AoA) & Tau (ToF) 3D Plot
%         figure('Name', 'Selective AoA & ToF MUSIC Peaks, with only peaked AoAs', 'NumberTitle', 'off')
%         mesh(tau, estimated_aoas, Pmusic(aoa_peak_indices, :))
%         xlabel('Time of Flight')
%         ylabel('Angle of Arrival in degrees')
%         zlabel('Spectrum Peaks')
%         title('AoA and ToF Estimation from Modified MUSIC Algorithm')
%         grid on
%     end
%     
%     if OUTPUT_TOF_MUSIC_PEAK_GRAPH && ~OUTPUT_SUPPRESSED && ~OUTPUT_FIGURES_SUPPRESSED
%         % Tau (ToF)
%         for ii = 1:1%length(estimated_aoas)
%             figure_name_string = sprintf('ToF Estimation as a Function of Tau w/ AoA: %f', ...
%                     estimated_aoas(ii));
%             figure('Name', figure_name_string, 'NumberTitle', 'off')
%             plot(tau, Pmusic(ii, :), '-k')
%             xlabel('Time of Flight \tau / degree')
%             ylabel('Spectrum function P(\theta, \tau)  / dB')
%             title(figure_name_string)
%             grid on
%         end
%     end
%     
%     % Set return values
%     % AoA is now a column vector
%     estimated_aoas = transpose(estimated_aoas);
%     % ToF is now a length(estimated_aoas) x length(tau) matrix, with -1 padding for unused cells
%     estimated_tofs = time_peak_indices;
%end

% %% Computes the steering vector for SpotFi. 
% % Each steering vector covers 2 antennas on 15 subcarriers each.
% % theta           -- the angle of arrival (AoA) in degrees
% % tau             -- the time of flight (ToF)
% % freq            -- the central frequency of the signal
% % sub_freq_delta  -- the frequency difference between subcarrier
% % ant_dist        -- the distance between each antenna
% % Return:
% % steering_vector -- the steering vector evaluated at theta and tau
% %
% % NOTE: All distance measurements are in meters
function steering_vector = compute_steering_vector11(theta, tau, freq, sub_freq_delta, ant_dist)
% nn =11时,smoothed_csi的shape: 40*22
    steering_vector = zeros(10, 1);
    k = 1;
    base_element = 1;
    for ii = 1:2
        for jj = 1:20 %steering vector, adjustable
            steering_vector(k, 1) = base_element * omega_tof_phase11(tau, sub_freq_delta)^(jj - 1);
            k = k + 1;
        end
        base_element = base_element * phi_aoa_phase11(theta, freq, ant_dist);
    end
end

% %% Compute the phase shifts across subcarriers as a function of ToF
% % tau             -- the time of flight (ToF)
% % frequency_delta -- the frequency difference between adjacent subcarriers
% % Return:
% % time_phase      -- complex exponential representing the phase shift from time of flight
function time_phase = omega_tof_phase11(tau, sub_freq_delta)
    time_phase = exp(-1i * 2 * pi * sub_freq_delta * tau);
end

% %% Compute the phase shifts across the antennas as a function of AoA
% % theta       -- the angle of arrival (AoA) in degrees
% % frequency   -- the frequency of the signal being used
% % d           -- the spacing between antenna elements
% % Return:
% % angle_phase -- complex exponential representing the phase shift from angle of arrival
function angle_phase = phi_aoa_phase11(theta, frequency, d)
    % Speed of light (in m/s)
    c = 3.0 * 10^8;
    % Convert to radians
    theta = theta / 180 * pi;
    angle_phase = exp(-1i * 2 * pi * d * sin(theta) * (frequency / c));
end

% %% Creates the smoothed CSI matrix by rearranging the various csi values in the default CSI matrix.
% % csi          -- the regular CSI matrix to use for creating the smoothed CSI matrix
% % Return:
% % smoothed_csi -- smoothed CSI matrix following the construction put forth in the SpotFi paper.
% %                   Each column in the matrix includes data from 2 antennas and 15 subcarriers each.
% %                   Has dimension 30x32. 
function smoothed_csi = smooth_csi11(csi) %nn=11时，smoothed_csi的形状为 40 *22

nn = 11; %related to number of sensors, adjustable
nnn = (30-nn+1)*2;

    smoothed_csi = zeros(nnn, 2*nn);
    
    for n=1:nn
        
    smoothed_csi(:,n)=reshape(csi(1:2,n:n+(30-nn)).',nnn,1);
    smoothed_csi(:,n+nn)=reshape(csi(2:3,n:n+(30-nn)).',nnn,1);
    
    end


%     smoothed_csi = zeros(size(csi, 2), size(csi, 2));
%     % Antenna 1 (values go in the upper left quadrant)
%     m = 1;
%     for ii = 1:1:15
%         n = 1;
%         for j = ii:1:(ii + 15)
%             smoothed_csi(m, n) = csi(1, j); % 1 + sqrt(-1) * j;
%             n = n + 1;
%         end
%         m = m + 1;
%     end
%     
%     % Antenna 2
%     % Antenna 2 has its values in the top right and bottom left
%     % quadrants, the first for loop handles the bottom left, the second for
%     % loop handles the top right
%     
%     % Bottom left of smoothed csi matrix
%     for ii = 1:1:15
%         n = 1;
%         for j = ii:1:(ii + 15)
%             smoothed_csi(m, n) = csi(2, j); % 2 + sqrt(-1) * j;
%             n = n + 1;
%         end
%         m = m + 1;
%     end
%     
%     % Top right of smoothed csi matrix
%     m = 1;
%     for ii = 1:1:15
%         n = 17;
%         for j = ii:1:(ii + 15)
%             smoothed_csi(m, n) = csi(2, j); %2 + sqrt(-1) * j;
%             n = n + 1;
%         end
%         m = m + 1;
%     end
%     
%     % Antenna 3 (values go in the lower right quadrant)
%     for ii = 1:1:15
%         n = 17;
%         for j = ii:1:(ii + 15)
%             smoothed_csi(m, n) = csi(3, j); %3 + sqrt(-1) * j;
%             n = n + 1;
%         end
%         m = m + 1;
%     end
end
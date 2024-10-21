
tic
SignalEndIdx = 150; %%%number of signals, adjustable
antenna_distance = 2.8e-2;


    frequency = 5.32 * 10^9;
    sub_freq_delta = (40 * 10^6) / 30;

x=squeeze(CSI(200,:,:));
% AOA=zeros(600,60,60,21);
% for index=1:600
% x = squeeze(CSI(index,:,:));
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
 
    
    % Data covarivance matrix
    R = x * x'/L(2); 
% 计算协方差矩阵的特征值和特征向量。    
[Utmp,D] = eig(R);
% 对特征值进行排序，并划分出信号子空间和噪声子空间。
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
    tau = (0 * 10^-9):(1 * 10^-9):(20 * 10^-9);%ToF adjustable
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

    PPmusic = zeros(length(theta),length(phi));
    PPstatic = sum(Pmusic,3);

     I=figure('Name', 'Azimuth & Elevation MUSIC Peaks', 'NumberTitle', 'off','visible','on');
        surf(phi,theta,PPstatic)
        ylabel('elevation')
        xlabel('azimuth')
        zlabel('Spectrum Peaks')
        title('Azimuth and Elevation Estimation from Modified MUSIC Algorithm')
        shading interp
        view([0,90])    
%         saveas(I,['C:\Users\yicha\Desktop\figure6\',int2str(index),'.jpg']);
% is=index;
%         AOA(is,:,:,:)=Pmusic;
% end
% 
%    
% toc
% disp(['Time: ',num2str(toc)]);

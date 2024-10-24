% [tau, theta, Pmusic] = estimate_aoa_tof('my_data/rb_01_01_01.dat', 'C:\Users\LibaoXing\Desktop\rb_01_01_01.mat');
process_dat_files('E:\WorkSpace\datasets\person_in_wifi_3d\CSI', 'E:\WorkSpace\datasets\person_in_wifi_3d\AOA-TOF', 'E:\WorkSpace\datasets\person_in_wifi_3d\AOA-TOF-PNG')

function process_dat_files(sourceFolder, targetFolder, targetPngFolder)
    % 创建或打开 log.txt 文件用于记录处理信息
    logFilePath = fullfile(targetFolder, 'log.txt');
    logFileID = fopen(logFilePath, 'a'); % 以追加模式打开文件

    % 找到源文件夹下所有 .dat 文件
    files = dir(fullfile(sourceFolder, '**', '*.mat'));
    
    for i = 1:length(files)
        % 获取 .dat 文件的完整路径
        sourceFilePath = fullfile(files(i).folder, files(i).name);
        
        % 计算目标文件的完整路径
        relativePath = strrep(files(i).folder, sourceFolder, '');
        targetFilePath = fullfile(targetFolder, relativePath, ...
                                   [files(i).name(1:end-4), '.mat']);

        targetPngPath = fullfile(targetPngFolder, relativePath, ...
                                   [files(i).name(1:end-4), '.png']);
        
        % 确保目标目录存在
        if ~exist(fileparts(targetFilePath), 'dir')
            mkdir(fileparts(targetFilePath));
        end

        if ~exist(fileparts(targetPngPath), 'dir')
            mkdir(fileparts(targetPngPath));
        end

        try
            % 调用 estimate_aoa_tof 函数处理文件并保存结果
            estimate_aoa_tof(sourceFilePath, targetFilePath, targetPngPath, false, true);
            % 输出结果并记录日志
            fprintf('结果已保存到: %s\n', targetFilePath);
            fprintf(logFileID, '成功处理文件: %s，结果已保存到: %s\n', sourceFilePath, targetFilePath);
        catch ME
            % 捕获错误并输出警告信息
            warning('处理文件 %s 时发生错误: %s', sourceFilePath, ME.message);
            fprintf(logFileID, '处理文件 %s 时发生错误: %s\n', sourceFilePath, ME.message); % 记录失败信息
            continue; % 继续处理下一个文件
        end
    end

    % 关闭 log.txt 文件
    fclose(logFileID);
end


function [tau, theta, Pmusic] = estimate_aoa_tof(source_path, save_path, png_save_path, visualize, save_png)
    % 参数设置
    SignalEndIdx = 25;
    sub_freq_delta = (40 * 10^6) / 114;
    frequency = 5e9;
    M = 3;
    antenna_distance = 2.6e-2;
    SubCarrInd = [-57:-1,1:57];
    N = length(SubCarrInd);
    T = 1;
    
    % 读取CSI数据
    sample_csi_traceTmp = load(source_path); 
    csi_phases = sample_csi_traceTmp.CSIphase;
    csi_amps = sample_csi_traceTmp.CSIlamp;
    sample_csi_trace0 = csi_amps .* exp(1j * csi_phases);

    % sample_csi_trace0 =sample_csi_trace0(:,42:71,:); % 从114个子载波里选取中心的30个子载波 shape：3*30*10

    size0=size(sample_csi_trace0,3);
    antenna1_card1(size0,N) = 0;
    antenna2_card1(size0,N) = 0;
    antenna3_card1(size0,N) = 0;

    
    % 提取CSI
    for k = 1:size0
        csi = sample_csi_trace0(:,:,k);
            antenna1_card1(k, :) = csi(1, :);
            antenna2_card1(k, :) = csi(2, :);
            antenna3_card1(k, :) = csi(3, :);
    end
    
    % 重塑CSI
    CSI = zeros(1, N*3, 10);
    for i = 1:10:10
        tempcsi = zeros(N*3, 10);
        for pack = i:(i + 9)
            sample_csi_trace = [antenna1_card1(pack, :)'; antenna2_card1(pack, :)'; antenna3_card1(pack, :)'];
            csi_plot = reshape(sample_csi_trace, N, M);
            [PhsSlope, PhsCons] = removePhsSlope(csi_plot, M, SubCarrInd, N);
            ToMult = exp(1i * (-PhsSlope * repmat(SubCarrInd(:), 1, M) - PhsCons * ones(N, M)));
            csi_plot = csi_plot .* ToMult;
            relChannel_noSlope = reshape(csi_plot, N, M, T);
            sanitized_csi0 = relChannel_noSlope(:);
            tempcsi(:, (pack - i + 1)) = sanitized_csi0;
        end
        CSI(round(i / 10) + 1, :, :) = tempcsi;
    end
    
    % 计算协方差矩阵
    ind = 1;
    x = squeeze(CSI(ind, :, :));
    L = size(x);
    R = x * x' / L(2);
    [Utmp, D] = eig(R);
    D = abs(D);
    [Dtmp, I] = sort(diag(D), 'descend');
    D = diag(Dtmp);
    U = Utmp(:, I);
    Qn = U(:, SignalEndIdx + 1:end);
    Qs = U(:, 1:SignalEndIdx);
    eigenvectors = Qn;

    % 计算MUSIC谱
    theta = -90:1:90;
    tau = -(5 * 10^-8):(0.05 * 10^-8):(5 * 10^-8);
    Pmusic = zeros(length(theta), length(tau));
    
    for ii = 1:length(theta)
        for jj = 1:length(tau)
            steering_vector = compute_steering_vector342(theta(ii), tau(jj), frequency, sub_freq_delta, antenna_distance);
            PP = steering_vector' * (eigenvectors * eigenvectors') * steering_vector;
            Pmusic(ii, jj) = abs(1 / PP);
        end
    end
    
    % 保存数据到.mat文件
    save(save_path, 'tau', 'theta', 'Pmusic');
    
    % 可视化（根据用户选择）
    if visualize
        figure('Name', 'AoA & ToF MUSIC Peaks', 'NumberTitle', 'off');
        surf(tau, theta, Pmusic);
        xlabel('Time of Flight');
        ylabel('Angle of Arrival in degrees');
        zlabel('Spectrum Peaks');
        title('AoA and ToF Estimation from Modified MUSIC Algorithm');
        shading interp;
    end
    
    if save_png
        I = figure('Name', 'AoA & ToF MUSIC Peaks', 'NumberTitle', 'off', 'Visible', 'off');
        surf(tau, theta, Pmusic);
        xlabel('Time of Flight');
        ylabel('Angle of Arrival in degrees');
        zlabel('Spectrum Peaks');
        title('AoA and ToF Estimation from Modified MUSIC Algorithm');
        shading interp;
        view([0,90]);
        saveas(I, png_save_path);
        close(I); % 关闭图形对象，释放资源
    end
end


function steering_vector = compute_steering_vector342(theta, tau, freq, sub_freq_delta, ant_dist)
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

% function steering_vector = compute_steering_vector90(theta, tau, freq, sub_freq_delta, ant_dist)
%     steering_vector = zeros(342, 1);
%     k = 1;
%     base_element = 1;
%     for ii = 1:3
%         for jj = 1:114
%             steering_vector(k, 1) = base_element * omega_tof_phase11(tau, sub_freq_delta)^(jj - 1);
%             k = k + 1;
%         end
%         base_element = base_element * phi_aoa_phase11(theta, freq, ant_dist);
%     end
% end

function time_phase = omega_tof_phase11(tau, sub_freq_delta)
    time_phase = exp(-1i * 2 * pi * sub_freq_delta * tau);
end

function angle_phase = phi_aoa_phase11(theta, frequency, d)
    c = 3.0 * 10^8;
    theta = theta / 180 * pi;
    angle_phase = exp(-1i * 2 * pi * d * sin(theta) * (frequency / c));
end
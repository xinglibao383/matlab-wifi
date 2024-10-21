%% Computes the steering vector for SpotFi. 
% Each steering vector covers 2 antennas on 15 subcarriers each.
% theta           -- the angle of arrival (AoA) in degrees
% tau             -- the time of flight (ToF)
% alpha -- AoD
% freq            -- the central frequency of the signal
% sub_freq_delta  -- the frequency difference between subcarrier
% ant_dist        -- the distance between each antenna
% Return:
% steering_vector -- the steering vector evaluated at theta and tau
%
% NOTE: All distance measurements are in meters
% ant_dist = 0.1;
% freq = 5.32 * 10^9;
% sub_freq_delta = (40 * 10^6) / 30;
% theta = 45;
% phi = 30;
% alpha = 60;
% tau = 100.0 * 10^-9;


function steering_vector = compute_steering_vector_AoD(theta,phi,tau, alpha, freq, sub_freq_delta, ant_dist)
    steering_vector = zeros(100, 1);
    k = 1;
    base_element = 1;
    base_element_1 = 1;
    for ii = 1:5
        for jj = 1:10
            steering_vector(k, 1) = base_element * omega_tof_phase(tau, sub_freq_delta)^(jj - 1);
            k = k + 1;
        end
        %base_element = base_element * phi_aoa_phase(theta, freq, ant_dist);
        base_element = base_element * phi_aoa_phase(theta, phi, freq, ant_dist);
    end
    
    for ff = 1:5
        for gg = 1:10
            steering_vector(k, 1) = base_element_1 * omega_tof_phase(tau, sub_freq_delta)^(gg - 1);
            k = k + 1;
        end
        base_element_1 = base_element_1 * phi_aoa_phase_2(theta, freq, ant_dist);
    end
    
%      base_element_2 = phi_aoa_phase_2(theta, phi, freq, ant_dist);
%     for ff = 1:15
%             steering_vector(k, 1) = base_element_2 * omega_tof_phase(tau, sub_freq_delta)^(ff - 1);
%             k = k + 1;
%     end
    
steering_vector = [steering_vector; steering_vector*phi_aoa_phase_3(alpha,freq, ant_dist); steering_vector*phi_aoa_phase_3(alpha,freq, ant_dist)*phi_aoa_phase_3(alpha,freq, ant_dist)];


end
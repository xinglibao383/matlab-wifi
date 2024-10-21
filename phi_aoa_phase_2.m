%% Compute the phase shifts across the antennas as a function of AoA
% theta       -- the elevation in degrees
% phi -- the azimuth in degrees
% frequency   -- the frequency of the signal being used
% d           -- the spacing between antenna elements
% Return:
% angle_phase -- complex exponential representing the phase shift from angle of arrival
function angle_phase_2 = phi_aoa_phase_2(theta,frequency, d)
    % Speed of light (in m/s)
    c = 3.0 * 10^8;
    % Convert to radians
    theta = theta / 180 * pi;
   % phi = phi / 180 * pi;
    angle_phase_2 = exp(-1i * 2 * pi * d * cos(theta) * (frequency / c));
end

% function angle_phase_2 = phi_aoa_phase_2(theta,phi, frequency, d)
%     % Speed of light (in m/s)
%     c = 3.0 * 10^8;
%     % Convert to radians
%     theta = theta / 180 * pi;
%     phi = phi / 180 * pi;
%     angle_phase_2 = exp(-1i * 2 * pi * d * cos(theta)*sin(phi) * (frequency / c));
% end
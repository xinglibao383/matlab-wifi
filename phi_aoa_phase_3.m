%% Compute the phase shifts across the antennas as a function of AoA
% theta       -- the elevation in degrees
% phi -- the azimuth in degrees
% alpha -- AoD
% frequency   -- the frequency of the signal being used
% d           -- the spacing between antenna elements
% Return:
% angle_phase -- complex exponential representing the phase shift from angle of arrival
function angle_phase_3 = phi_aoa_phase_3(alpha, frequency, d)
    % Speed of light (in m/s)
    c = 3.0 * 10^8;
    % Convert to radians
    alpha = alpha / 180 * pi;
    angle_phase_3 = exp(-1i * 2 * pi * d * sin(alpha) * (frequency / c));
end

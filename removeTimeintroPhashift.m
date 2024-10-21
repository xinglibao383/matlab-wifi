function [csi_matrix2, phaseshiftA, phaseshiftB] = removeTimeintroPhashift(csi_matrix1,csi_matrix2)
% 返回去除时间引入相移后的CSI矩阵 csi_matrix2 和相移量 phaseshiftA、phaseshiftB。
% csi_matrix=sample_csi_trace_sanitized(1:30,:);
% csi_matrix2=sample_csi_trace_sanitized2(1:30,:);
SubCarrInd = [-58,-54,-50,-46,-42,-38,-34,-30,-26,-22,-18,-14,-10,-6,-2,2,6,10,14,18,22,26,30,34,38,42,46,50,54,58]; % WiFi subcarrier indices at which CSI is available
% N = length(SubCarrInd); % number of subcarriers
% M = 3;
    
    R = abs(csi_matrix2);
    phase_matrix1 = unwrap(angle(csi_matrix1));
    phase_matrix2 = unwrap(angle(csi_matrix2));
    packet_one_phase_matrix=unwrap(angle(csi_matrix2));
    phasediff = unwrap(phase_matrix1-phase_matrix2);
    
    %plot(phasediff)
fit_X = SubCarrInd(:);
%fit_X = 1:1:30;
fit_Y = phasediff;
result = polyfit(fit_X, fit_Y, 1);

% x1 = 1:1:30;
% f1 = polyval(result,x1);
% plot(x1,f1)
% hold on
% plot(phasediff)
% hold on
% plot(fit_X*result(1)-result(2))


for n=1:30
   phase_matrix2(n) = packet_one_phase_matrix(n) + SubCarrInd(n)*result(1)+result(2); 
end

csi_matrix2 = R .* exp(1i * phase_matrix2);
phaseshiftA = result(1);
phaseshiftB = result(2);

% plot(phase_matrix1)
% hold on
% plot(phase_matrix2)

end
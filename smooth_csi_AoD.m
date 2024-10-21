%% Creates the smoothed CSI matrix by rearranging the various csi values in the default CSI matrix.
% csi          -- the regular CSI matrix to use for creating the smoothed CSI matrix
% Return:
% smoothed_csi -- smoothed CSI matrix following the construction put forth in the SpotFi paper.
%                   Each column in the matrix includes data from 2 antennas and 15 subcarriers each.
%                   Has dimension 30x32. 


function smoothed_csi = smooth_csi_AoD(csi1, csi2, csi3, abc)

if abc == 1
 A = 1; B = 2; C = 3;
elseif abc == 2
 A = 1; B = 3; C = 2;
elseif abc == 3
 A = 2; B = 1; C = 3;
elseif abc == 4
 A = 2; B = 3; C = 1;
elseif abc == 5
 A = 3; B = 1; C = 2;
else
 A = 3; B = 2; C = 1;
end
    %smoothed_csi = zeros(180,16);
    %smoothed_csi = zeros(size(csi, 2), size(csi, 2));
    
    for ii = 1:16
       % ii = 1
        csitmp1 = [];
        csitmp1 = [csi1(A,(1+ii-1):(15+ii-1)) csi1(B,(1+ii-1):(15+ii-1)) csi1(A,(1+ii-1):(15+ii-1)) csi1(C,(1+ii-1):(15+ii-1))];
        csitmp1 = csitmp1.';
        csitmp11(:,ii) = csitmp1;
        
        
        csitmp2 = [];
        csitmp2 = [csi2(A,(1+ii-1):(15+ii-1)) csi2(B,(1+ii-1):(15+ii-1)) csi2(A,(1+ii-1):(15+ii-1)) csi2(C,(1+ii-1):(15+ii-1))];
        csitmp2 = csitmp2.';
        csitmp22(:,ii) = csitmp2;
        
        csitmp3 = [];
        csitmp3 = [csi3(A,(1+ii-1):(15+ii-1)) csi3(B,(1+ii-1):(15+ii-1)) csi3(A,(1+ii-1):(15+ii-1)) csi3(C,(1+ii-1):(15+ii-1))];
        csitmp3 = csitmp3.';
        csitmp33(:,ii) = csitmp3;
        
%         smoothed_csi(:,ii) = [csitmp1; csitmp2; csitmp3]; 
%         smoothed_csi(:,ii) = [csitmp1; csitmp3; csitmp2];
%         smoothed_csi(:,ii) = [csitmp2; csitmp1; csitmp3];
%         smoothed_csi(:,ii) = [csitmp2; csitmp3; csitmp1];
%         smoothed_csi(:,ii) = [csitmp3; csitmp2; csitmp1];
         %smoothed_csi(:,ii) = [csitmp1; csitmp2; csitmp3];
    end
    smoothed_csi = [csitmp11 csitmp22 csitmp33];
end





% function smoothed_csi = smooth_csi(csi)
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
% end
%% Particle viewer plots
% This code will plot the plasma wave generation and both AC & DC velocity 
% oscillations of the particles
% Copyright (c) 2026 Md Faiyaz Bin Hassan and Dr. Shubhendu Bhardwaj
% 
% This file is part of the ADI-FDTD Particle-In-Cell modeling software.
% It is released under the MIT License. See the LICENSE file for details.


figure;
plot(XX(gate_y_start+2, POS1_2deg_x(1):POS1_2deg_x(2)-1)-XX(gate_y_start+2, POS1_2deg_x(1)),...
    diag_figure(gate_y_start+2, POS1_2deg_x(1):POS1_2deg_x(2)-1)*(40)-5000, 'k')
hold on;
scatter(part_x, part_v_ac, '.r');
title('AC velocity of particles');

figure;
plot(XX(gate_y_start+2, POS1_2deg_x(1):POS1_2deg_x(2)-1)-XX(gate_y_start+2, POS1_2deg_x(1)),...
    diag_figure(gate_y_start+2, POS1_2deg_x(1):POS1_2deg_x(2)-1)*(40)+1e5, 'k')
hold on;
scatter(part_x, part_v, '.r');
title('velocity of particles');

figure;
subplot(2,1,1); plot(PIC_per_cell*MF/dx_2DEG);
subplot(2,1,2); plot(J);

barr_y_start =barr_y_start-10
barr_y_end=barr_y_start+30;
% figure;
% pcolor(XX(barr_y_start:barr_y_end,POS1_2deg_x(1):POS1_2deg_x(2)), ...
%           YY(barr_y_start:barr_y_end,POS1_2deg_x(1):POS1_2deg_x(2)),...
%           log10(abs(Ex(barr_y_start:barr_y_end,POS1_2deg_x(1):POS1_2deg_x(2)))));
% 
%       caxis([-2 2]*1e4); colorbar

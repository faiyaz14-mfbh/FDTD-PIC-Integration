%% =========================================================================
% ADI-FDTD-PIC: Coupled Electromagnetic / Monte-Carlo Particle Solver
% for 2DEG Channel Current in GaN/AlGaN HEMT Structures
%
% Author: Md Faiyaz Bin Hassan and Dr. Shubhendu Bhardwaj
% Advanced Wireless and EM Research Lab, University of Nebraska-Lincoln
% Advisor: Prof. Shubhendu Bhardwaj
%
% If you use or adapt any version of this code, please cite:
% M. F. B. Hassan and S. Bhardwaj, "Efficiency and Error Optimizations in
% ADI-FDTD Particle-In-Cell Model for Wave-Particle Interactions and 
% Slow-wave Effects in THz Gap," in Proc. IEEE International Symposium on 
% Antennas and Propagation and USNC-URSI Radio Science Meeting (AP-S/URSI),
% Detroit, MI, USA, 2026.
%
%% -------------------------------------------------------------------------
% DESCRIPTION
% This script couples an Alternating-Direction-Implicit FDTD (ADI-FDTD)
% electromagnetic solver with a Particle-In-Cell (PIC)
% particle mover to model plasma generation and wave-electron interaction
% in a HEMT two-dimensional electron gas (2DEG) channel. 
% 
% It checks that the coupled EM-MC model reproduces the expected current 
% density in the channel.
%
% Notes on this version:
%   - The 2DEG is presented with a interdigitated gate. The gate is
%   designed for 300 GHz frequency. The period needs to be modified once
%   the frequency is changed.
%   - The macro-particle scaling factor (MF) definition was revised to
%     be based on particle behavior at the moment of injection.
%
% REQUIREMENTS
%   - MATLAB (no additional toolboxes required for now)
%   - Requires precomputed Monte-Carlo current data file:
%     'MC_Current_1p0e5_dx_10nm_L_5um_delV_200_ttdynamic_CN20.mat'
%
% USAGE
%   Run directly in MATLAB. Key simulation parameters
% =========================================================================

clear;
f =300e9;
close all;

%% 2DEG part
timestamp = fix(clock);
uniq_identifier = '_';
for i=1:length(timestamp)
    uniq_identifier = [uniq_identifier num2str(timestamp(i))];
end

filename = ['Output_' uniq_identifier];


%% sim mode
% Flags controlling output/debug behavior for this run.
% recording_on: enable field/data recording during the run
% mesh_only: if 1, generate and check the mesh without running the solver

recording_on   = 0;
mesh_only      = 0; %keep it 0 for a full simulation run

periodic_file_dump_on    = 1;
period_file_dump         =  0.5e-12;


%% FDTD init
% constants 

epsi_o    =  8.854e-12;       % F/m
mu_o      =  4*pi*1e-7;
c         =  1/sqrt(mu_o*epsi_o);
lambda    =  (c/f);                 % m N_HD  1 THz
omega     =  2*pi*c/lambda;
eta       =  sqrt(mu_o/epsi_o);

epsi_sub  = 9.5;
epsi_barr = 9.5;


%%  non-uniform mesh
% Builds a graded, non-uniform grid in x and y. Fine resolution (nm-scale)
% is used across the 2DEG/gate/barrier layers where field gradients and
% particle transport matter most; coarser resolution is used in the bulk
% metal/substrate/air regions to keep the overall grid size manageable.

% y-directed mesh

len_y_metal=1e-6;  dy_metal=100e-9;           M_metal=round(len_y_metal/dy_metal);
len_y_air2=1e-6;   dy_air2=50e-9;             M_air2=round(len_y_air2/dy_air2);
len_y_gate=30e-9;  dy_gate=10e-9;             M_gate=round(len_y_gate/dy_gate);
len_y_barr=30e-9;  dy_barr=10e-9;             M_barr=round(len_y_barr/dy_barr);
len_y_2DEG=10e-9;  dy_2DEG=len_y_2DEG;        M_2DEG=round(len_y_2DEG/dy_2DEG);        
len_y_subs1=1e-6;  dy_subs1=10e-9;            M_subs1=round(len_y_subs1/dy_subs1);              
len_y_subs2=1e-6;  dy_subs2= 50e-9;           M_subs2=round(len_y_subs2/dy_subs2);
len_y_metal2=1e-6; dy_metal2=100e-9;          M_metal2=round(len_y_metal2/dy_metal2);


t_2DEG   = dy_2DEG;

M        =  M_metal+M_air2+M_gate+M_barr+M_2DEG+M_subs1+M_subs2+M_metal2;

% x-directed mesh
len_PML     = 80e-6;       
len_PML_x   = 80e-6;

len_x_1=len_PML_x;    dx_1=10e-6;      N_1=round(len_x_1/dx_1);
len_x_2=80e-6;         dx_2=5e-6;       N_2=round(len_x_2/dx_2);
len_x_3=1e-6;         dx_3=100e-9;     N_3=round(len_x_3/dx_3);
len_x_2DEG=5e-6;      dx_2DEG=10e-9;   N_2DEG=round(len_x_2DEG/dx_2DEG);
len_x_4=1e-6;         dx_4=100e-9;     N_4=round(len_x_4/dx_4);
len_x_5=80e-6;         dx_5=5e-6;       N_5=round(len_x_5/dx_5);
len_x_6=len_PML_x;    dx_6=10e-6;      N_6=round(len_x_6/dx_6);

N            =          N_1+N_2+N_3+N_2DEG+N_4+N_5+N_6;

%% mesh sizes along x direction
DX      = cat(2,ones(M,N_1)*dx_1,        ...
                ones(M,N_2)*dx_2,        ...
                ones(M,N_3)*dx_3,        ...
                ones(M,N_2DEG)*dx_2DEG,  ...
                ones(M,N_4)*dx_4,        ... 
                ones(M,N_5)*dx_5,        ...
                ones(M,N_6)*dx_6); 

DY      = cat(1,   ...
                ones(M_metal,N)*dy_metal,   ...
                ones(M_air2,N)*dy_air2,   ...
                ones(M_gate,N)*dy_gate,   ...
                ones(M_barr,N)*dy_barr,   ...
                ones(M_2DEG,N)*dy_2DEG,   ...
                ones(M_subs1,N)*dy_subs1, ...
                ones(M_subs2,N)*dy_subs2, ...
                ones(M_metal2,N)*dy_metal2   ...
                );

it_ADI = 1; % keep it 1 for ADI-FDTD
CN     = 20; % Courant Number
%% delta t for FDTD
dt      = CN/(c*sqrt(1/min(DX(1,1:N))^2 + 1/min(DY(1:M,1))^2));
%% plot grid for diganosis
    X           =  cumsum(DX(1,1:N));
    Y           =  cumsum(DY(1:M,1));
    [XX YY]     =  meshgrid(X,Y);
    diag_figure =  ones(M,N);

%% material properties
    mu       =  ones(M,N)*mu_o;
    epsi     =  ones(M,N)*epsi_o;

%% MC init
% Loads precomputed Monte-Carlo transport data (particle positions,
% velocities, drift velocity, sheet density) used to initialize the
% particle population for the coupled EM-MC run.
load('MC_Current_1p0e5_dx_10nm_L_5um_delV_200_ttdynamic_CN20.mat');

% constants
   qe = 1.602e-19;
   m_o = 9.1093e-31;
   me  = 0.2*m_o;
 
% 2DEG Mesh Edge points
    X_2DEG_edge                   = dx_2DEG*(0:N_2DEG);
    
% 2DEG  cell centers, very last one is extra
     X_2DEG_CellCenter              = dx_2DEG*(0:N_2DEG)+dx_2DEG/2;
     
% j and e at the edges 
     j_sh_density_macro_1by2  = zeros(1, N_2DEG+1);
     j_sh_density_macro       = zeros(1, N_2DEG+1);
     e                        = zeros(1, N_2DEG+1);
     J                        = zeros(1, N_2DEG+1);   

    part_v_dc  =  part_v;
    part_v_ac  =  zeros(size(part_v));
    part_E     =  zeros(part_max,1);
    
% Macro-particle scaling: each simulated "particle" represents many
% real electrons. PIC sets the initial particle count per cell from
% injection rate and time step; MF is the ratio of the physical sheet
% density to the simulated macro-particle sheet density, used to
% scale simulated current back to physical current.
    PIC                    = np_i*dx_2DEG/(v_drift*dt);
    macro_sheet_density    = PIC / (dx_2DEG);  % m^-2
    MF                     = n_sh_m2 /macro_sheet_density;

%% intialize field grids TE mode

Ex         = zeros(M,N);
Ex_k_half  = zeros(M,N); 
Ey    = zeros(M,N);
Ey_k_half  = zeros(M,N);
Ey_n_full  = zeros(M,N);
Hzy   = zeros(M,N);
Hzy_k_half  = zeros(M,N);
Hzx   = zeros(M,N);
Hzx_k_half  = zeros(M,N);
Hz    = zeros(M,N);
Hz_k_half    = zeros(M,N);

Ex_k         = zeros(M,N);
Ey_k    = zeros(M,N);
Hzy_k   = zeros(M,N);
Hzx_k   = zeros(M,N);
Hz_k    = zeros(M,N);

sigma_ey  = zeros(M,N);
sigma_ex  = zeros(M,N);

% converting to grid points
PML_x   =  round(len_PML_x/dx_1);
PML_y   =  round(len_PML/dy_metal);
PML_y_t =  round(len_y_subs2/dy_subs2 + len_y_barr/dy_barr + (len_PML-len_y_subs2-len_y_barr)/dy_metal);

%% PML conditions 
% electrical conductivity
% anisotropic conductivity is implemented for the same of PML boundary
% conditions on the border, but they are used even for the main region
% using sigma_ex = sigma_ey , sigma_mx = sigma_my 
% for wave propagting in   y- direction

m=2;
gradPMLx= ((1:PML_x)' / PML_x).^m * 100;
for u=N-PML_x+1:N
    sigma_ex(:,u-1)        = (sigma_ex(:,N-PML_x) +1) * (gradPMLx(u-(N-PML_x)));
    sigma_ex(:,N-u+1)      =  (sigma_ex(:,N-PML_x) +1) * (gradPMLx(u-(N-PML_x)));
    diag_figure(:,u-1)     = 8*ones(size(diag_figure(:,u)));
    diag_figure(:,N-u+1)   = 8*ones(size(diag_figure(:,N-u+1)));
end

% magnetic cond
sigma_mx   = zeros(M,N); % analogous to Hz
sigma_my   = zeros(M,N); % analogous to Hz

% PML conditions - magnetic conductivity
%sigma_my(M-PML_y:M,:)  =  sigma_ey(M-PML_y:M,:) .* mu(M-PML_y:M,:)./epsi(M-PML_y:M,:);
%sigma_my(1:PML_y,:)    =  sigma_ey(1:PML_y,:)   .* mu(1:PML_y,:)./epsi(1:PML_y,:);
sigma_mx(:,N-PML_x:N)  =  sigma_ex(:,N-PML_x:N) .* mu(:,N-PML_x:N)./epsi(:,N-PML_x:N);
sigma_mx(:,1:PML_x)    =  sigma_ex(:,1:PML_x)   .* mu(:,1:PML_x)./epsi(:,1:PML_x);

% Coeeficients for one D wave for SF/TF reference
% this should be calculated on bare skeleton with PML- without any features
% introduced in the grid, meaning conductivity 


%% Geometry

% 2DEG
    POS1_2deg_x   = [N_1+N_2+N_3   N_1+N_2+N_3+N_2DEG];
    POS1_2deg_y   = M_metal+M_air2+M_gate+M_barr+1;
    diag_figure(POS1_2deg_y,POS1_2deg_x(1):POS1_2deg_x(2)) = 60*ones(size(diag_figure(POS1_2deg_y,POS1_2deg_x(1):POS1_2deg_x(2))));
    
% waveguide (This was added based on parallel plate waveguide type feeding)
    sigma_ey(1:M_metal,:)         = 1*1e20;                   % pec TOP
    sigma_ex(1:M_metal,:)         = 1*1e20;                   % pec TOP 
    sigma_ey(M-M_metal2+1:M, :)   = 1*1e20;                   % pec BOTTOM 
    sigma_ex(M-M_metal2+1:M, :)   = 1*1e20;                   % pec BOTTOM 
    
    diag_figure(1:M_metal,:)      = 200*ones(size(diag_figure(1:M_metal,:)));
    diag_figure(M-M_metal2+1:M,:) = 100*ones(size(diag_figure(M-M_metal2+1:M, :)));
    
%  introducing periodc grating regions
    gate_y_start  = M_metal+M_air2+1;
    gate_y_end    = M_metal+M_air2+M_gate;

    grating_x_start= POS1_2deg_x(1)+6;
    num_of_grat  =  29;
    gate_gap     =  80e-9;
    gate_w       =  90e-9;
    grating_x_start_tmp = grating_x_start;

    for gr=1:num_of_grat
        sigma_ey(gate_y_start:gate_y_end, grating_x_start_tmp:grating_x_start_tmp+round(gate_w/dx_2DEG))   = 1*1e20;        % pec strip 
        sigma_ex(gate_y_start:gate_y_end, grating_x_start_tmp:grating_x_start_tmp+round(gate_w/dx_2DEG))   = 1*1e20;  
        diag_figure(gate_y_start:gate_y_end, grating_x_start_tmp:grating_x_start_tmp+round(gate_w/dx_2DEG)) = 200*ones(size(diag_figure(gate_y_start:gate_y_end,grating_x_start_tmp:grating_x_start_tmp+round(gate_w/dx_2DEG))));   % end
        grating_x_start_tmp = grating_x_start_tmp + (gate_gap+gate_w)/dx_2DEG;
    end
    
% dielectric substrate
    sub_y_start  = M_metal+M_air2+M_gate+M_barr+1;  sub_y_end   = sub_y_start+M_2DEG+M_subs1+M_subs2-1;
    epsi(sub_y_start:sub_y_end,POS1_2deg_x(1):POS1_2deg_x(2))        = epsi(sub_y_start:sub_y_end,POS1_2deg_x(1):POS1_2deg_x(2))*epsi_sub;
    diag_figure(sub_y_start:sub_y_end,POS1_2deg_x(1):POS1_2deg_x(2)) = 30*ones(size(diag_figure(sub_y_start:sub_y_end,POS1_2deg_x(1):POS1_2deg_x(2))));

% dielectric barrier
    barr_y_start = M_metal+M_air2+M_gate+1;         barr_y_end  = M_metal+M_air2+M_gate+M_barr;
    epsi(barr_y_start:barr_y_end,POS1_2deg_x(1):POS1_2deg_x(2))        = epsi(barr_y_start:barr_y_end,POS1_2deg_x(1):POS1_2deg_x(2))*epsi_sub;
    diag_figure(barr_y_start:barr_y_end,POS1_2deg_x(1):POS1_2deg_x(2)) = 50*ones(size(diag_figure(barr_y_start:barr_y_end,POS1_2deg_x(1):POS1_2deg_x(2))));

%% measure
% source

% introducing source in the periodc grating regions
    grating_x_start_tmp = grating_x_start;
    source_x = zeros(1, size(diag_figure,2));
    for gr=1:num_of_grat-1    
        if(mod(gr,2)==0)
          source_profile(grating_x_start_tmp-grating_x_start+round(gate_w/dx_2DEG) ...
                        :grating_x_start_tmp-grating_x_start+round(gate_w/dx_2DEG)+round(gate_gap/dx_2DEG))=-1;
                    
        else
         source_profile(grating_x_start_tmp-grating_x_start+round(gate_w/dx_2DEG) ...
                        :grating_x_start_tmp-grating_x_start+round(gate_w/dx_2DEG)+round(gate_gap/dx_2DEG))=1;
                      
        end
        grating_x_start_tmp = grating_x_start_tmp + (gate_gap+gate_w)/dx_2DEG;
    end

    source_x = grating_x_start+1:grating_x_start+size(source_profile,2); % N_1+3;
    source_y =  M_metal+M_air2+round(M_gate/2);
    diag_figure(source_y, source_x) = 180*source_profile;
    
    source_present_x=zeros(1,size(diag_figure,2));
     source_profile_x=zeros(1,size(diag_figure,2));
    source_present_x(source_x) =1;
    source_profile_x(source_x)=source_profile;
%% These are probing coordinates of the field at the left and right side of the devices (Don't need them for now)
    % meas_IP_x = POS1_2deg_x(1)-3;
    % meas_IP_y = gate_y_start:sub_y_end;
    % diag_figure(meas_IP_y, meas_IP_x) = 3000;
    % 
    % meas_OP_x = POS1_2deg_x(2)+3;
    % meas_OP_y = gate_y_start:sub_y_end;
    % diag_figure(meas_OP_y, meas_OP_x) = 190;
    % 
    % meas_OP2_x = N-N_6-4;
    % meas_OP2_y = M_metal+2:M-M_metal2-2;
    % diag_figure(meas_OP2_y, meas_OP2_x) = 200;
    % 
    % meas_barr_x = POS1_2deg_x(1):POS1_2deg_x(2);
    % meas_barr_y = M_metal+M_air2+M_gate+round(M_barr/2);
    % diag_figure(meas_barr_y, meas_barr_x) = 200;

%% coeeficients for update equations

c_exe    = (1-sigma_ey*dt./(4*epsi)) ./ (1+sigma_ey*dt./(4*epsi));
c_exh    = 1./(1+sigma_ey*dt./(4*epsi));
c_eye    = (1-sigma_ex*dt./(4*epsi)) ./ (1+sigma_ex*dt./(4*epsi));
c_eyh    = 1./(1+sigma_ex*dt./(4*epsi));

%calculate H field update coefficients
c_hzhy    = (1-sigma_my*dt./(4*mu)) ./ (1+sigma_my*dt./(4*mu));
c_hzey    = 1./(1+sigma_my*dt./(4*mu));
c_hzhx    = (1-sigma_mx*dt./(4*mu)) ./ (1+sigma_mx*dt./(4*mu));
c_hzex    = 1./(1+sigma_mx*dt./(4*mu));

c_itex   = 1 ./ (1+sigma_ey*dt./(4*epsi));
c_itey   = 1 ./ (1+sigma_ex*dt./(4*epsi));
c_ithx   = 1 ./ (1+sigma_mx*dt./(4*mu));
c_ithz   = 1 ./ (1+sigma_my*dt./(4*mu));

%% grid: illustration
%  Ey  Hz  Ey  Hz  Ey  Hz  Ey  Hz
%      Ex      Ex      Ex      Ex
%  Ey  Hz  Ey  Hz  Ey  Hz  Ey  Hz
%      Ex      Ex      Ex      Ex
%  Ey  Hz  Ey  Hz  Ey  Hz  Ey  Hz
%      Ex      Ex      Ex      Ex

% ^  this is y direction
% >  this is x direction

% difference between E and H calculations in 0.5 dleta_t, however each H
% calculation is delta_t apart- as per leapfrog algorithm
% column 1 in each matrix is along y-direction

%% in this code we are talking about the TEz polarized wave:
%  Hz, Ex, Ey components are being considered.

%% diag
figure;
h= pcolor(diag_figure);
set(h, 'EdgeColor', 'none');
% set(gca, 'XLim',[N_1 N-N_1], 'YLim',[M_1 M-M_1]);
colorbar;

%% diag
        figure;
        h= pcolor(XX, YY,diag_figure);
        set(h, 'EdgeColor', 'none');
%         set(gca, 'XLim',[len_x_1 X(end)-len_x_5], 'YLim',[len_y_1 Y(end)-len_y_6]);
        colorbar;

%%

%  sigma_ex(M_1:M_1+M_2-1, N_1:N_1+N_2-1) = 1e20;
%  sigma_ey(POS_2deg_y,POS_2deg_x(1):POS_2deg_x(2)-1) = 1e10/2;
%  sigma_ex(PEC_y_start:M-2,PML_x+10:N-PML_x-11)      = 1e10/5;
%  sigma_ex(POS_2deg_y,POS_2deg_x(1):POS_2deg_x(2)-1) = 1e10/2;
%  figure;    pcolor(XX,YY,sigma_ex);
  
T_toal = 20000;
k      = 10;


%% time iteration loop
t_total = 1e-10;
tt_max  = fix(t_total/(dt));
tic

tt_rec  = 10;
t_record_step = dt*tt_rec;

if (~mesh_only)
figure;
for tt =1:tt_max  % tt = 1:150 K
    
    tt    
 f=300e9;
   Ex(source_y, source_x)    =  1e4*sin(2*pi*f*dt*tt)*source_profile ;
    Hzy_1D(1)           = sin(2*pi*f*(tt)*dt)/377;
%% making matrix LU de-composition n+1/2 update
    delEy    = zeros(M,N);
    delHz    = zeros(M,N);

for k=1:it_ADI    
    if k~=1
        delEy = Ey_k - Ey;
        delHz = Hz_k - Hz;
    end
%% making matrix LU de-composition n+1/2 update
for nn=1:N-1

    % LU decomposition method for matrix inversion
  clear a b c beta alpha g;
  a(1:M-1)  = 1 + c_exh(1:M-1,nn) * dt^2 ./ (4* epsi(1:M-1,nn).* DY(1:M-1, nn)) ...
                            .*(c_hzey(2:M,nn) ./ mu(2:M,nn) ./ DY(2:M, nn) + c_hzey(1:M-1,nn) ./ mu(1:M-1,nn) ./DY(1:M-1,nn));

  b(1)      =   0;
  b(2:M-1)  =  -0.25*c_exh(2:M-1,nn).* c_hzey(2:M-1,nn) *dt^2 ./ ( epsi(2:M-1,nn) .* mu(2:M-1,nn)) ./ (DY(2:M-1, nn).^2);
  
  c(1:M-2)  =   -0.25*c_exh(1:M-2,nn).* c_hzey(2:M-1,nn) *dt^2  ./ ( epsi(1:M-2,nn) .* mu(2:M-1,nn)) ./ (DY(1:M-2, nn).*DY(2:M-1, nn));
  c(M-1)    =   0;
  
  B(1:M-1)  = c_exe(1:M-1,nn) .* Ex(1:M-1,nn) ...
                        + c_exh(1:M-1,nn) * dt ./(2*epsi(1:M-1,nn).*DY(1:M-1,nn)).* ...
                                         (c_hzhy(2:M,nn) .*Hzy(2:M,nn) + c_hzhx(2:M,nn) .*Hzx(2:M,nn)- ...
                                          c_hzhy(1:M-1,nn) .* Hzy(1:M-1,nn) - c_hzhx(1:M-1,nn) .* Hzx(1:M-1,nn)) ...
                        - c_exh(1:M-1,nn) * dt^2 ./(4*epsi(1:M-1,nn).*DY(1:M-1,nn)).*...
                                      ((c_hzex(2:M,nn)./ (mu(2:M,nn).*DX(2:M,nn))).*(Ey(2:M,nn+1)- Ey(2:M,nn))...
                                       -(c_hzex(1:M-1,nn)./ (mu(1:M-1,nn).*DX(1:M-1,nn))).*(Ey(1:M-1,nn+1)- Ey(1:M-1,nn)));

  % B(1:M-1)  =   B(1:M-1) - dt^2 .* ( c_itex(1:M-1,nn) ./ (8*epsi(1:M-1,nn).*mu(1:M-1,nn).*DY(1:M-1,nn)) ...
  %                     .* ((1./DX(2:M,nn)).*(delEy(2:M,nn+1) - delEy(2:M,nn)) - (1./DX(1:M-1,nn)).*(delEy(1:M-1,nn+1) - delEy(1:M-1,nn))))';

%% SF/TF correction for Ex           %% woprking on this:
if (source_present_x)
   B(source_y)   =   B(source_y)  +  0.5*c_exh(source_y,nn).*  ...
                                    (dt./DY(source_y,nn)./epsi(source_y,nn)).*Hzy_1D(1)*source_profile_x(nn);
end
 
                      
 alpha(1) = a(1);
       for i=2:M-1
           beta(i)  = b(i)/alpha(i-1);
           alpha(i) = a(i)-beta(i)*c(i-1);
       end

       g(1) = B(1);
       for i = 2:M-1
           g(i) = B(i) - beta(i)*g(i-1);
       end
       
       Ex_k_half(M-1,nn) = g(M-1)/alpha(M-1);
       for i = M-2:-1:1
           Ex_k_half(i,nn) = (g(i) - c(i)*Ex_k_half(i+1,nn))/alpha(i);
       end

%    above LU method is correct
%    Matrix method-   AA Ex = B
%    AA = gallery('tridiag',b(2:M-1),a,c(1:M-2));
%    test2  = inv(AA) * B';

end

% Ey update (x-propagating)          (Ex_k_half)
  Ey_k_half(1:M,2:N)    =    c_eye(1:M,2:N).*Ey(1:M,2:N)   ...
                     - 0.5*c_eyh(1:M,2:N).*(dt./DX(1:M,2:N)./epsi(1:M,2:N)).*(Hz(1:M,2:N) - Hz(1:M,1:N-1));
% Hzy update  (y-propagating)
  Hzy_k_half(2:M,1:N-1)  = c_hzhy(2:M,1:N-1).*Hzy(2:M,1:N-1)    ...
                       +0.5*c_hzey(2:M,1:N-1).*(dt./DY(2:M,1:N-1)./mu(2:M,1:N-1)).*(Ex_k_half(2:M,1:N-1) - Ex_k_half(1:M-1,1:N-1));
         
%Hzx update  (x-propagating)
  Hzx_k_half(2:M,1:N-1)   =  c_hzhx(2:M,1:N-1).*Hzx(2:M,1:N-1)    ...
                           - 0.5*c_hzex(2:M,1:N-1).*(dt./DX(2:M,1:N-1)./mu(2:M,1:N-1)).*(Ey(2:M,2:N) - Ey(2:M,1:N-1));
            
%  Hzy_n_half(y_sr,x_sr(1:end-1))  =  Hzy_n_half(y_sr,x_sr(1:end-1))  ...
%                                +   0.5*c_hzey(y_sr,x_sr(1:end-1)).*(dt./DY(y_sr,x_sr(1:end-1))./mu(y_sr,x_sr(1:end-1))).*Ex_1D(1);

 Hz_k_half            = Hzy_k_half+Hzx_k_half;
 
%     Hz_n_half(1,1:N-1)   = Hz_n_half(2,1:N-1);
%     Hz_n_half(1:M,N)     = Hz_n_half(1:M,N-1);
%     Ex_n_half(:,N)       = Ex_n_half(:,N-1);
%     Ey_n_half(:,N)       = Ex_n_half(:,N-1);
                   
%%

for mm=2:M
    
   clear a b c beta alpha g B;
   a(1:N-1)  = 1 + c_eyh(mm,2:N) * dt^2 ./ (4* epsi(mm, 2:N) .* DX(mm, 2:N)) ...
                             .*(c_hzex(mm,1:N-1) ./ (mu(mm,1:N-1).* DX(mm,1:N-1)) + c_hzex(mm,2:N) ./ ( mu(mm,2:N).*DX(mm,2:N)) );
     
  b(1)      =   0;
  b(2:N-1)  =   -c_eyh(mm,3:N).* c_hzex(mm,2:N-1) * dt^2 ./ ( 4* epsi(mm,3:N) .* mu(mm,2:N-1) .* DX(mm,3:N).*DX(mm,2:N-1));
  c(1:N-2)  =   -c_eyh(mm, 2:N-1).* c_hzex(mm,2:N-1) *dt^2  ./ (4*epsi(mm, 2:N-1) .* mu(mm,2:N-1) .* DX(mm,2:N-1).*DX(mm,2:N-1));
  c(N-1)    =   0; 

  B(1:N-1)  = c_eye(mm,2:N) .* Ey_k_half(mm,2:N) ...
                        -  c_eyh(mm,2:N) * dt./(2*epsi(mm,2:N).*DX(mm,2:N)) .* ...
                                         (c_hzhx(mm,2:N).*Hzx_k_half(mm,2:N) + c_hzhy(mm,2:N) .* Hzy_k_half(mm,2:N)   ...
                                       - c_hzhx(mm,1:N-1).*Hzx_k_half(mm,1:N-1) - c_hzhy(mm,1:N-1) .* Hzy_k_half(mm,1:N-1)) ...
                        -  c_eyh(mm,2:N) * dt^2 ./(2*epsi(mm,2:N).*DX(mm,2:N)) .* ...
                                     ((c_hzey(mm,1:N-1)./ (2*mu(mm,1:N-1).*DY(mm,1:N-1))).*(-Ex_k_half(mm,1:N-1) + Ex_k_half(mm-1,1:N-1))...
                                     +(c_hzey(mm,2:N)./ (2*mu(mm,2:N).*DY(mm,2:N))).*(Ex_k_half(mm,2:N)- Ex_k_half(mm-1,2:N)));
  
    alpha(1) = a(1);
       
       for i=2:N-1
           beta(i)  = b(i)/alpha(i-1);
           alpha(i) = a(i)-beta(i)*c(i-1);
       end

       g(1) = B(1);
       for i = 2:N-1
           g(i) = B(i) - beta(i)*g(i-1);
       end
       
       Ey_k(mm, N) = g(N-1)/alpha(N-1);
       for i = N-2:-1:1
           Ey_k(mm,i+1) = (g(i) - c(i)*Ey_k(mm,i+2))/alpha(i);
       end

% above LU method is correct as verfieid by following       
% Matrix method-   AA Ex = B
%    AA = gallery('tridiag',b(2:N-1),a,c(1:N-2));
%    test2  = inv(AA) * B';
%    Ey_n_full(mm,2:N) = test2;
      
end

% Ex update (x-propagating)          
  Ex_k(1:M-1,1:N)  =    c_exe(1:M-1,1:N).*Ex_k_half(1:M-1,1:N)  ...
                      + 0.5*c_exh(1:M-1,1:N) .* dt./ (DY(1:M-1,1:N) .* epsi(1:M-1,1:N)).*(Hz_k_half(2:M,1:N) - Hz_k_half(1:M-1,1:N));
          
 % Ex_k(1:M-1,1:N-1) = Ex_k(1:M-1,1:N-1)  ...
 %                - dt^2 * c_itex(1:M-1,1:N-1) ./ (8*epsi(1:M-1,1:N-1).*mu(1:M-1,1:N-1).*DY(1:M-1,1:N-1)) .* ...
 %           ((1./DX(2:M,1:N-1)).*(delEy(2:M,2:N)-delEy(2:M,1:N-1))-(1./DX(1:M-1,1:N-1)).*(delEy(1:M-1,2:N)-delEy(1:M-1,1:N-1)));
 % 
    % SF/TF correction for Ex       
     Ex_k(source_y,1:N)   =   Ex_k(source_y,1:N)  + 0.5*c_exh(source_y,1:N).*   ...
                                      (dt./DY(source_y,1:N)./epsi(source_y,1:N)).*Hzy_1D(1).*source_profile_x;   

%     Ex_k(POS1_2deg_y,POS1_2deg_x(1):POS1_2deg_x(2)) = Ex_k(POS1_2deg_y,POS1_2deg_x(1):POS1_2deg_x(2))...
%                                                     - 0.5*J1_e*dt./epsi(POS1_2deg_y,POS1_2deg_x(1):POS1_2deg_x(2));
                 
% Hzy update  (y-propagating)
  Hzy_k(2:M,1:N-1)  = c_hzhy(2:M,1:N-1).*Hzy_k_half(2:M,1:N-1)    ...
                  + 0.5*c_hzey(2:M,1:N-1).*(dt./DY(2:M,1:N-1)./mu(2:M,1:N-1)).*(Ex_k_half(2:M,1:N-1) - Ex_k_half(1:M-1,1:N-1));
                   
 % Hzx update  (x-propagating)
  Hzx_k(2:M,1:N-1)   =  c_hzhx(2:M,1:N-1).*Hzx_k_half(2:M,1:N-1)    ...
                   - 0.5*c_hzex(2:M,1:N-1).*(dt./DX(2:M,1:N-1)./mu(2:M,1:N-1)).*(Ey_k(2:M,2:N) - Ey_k(2:M,1:N-1));             
             
%  Ey =  Ey_n_full;
  
%   Hzy(y_sr,x_sr(1:end-1))  =  Hzy(y_sr,x_sr(1:end-1))  ...
%                                +  0.5* c_hzey(y_sr,x_sr(1:end-1)).*(dt./DY(y_sr,x_sr(1:end-1))./mu(y_sr,x_sr(1:end-1))).*Ex_1D(1);
                      
%  PEC boundary on top and far end
   Hz_k(2:M,1:N-1) = Hzy_k(2:M,1:N-1) +  Hzx_k(2:M,1:N-1);
    
Ex(POS1_2deg_y,POS1_2deg_x(1):POS1_2deg_x(2)) = Ex(POS1_2deg_y,POS1_2deg_x(1):POS1_2deg_x(2))...
                                                - 0.5*J*dt./epsi(POS1_2deg_y,POS1_2deg_x(1):POS1_2deg_x(2));


end

    Ex  = Ex_k;
    Ey  = Ey_k;
    Hz  = Hz_k;
    Hzx = Hzx_k;
    Hzy = Hzy_k;
    

%% Runtime plots (can change them to whatever we want to look at during the simulation)
    if mod(tt,tt_rec)==0
          subplot(4,1,1); pcolor(Ey);  caxis([-100 100]*10); drawnow; colorbar;
          subplot(4,1,2); plot(X_2DEG_edge, PIC_per_cell, '-R'); title('No. of particles'); drawnow;         
          subplot(4,1,3); plot(e);%, vel_times_PICPC(1:N_2DEG)./PIC_per_cell(1:N_2DEG), '-B');    
          title('e'); drawnow;
          subplot(4,1,4); plot(j_sh_density_macro, '-B');         title('J'); drawnow;

        fprintf('FDTD sim Time = %f ps;  Sim. time %f s\n', tt*dt/1e-12, toc);
     
   end
    

%% 2DEG channel MC (Wave-particle interaction part)

% time iteration for the MC mover
  
%====== 1.  field coupled to channel
% electric field at the edges of the cells  (N_2DEG-1)

e(1:N_2DEG+1) = Ex(POS1_2deg_y,POS1_2deg_x(1):POS1_2deg_x(2));

%====== 2. GENERATE NEW PARTICLES at every tt_insert

  if mod(tt,tt_insert)==0
        if (part_max - np <= np_insert)
             np_insert = part_max - np;  % random velocity and n matrices may have issue with this
        end
% particle position initialization in the very first cell.
            part_x(np+1:np+np_insert,1) = rand(np_insert, 1 ) * dx_2DEG;    
% normal velocty distribution initialization
            part_v(np+1:np+np_insert,1) = normrnd(v_drift, v_drift/200, [np_insert, 1]); 
            part_v_dc(np+1:np+np_insert,1) = part_v(np+1:np+np_insert,1);
            np = np + np_insert;
            disp('inserted')
  end

  
% e  e  e  e  e  e ... e  e 
%(N_2DEG+1 elements representing edge field values for first edge 
% to N_2DEG's cell's edge)
%   n  n  n  n  n  n ... n  (N_2DEG+1 elements representing PICs,
%  last cell is extra, so first N_2DEG elements are PIC CENTERS)


part_E = interp1(X_2DEG_edge, e(1:N_2DEG+1), part_x);
%e_center = interp1(X_2DEG_edge, e(1:N_2DEG+1), X_2DEG_CellCenter);

%====== 3. update each particle's position and velocity
% Simple explicit (Euler) push: force from the interpolated local field
% updates velocity, velocity updates position. Particles that exit the
% channel (x <= 0 or x >= channel length) are removed by swapping in the
% last active particle and decrementing the count (fast removal without
% shifting the whole array).
p = 1;
    
while ( p <= np )
       xx       = 1.0 + part_x(p) / dx_2DEG;
       part_pos = floor ( xx );

%  E field and Particle position arrays

      f    = -qe * part_E(p);
      %f         = -qe *e_center(part_pos);
      acc       = f / me;
      part_v(p)    = part_v(p) + acc * (dt);
      part_v_ac(p) = part_v_ac(p)+acc * (dt);
      part_x(p)    = part_x(p) + part_v(p) * (dt);
        
% Absorption:
% Kill the particle by replacing it with the last particle.
     if (part_x(p) >= len_x_2DEG)||(part_x(p) <= 0)
        part_x(p)= part_x(np);
        part_v(p) = part_v(np);
        part_v_dc(p) = part_v_dc(np);
        part_v_ac(p) = part_v_ac(np);
        part_E(p)    =part_E(np);
        np = np - 1;
        p = p - 1;
     end
      p = p + 1;
      
end
    
% =========  4. CALCULATE THE CURRENT DENSITY.
% Current density is built by depositing each particle's velocity onto
% the two nearest grid edges, summing over all particles per cell, then 
% scaling by the macro-particle factor(MF) and electron charge to get 
% physical sheet current density.

PIC_per_cell                = zeros(1, N_2DEG+1);
vel_times_PICPC             = zeros(1, N_2DEG+1);

    for p = 1 : np
        
      fi = 1.0 + part_x(p) / dx_2DEG;
      i  = floor ( fi );
      hx = fi - i;

      % particles in each cell
      PIC_per_cell(i)     = PIC_per_cell(i)   + (1.0 - hx );
      vel_times_PICPC(i)  = vel_times_PICPC(i)  ...     
                               + (1.0-hx)*(part_v_ac(p));
      
      PIC_per_cell(i+1)      = PIC_per_cell(i+1)  +  hx;
      vel_times_PICPC(i+1)   = vel_times_PICPC(i+1)     ...
                              +  hx*(part_v_ac(p));
     
    end

% There are N_2DEG+1 cells, we ignore the first and the last cell, 
% so we have N_2DEG-1 cells
    
% current at the center of N_2DEG+1 cells
 j_sh_density_macro                = vel_times_PICPC/dx_2DEG;
 j_sh_density_macro_1by2(2:N_2DEG) =  ...
                                    interp1(X_2DEG_CellCenter(1:N_2DEG),...
                                            j_sh_density_macro(1:N_2DEG), ...
                                            X_2DEG_edge(2:N_2DEG));
 
 j_sh_density_macro_1by2(1)        = j_sh_density_macro_1by2(2);
 j_sh_density_macro_1by2(N_2DEG+1) = j_sh_density_macro_1by2(N_2DEG);
  
% % Ignore 1st and last cell, now we have N_2DEG-1 cells with N_2DEG edges
% % current at theses N_2DEG edges is,

 J        =  -j_sh_density_macro_1by2*MF*qe/t_2DEG;


%% Matlab workspace saving perdically for additional data storage
        if(periodic_file_dump_on==1)
           if mod(tt-1,round(period_file_dump/dt))==0
               filNameWorkspace = [filename num2str(tt) '.mat'];
               save(filNameWorkspace);
           end
        end
end
end

        save(filename);
toc
%

disp('end of calculations')
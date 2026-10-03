%% ==============================================================================
%% AUTOMOTIVE DUAL-CLUTCH TRANSMISSION (DCT) SHAFT FATIGUE DESIGN
%% Machine Element Design I - Sharif University of Technology
%% Authors: Zahra Khazaei & Arshia Edalatpajouh | Instructor: Dr. Keyvan Mohammadi
%% ==============================================================================

clc;
clear;
close all;

%% 1. Engine Input Parameters and Material Properties
T_engine = 390.0;       % Maximum engine torque (N*m)
phi_deg = 20.0;         % Pressure angle (degrees)
phi = deg2rad(phi_deg);
FS = 3.0;               % Design Factor of Safety

% Material: AISI 1045 Hot Rolled (Shafts)
S_ut = 570.0;           % Ultimate tensile strength (MPa)
S_y = 310.0;            % Yield strength (MPa)

% Material: AISI 1020 Cold Drawn (Keys)
S_ut_key = 1020.0;      % Ultimate tensile strength (MPa)
S_y_key = 848.0;        % Yield strength (MPa)

% Gear Dimensions (m)
r_g1 = 0.150 / 2.0;     % Gear G1 pitch radius (0.075 m)
r_g3 = 0.250 / 2.0;     % Gear G3 pitch radius (0.125 m)

%% 2. Marin Modified Endurance Limit Calculation
S_e_prime = 0.504 * S_ut; % 285.0 MPa (Rotary beam endurance limit)
k_a = 1.58 * (S_ut ^ -0.085); % 0.9213 (Ground finish)
k_b = 1.24 * (50.0 ^ -0.107); % 0.8159 (Size factor for d ~ 50 mm)
k_c = 1.0;                    % Bending load
k_d = 1.0;                    % T < 200 °C
k_e = 0.702;                  % 99.99% Reliability (z = 3.719)
k_f = 1.0;                    % Miscellaneous factor

S_e = k_a * k_b * k_c * k_d * k_e * k_f * S_e_prime; % ~ 150.39 MPa

fprintf('====================================================================\n');
fprintf('AUTOMOTIVE DUAL-CLUTCH TRANSMISSION SHAFT DESIGN REPORT\n');
fprintf('====================================================================\n');
fprintf('Engine Torque: %.1f N*m | Power: 254 hp | Safety Factor: %.1f\n', T_engine, FS);
fprintf('AISI 1045 HR: Sut = %.1f MPa, Sy = %.1f MPa\n', S_ut, S_y);
fprintf('Marin Factors: ka = %.3f, kb = %.3f, ke = %.3f (99.99%%)\n', k_a, k_b, k_e);
fprintf('Modified Endurance Limit (Se): %.2f MPa\n', S_e);
fprintf('--------------------------------------------------------------------\n');

%% 3. Static Equilibrium and Gear Loads (1st Gear Mode)
% Gear 1 on Shaft 1 (Red Solid Shaft)
W_t1 = T_engine / r_g1;         % 5200 N (Tangential force)
W_r1 = W_t1 * tan(phi);         % 1892.6 N (Radial force)

% Span between Bearings B4 (x = 0) and B1 (x = 0.675 m); G1 at x = 0.630 m
L_span = 0.675;
x_g1 = 0.630;

% Vertical plane (xy):
B1_xy = (W_r1 * x_g1) / L_span; % 1766.46 N
B4_xy = W_r1 - B1_xy;           % 126.17 N

% Horizontal plane (xz):
B1_xz = (W_t1 * x_g1) / L_span; % 4853.33 N
B4_xz = W_t1 - B1_xz;           % 346.67 N

fprintf('Shaft 1 Gear 1 Loading: Wt = %.1f N, Wr = %.1f N\n', W_t1, W_r1);
fprintf('Bearing Reactions B4: xy = %.1f N, xz = %.1f N\n', B4_xy, B4_xz);
fprintf('Bearing Reactions B1: xy = %.1f N, xz = %.1f N\n', B1_xy, B1_xz);
fprintf('--------------------------------------------------------------------\n');

%% 4. Westinghouse Fatigue Sizing for Critical Sections
% Critical section locations along Shaft 1 (mm):
sec_x = [10, 40, 170, 410, 615, 630, 665];
sec_names = {'B4 Seat', 'Spline Entry', 'B3 Seat', 'B2 Seat', 'G1 Shoulder', 'G1 Keyway', 'B1 Seat'};
sec_Kf = [1.975, 1.975, 1.975, 1.975, 1.975, 1.300, 1.975];
sec_Kfs = [1.595, 1.595, 1.595, 1.595, 1.595, 1.300, 1.595];
sec_has_T = [false, true, true, true, true, false, false];

d_calc = zeros(size(sec_x));

for idx = 1:length(sec_x)
    x_m = sec_x(idx) / 1000.0;
    if x_m <= x_g1
        m_xy = B4_xy * x_m;
        m_xz = B4_xz * x_m;
    else
        m_xy = B4_xy * x_m - W_r1 * (x_m - x_g1);
        m_xz = B4_xz * x_m - W_t1 * (x_m - x_g1);
    end
    M_a = sqrt(m_xy^2 + m_xz^2);
    M_m = 0.0; % Rotating bending -> fully reversed
    
    if sec_has_T(idx)
        T_m = T_engine / 2.0; % 195 N*m
        T_a = T_engine / 2.0; % 195 N*m
    else
        T_m = 0.0;
        T_a = 0.0;
    end
    
    % Westinghouse Solid Shaft Diameter Formula:
    bending_term = (M_m / S_y) + (M_a / S_e);
    term_1 = (sec_Kf(idx)^2) * (bending_term^2);
    torsion_term = (T_m / S_y) + (T_a / S_e);
    term_2 = (sec_Kfs(idx)^2) * (torsion_term^2);
    inner_sqrt = sqrt(term_1 + term_2);
    d_calc(idx) = ((32.0 * 1000.0 * FS / pi) * inner_sqrt)^(1.0 / 3.0);
    
    fprintf('Section %d (%-14s @ x = %3d mm): Ma = %6.2f N*m -> d_min = %5.2f mm\n', ...
        idx, sec_names{idx}, sec_x(idx), M_a, d_calc(idx));
end

%% 5. Rectangular Key Sizing for Gear G1
% Key dimensions: width b = 16 mm, height h = 10 mm on d = 52 mm shaft
b_key = 16.0; % mm
h_key = 10.0; % mm
d_shaft_key = 52.0; % mm

% Required key length against crushing/bearing stress (Sy_key = 848 MPa, FS = 3)
L_key_min = (4.0 * T_engine * 1000.0 * FS) / (d_shaft_key * h_key * S_y_key);
fprintf('--------------------------------------------------------------------\n');
fprintf('Rectangular Key for Gear G1: 16 x 10 mm on d = 52 mm\n');
fprintf('Calculated Minimum Key Length (Crushing): %.2f mm\n', L_key_min);
fprintf('Selected Standard Key Length: 40.0 mm (with 40 mm bushing spacer)\n');
fprintf('====================================================================\n');

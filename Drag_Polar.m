clc
clear
close all;
%% Initial
syms CL
rho_cruise = 0.000632184;
mu = 5.008e-4;
U = 600.618; % ft/s
M = 0.6;
W_TO = 1.444106e+05;
q_bar = 0.5 * rho_cruise * (U ^ 2);
sweep_LE_w = 28;
sweep_LE_2 = 20.59; % Figure 8.45 % 8.46
sweep_LE_4 = 24.4;
d_f = 12.336;
l_f = 98.03;
t_over_c_w = 0.108; % Airfoil
t_over_c_max = 0.14459; % Airfoil 
t_over_c_h = 0.1; % Airfoil % Caution
t_over_c_v = 0.1; % Airfoil % Caution
b = 116.076; % GEO
i_h = deg2rad(-2);
i_w = 0;
i_v = 0;
eps_w = 2;
taper_w = 0.16;
z_h = 16.914; % GEO % Caution
z_w = 5.9055; % GEO % Caution
x_h = 86.03; % GEO % Caution
S_plf_fus = 1089.9711; % GEO
S_w = 1425.81;
S_h = 284.9906; %  GEO % Caution
S_v = 203.7801; %  GEO % Caution
S_fus = 127.4653;
S_b_fus = 7.11494784;
S_wet_w = 2.4824e+03; % Part1
S_wet_fus = 3433.7505; %  GEO
S_wet_h =  S_h * 2; %  GEO % Caution
S_wet_v =  S_v * 2; %  GEO % Caution
AR = 9.45;
dihedral_w = 6;
C_bar_w = 14.45;
AR_h = 6.26; % GEO % Caution
AR_v = 1.91; % GEO % Caution
z_wh = 6.264; % GEO % Caution

%% Lift Coefficient

%%% CL_alpha_Wing
m = (z_wh * 2) / b;
K_AR = (1/AR) - (1 / (1 + power(AR, 1.7))); 
K_landa = (10 - (3 * taper_w))/7;
K_mr = (1 - (m/2))/power(deg2rad(dihedral_w), 0.333);
CL_w = 1.05 * CL ;
Cl_alpha_Mzero_w = 5.16; % CFD Data % Caution
Cl_alpha_M_w = Cl_alpha_Mzero_w; % CFD Data % Caution
Cl_alpha_w = Cl_alpha_Mzero_w/sqrt(1 - power(M, 2));
k = Cl_alpha_M_w/(2* pi);
betha = sqrt(1 - power(M, 2));
CL_alpha_w = (2 * pi * AR) / (2 + sqrt(((power(AR,2) * power(betha, 2)) / (power(k, 2) * (1 + (power(tand(sweep_LE_2), 2) / power(betha, 2))))) + 4));
%%% Wing Zero-Lift Drag Coefficient
R_wf = 1;
R_Ls = 1.25;  % Figure 4.2
Cf_w = 0.0027; % Figure 4.3
L_prime = 1.2;
CD_0_w = (R_wf * Cf_w *  R_Ls * (1 + (L_prime * t_over_c_w) + (100 * power(t_over_c_w, 4))) * (S_wet_w / S_w));
%%% CL_alpha_horizental
etha_h = 1 - ((power(cosd((pi * z_h)/(2*z_w)), 2)) * (2.42 * sqrt(CD_0_w))/((x_h/C_bar_w)+0.3));
Cl_alpha_Mzero_h = 6.2452; % CFD Data % Caution
Cl_alpha_h = Cl_alpha_Mzero_h/sqrt(1 - power(M, 2));
k_h = Cl_alpha_h/(2* pi);
betha_h =  sqrt(1 - power(M, 2));
CL_alpha_h = (2 * pi * AR_h) / (2 + sqrt((power(AR_h,2) * power(betha_h, 2)) / ...
    (power(k_h, 2) * (1 + (power(tand(sweep_LE_2), 2) / power(betha_h, 2)))) + 4));  % Caution
%%% CL_alpha
deps_over_dalpha = 4.44 * power(K_AR * K_landa * K_mr * sqrt(cosd(sweep_LE_4)), 1.19); % Caution
K_wf = 1 + (0.025 * (d_f/b)) - (0.25 * power(d_f/b, 2));
CL_alpha_wf = K_wf * CL_alpha_w;
CL_alpha = CL_alpha_wf + (CL_alpha_h * etha_h * (S_h/S_w) * (1 - deps_over_dalpha));
%%% CL_zero
delta_alpha_0_over_eps_w = -0.5; % Figure8.41 % Caution
alpha_0_l = deg2rad(-4); % Airfoil % Caution
alpha_0_l_M_over_alpha_0_l_M0dot3 = 0.019; % Figure8.42 % Caution
alpha_zero_L_w = (alpha_0_l + ((delta_alpha_0_over_eps_w) * eps_w)) * (alpha_0_l_M_over_alpha_0_l_M0dot3);
CL_0_wf = (i_w - alpha_zero_L_w) * CL_alpha_wf;
CL_0 = CL_0_wf;
%% Wing
Re_w_c_bar = rho_cruise * U * C_bar_w / mu;
R_l_LER = (AR * 0.2) / cos(sweep_LE_w); % Read from Fig4.7

R = 0.94;
v = -0.0001; % Caution
w = -0.00015; % Caution

e = 1.1 * (Cl_alpha_w / AR) / ((R * (CL_alpha_w / AR)) + ((1 - R) * pi)); 

CD_L_w = (power(CL_w, 2) / (pi * AR * e)) + (2 * pi * CL_w * eps_w * v) + (4 * power(pi, 2) ...
    * power(eps_w, 2)*w);

CD_wing = CD_0_w + CD_L_w;
%% Fuslage
% Fuslage Zero-Lift Drag Coefficient
R_N_fus = 156740130; % fus Re
Cf_fus = 0.0016; % Figure 4.3

CD_0_fus_base = 1; % AAA % Caution
d_b = sqrt((4/pi) * S_b_fus); 
CD_b_fus = (0.029 * power((d_b/d_f), 3) / ...
    sqrt(CD_0_fus_base * (S_w/S_fus))) * (S_fus/S_w);

CD_0_fus = (R_wf * Cf_fus * (1 + (60/power(l_f/d_f,3 )) ...
    + 0.0025 * (l_f/d_f)) * S_wet_fus / S_w) + CD_b_fus;
%

etha = 0.65; % Figure 4.19

alpha = ((W_TO/(q_bar * S_w)) - CL_0) / CL_alpha;

C_d_c = 0.12; % Figure 4.20 % Caution

CD_L_fus = ((2 * power(alpha, 2) * S_b_fus) ...
    / S_w) + (etha * C_d_c * power(alpha, 3) ...
    * S_plf_fus) / S_w;


CD_fus = CD_0_fus + CD_L_fus;

%% Empenage

% Horizontal
R_Ls_h = 0.118; % Figure 4.2 % Caution
L_prime_h = 1.2; % Figure 4.4 % Caution
Cf_h = 0.005; % Figure 4,3 % Caution
CD_0_h = (R_wf* R_Ls_h * Cf_h * (1 + (L_prime_h ...
    * t_over_c_h) + (100 * power(t_over_c_h, 4))) ...
    * S_wet_h) / S_h; % S_h might be replace with S_w % Caution
%
e_h = 0.5;
alpha_0_L_h = 0;
alpha_h = alpha * (1 - deps_over_dalpha) + i_h;
CL_h = CL_alpha_h * (alpha_h - alpha_0_L_h);
CD_L_h = ((power(CL_h, 2) / (pi * AR_h * e_h)) * S_h / S_w);

% Vertical
R_Ls_v = 0.15; % Figure 4.2 % Caution
L_prime_v = 1.2; % Figure 4.4 % Caution
Cf_v = 0.005; % Figure 4,3 % Caution
CD_0_v = (R_wf* R_Ls_v * Cf_v * (1 + (L_prime_v ...
    * t_over_c_v) + (100 * power(t_over_c_v, 4))) ...
    * S_wet_v) / S_v; % S_v might be replace with S_w % Caution
%
e_v = 0.5;
CL_v = 0; % Page 69 Part 6 says this
CD_L_v = ((power(CL_v, 2) / (pi * AR_v * e_v)) * S_v / S_w);

CD_emp = CD_L_v + CD_0_v + CD_L_h + CD_0_h;

%% Gear
CD_gear = 4.300013e-02+((4.210448e-02)*CL^2); 

%% Naccele & Pylons
CD_np = 0.00049384134448129;
%%
CD_store = 0;
%%
CD_cw = 0.00076578;
%% Nozzele
CD_Nozzele = 0.006225520;
%%
CD_trim = 0.0154599245351633;
%%
CD_int = 0.0067433;
%%
CD_misc = 0.0067876;
%% Landing Flap (deg=40) Single Slotted flap
CD_flap_L = 0.02996719 + 0.0567252;
%% Take-off Flap (deg=15) Single Slotted flap
CD_flap_TO = 0.0049268 + 0.00952896;
%%
CD_Airplane = CD_wing + CD_fus + CD_emp + CD_np + CD_cw ...
    + CD_store + CD_trim + CD_int + CD_misc + CD_Nozzele;

Cl = linspace(0,5,1000);

CD_Clean = double(subs(CD_Airplane,CL,Cl));
CD_TO_gear = double(subs(4.300013e-02+(4.210448e-02)* (CL^2),CL,Cl));
CD_L_gear = double(subs(4.300013e-02+(4.491145e-02)* (CL^2),CL,Cl));
CD_TO = double(subs(CD_Airplane + CD_flap_TO,CL,Cl));
CD_L = double(subs(CD_Airplane + CD_flap_L,CL,Cl));


figure(1), hold on;
plot(CD_Clean,Cl, 'LineWidth',3)
plot(CD_TO,Cl, 'LineWidth', 3, 'LineStyle', ':');
plot(CD_L,Cl, 'LineWidth', 3, 'LineStyle', '-.');
plot(CD_TO_gear,Cl, 'LineWidth', 3, 'LineStyle', ':');
plot(CD_L_gear,Cl, 'LineWidth', 3, 'LineStyle', '-.');
title("\fontsize{18}\fontname{Times}Drag Polar");
xlabel("\fontsize{14}\fontname{Times}Drag Coefficient");
ylabel("\fontsize{14}\fontname{Times}Lift Coefficient");
grid on;
xlim([0,0.3]);
ylim([0, 2]);

legend("\fontsize{12}\fontname{Times}Clean", ... 
    "\fontsize{12}\fontname{Times}Take-off Gear-Down", ...
    "\fontsize{12}\fontname{Times}Landing Gear-Down", ...
    "\fontsize{12}\fontname{Times}Take-off", ...
    "\fontsize{12}\fontname{Times}Landing")

Cl = W_TO / (0.5 * S_w * (U ^ 2) * rho_cruise);
CD_Airplane = double(subs(CD_Airplane,CL,Cl));
L_over_D = Cl/CD_Airplane;
fprintf("Lift to Drag Ratio: %.2d\n", L_over_D);
fprintf("Lift Coefficient in Cruise: %.2d\n", Cl);
fprintf("Total Drag Coefficient in Clean phase: %.2d\n", CD_Airplane);
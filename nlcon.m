function [c, ceq] = nlcon(x)
%NLCON Summary of this function goes here
%   Detailed explanation goes here

% Rcb_l = x(1);
% Tst_l = x(2);
% Hst_l = x(3);
% Tsp_l = x(4);
% Hsp_l = x(5);
% Tl_l = x(6);
% Hl_l = x(7);
% Rcb_k = x(8);
% Tst_k = x(9);
% Hst_k = x(10);
% Tsp_k = x(11);
% Hsp_k = x(12);
% Tl_k = x(13);
% Hl_k = x(14);

% No nonlinear equalitird
ceq = [];

%% Laptop and Keyboard Stands

% Material properties
rho = 350; % (x(1)) 350 kg / m^3
Est = 8500 * 1000 * 1000; % Pa
Ecb = 8500 * 1000 * 1000; % Pa

% Physics constants
g = 9.81;

% Constant dimensions
theta_max = 45;
theta_min = 5;

% Bearings
R_b = .0254 / 2 / 2; % half inch diameter
sigma_b = 6200 * 1000; % Pa

% Locking Pins
R_l = .0254 / 2 / 2; % half inch diameter
sigma_l = 6200 * 1000; % Pa

% Factor of safety
K = 2;

%% Laptop Stand

% Constant Dimensions
% From requirements (and theta_max)
Lst_l = 0.384;
Lcb_l = 0.35;

% Based on the above 2
Lsp_l = Lst_l;
Ll_l = Lst_l / 2;

% To be conservative on strength requirements
dsp_l = 0;

% Parameterized Dimensions
% Crossbar Dimensions - Laptop Stand
Rcb_l = x(1);

% Strut Dimensions - Laptop Stand
Tst_l = x(2);
Hst_l = x(3);

% Support Dimensions - Laptop Stand
Tsp_l = x(4);
Hsp_l = x(5);

% Locking Bar Dimensions - Laptop Stand
Tl_l = x(6);
Hl_l = x(7);
dl_l = (Lcb_l / 2) - Tst_l;

% Specification of constraints
% Should be conservative in most situations - need to think about this

Wcb_l = rho * (pi * Rcb_l^2 * Lcb_l) * g;
Wst_l = rho * (Lst_l * Tst_l * Hst_l) * g;
Wsp_l = rho * (Lsp_l * Tsp_l * Hsp_l) * g;
Wl_l = rho * (Tl_l * Ll_l * Hl_l) * g;
Flz = 24.45 + Wsp_l;

c(34) = (5 * 2 * Lsp_l^3 / 6 / Est / (Hsp_l * Tsp_l^3 / 12)) - (0.001 / K);

c(18) = (Flz * (Lst_l * cos(theta_min))^3 / 48 / Est / (Hsp_l^3 * Tsp_l / 12)) - (0.001 / K);

Icb = (pi * Rcb_l^4) / 4;

Fna = (Flz / 4) + (2 * Wcb_l) + (3*Wst_l) + Wl_l;

dmax_cb = Flz*dl_l*(3*Lcb_l^2 - 4*dl_l^2) / (48 * Ecb * (pi * Rcb_l^4 / 4));
c(1) = dmax_cb - (0.001 / K);

dmax_st = @(theta) Flz*(cos(theta) + sin(theta)) * (Lst_l^3) / (48 * Est * (Tst_l * Hst_l^3 / 12));
c(2) = dmax_st(theta_max) - (0.003 / K);
c(20) = dmax_st(theta_min) - (0.003 / K);
d_y_st_l = 5 * 2 * (Lst_l / 2)^3 / 6 / Est / (Hst_l * Tst_l^3 / 12);
c(16) = d_y_st_l - (0.001 / K);

% shear force on s pin
Sz = (Flz/2) - Wcb_l - Wst_l + (dsp_l * Flz / (2*Lst_l));
Pz = Flz/2;
Sx = @(theta) (Sz/(tan(theta))) - Pz;

tau_st_l = @(theta) sqrt(Sx(theta)^2 + Sz^2) / (pi*Rcb_l^2);
c(3) = tau_st_l(theta_max) - (sigma_l / K);
c(21) = tau_st_l(theta_min) - (sigma_l / K);

% Shear force on pin V in lockbar
Gz = Fna + Wcb_l/2;
Gx = @(theta) -Sx(theta);
Uz = Gz - Wl_l;
Ux = @(theta) Gx(theta);

tau_l_l = @(theta) sqrt(Uz^2 + Ux(theta)^2) / (pi*R_l^2);
c(4) = tau_l_l(theta_max) - (sigma_l / K);
c(22) = tau_l_l(theta_min) - (sigma_l / K);
c(36) = (5 * 2 * Ll_l^3 / 6 / Est / (Hl_l * Tl_l^3 / 12)) - (0.001 / K);

dmax = @(theta) (dl_l * (3*(Lcb_l)^2 - 4*(dl_l)^2) * sqrt((Gx(theta)^2) + (Gz^2))) / (24*Ecb*Icb);
c(5) = dmax(theta_max) - (0.001 / K);
c(23) = dmax(theta_min) - (0.001 / K);

%% Keyboard Stand

% Constant dimensions
Lst_k = 0.318;
Lcb_k = 0.4;

% Based on the above 2
Lsp_k = Lst_k;
Ll_k = Lst_k / 2;

% To be conservative on strength requirements
dsp_k = 0;

% Crossbar Bar Dimensions - Keyboard Stand
Rcb_k = x(8);

% Strut Dimensions - Keyboard Stand
Tst_k = x(9);
Hst_k = x(10);

% Support Dimensions - Keyboard Stand
Tsp_k = x(11);
Hsp_k = x(12);

% Locking Bar Dimensions - Keyboard Stand
Tl_k = x(13);
Hl_k = x(14);
dl_k = (Lcb_k / 2) - Tst_k;

% Input Params - Keyboard Stand
Wk = 0.5 * g; % kg * m / s^2
Fkx = 5.2; % N, 25sin12
Wsp_k = (Lsp_k * Tsp_k * Hsp_k) * rho * g;
Fkz = 24.45 + Wsp_k; % N, 25cos12
W_kz = Wk + Fkz;
theta_max = 45;
theta_min = 5;

% Objective Function Params - Keyboard Stand
Wcb_k = (pi * Rcb_k^2 * Lcb_k) * rho * g;
Wst_k = (Lst_k * Tst_k * Hst_k) * rho * g;
Wl_k = (Tl_k * Ll_k * Hl_k) * rho * g;

% External Reaction Forces - Keyboard Stand
Ffx = 0.5 * (Fkx);
Fng = @(theta) 0.5*(W_kz) + 2*Wst_k + 1.5*Wcb_k + Wl_k + 2*Fkx*tan(theta);
% Fnf = 0.5*(W_kz) + 2*Wst_k + 1.5*Wcb_k + Wl_k - 2*Fkx*tan(theta);

% Top Crossbar Reaction Forces - Keyboard Stand
Rcb_az_k = @(theta) W_kz * (Lsp_k / (2 * Lst_k * cos(theta)));
Rcb_bz_k = @(theta) W_kz * (1 - (Lsp_k / (2 * Lst_k * cos(theta))));
Rcb_bx_k = -Fkx;

Acb_b_k = Tsp_k * (2 * Rcb_k);
c(6) = (sqrt(Rcb_bx_k^2 + Rcb_bz_k(theta_max)^2) / Acb_b_k) - (sigma_b / K);
c(24) = (sqrt(Rcb_bx_k^2 + Rcb_bz_k(theta_min)^2) / Acb_b_k) - (sigma_b / K);

c(19) = W_kz * (Lst_k * cos(theta_min))^3 / 48 / Est / (Hsp_k^3 * Tsp_k / 12) - (0.001 / K);

c(35) = (5 * 2 * Lsp_k^3 / 6 / Est / (Hsp_k * Tsp_k^3 / 12)) - (0.001 / K);

% Top Crossbar Bending Deflection - Keyboard Stand
Ist = (Wst_k^3 * Tst_k) / 12;
Icb = (pi * Rcb_k^4) / 4;

Jcb_az_k = @(theta) ((Wcb_k * Lcb_k^3) + (2 * Rcb_az_k(theta) * ((Lcb_k/2) - dsp_k) * ((3 * Lcb_k^2) - (4 * ((Lcb_k/2) - dsp_k)^2))))/ (48 * Ecb * Icb);
Jcb_bz_k = @(theta) ((Wcb_k * Lcb_k^3) + (2 * Rcb_bz_k(theta) * ((Lcb_k/2) - dsp_k) * ((3 * Lcb_k^2) - (4 * ((Lcb_k/2) - dsp_k)^2))))/ (48 * Ecb * Icb);
Jcb_bx_k = (2 * Rcb_bx_k * ((Lcb_k/2) - dsp_k) * ((3 * Lcb_k^2) - (4 * ((Lcb_k/2) - dsp_k)^2))) / (48 * Ecb * Icb);

% Top crossbar bending deflection
c(7) = sqrt(Jcb_bz_k(theta_max)^2 + Jcb_bx_k^2) - (0.001 / K);
c(25) = sqrt(Jcb_bz_k(theta_min)^2 + Jcb_bx_k^2) - (0.001 / K);

c(8) = Jcb_az_k(theta_max) - (0.001 / K);
c(25) = Jcb_az_k(theta_min) - (0.001 / K);

% Strut Deflection - Keyboard Stand
Kx = Rcb_bx_k;
Kz = @(theta) Rcb_bz_k(theta);
Lz = @(theta) Rcb_az_k(theta);

Jkm = @(theta) (((Kx*sin(theta) + Kz(theta)*cos(theta)) * (Lst_k/2)^3) / (3*Est*Ist)) + ((Wst_k/2) * cos(theta) * (Lst_k/4)^2 - ((3 * Lst_k/2) - (Lst_k/4))) / (6*Est*Ist);
Jln = @(theta) ((Lz(theta)*cos(theta) * (Lst_k/2)^3) / (3*Est*Ist)) + (((Wst_k/2) * cos(theta) * (Lst_k/4)^2 * ((3*(Lst_k/2)) - Lst_k/4)) / (6*Est*Ist));

c(9) = Jkm(theta_max) - (0.005 / K);
c(26) = Jkm(theta_min) - (0.005 / K);
c(10) = Jln(theta_max) - (0.005 / K);
c(27) = Jln(theta_min) - (0.005 / K);

d_y_st_k = 5 * 2 * (Lst_k / 2)^3 / 6 / Est / (Hst_k * Tst_k^3 / 12);
c(17) = d_y_st_k - (0.001 / K);

% Reaction forces between Xs - Keyboard Stand
% Dz = (0.5*Wk) + Wst_k + (0.75*Wcb_k) + (Fkx*tan(theta));
Cz = @(theta) (0.5*W_kz) + Wst_k + (0.75*Wcb_k) + (Fkx*tan(theta));
Cx = @(theta) (Cz(theta) + Lz(theta)) / tan(theta);
% Dx = -Fkx - Cx;

% Bottom Strut Reaction Forces - Keyboard Stand
Az = Wl_k + (0.5*Wcb_k);
Bz = Az;
Ez = @(theta) Bz + Cz(theta) - Fng(theta);
Ex = @(theta) ((1/tan(theta)) * ((2*Cz(theta)) - Ez(theta))) + (2*Cx(theta)) - (Hst_k * Ffx / (Lst_k * sin(theta)));

% Shear Failure - Keyboard Stand
c(11) = (sqrt(Ex(theta_max)^2 + Ez(theta_max)^2) / pi*R_b^2) - (sigma_b / K);
c(28) = (sqrt(Ex(theta_min)^2 + Ez(theta_min)^2) / pi*R_b^2) - (sigma_b / K);

% Bottom Crossbar Reaction Forces - Keyboard Stand
Ax = @(theta) Ex(theta) - Cx(theta) - Ffx;
% Bx = Ax;
Qx = @(theta) Ax(theta);
Qz = Az - (0.5 * Wcb_k);

c(12) = (sqrt(Qx(theta_max)^2 + Qz^2) / (2 * Rcb_k * Tl_k)) - (sigma_b / K);
c(29) = (sqrt(Qx(theta_min)^2 + Qz^2) / (2 * Rcb_k * Tl_k)) - (sigma_b / K);

% Locking Pin Shear Failure - Keyboard Stand
Ux = @(theta) Qx(theta);
Uz = Qz - Wl_k;
sigma_u = @(theta) sqrt(Ux(theta)^2 + Uz^2) / (pi * R_l^2);
c(13) = sigma_u(theta_max) - (sigma_l / K);
c(30) = sigma_u(theta_min) - (sigma_l / K);

c(14) = (sqrt(Ux(theta_max)^2 + Uz^2) / (2 * R_l * Tl_k)) - (sigma_l / K);
c(31) = (sqrt(Ux(theta_min)^2 + Uz^2) / (2 * R_l * Tl_k)) - (sigma_l / K);

c(33) = (5 * 2 * Ll_k^3 / 6 / Est / (Hl_k * Tl_k^3 / 12)) - (0.001 / K);

% Bottom Crossbar Bending Deflection - Keyboard Stand
Uaf_z = ((Wcb_k * Lcb_k^3) + (2 * Qz * ((Lcb_k/2) - dl_k) * ((3*Lcb_k^2) - (4*((Lcb_k/2) - dl_k)^2)))) / (48*Ecb*Icb);
Uaf_x = @(theta) (2 * Qx(theta) * ((Lcb_k/2) - dl_k) * ((3*Lcb_k^2) - (4 * ((Lcb_k/2) - dl_k)^2))) / (48*Ecb*Icb);

c(15) = sqrt(Uaf_z^2 + Uaf_x(theta_max)^2) - (0.001 / K);
c(32) = sqrt(Uaf_z^2 + Uaf_x(theta_min)^2) - (0.001 / K);

end


%==========================================================================
% AUTHOR: David L. Tran
%
% DESCRIPTION: Numerically solves the conservative 1D Euler equation in a
% finite volume formulation using the (1) MacCormack scheme with artificial
% dissipation, (2) Lax-Wendroff method with artificial dissipation, and (3)
% the Steger-Warming flux-vector splitting scheme. Initial conditions
% correspond to those of Sod's shock tube. This is a Riemann problem
% wherein sharp discontinuities are present in the exact solution.
%
%==========================================================================

%% Clear Cache
clc; close all; clearvars;

%% Variables

%left of diaphragm
p_L = 1e5;          %pressure in [N/m^2]
rho_L = 1;          %density in [kg/m^3]
u_L = 0;            %velocity in [m/s]

%right of diaphragm
p_R = 1e4;          %pressure in [N/m^2]
rho_R = 0.125;      %density in [kg/m^3]
u_R = 0;            %velocity in [m/s]

%gas properties
gamma = 1.4;        %specific heat ratio of air
R = 287;            %specific gas constant in [J/(kg*K)]

%grid properties in space and time
IL_1 = 300;         %number of grid points for sim 1
IL_2 = 600;         %number of grid points for sim 2
IL = [IL_1; IL_2];

t_steps_1 = 70;     %number of time steps for sim 1
t_steps_2 = 140;    %number of time steps for sim 2
t_steps = [t_steps_1; t_steps_2];

epsilon = 0.6;      %artificial dissipation parameter for LW and MacCormack

nu_1 = 0.8;         %first Courant number
nu_2 = 1.3;         %second Courant number
nu = [nu_1; nu_2];

%Analytical calculation
P = 3.031;          %pressure ratio p_2/p_R
alpha = (gamma + 1) / (gamma - 1);

c_L = sqrt(gamma * p_L / rho_L); %speed of sound on left side [m/s]

c_R = sqrt(gamma * p_R / rho_R); %speed of sound on right side [m/s]

p_2 = p_R * P;
rho_2 = rho_R * (1 + alpha * P) / (alpha + P);
u_S = c_R * sqrt((gamma - 1) / (2 * gamma) * (1 + alpha * P)) + u_R;
u_2 = (P - 1) / (sqrt(1 + alpha * P))*(1 / sqrt(gamma * (gamma - 1) / 2)) * c_R + u_R;

p_3 = p_2;
u_3 = u_2;
rho_3 = (p_3 / p_L)^(1 / gamma) * rho_L;
c_3 = sqrt(gamma * p_3 / rho_3);

%GRID 1 ARRAYS (time steps = 70, IL = 300)
delta_x_1 = 6.25e-2;        %grid 1 spacing in [m]

x_arr_1_left = linspace(0-delta_x_1/2, (0-delta_x_1/2 * IL(1)), IL(1)/2);
x_arr_1_right = linspace(0+delta_x_1/2, (0+delta_x_1/2 * IL(1)), IL(1)/2);
x_arr_1 = [flip(x_arr_1_left), x_arr_1_right];

%GRID 2 ARRAYS (time steps = 140, IL = 600)
x_arr_2_left = linspace((0-delta_x_1/2)/2, (0-delta_x_1/2 * IL(1)), IL(2)/2);
x_arr_2_right = linspace((0+delta_x_1/2)/2, (0+delta_x_1/2 * IL(1)), IL(2)/2);
x_arr_2 = [flip(x_arr_2_left), x_arr_2_right];
delta_x_2 = x_arr_2(2) - x_arr_2(1);

%MacCormack with Artificial Dissipation Numerical Arrays
U_MC_1 = zeros(IL(1), 3, t_steps(1)+1, length(nu));
U_MC_1_p = zeros(IL(1), 3, t_steps(1)+1, length(nu));
E_MC_1 = zeros(IL(1), 3, t_steps(1)+1, length(nu));
E_MC_1_p = zeros(IL(1), 3, t_steps(1)+1, length(nu));
D_1 = zeros(IL(1), 3, t_steps(1)+1, length(nu));
D_1_p = zeros(IL(1), 3, t_steps(1)+1, length(nu));
p_MC_1 = zeros(IL(1), t_steps(1)+1, length(nu));
p_MC_1_p = zeros(IL(1), t_steps(1)+1, length(nu));
c_MC = zeros(IL(1), 1);
u_plus_c_MC = zeros(IL(1), 1);

U_MC_2 = zeros(IL(2), 3, t_steps(2)+1, length(nu));
U_MC_2_p = zeros(IL(2), 3, t_steps(2)+1, length(nu));
E_MC_2 = zeros(IL(2), 3, t_steps(2)+1, length(nu));
E_MC_2_p = zeros(IL(2), 3, t_steps(2)+1, length(nu));
D_2 = zeros(IL(2), 3, t_steps(2)+1, length(nu));
D_2_p = zeros(IL(2), 3, t_steps(2)+1, length(nu));
p_MC_2 = zeros(IL(2), t_steps(2)+1, length(nu));
p_MC_2_p = zeros(IL(2), t_steps(2)+1, length(nu));
c_MC_2 = zeros(IL(2), 1);
u_plus_c_MC_2 = zeros(IL(2), 1);

%MacCormack Analytical Arrays
U_MC_1_anlyt = zeros(IL(1), 3, length(nu));
E_MC_1_anlyt = zeros(IL(1), 3, length(nu));

U_MC_2_anlyt = zeros(IL(2), 3, length(nu));
E_MC_2_anlyt = zeros(IL(2), 3, length(nu));


%MacCormack with Artificial Dissipation Time Arrays
t_total_MC_1 = zeros(length(nu), 1);
delta_t_MC = zeros(length(nu), 1);
delta_t_MC(1) = sqrt(gamma * p_L / rho_L);
delta_t_MC(end) = sqrt(gamma * p_R / rho_R);

t_total_MC_2 = zeros(length(nu), 1);
delta_t_MC_2 = zeros(length(nu), 1);


%Lax-Wendroff with Artificial Dissipation arrays
U_LW_1 = zeros(IL(1), 3, t_steps(1)+1, length(nu));
E_LW_1 = zeros(IL(1), 3, t_steps(1)+1, length(nu));
p_LW_1 = zeros(IL(1), t_steps(1)+1, length(nu));
c_LW_1 = zeros(IL(1),1);
u_plus_c_LW = zeros(IL(1), 1);

t_total_LW_1 = zeros(length(nu), 1);
delta_t_LW = zeros(length(nu), 1);

U_LW_1_anlyt = zeros(IL(1), 3, length(nu));
E_LW_1_anlyt = zeros(IL(1), 3, length(nu));

U_LW_2 = zeros(IL(2), 3, t_steps(2)+1, length(nu));
E_LW_2 = zeros(IL(2), 3, t_steps(2)+1, length(nu));
p_LW_2 = zeros(IL(2), t_steps(2)+1, length(nu));
c_LW_2 = zeros(IL(2),1);
u_plus_c_LW_2 = zeros(IL(2), 1);

t_total_LW_2 = zeros(length(nu), 1);
delta_t_LW_2 = zeros(length(nu), 1);

U_LW_2_anlyt = zeros(IL(2), 3, length(nu));
E_LW_2_anlyt = zeros(IL(2), 3, length(nu));


%Steger-Warming scheme arrays
U_SW_1 = zeros(IL(1), 3, t_steps(1)+1, length(nu));
E_SW_1 = zeros(IL(1), 3, t_steps(1)+1, length(nu));
c_SW_1 = zeros(IL(1),1);
u_plus_c_SW = zeros(IL(1), 1);

t_total_SW_1 = zeros(length(nu), 1);
delta_t_SW = zeros(length(nu), 1);

U_SW_1_anlyt = zeros(IL(1), 3, length(nu));
E_SW_1_anlyt = zeros(IL(1), 3, length(nu));

t_total_SW_2 = zeros(length(nu), 1);
delta_t_SW_2 = zeros(length(nu), 1);

U_SW_2_anlyt = zeros(IL(2), 3, length(nu));
E_SW_2_anlyt = zeros(IL(2), 3, length(nu));


%% Loops

%%%%%%%%%%%
%U and E terms
%U1 = rho
%U2 = rho*u
%U3 = e
%E1 = rho*u = U2
%E2 = rho*(u^2)+P = (U2^2/U1)+(gamma-1)*(U3-(U2^2)/(2*U1))
%E3 = (e+P)*u = U2*U3/U1+(U2/U1)*(gamma-1)*(U3-(U2^2)/(2*U1))
%%%%%%%%%%%%%%

%Grid 1
for k = 1:length(nu)
%Loop through all Courant numbers
    for i = 1:IL(1)
    %Set initial values for all of space
        if x_arr_1(i) <= 0
            %Left half
            U_MC_1(i,1,1,k) = rho_L;
            U_MC_1(i,2,1,k) = rho_L * u_L;
            U_MC_1(i,3,1,k) = p_L / (gamma - 1) + (rho_L / 2) * u_L^2;

            p_MC_1(i,1,k) = p_L;
            
            E_MC_1(i,1,1,k) = U_MC_1(i,2,1,k);
            E_MC_1(i,2,1,k) = U_MC_1(i,2,1,k)^2 / U_MC_1(i,1,1,k) + (gamma - 1) * (U_MC_1(i,3,1,k) - (U_MC_1(i,2,1,k))^2 / (2*U_MC_1(i,1,1,k)));
            E_MC_1(i,3,1,k) = U_MC_1(i,2,1,k) * U_MC_1(i,3,1,k) / U_MC_1(i,1,1,k) + (U_MC_1(i,2,1,k) / U_MC_1(i,1,1,k)) * (gamma - 1) * (U_MC_1(i,3,1,k) - (U_MC_1(i,2,1,k))^2 / (2*U_MC_1(i,1,1,k)));
            
            c_MC(i) = sqrt(gamma * ((gamma - 1) * (U_MC_1(i,3,1,k) - (U_MC_1(i,2,1,k))^2 / (2*U_MC_1(i,1,1,k)))) / (U_MC_1(i,1,1,k)));
            u_plus_c_MC(i,1) = abs(U_MC_1(i, 2, 1, k) / U_MC_1(i, 1, 1, k)) + c_MC(i);

        else
            %Right half
            U_MC_1(i,1,1,k) = rho_R;
            U_MC_1(i,2,1,k) = rho_R * u_R;
            U_MC_1(i,3,1,k) = p_R / (gamma - 1) + (rho_R / 2) * u_R^2;

            p_MC_1(i,1,k) = p_R;

            E_MC_1(i,1,1,k) = U_MC_1(i,2,1,k);
            E_MC_1(i,2,1,k) = U_MC_1(i,2,1,k)^2 / U_MC_1(i,1,1,k) + (gamma - 1) * (U_MC_1(i,3,1,k) - (U_MC_1(i,2,1,k))^2 / (2*U_MC_1(i,1,1,k)));
            E_MC_1(i,3,1,k) = U_MC_1(i,2,1,k) * U_MC_1(i,3,1,k) / U_MC_1(i,1,1,k) + (U_MC_1(i,2,1,k) / U_MC_1(i,1,1,k)) * (gamma - 1) * (U_MC_1(i,3,1,k) - (U_MC_1(i,2,1,k))^2 / (2*U_MC_1(i,1,1,k)));
            
            c_MC(i,1) = sqrt(gamma * ((gamma - 1) * (U_MC_1(i,3,1,k) - (U_MC_1(i,2,1,k))^2 / (2*U_MC_1(i,1,1,k)))) / (U_MC_1(i,1,1,k)));
            u_plus_c_MC(i,1) = abs(U_MC_1(i, 2, 1, k) / U_MC_1(i, 1, 1, k)) + c_MC(i);

        end

    end

    delta_t_MC(k, 1) = nu(k) * delta_x_1 / max(u_plus_c_MC);
    t_total_MC_1(k, 1) = delta_t_MC(k,1);

    %Set Dirichlet BCs for all of time
    U_MC_1(1, 1, :, k) = rho_L;
    U_MC_1(1, 2, :, k) = rho_L * u_L;
    U_MC_1(1, 3, :, k) = p_L / (gamma - 1) + (rho_L / 2) * u_L^2;
    U_MC_1(end, 1, :, k) = rho_R;
    U_MC_1(end, 2, :, k) = rho_R * u_R;
    U_MC_1(end, 3, :, k) = p_R / (gamma - 1) + (rho_R / 2) * u_R^2;

    E_MC_1(1, 1, :, k) = U_MC_1(1, 2, 1, k);
    E_MC_1(1, 2, :, k) = U_MC_1(1,2,1,k)^2 / U_MC_1(1,1,1,k) + (gamma - 1) * (U_MC_1(1,3,1,k) - (U_MC_1(1,2,1,k))^2 / (2*U_MC_1(1,1,1,k)));
    E_MC_1(1, 3, :, k) = U_MC_1(1,2,1,k) * U_MC_1(1,3,1,k) / U_MC_1(1,1,1,k) + (U_MC_1(1,2,1,k) / U_MC_1(1,1,1,k)) * (gamma - 1) * (U_MC_1(1,3,1,k) - (U_MC_1(1,2,1,k))^2 / (2*U_MC_1(1,1,1,k)));
    E_MC_1(end, 1, :, k) = U_MC_1(end, 2, 1, k);
    E_MC_1(end, 2, :, k) = U_MC_1(end,2,1,k)^2 / U_MC_1(end,1,1,k) + (gamma - 1) * (U_MC_1(end,3,1,k) - (U_MC_1(end,2,1,k))^2 / (2*U_MC_1(end,1,1,k)));
    E_MC_1(end, 3, :, k) = U_MC_1(end,2,1,k) * U_MC_1(end,3,1,k) / U_MC_1(end,1,1,k) + (U_MC_1(end,2,1,k) / U_MC_1(end,1,1,k)) * (gamma - 1) * (U_MC_1(end,3,1,k) - (U_MC_1(end,2,1,k))^2 / (2*U_MC_1(end,1,1,k)));

    p_MC_1(1,:,k) = p_L;
    p_MC_1(end,:,k) = p_R;

    U_MC_1_p(1, 1, :, k) = U_MC_1(1, 1, 1, k);
    U_MC_1_p(1, 2, :, k) = U_MC_1(1, 2, 1, k);
    U_MC_1_p(1, 3, :, k) = U_MC_1(1, 3, 1, k);
    U_MC_1_p(end, 1, :, k) = U_MC_1(end, 1, 1, k);
    U_MC_1_p(end, 2, :, k) = U_MC_1(end, 2, 1, k);
    U_MC_1_p(end, 3, :, k) = U_MC_1(end, 3, 1, k);

    u_plus_c_MC(1,1) = sqrt(gamma * p_L / rho_L);
    u_plus_c_MC(end,1) = sqrt(gamma * p_R / rho_R);

    for n = 2:t_steps(1)+1
    %Loop through all of time
        for i = 2:IL(1)-1
        %Loop through all of space for predictors
            %Calculate predictors for U
            U_MC_1_p(i, 1, n, k) = U_MC_1(i, 1, n-1, k) - delta_t_MC(k, 1) / delta_x_1 * ( (E_MC_1(i+1, 1, n-1, k) - epsilon * ((u_plus_c_MC(i) + u_plus_c_MC(i+1)) / 2) * (abs(p_MC_1(i+1,n-1,k) - 2*p_MC_1(i,n-1,k) + p_MC_1(i-1,n-1,k)) / (p_MC_1(i+1,n-1,k) + 2*p_MC_1(i,n-1,k) + p_MC_1(i-1,n-1,k))  ) * (U_MC_1(i+1, 1, n-1, k) - U_MC_1(i, 1, n-1, k)) ) - (E_MC_1(i, 1, n-1, k) - epsilon * ((u_plus_c_MC(i) + u_plus_c_MC(i-1)) / 2) * (abs(p_MC_1(i+1,n-1,k) - 2*p_MC_1(i,n-1,k) + p_MC_1(i-1,n-1,k)) / (p_MC_1(i+1,n-1,k) + 2*p_MC_1(i,n-1,k) + p_MC_1(i-1,n-1,k))) * (U_MC_1(i, 1, n-1, k) - U_MC_1(i-1, 1, n-1, k)) ) );
            U_MC_1_p(i, 2, n, k) = U_MC_1(i, 2, n-1, k) - delta_t_MC(k, 1) / delta_x_1 * ( (E_MC_1(i+1, 2, n-1, k) - epsilon * ((u_plus_c_MC(i) + u_plus_c_MC(i+1)) / 2) * (abs(p_MC_1(i+1,n-1,k) - 2*p_MC_1(i,n-1,k) + p_MC_1(i-1,n-1,k)) / (p_MC_1(i+1,n-1,k) + 2*p_MC_1(i,n-1,k) + p_MC_1(i-1,n-1,k))  ) * (U_MC_1(i+1, 2, n-1, k) - U_MC_1(i, 2, n-1, k)) ) - (E_MC_1(i, 2, n-1, k) - epsilon * ((u_plus_c_MC(i) + u_plus_c_MC(i-1)) / 2) * (abs(p_MC_1(i+1,n-1,k) - 2*p_MC_1(i,n-1,k) + p_MC_1(i-1,n-1,k)) / (p_MC_1(i+1,n-1,k) + 2*p_MC_1(i,n-1,k) + p_MC_1(i-1,n-1,k))) * (U_MC_1(i, 2, n-1, k) - U_MC_1(i-1, 2, n-1, k)) ) );
            U_MC_1_p(i, 3, n, k) = U_MC_1(i, 3, n-1, k) - delta_t_MC(k, 1) / delta_x_1 * ( (E_MC_1(i+1, 3, n-1, k) - epsilon * ((u_plus_c_MC(i) + u_plus_c_MC(i+1)) / 2) * (abs(p_MC_1(i+1,n-1,k) - 2*p_MC_1(i,n-1,k) + p_MC_1(i-1,n-1,k)) / (p_MC_1(i+1,n-1,k) + 2*p_MC_1(i,n-1,k) + p_MC_1(i-1,n-1,k))  ) * (U_MC_1(i+1, 3, n-1, k) - U_MC_1(i, 3, n-1, k)) ) - (E_MC_1(i, 3, n-1, k) - epsilon * ((u_plus_c_MC(i) + u_plus_c_MC(i-1)) / 2) * (abs(p_MC_1(i+1,n-1,k) - 2*p_MC_1(i,n-1,k) + p_MC_1(i-1,n-1,k)) / (p_MC_1(i+1,n-1,k) + 2*p_MC_1(i,n-1,k) + p_MC_1(i-1,n-1,k))) * (U_MC_1(i, 3, n-1, k) - U_MC_1(i-1, 3, n-1, k)) ) );

            %Calculate predictors for E based on predictor U
            E_MC_1_p(i, 1, n, k) = U_MC_1_p(i,2,n,k);
            E_MC_1_p(i, 2, n, k) = U_MC_1_p(i,2,n,k)^2 / U_MC_1_p(i,1,n,k) + (gamma - 1) * (U_MC_1_p(i,3,n,k) - (U_MC_1_p(i,2,n,k))^2 / (2*U_MC_1_p(i,1,n,k)));
            E_MC_1_p(i, 3, n, k) = U_MC_1_p(i,2,n,k) * U_MC_1_p(i,3,n,k) / U_MC_1_p(i,1,n,k) + (U_MC_1_p(i,2,n,k) / U_MC_1_p(i,1,n,k)) * (gamma - 1) * (U_MC_1_p(i,3,n,k) - (U_MC_1_p(i,2,n,k))^2 / (2*U_MC_1_p(i,1,n,k)));

            %Calculate predictor for p for use in corrector step
            p_MC_1_p(i,n,k) = (gamma - 1) * (U_MC_1_p(i,3,n,k) - (U_MC_1_p(i,2,n,k))^2 / (2*U_MC_1_p(i,1,n,k))  );          

            %Calculate predictor for |u|+c (this array is overwritten
            %everytime anyways)
            c_MC(i) = real(sqrt(gamma * ((gamma - 1) * (U_MC_1_p(i,3,n,k) - (U_MC_1_p(i,2,n,k))^2 / (2*U_MC_1_p(i,1,n,k)))) / (U_MC_1_p(i,1,n,k))));
            u_plus_c_MC(i,1) = abs(U_MC_1_p(i, 2, n, k) / U_MC_1_p(i, 1, n, k)) + c_MC(i);

        end

        for i = 2:IL(1)-1
        %Loop through all of space for correctors
            %Calculate correctors for U
            U_MC_1(i, 1, n, k) = 0.5 * (U_MC_1(i,1,n-1,k) + U_MC_1_p(i,1,n,k) - delta_t_MC(k,1)/delta_x_1*( (E_MC_1_p(i,1,n,k) - epsilon*((u_plus_c_MC(i)+u_plus_c_MC(i+1)) / 2)*(abs(p_MC_1_p(i+1,n,k)-2*p_MC_1_p(i,n,k)+p_MC_1_p(i-1,n,k)) / (p_MC_1_p(i+1,n,k)+2*p_MC_1_p(i,n,k)+p_MC_1_p(i-1,n,k)))*(U_MC_1_p(i+1,1,n,k) - U_MC_1_p(i,1,n,k))) - (E_MC_1_p(i-1,1,n,k) - epsilon*((u_plus_c_MC(i)+u_plus_c_MC(i-1)) / 2)*(abs(p_MC_1_p(i+1,n,k)-2*p_MC_1_p(i,n,k)+p_MC_1_p(i-1,n,k)) / (p_MC_1_p(i+1,n,k)+2*p_MC_1_p(i,n,k)+p_MC_1_p(i-1,n,k)))*(U_MC_1_p(i,1,n,k) - U_MC_1_p(i-1,1,n,k))) )  );
            U_MC_1(i, 2, n, k) = 0.5 * (U_MC_1(i,2,n-1,k) + U_MC_1_p(i,2,n,k) - delta_t_MC(k,1)/delta_x_1*( (E_MC_1_p(i,2,n,k) - epsilon*((u_plus_c_MC(i)+u_plus_c_MC(i+1)) / 2)*(abs(p_MC_1_p(i+1,n,k)-2*p_MC_1_p(i,n,k)+p_MC_1_p(i-1,n,k)) / (p_MC_1_p(i+1,n,k)+2*p_MC_1_p(i,n,k)+p_MC_1_p(i-1,n,k)))*(U_MC_1_p(i+1,2,n,k) - U_MC_1_p(i,2,n,k))) - (E_MC_1_p(i-1,2,n,k) - epsilon*((u_plus_c_MC(i)+u_plus_c_MC(i-1)) / 2)*(abs(p_MC_1_p(i+1,n,k)-2*p_MC_1_p(i,n,k)+p_MC_1_p(i-1,n,k)) / (p_MC_1_p(i+1,n,k)+2*p_MC_1_p(i,n,k)+p_MC_1_p(i-1,n,k)))*(U_MC_1_p(i,2,n,k) - U_MC_1_p(i-1,2,n,k))) )  );
            U_MC_1(i, 3, n, k) = 0.5 * (U_MC_1(i,3,n-1,k) + U_MC_1_p(i,3,n,k) - delta_t_MC(k,1)/delta_x_1*( (E_MC_1_p(i,3,n,k) - epsilon*((u_plus_c_MC(i)+u_plus_c_MC(i+1)) / 2)*(abs(p_MC_1_p(i+1,n,k)-2*p_MC_1_p(i,n,k)+p_MC_1_p(i-1,n,k)) / (p_MC_1_p(i+1,n,k)+2*p_MC_1_p(i,n,k)+p_MC_1_p(i-1,n,k)))*(U_MC_1_p(i+1,3,n,k) - U_MC_1_p(i,3,n,k))) - (E_MC_1_p(i-1,3,n,k) - epsilon*((u_plus_c_MC(i)+u_plus_c_MC(i-1)) / 2)*(abs(p_MC_1_p(i+1,n,k)-2*p_MC_1_p(i,n,k)+p_MC_1_p(i-1,n,k)) / (p_MC_1_p(i+1,n,k)+2*p_MC_1_p(i,n,k)+p_MC_1_p(i-1,n,k)))*(U_MC_1_p(i,3,n,k) - U_MC_1_p(i-1,3,n,k))) )  );
            
            if U_MC_1(i, 2, n, k) < 0
                U_MC_1(i, 1, n, k) = rho_L;
                U_MC_1(i, 2, n, k) = rho_L * u_L;
                U_MC_1(i, 3, n, k) = U_MC_1(i, 3, n-1, k);
            end
            

            %Calculate new flux vectors E based on corrected U
            E_MC_1(i, 1, n, k) = U_MC_1(i,2,n,k);
            E_MC_1(i, 2, n, k) = U_MC_1(i,2,n,k)^2 / U_MC_1(i,1,n,k) + (gamma - 1) * (U_MC_1(i,3,n,k) - (U_MC_1(i,2,n,k))^2 / (2*U_MC_1(i,1,n,k)));
            E_MC_1(i, 3, n, k) = U_MC_1(i,2,n,k) * U_MC_1(i,3,n,k) / U_MC_1(i,1,n,k) + (U_MC_1(i,2,n,k) / U_MC_1(i,1,n,k)) * (gamma - 1) * (U_MC_1(i,3,n,k) - (U_MC_1(i,2,n,k))^2 / (2*U_MC_1(i,1,n,k)));

            %Calculate new speed of sound at each point based on correctors
            c_MC(i) = real(sqrt(gamma * ((gamma - 1) * (U_MC_1(i,3,n,k) - (U_MC_1(i,2,n,k))^2 / (2*U_MC_1(i,1,n,k)))) / (U_MC_1(i,1,n,k))));
            u_plus_c_MC(i,1) = abs(U_MC_1(i, 2, n, k) / U_MC_1(i, 1, n, k)) + c_MC(i);

            %Calculate new pressure
            p_MC_1(i,n,k) = ((gamma - 1) * (U_MC_1(i,3,n,k) - (U_MC_1(i,2,n,k))^2 / (2*U_MC_1(i,1,n,k))));

        end
        %Set new time step based on max(|u|+c) and calculate total time
        delta_t_MC(k,1) = nu(k) * delta_x_1 / max(u_plus_c_MC);
        t_total_MC_1(k, 1) = t_total_MC_1(k,1) + delta_t_MC(k,1);

    end

end


%Grid 2
for k = 1:length(nu)
%Loop through all Courant numbers
    for i = 1:IL(2) 
    %Set initial values for all of space
        if x_arr_2(i) <= 0
            %Left half
            U_MC_2(i,1,1,k) = rho_L;
            U_MC_2(i,2,1,k) = rho_L * u_L;
            U_MC_2(i,3,1,k) = p_L / (gamma - 1) + (rho_L / 2) * u_L^2;

            p_MC_2(i,1,k) = p_L;
            
            E_MC_2(i,1,1,k) = U_MC_2(i,2,1,k);
            E_MC_2(i,2,1,k) = U_MC_2(i,2,1,k)^2 / U_MC_2(i,1,1,k) + (gamma - 1) * (U_MC_2(i,3,1,k) - (U_MC_2(i,2,1,k))^2 / (2*U_MC_2(i,1,1,k)));
            E_MC_2(i,3,1,k) = U_MC_2(i,2,1,k) * U_MC_2(i,3,1,k) / U_MC_2(i,1,1,k) + (U_MC_2(i,2,1,k) / U_MC_2(i,1,1,k)) * (gamma - 1) * (U_MC_2(i,3,1,k) - (U_MC_2(i,2,1,k))^2 / (2*U_MC_2(i,1,1,k)));
            
            c_MC_2(i) = sqrt(gamma * ((gamma - 1) * (U_MC_2(i,3,1,k) - (U_MC_2(i,2,1,k))^2 / (2*U_MC_2(i,1,1,k)))) / (U_MC_2(i,1,1,k)));
            u_plus_c_MC_2(i,1) = abs(U_MC_2(i, 2, 1, k) / U_MC_2(i, 1, 1, k)) + c_MC_2(i);

        else
            %Right half
            U_MC_2(i,1,1,k) = rho_R;
            U_MC_2(i,2,1,k) = rho_R * u_R;
            U_MC_2(i,3,1,k) = p_R / (gamma - 1) + (rho_R / 2) * u_R^2;

            p_MC_2(i,1,k) = p_R;

            E_MC_2(i,1,1,k) = U_MC_2(i,2,1,k);
            E_MC_2(i,2,1,k) = U_MC_2(i,2,1,k)^2 / U_MC_2(i,1,1,k) + (gamma - 1) * (U_MC_2(i,3,1,k) - (U_MC_2(i,2,1,k))^2 / (2*U_MC_2(i,1,1,k)));
            E_MC_2(i,3,1,k) = U_MC_2(i,2,1,k) * U_MC_2(i,3,1,k) / U_MC_2(i,1,1,k) + (U_MC_2(i,2,1,k) / U_MC_2(i,1,1,k)) * (gamma - 1) * (U_MC_2(i,3,1,k) - (U_MC_2(i,2,1,k))^2 / (2*U_MC_2(i,1,1,k)));
            
            c_MC_2(i) = sqrt(gamma * ((gamma - 1) * (U_MC_2(i,3,1,k) - (U_MC_2(i,2,1,k))^2 / (2*U_MC_2(i,1,1,k)))) / (U_MC_2(i,1,1,k)));
            u_plus_c_MC_2(i,1) = abs(U_MC_2(i, 2, 1, k) / U_MC_2(i, 1, 1, k)) + c_MC_2(i);

        end

    end

    delta_t_MC_2(k, 1) = nu(k) * delta_x_2 / max(u_plus_c_MC_2);
    t_total_MC_2(k, 1) = delta_t_MC_2(k,1);

    %Set Dirichlet BCs for all of time
    U_MC_2(1, 1, :, k) = rho_L;
    U_MC_2(1, 2, :, k) = rho_L * u_L;
    U_MC_2(1, 3, :, k) = p_L / (gamma - 1) + (rho_L / 2) * u_L^2;
    U_MC_2(end, 1, :, k) = rho_R;
    U_MC_2(end, 2, :, k) = rho_R * u_R;
    U_MC_2(end, 3, :, k) = p_R / (gamma - 1) + (rho_R / 2) * u_R^2;

    E_MC_2(1, 1, :, k) = U_MC_2(1, 2, 1, k);
    E_MC_2(1, 2, :, k) = U_MC_2(1,2,1,k)^2 / U_MC_2(1,1,1,k) + (gamma - 1) * (U_MC_2(1,3,1,k) - (U_MC_2(1,2,1,k))^2 / (2*U_MC_2(1,1,1,k)));
    E_MC_2(1, 3, :, k) = U_MC_2(1,2,1,k) * U_MC_2(1,3,1,k) / U_MC_2(1,1,1,k) + (U_MC_2(1,2,1,k) / U_MC_2(1,1,1,k)) * (gamma - 1) * (U_MC_2(1,3,1,k) - (U_MC_2(1,2,1,k))^2 / (2*U_MC_2(1,1,1,k)));
    E_MC_2(end, 1, :, k) = U_MC_2(end, 2, 1, k);
    E_MC_2(end, 2, :, k) = U_MC_2(end,2,1,k)^2 / U_MC_2(end,1,1,k) + (gamma - 1) * (U_MC_2(end,3,1,k) - (U_MC_2(end,2,1,k))^2 / (2*U_MC_2(end,1,1,k)));
    E_MC_2(end, 3, :, k) = U_MC_2(end,2,1,k) * U_MC_2(end,3,1,k) / U_MC_2(end,1,1,k) + (U_MC_2(end,2,1,k) / U_MC_2(end,1,1,k)) * (gamma - 1) * (U_MC_2(end,3,1,k) - (U_MC_2(end,2,1,k))^2 / (2*U_MC_2(end,1,1,k)));

    p_MC_2(1,:,k) = p_L;
    p_MC_2(end,:,k) = p_R;

    U_MC_2_p(1, 1, :, k) = U_MC_2(1, 1, 1, k);
    U_MC_2_p(1, 2, :, k) = U_MC_2(1, 2, 1, k);
    U_MC_2_p(1, 3, :, k) = U_MC_2(1, 3, 1, k);
    U_MC_2_p(end, 1, :, k) = U_MC_2(end, 1, 1, k);
    U_MC_2_p(end, 2, :, k) = U_MC_2(end, 2, 1, k);
    U_MC_2_p(end, 3, :, k) = U_MC_2(end, 3, 1, k);

    u_plus_c_MC_2(1,1) = sqrt(gamma * p_L / rho_L);
    u_plus_c_MC_2(1,end) = sqrt(gamma * p_R / rho_R);

    for n = 2:t_steps(2)+1
    %Loop through all of time
        for i = 2:IL(2)-1
        %Loop through all of space for predictors
            %Calculate predictors for U
            U_MC_2_p(i, 1, n, k) = U_MC_2(i, 1, n-1, k) - delta_t_MC_2(k, 1) / delta_x_2 * ( (E_MC_2(i+1, 1, n-1, k) - epsilon * ((u_plus_c_MC_2(i) + u_plus_c_MC_2(i+1)) / 2) * (abs(p_MC_2(i+1,n-1,k) - 2*p_MC_2(i,n-1,k) + p_MC_2(i-1,n-1,k)) / (p_MC_2(i+1,n-1,k) + 2*p_MC_2(i,n-1,k) + p_MC_2(i-1,n-1,k))  ) * (U_MC_2(i+1, 1, n-1, k) - U_MC_2(i, 1, n-1, k)) ) - (E_MC_2(i, 1, n-1, k) - epsilon * ((u_plus_c_MC_2(i) + u_plus_c_MC_2(i-1)) / 2) * (abs(p_MC_2(i+1,n-1,k) - 2*p_MC_2(i,n-1,k) + p_MC_2(i-1,n-1,k)) / (p_MC_2(i+1,n-1,k) + 2*p_MC_2(i,n-1,k) + p_MC_2(i-1,n-1,k))) * (U_MC_2(i, 1, n-1, k) - U_MC_2(i-1, 1, n-1, k)) ) );
            U_MC_2_p(i, 2, n, k) = U_MC_2(i, 2, n-1, k) - delta_t_MC_2(k, 1) / delta_x_2 * ( (E_MC_2(i+1, 2, n-1, k) - epsilon * ((u_plus_c_MC_2(i) + u_plus_c_MC_2(i+1)) / 2) * (abs(p_MC_2(i+1,n-1,k) - 2*p_MC_2(i,n-1,k) + p_MC_2(i-1,n-1,k)) / (p_MC_2(i+1,n-1,k) + 2*p_MC_2(i,n-1,k) + p_MC_2(i-1,n-1,k))  ) * (U_MC_2(i+1, 2, n-1, k) - U_MC_2(i, 2, n-1, k)) ) - (E_MC_2(i, 2, n-1, k) - epsilon * ((u_plus_c_MC_2(i) + u_plus_c_MC_2(i-1)) / 2) * (abs(p_MC_2(i+1,n-1,k) - 2*p_MC_2(i,n-1,k) + p_MC_2(i-1,n-1,k)) / (p_MC_2(i+1,n-1,k) + 2*p_MC_2(i,n-1,k) + p_MC_2(i-1,n-1,k))) * (U_MC_2(i, 2, n-1, k) - U_MC_2(i-1, 2, n-1, k)) ) );
            U_MC_2_p(i, 3, n, k) = U_MC_2(i, 3, n-1, k) - delta_t_MC_2(k, 1) / delta_x_2 * ( (E_MC_2(i+1, 3, n-1, k) - epsilon * ((u_plus_c_MC_2(i) + u_plus_c_MC_2(i+1)) / 2) * (abs(p_MC_2(i+1,n-1,k) - 2*p_MC_2(i,n-1,k) + p_MC_2(i-1,n-1,k)) / (p_MC_2(i+1,n-1,k) + 2*p_MC_2(i,n-1,k) + p_MC_2(i-1,n-1,k))  ) * (U_MC_2(i+1, 3, n-1, k) - U_MC_2(i, 3, n-1, k)) ) - (E_MC_2(i, 3, n-1, k) - epsilon * ((u_plus_c_MC_2(i) + u_plus_c_MC_2(i-1)) / 2) * (abs(p_MC_2(i+1,n-1,k) - 2*p_MC_2(i,n-1,k) + p_MC_2(i-1,n-1,k)) / (p_MC_2(i+1,n-1,k) + 2*p_MC_2(i,n-1,k) + p_MC_2(i-1,n-1,k))) * (U_MC_2(i, 3, n-1, k) - U_MC_2(i-1, 3, n-1, k)) ) );

            %Calculate predictors for E based on predictor U
            E_MC_2_p(i, 1, n, k) = U_MC_2_p(i,2,n,k);
            E_MC_2_p(i, 2, n, k) = U_MC_2_p(i,2,n,k)^2 / U_MC_2_p(i,1,n,k) + (gamma - 1) * (U_MC_2_p(i,3,n,k) - (U_MC_2_p(i,2,n,k))^2 / (2*U_MC_2_p(i,1,n,k)));
            E_MC_2_p(i, 3, n, k) = U_MC_2_p(i,2,n,k) * U_MC_2_p(i,3,n,k) / U_MC_2_p(i,1,n,k) + (U_MC_2_p(i,2,n,k) / U_MC_2_p(i,1,n,k)) * (gamma - 1) * (U_MC_2_p(i,3,n,k) - (U_MC_2_p(i,2,n,k))^2 / (2*U_MC_2_p(i,1,n,k)));

            %Calculate predictor for p for use in corrector step
            p_MC_2_p(i,n,k) = (gamma - 1) * (U_MC_2_p(i,3,n,k) - (U_MC_2_p(i,2,n,k))^2 / (2*U_MC_2_p(i,1,n,k))  );          

            %Calculate predictor for |u|+c (this array is overwritten
            %everytime anyways)
            c_MC_2(i) = real(sqrt(gamma * ((gamma - 1) * (U_MC_2_p(i,3,n,k) - (U_MC_2_p(i,2,n,k))^2 / (2*U_MC_2_p(i,1,n,k)))) / (U_MC_2_p(i,1,n,k))));
            u_plus_c_MC_2(i,1) = abs(U_MC_2_p(i, 2, n, k) / U_MC_2_p(i, 1, n, k)) + c_MC_2(i);

        end

        for i = 2:IL(2)-1
        %Loop through all of space for correctors
            %Calculate correctors for U
            U_MC_2(i, 1, n, k) = 0.5 * (U_MC_2(i,1,n-1,k) + U_MC_2_p(i,1,n,k) - delta_t_MC_2(k,1)/delta_x_2*( (E_MC_2_p(i,1,n,k) - epsilon*((u_plus_c_MC_2(i)+u_plus_c_MC_2(i+1)) / 2)*(abs(p_MC_2_p(i+1,n,k)-2*p_MC_2_p(i,n,k)+p_MC_2_p(i-1,n,k)) / (p_MC_2_p(i+1,n,k)+2*p_MC_2_p(i,n,k)+p_MC_2_p(i-1,n,k)))*(U_MC_2_p(i+1,1,n,k) - U_MC_2_p(i,1,n,k))) - (E_MC_2_p(i-1,1,n,k) - epsilon*((u_plus_c_MC_2(i)+u_plus_c_MC_2(i-1)) / 2)*(abs(p_MC_2_p(i+1,n,k)-2*p_MC_2_p(i,n,k)+p_MC_2_p(i-1,n,k)) / (p_MC_2_p(i+1,n,k)+2*p_MC_2_p(i,n,k)+p_MC_2_p(i-1,n,k)))*(U_MC_2_p(i,1,n,k) - U_MC_2_p(i-1,1,n,k))) )  );
            U_MC_2(i, 2, n, k) = 0.5 * (U_MC_2(i,2,n-1,k) + U_MC_2_p(i,2,n,k) - delta_t_MC_2(k,1)/delta_x_2*( (E_MC_2_p(i,2,n,k) - epsilon*((u_plus_c_MC_2(i)+u_plus_c_MC_2(i+1)) / 2)*(abs(p_MC_2_p(i+1,n,k)-2*p_MC_2_p(i,n,k)+p_MC_2_p(i-1,n,k)) / (p_MC_2_p(i+1,n,k)+2*p_MC_2_p(i,n,k)+p_MC_2_p(i-1,n,k)))*(U_MC_2_p(i+1,2,n,k) - U_MC_2_p(i,2,n,k))) - (E_MC_2_p(i-1,2,n,k) - epsilon*((u_plus_c_MC_2(i)+u_plus_c_MC_2(i-1)) / 2)*(abs(p_MC_2_p(i+1,n,k)-2*p_MC_2_p(i,n,k)+p_MC_2_p(i-1,n,k)) / (p_MC_2_p(i+1,n,k)+2*p_MC_2_p(i,n,k)+p_MC_2_p(i-1,n,k)))*(U_MC_2_p(i,2,n,k) - U_MC_2_p(i-1,2,n,k))) )  );
            U_MC_2(i, 3, n, k) = 0.5 * (U_MC_2(i,3,n-1,k) + U_MC_2_p(i,3,n,k) - delta_t_MC_2(k,1)/delta_x_2*( (E_MC_2_p(i,3,n,k) - epsilon*((u_plus_c_MC_2(i)+u_plus_c_MC_2(i+1)) / 2)*(abs(p_MC_2_p(i+1,n,k)-2*p_MC_2_p(i,n,k)+p_MC_2_p(i-1,n,k)) / (p_MC_2_p(i+1,n,k)+2*p_MC_2_p(i,n,k)+p_MC_2_p(i-1,n,k)))*(U_MC_2_p(i+1,3,n,k) - U_MC_2_p(i,3,n,k))) - (E_MC_2_p(i-1,3,n,k) - epsilon*((u_plus_c_MC_2(i)+u_plus_c_MC_2(i-1)) / 2)*(abs(p_MC_2_p(i+1,n,k)-2*p_MC_2_p(i,n,k)+p_MC_2_p(i-1,n,k)) / (p_MC_2_p(i+1,n,k)+2*p_MC_2_p(i,n,k)+p_MC_2_p(i-1,n,k)))*(U_MC_2_p(i,3,n,k) - U_MC_2_p(i-1,3,n,k))) )  );
            
            if U_MC_2(i, 2, n, k) < 0
                U_MC_2(i, 1, n, k) = rho_L;
                U_MC_2(i, 2, n, k) = rho_L * u_L;
                U_MC_2(i, 3, n, k) = U_MC_2(i, 3, n-1, k);
            end
            

            %Calculate new flux vectors E based on corrected U
            E_MC_2(i, 1, n, k) = U_MC_2(i,2,n,k);
            E_MC_2(i, 2, n, k) = U_MC_2(i,2,n,k)^2 / U_MC_2(i,1,n,k) + (gamma - 1) * (U_MC_2(i,3,n,k) - (U_MC_2(i,2,n,k))^2 / (2*U_MC_2(i,1,n,k)));
            E_MC_2(i, 3, n, k) = U_MC_2(i,2,n,k) * U_MC_2(i,3,n,k) / U_MC_2(i,1,n,k) + (U_MC_2(i,2,n,k) / U_MC_2(i,1,n,k)) * (gamma - 1) * (U_MC_2(i,3,n,k) - (U_MC_2(i,2,n,k))^2 / (2*U_MC_2(i,1,n,k)));

            %Calculate new speed of sound at each point based on correctors
            c_MC_2(i) = real(sqrt(gamma * ((gamma - 1) * (U_MC_2(i,3,n,k) - (U_MC_2(i,2,n,k))^2 / (2*U_MC_2(i,1,n,k)))) / (U_MC_2(i,1,n,k))));
            u_plus_c_MC_2(i,1) = abs(U_MC_2(i, 2, n, k) / U_MC_2(i, 1, n, k)) + c_MC_2(i);

            %Calculate new pressure
            p_MC_2(i,n,k) = ((gamma - 1) * (U_MC_2(i,3,n,k) - (U_MC_2(i,2,n,k))^2 / (2*U_MC_2(i,1,n,k))));

        end
        %Set new time step based on max(|u|+c) and calculate total time
        delta_t_MC_2(k,1) = nu(k) * delta_x_2 / max(u_plus_c_MC_2);
        t_total_MC_2(k, 1) = t_total_MC_2(k,1) + delta_t_MC_2(k,1);

    end

end


%Analytical calculation for MacCormack
c_L = sqrt(gamma * p_L / rho_L);
x_L = (u_L - c_L) .* t_total_MC_1;

c_R = sqrt(gamma * p_R / rho_R);

p_2 = p_R * P;
rho_2 = rho_R * (1 + alpha * P) / (alpha + P);
u_S = c_R * sqrt((gamma - 1) / (2 * gamma) * (1 + alpha * P)) + u_R;
u_2 = (P - 1) / (sqrt(1 + alpha * P))*(1 / sqrt(gamma * (gamma - 1) / 2)) * c_R + u_R;

p_3 = p_2;
u_3 = u_2;
rho_3 = (p_3 / p_L)^(1 / gamma) * rho_L;
c_3 = sqrt(gamma * p_3 / rho_3);

%analytical for MC grid 1
for k = 1:length(nu)
    %Note: U_anlyt stores rho, u, and p (not the normal U parameters of
    %rho, rho*u, and e)
    for i = 1:IL(1)
        if x_arr_1(i) >= u_S * t_total_MC_1(k)
            U_MC_1_anlyt(i, 1, k) = rho_R;
            U_MC_1_anlyt(i, 2, k) = u_R;
            U_MC_1_anlyt(i, 3, k) = p_R;
            
        elseif (x_arr_1(i) > u_2 * t_total_MC_1(k)) && (x_arr_1(i) < u_S * t_total_MC_1(k)) 
            U_MC_1_anlyt(i, 1, k) = rho_2;
            U_MC_1_anlyt(i, 2, k) = u_2;
            U_MC_1_anlyt(i, 3, k) = p_2;

        elseif (x_arr_1(i) > (u_3 - c_3) * t_total_MC_1(k)) && (x_arr_1(i) < u_3 * t_total_MC_1(k))
            U_MC_1_anlyt(i, 1, k) = rho_3;
            U_MC_1_anlyt(i, 2, k) = u_3;
            U_MC_1_anlyt(i, 3, k) = p_3;

        elseif (x_arr_1(i) > (u_L - c_L) * t_total_MC_1(k)) && (x_arr_1(i) < (u_3 - c_3) * t_total_MC_1(k))
            u_5 = 2 / (gamma + 1) * ( (x_arr_1(i) / t_total_MC_1(k)) + c_L + (gamma - 1) / 2 * u_L);
            c_5 = u_5 - x_arr_1(i) / t_total_MC_1(k);
            p_5 = p_L * (c_5 / c_L)^((2 * gamma) / (gamma - 1));
            rho_5 = (p_5 / p_L)^(1 / gamma) * rho_L;
            U_MC_1_anlyt(i, 1, k) = rho_5;
            U_MC_1_anlyt(i, 2, k) = u_5;
            U_MC_1_anlyt(i, 3, k) = p_5;

        else
            U_MC_1_anlyt(i, 1, k) = rho_L;
            U_MC_1_anlyt(i, 2, k) = u_L;
            U_MC_1_anlyt(i, 3, k) = p_L;
        end
    end
end

%analytical for MC grid 2
for k = 1:length(nu)
    %Note: U_anlyt stores rho, u, and p (not the normal U parameters of
    %rho, rho*u, and e)
    for i = 1:IL(2)
        if x_arr_2(i) >= u_S * t_total_MC_2(k)
            U_MC_2_anlyt(i, 1, k) = rho_R;
            U_MC_2_anlyt(i, 2, k) = u_R;
            U_MC_2_anlyt(i, 3, k) = p_R;
            
        elseif (x_arr_2(i) > u_2 * t_total_MC_2(k)) && (x_arr_2(i) < u_S * t_total_MC_2(k)) 
            U_MC_2_anlyt(i, 1, k) = rho_2;
            U_MC_2_anlyt(i, 2, k) = u_2;
            U_MC_2_anlyt(i, 3, k) = p_2;

        elseif (x_arr_2(i) > (u_3 - c_3) * t_total_MC_2(k)) && (x_arr_2(i) < u_3 * t_total_MC_2(k))
            U_MC_2_anlyt(i, 1, k) = rho_3;
            U_MC_2_anlyt(i, 2, k) = u_3;
            U_MC_2_anlyt(i, 3, k) = p_3;

        elseif (x_arr_2(i) > (u_L - c_L) * t_total_MC_2(k)) && (x_arr_2(i) < (u_3 - c_3) * t_total_MC_2(k))
            u_5 = 2 / (gamma + 1) * ( (x_arr_2(i) / t_total_MC_2(k)) + c_L + (gamma - 1) / 2 * u_L);
            c_5 = u_5 - x_arr_2(i) / t_total_MC_2(k);
            p_5 = p_L * (c_5 / c_L)^((2 * gamma) / (gamma - 1));
            rho_5 = (p_5 / p_L)^(1 / gamma) * rho_L;
            U_MC_2_anlyt(i, 1, k) = rho_5;
            U_MC_2_anlyt(i, 2, k) = u_5;
            U_MC_2_anlyt(i, 3, k) = p_5;

        else
            U_MC_2_anlyt(i, 1, k) = rho_L;
            U_MC_2_anlyt(i, 2, k) = u_L;
            U_MC_2_anlyt(i, 3, k) = p_L;
        end
    end
end


%Grid 1 Lax-Wendroff with Artificial Dissipation
for k = 1:length(nu)
    %Loop through all of time    
        for i = 1:IL(1)
        %Set initial conditions
            if x_arr_1(i) <= 0
            %Left half
                U_LW_1(i,1,1,k) = rho_L;
                U_LW_1(i,2,1,k) = rho_L * u_L;
                U_LW_1(i,3,1,k) = p_L / (gamma - 1) + (rho_L / 2) * u_L^2;
    
                E_LW_1(i,1,1,k) = U_LW_1(i,2,1,k);
                E_LW_1(i,2,1,k) = (U_LW_1(i,2,1,k))^2 / U_LW_1(i,1,1,k) + (gamma - 1) * (U_LW_1(i,3,1,k) - (U_LW_1(i,2,1,k))^2 / (2*U_LW_1(i,1,1,k)));
                E_LW_1(i,3,1,k) = (U_LW_1(i,2,1,k)*U_LW_1(i,3,1,k)) / U_LW_1(i,1,1,k) + (U_LW_1(i,2,1,k) / U_LW_1(i,1,1,k)) * ((gamma - 1) * (U_LW_1(i,3,1,k) - (U_LW_1(i,2,1,k))^2 / (2*U_LW_1(i,1,1,k))) );
    
                p_LW_1(i,1,k) = p_L;
    
                c_LW_1(i,1) = sqrt(gamma * ((gamma - 1) * (U_LW_1(i,3,1,k) - (U_LW_1(i,2,1,k))^2 / (2*U_LW_1(i,1,1,k)))) / (U_LW_1(i,1,1,k)));
                u_plus_c_LW(i,1) = abs(U_LW_1(i, 2, 1, k) / U_LW_1(i, 1, 1, k)) + c_LW_1(i);
    
            else
            %Right half
                U_LW_1(i,1,1,k) = rho_R;
                U_LW_1(i,2,1,k) = rho_R * u_R;
                U_LW_1(i,3,1,k) = p_R / (gamma - 1) + (rho_R / 2) * u_R^2;
    
                E_LW_1(i,1,1,k) = U_LW_1(i,2,1,k);
                E_LW_1(i,2,1,k) = (U_LW_1(i,2,1,k))^2 / U_LW_1(i,1,1,k) + (gamma - 1) * (U_LW_1(i,3,1,k) - (U_LW_1(i,2,1,k))^2 / (2*U_LW_1(i,1,1,k)));
                E_LW_1(i,3,1,k) = (U_LW_1(i,2,1,k)*U_LW_1(i,3,1,k)) / U_LW_1(i,1,1,k) + (U_LW_1(i,2,1,k) / U_LW_1(i,1,1,k)) * ((gamma - 1) * (U_LW_1(i,3,1,k) - (U_LW_1(i,2,1,k))^2 / (2*U_LW_1(i,1,1,k))) );
    
                p_LW_1(i,1,k) = p_R;
    
                c_LW_1(i,1) = sqrt(gamma * ((gamma - 1) * (U_LW_1(i,3,1,k) - (U_LW_1(i,2,1,k))^2 / (2*U_LW_1(i,1,1,k)))) / (U_LW_1(i,1,1,k)));
                u_plus_c_LW(i,1) = abs(U_LW_1(i, 2, 1, k) / U_LW_1(i, 1, 1, k)) + c_LW_1(i);
    
            end
        end  
        %set initial time step
        delta_t_LW(k, 1) = nu(k) * delta_x_1 / max(u_plus_c_LW);
        t_total_LW_1(k, 1) = delta_t_LW(k,1);
    
        %Set Dirichlet BCs for all of time
        U_LW_1(1,1,:,k) = rho_L;
        U_LW_1(1,2,:,k) = rho_L * u_L;
        U_LW_1(1,3,:,k) = p_L / (gamma - 1) + (rho_L / 2) * u_L^2;
        U_LW_1(end,1,:,k) = rho_R;
        U_LW_1(end,2,:,k) = rho_R * u_R;
        U_LW_1(end,3,:,k) = p_R / (gamma - 1) + (rho_R / 2) * u_R^2;
    
        E_LW_1(1,1,:,k) = U_LW_1(1,2,1,k);
        E_LW_1(1,2,:,k) = (U_LW_1(1,2,1,k))^2 / U_LW_1(1,1,1,k) + (gamma - 1) * (U_LW_1(1,3,1,k) - (U_LW_1(1,2,1,k))^2 / (2*U_LW_1(1,1,1,k)));
        E_LW_1(1,3,:,k) = (U_LW_1(1,2,1,k)*U_LW_1(1,3,1,k)) / U_LW_1(1,1,1,k) + (U_LW_1(1,2,1,k) / U_LW_1(1,1,1,k)) * ((gamma - 1) * (U_LW_1(1,3,1,k) - (U_LW_1(1,2,1,k))^2 / (2*U_LW_1(1,1,1,k))) );
        E_LW_1(end,1,:,k) = U_LW_1(end,2,1,k);
        E_LW_1(end,2,:,k) = (U_LW_1(end,2,1,k))^2 / U_LW_1(end,1,1,k) + (gamma - 1) * (U_LW_1(end,3,1,k) - (U_LW_1(end,2,1,k))^2 / (2*U_LW_1(end,1,1,k)));
        E_LW_1(end,3,:,k) = (U_LW_1(end,2,1,k)*U_LW_1(end,3,1,k)) / U_LW_1(end,1,1,k) + (U_LW_1(end,2,1,k) / U_LW_1(end,1,1,k)) * ((gamma - 1) * (U_LW_1(end,3,1,k) - (U_LW_1(end,2,1,k))^2 / (2*U_LW_1(end,1,1,k))) );
        
        p_LW_1(1,:,k) = p_L;
        p_LW_1(end,:,k) = p_R;

        c_LW_1(1,1) = sqrt(gamma * p_L / rho_L);
        u_plus_c_LW(1,1) = abs(U_LW_1(1, 2, 1, k) / U_LW_1(1, 1, 1, k)) + c_LW_1(1);
        c_LW_1(end,1) = sqrt(gamma * p_R / rho_R);
        u_plus_c_LW(end,1) = abs(U_LW_1(end, 2, 1, k) / U_LW_1(end, 1, 1, k)) + c_LW_1(end);
        for n = 2:t_steps(1)+1
            for i = 2:IL(1)-1
            %Perform iterative procedure for all of interior grid points
                %assign jacobians at i = 1 + 1/2 and i = 1 - 1/2
                A_i_plus_half_1 = [0, 1, 0];
                A_i_plus_half_2 = [-(3- gamma) * ((U_LW_1(i,2,n-1,k) + U_LW_1(i+1,2,n-1,k))/2)^2 / (2 * ((U_LW_1(i,1,n-1,k) + U_LW_1(i+1,1,n-1,k)) / 2)^2), (3 - gamma) * ((U_LW_1(i,2,n-1,k) + U_LW_1(i+1,2,n-1,k)) / 2) / ((U_LW_1(i,1,n-1,k) + U_LW_1(i+1,1,n-1,k)) / 2), (gamma - 1)];
                A_i_plus_half_3 = [(gamma - 1) * (((U_LW_1(i,2,n-1,k) + U_LW_1(i+1,2,n-1,k)) / 2)^3 / ((U_LW_1(i,1,n-1,k) + U_LW_1(i+1,1,n-1,k)) / 2)^3) - gamma * ( (U_LW_1(i,2,n-1,k) + U_LW_1(i+1,2,n-1,k)) / 2) / ((U_LW_1(i,1,n-1,k) + U_LW_1(i+1,1,n-1,k)) / 2)^2 * ((U_LW_1(i,3,n-1,k) + U_LW_1(i+1,3,n-1,k)) / 2), gamma * ((U_LW_1(i,3,n-1,k) + U_LW_1(i+1,3,n-1,k)) / 2) / ((U_LW_1(i,1,n-1,k) + U_LW_1(i+1,1,n-1,k)) / 2) - 1.5 * (gamma - 1) * ((U_LW_1(i,2,n-1,k) + U_LW_1(i+1,2,n-1,k)) / 2)^2 / ((U_LW_1(i,1,n-1,k) + U_LW_1(i+1,1,n-1,k)) / 2)^2, gamma * ((U_LW_1(i,2,n-1,k) + U_LW_1(i+1,2,n-1,k)) / 2) / ((U_LW_1(i,1,n-1,k) + U_LW_1(i+1,1,n-1,k)) / 2)];

                A_i_minus_half_1 = [0, 1, 0];
                A_i_minus_half_2 = [-(3- gamma) * ((U_LW_1(i,2,n-1,k) + U_LW_1(i-1,2,n-1,k))/2)^2 / (2 * ((U_LW_1(i,1,n-1,k) + U_LW_1(i-1,1,n-1,k)) / 2)^2), (3 - gamma) * ((U_LW_1(i,2,n-1,k) + U_LW_1(i-1,2,n-1,k)) / 2) / ((U_LW_1(i,1,n-1,k) + U_LW_1(i-1,1,n-1,k)) / 2), (gamma - 1)];
                A_i_minus_half_3 = [(gamma - 1) * (((U_LW_1(i,2,n-1,k) + U_LW_1(i-1,2,n-1,k)) / 2)^3 / ((U_LW_1(i,1,n-1,k) + U_LW_1(i-1,1,n-1,k)) / 2)^3) - gamma * ( (U_LW_1(i,2,n-1,k) + U_LW_1(i-1,2,n-1,k)) / 2) / ((U_LW_1(i,1,n-1,k) + U_LW_1(i-1,1,n-1,k)) / 2)^2 * ((U_LW_1(i,3,n-1,k) + U_LW_1(i-1,3,n-1,k)) / 2), gamma * ((U_LW_1(i,3,n-1,k) + U_LW_1(i-1,3,n-1,k)) / 2) / ((U_LW_1(i,1,n-1,k) + U_LW_1(i-1,1,n-1,k)) / 2) - 1.5 * (gamma - 1) * ((U_LW_1(i,2,n-1,k) + U_LW_1(i-1,2,n-1,k)) / 2)^2 / ((U_LW_1(i,1,n-1,k) + U_LW_1(i-1,1,n-1,k)) / 2)^2, gamma * ((U_LW_1(i,2,n-1,k) + U_LW_1(i-1,2,n-1,k)) / 2) / ((U_LW_1(i,1,n-1,k) + U_LW_1(i-1,1,n-1,k)) / 2)];
                
                A_i_plus_half = [A_i_plus_half_1; A_i_plus_half_2; A_i_plus_half_3];
                A_i_minus_half = [A_i_minus_half_1; A_i_minus_half_2; A_i_minus_half_3];

                %calculate new U
                U_LW_1(i,:,n,k) = (U_LW_1(i,:,n-1,k).' - delta_t_LW(k,1) / (2 * delta_x_1) .* ((E_LW_1(i+1,:,n-1,k).' - epsilon .* (U_LW_1(i+1,:,n-1,k).' - U_LW_1(i,:,n-1,k).') .* ((u_plus_c_LW(i) + u_plus_c_LW(i+1)) / 2) .* (abs(p_LW_1(i+1,n-1,k) - 2*p_LW_1(i,n-1,k) + p_LW_1(i-1,n-1,k)) / (p_LW_1(i+1,n-1,k) + 2*p_LW_1(i,n-1,k) + p_LW_1(i-1,n-1,k)))) - (E_LW_1(i-1,:,n-1,k).' - epsilon .* (U_LW_1(i,:,n-1,k).' - U_LW_1(i-1,:,n-1,k).') .* ((u_plus_c_LW(i) + u_plus_c_LW(i-1)) / 2) .* (abs(p_LW_1(i+1,n-1,k) - 2*p_LW_1(i,n-1,k) + p_LW_1(i-1,n-1,k)) / (p_LW_1(i+1,n-1,k) + 2*p_LW_1(i,n-1,k) + p_LW_1(i-1,n-1,k))))) + delta_t_LW(k,1)^2 / (2 * delta_x_1^2) .* (A_i_plus_half * ( (E_LW_1(i+1,:,n-1,k).' - epsilon .* (U_LW_1(i+1,:,n-1,k).' - U_LW_1(i,:,n-1,k).') .* ((u_plus_c_LW(i) + u_plus_c_LW(i+1)) / 2) .* (abs(p_LW_1(i+1,n-1,k) - 2*p_LW_1(i,n-1,k) + p_LW_1(i-1,n-1,k)) / (p_LW_1(i+1,n-1,k) + 2*p_LW_1(i,n-1,k) + p_LW_1(i-1,n-1,k)))) - (E_LW_1(i,:,n-1,k).' - epsilon .* (U_LW_1(i,:,n-1,k).' - U_LW_1(i-1,:,n-1,k).') .* ((u_plus_c_LW(i) + u_plus_c_LW(i-1))/ 2) .* (abs(p_LW_1(i+1,n-1,k) - 2*p_LW_1(i,n-1,k) + p_LW_1(i-1,n-1,k)) / (p_LW_1(i+1,n-1,k) + 2*p_LW_1(i,n-1,k) + p_LW_1(i-1,n-1,k))) ) ) - A_i_minus_half * ( (E_LW_1(i,:,n-1,k).' - epsilon .* (U_LW_1(i+1,:,n-1,k).' - U_LW_1(i,:,n-1,k).') .* ((u_plus_c_LW(i) + u_plus_c_LW(i+1)) / 2) .* (abs(p_LW_1(i+1,n-1,k) - 2*p_LW_1(i,n-1,k) + p_LW_1(i-1,n-1,k)) / (p_LW_1(i+1,n-1,k) + 2*p_LW_1(i,n-1,k) + p_LW_1(i-1,n-1,k)))) - (E_LW_1(i-1,:,n-1,k).' - epsilon .* (U_LW_1(i,:,n-1,k).' - U_LW_1(i-1,:,n-1,k).') .* ((u_plus_c_LW(i)+u_plus_c_LW(i-1)) / 2) .* (abs(p_LW_1(i+1,n-1,k) - 2*p_LW_1(i,n-1,k) + p_LW_1(i-1,n-1,k)) / (p_LW_1(i+1,n-1,k) + 2*p_LW_1(i,n-1,k) + p_LW_1(i-1,n-1,k))) ) ) )).';

                %update E based on new U
                E_LW_1(i,1,n,k) = U_LW_1(i,2,n,k);
                E_LW_1(i,2,n,k) = (U_LW_1(i,2,n,k))^2 / U_LW_1(i,1,n,k) + (gamma - 1) * (U_LW_1(i,3,n,k) - (U_LW_1(i,2,n,k))^2 / (2*U_LW_1(i,1,n,k)));
                E_LW_1(i,3,n,k) = (U_LW_1(i,2,n,k)*U_LW_1(i,3,n,k)) / U_LW_1(i,1,n,k) + (U_LW_1(i,2,n,k) / U_LW_1(i,1,n,k)) * ((gamma - 1) * (U_LW_1(i,3,n,k) - (U_LW_1(i,2,n,k))^2 / (2*U_LW_1(i,1,n,k))) );
        
                %calculate new pressure based on U
                p_LW_1(i,n,k) = (gamma - 1) * (U_LW_1(i,3,n,k) - (U_LW_1(i,2,n,k))^2 / (2*U_LW_1(i,1,n,k)));
    
                %Calculate new speed of sound and max(|u|+c)
                c_LW_1(i,1) = real(sqrt(gamma * ((gamma - 1) * (U_LW_1(i,3,n,k) - (U_LW_1(i,2,n,k))^2 / (2*U_LW_1(i,1,n,k)))) / (U_LW_1(i,1,n,k))));
                
                
            end
            u_plus_c_LW = abs(U_LW_1(:, 2, n, k) ./ U_LW_1(:, 1, n, k)) + c_LW_1;
    
            %Calculate new time step and total time elapsed
            delta_t_LW(k,1) = nu(k) * delta_x_1 / max(u_plus_c_LW);
            t_total_LW_1(k,1) = t_total_LW_1(k,1) + delta_t_LW(k,1);
        end
    
end

%analytical for LW grid 1
for k = 1:length(nu)
    %Note: U_anlyt stores rho, u, and p (not the normal U parameters of
    %rho, rho*u, and e)
    for i = 1:IL(1)
        if x_arr_1(i) >= u_S * t_total_LW_1(k)
            U_LW_1_anlyt(i, 1, k) = rho_R;
            U_LW_1_anlyt(i, 2, k) = u_R;
            U_LW_1_anlyt(i, 3, k) = p_R;
            
        elseif (x_arr_1(i) > u_2 * t_total_LW_1(k)) && (x_arr_1(i) < u_S * t_total_LW_1(k)) 
            U_LW_1_anlyt(i, 1, k) = rho_2;
            U_LW_1_anlyt(i, 2, k) = u_2;
            U_LW_1_anlyt(i, 3, k) = p_2;

        elseif (x_arr_1(i) > (u_3 - c_3) * t_total_LW_1(k)) && (x_arr_1(i) < u_3 * t_total_LW_1(k))
            U_LW_1_anlyt(i, 1, k) = rho_3;
            U_LW_1_anlyt(i, 2, k) = u_3;
            U_LW_1_anlyt(i, 3, k) = p_3;

        elseif (x_arr_1(i) > (u_L - c_L) * t_total_LW_1(k)) && (x_arr_1(i) < (u_3 - c_3) * t_total_LW_1(k))
            u_5 = 2 / (gamma + 1) * ( (x_arr_1(i) / t_total_LW_1(k)) + c_L + (gamma - 1) / 2 * u_L);
            c_5 = u_5 - x_arr_1(i) / t_total_LW_1(k);
            p_5 = p_L * (c_5 / c_L)^((2 * gamma) / (gamma - 1));
            rho_5 = (p_5 / p_L)^(1 / gamma) * rho_L;
            U_LW_1_anlyt(i, 1, k) = rho_5;
            U_LW_1_anlyt(i, 2, k) = u_5;
            U_LW_1_anlyt(i, 3, k) = p_5;

        else
            U_LW_1_anlyt(i, 1, k) = rho_L;
            U_LW_1_anlyt(i, 2, k) = u_L;
            U_LW_1_anlyt(i, 3, k) = p_L;
        end
    end
end

%Grid 2 Lax-Wendroff with Artificial Dissipation
for k = 1:length(nu)
    %Loop through all of time    
        for i = 1:IL(2)
        %Set initial conditions
            if x_arr_2(i) <= 0
            %Left half
                U_LW_2(i,1,1,k) = rho_L;
                U_LW_2(i,2,1,k) = rho_L * u_L;
                U_LW_2(i,3,1,k) = p_L / (gamma - 1) + (rho_L / 2) * u_L^2;
    
                E_LW_2(i,1,1,k) = U_LW_2(i,2,1,k);
                E_LW_2(i,2,1,k) = (U_LW_2(i,2,1,k))^2 / U_LW_2(i,1,1,k) + (gamma - 1) * (U_LW_2(i,3,1,k) - (U_LW_2(i,2,1,k))^2 / (2*U_LW_2(i,1,1,k)));
                E_LW_2(i,3,1,k) = (U_LW_2(i,2,1,k)*U_LW_2(i,3,1,k)) / U_LW_2(i,1,1,k) + (U_LW_2(i,2,1,k) / U_LW_2(i,1,1,k)) * ((gamma - 1) * (U_LW_2(i,3,1,k) - (U_LW_2(i,2,1,k))^2 / (2*U_LW_2(i,1,1,k))) );
    
                p_LW_2(i,1,k) = p_L;
    
                c_LW_2(i,1) = sqrt(gamma * ((gamma - 1) * (U_LW_2(i,3,1,k) - (U_LW_2(i,2,1,k))^2 / (2*U_LW_2(i,1,1,k)))) / (U_LW_2(i,1,1,k)));
                u_plus_c_LW_2(i,1) = abs(U_LW_2(i, 2, 1, k) / U_LW_2(i, 1, 1, k)) + c_LW_2(i);
    
            else
            %Right half
                U_LW_2(i,1,1,k) = rho_R;
                U_LW_2(i,2,1,k) = rho_R * u_R;
                U_LW_2(i,3,1,k) = p_R / (gamma - 1) + (rho_R / 2) * u_R^2;
    
                E_LW_2(i,1,1,k) = U_LW_2(i,2,1,k);
                E_LW_2(i,2,1,k) = (U_LW_2(i,2,1,k))^2 / U_LW_2(i,1,1,k) + (gamma - 1) * (U_LW_2(i,3,1,k) - (U_LW_2(i,2,1,k))^2 / (2*U_LW_2(i,1,1,k)));
                E_LW_2(i,3,1,k) = (U_LW_2(i,2,1,k)*U_LW_2(i,3,1,k)) / U_LW_2(i,1,1,k) + (U_LW_2(i,2,1,k) / U_LW_2(i,1,1,k)) * ((gamma - 1) * (U_LW_2(i,3,1,k) - (U_LW_2(i,2,1,k))^2 / (2*U_LW_2(i,1,1,k))) );
    
                p_LW_2(i,1,k) = p_R;
    
                c_LW_2(i,1) = sqrt(gamma * ((gamma - 1) * (U_LW_2(i,3,1,k) - (U_LW_2(i,2,1,k))^2 / (2*U_LW_2(i,1,1,k)))) / (U_LW_2(i,1,1,k)));
                u_plus_c_LW_2(i,1) = abs(U_LW_2(i, 2, 1, k) / U_LW_2(i, 1, 1, k)) + c_LW_2(i);
    
            end
        end  
        %set initial time step
        delta_t_LW_2(k, 1) = nu(k) * delta_x_2 / max(u_plus_c_LW_2);
        t_total_LW_2(k, 1) = delta_t_LW(k,1);
    
        %Set Dirichlet BCs for all of time
        U_LW_2(1,1,:,k) = rho_L;
        U_LW_2(1,2,:,k) = rho_L * u_L;
        U_LW_2(1,3,:,k) = p_L / (gamma - 1) + (rho_L / 2) * u_L^2;
        U_LW_2(end,1,:,k) = rho_R;
        U_LW_2(end,2,:,k) = rho_R * u_R;
        U_LW_2(end,3,:,k) = p_R / (gamma - 1) + (rho_R / 2) * u_R^2;
    
        E_LW_2(1,1,:,k) = U_LW_2(1,2,1,k);
        E_LW_2(1,2,:,k) = (U_LW_2(1,2,1,k))^2 / U_LW_2(1,1,1,k) + (gamma - 1) * (U_LW_2(1,3,1,k) - (U_LW_2(1,2,1,k))^2 / (2*U_LW_2(1,1,1,k)));
        E_LW_2(1,3,:,k) = (U_LW_2(1,2,1,k)*U_LW_2(1,3,1,k)) / U_LW_2(1,1,1,k) + (U_LW_2(1,2,1,k) / U_LW_2(1,1,1,k)) * ((gamma - 1) * (U_LW_2(1,3,1,k) - (U_LW_2(1,2,1,k))^2 / (2*U_LW_2(1,1,1,k))) );
        E_LW_2(end,1,:,k) = U_LW_2(end,2,1,k);
        E_LW_2(end,2,:,k) = (U_LW_2(end,2,1,k))^2 / U_LW_2(end,1,1,k) + (gamma - 1) * (U_LW_2(end,3,1,k) - (U_LW_2(end,2,1,k))^2 / (2*U_LW_2(end,1,1,k)));
        E_LW_2(end,3,:,k) = (U_LW_2(end,2,1,k)*U_LW_2(end,3,1,k)) / U_LW_2(end,1,1,k) + (U_LW_2(end,2,1,k) / U_LW_2(end,1,1,k)) * ((gamma - 1) * (U_LW_2(end,3,1,k) - (U_LW_2(end,2,1,k))^2 / (2*U_LW_2(end,1,1,k))) );
        
        p_LW_2(1,:,k) = p_L;
        p_LW_2(end,:,k) = p_R;

        c_LW_2(1,1) = sqrt(gamma * p_L / rho_L);
        u_plus_c_LW_2(1,1) = abs(U_LW_2(1, 2, 1, k) / U_LW_2(1, 1, 1, k)) + c_LW_2(1);
        c_LW_2(end,1) = sqrt(gamma * p_R / rho_R);
        u_plus_c_LW_2(end,1) = abs(U_LW_2(end, 2, 1, k) / U_LW_2(end, 1, 1, k)) + c_LW_2(end);
        for n = 2:t_steps(2)+1
            for i = 2:IL(2)-1
            %Perform iterative procedure for all of interior grid points
                %assign jacobians at i = 1 + 1/2 and i = 1 - 1/2
                A_i_plus_half_1 = [0, 1, 0];
                A_i_plus_half_2 = [-(3- gamma) * ((U_LW_2(i,2,n-1,k) + U_LW_2(i+1,2,n-1,k))/2)^2 / (2 * ((U_LW_2(i,1,n-1,k) + U_LW_2(i+1,1,n-1,k)) / 2)^2), (3 - gamma) * ((U_LW_2(i,2,n-1,k) + U_LW_2(i+1,2,n-1,k)) / 2) / ((U_LW_2(i,1,n-1,k) + U_LW_2(i+1,1,n-1,k)) / 2), (gamma - 1)];
                A_i_plus_half_3 = [(gamma - 1) * (((U_LW_2(i,2,n-1,k) + U_LW_2(i+1,2,n-1,k)) / 2)^3 / ((U_LW_2(i,1,n-1,k) + U_LW_2(i+1,1,n-1,k)) / 2)^3) - gamma * ( (U_LW_2(i,2,n-1,k) + U_LW_2(i+1,2,n-1,k)) / 2) / ((U_LW_2(i,1,n-1,k) + U_LW_2(i+1,1,n-1,k)) / 2)^2 * ((U_LW_2(i,3,n-1,k) + U_LW_2(i+1,3,n-1,k)) / 2), gamma * ((U_LW_2(i,3,n-1,k) + U_LW_2(i+1,3,n-1,k)) / 2) / ((U_LW_2(i,1,n-1,k) + U_LW_2(i+1,1,n-1,k)) / 2) - 1.5 * (gamma - 1) * ((U_LW_2(i,2,n-1,k) + U_LW_2(i+1,2,n-1,k)) / 2)^2 / ((U_LW_2(i,1,n-1,k) + U_LW_2(i+1,1,n-1,k)) / 2)^2, gamma * ((U_LW_2(i,2,n-1,k) + U_LW_2(i+1,2,n-1,k)) / 2) / ((U_LW_2(i,1,n-1,k) + U_LW_2(i+1,1,n-1,k)) / 2)];

                A_i_minus_half_1 = [0, 1, 0];
                A_i_minus_half_2 = [-(3- gamma) * ((U_LW_2(i,2,n-1,k) + U_LW_2(i-1,2,n-1,k))/2)^2 / (2 * ((U_LW_2(i,1,n-1,k) + U_LW_2(i-1,1,n-1,k)) / 2)^2), (3 - gamma) * ((U_LW_2(i,2,n-1,k) + U_LW_2(i-1,2,n-1,k)) / 2) / ((U_LW_2(i,1,n-1,k) + U_LW_2(i-1,1,n-1,k)) / 2), (gamma - 1)];
                A_i_minus_half_3 = [(gamma - 1) * (((U_LW_2(i,2,n-1,k) + U_LW_2(i-1,2,n-1,k)) / 2)^3 / ((U_LW_2(i,1,n-1,k) + U_LW_2(i-1,1,n-1,k)) / 2)^3) - gamma * ( (U_LW_2(i,2,n-1,k) + U_LW_2(i-1,2,n-1,k)) / 2) / ((U_LW_2(i,1,n-1,k) + U_LW_2(i-1,1,n-1,k)) / 2)^2 * ((U_LW_2(i,3,n-1,k) + U_LW_2(i-1,3,n-1,k)) / 2), gamma * ((U_LW_2(i,3,n-1,k) + U_LW_2(i-1,3,n-1,k)) / 2) / ((U_LW_2(i,1,n-1,k) + U_LW_2(i-1,1,n-1,k)) / 2) - 1.5 * (gamma - 1) * ((U_LW_2(i,2,n-1,k) + U_LW_2(i-1,2,n-1,k)) / 2)^2 / ((U_LW_2(i,1,n-1,k) + U_LW_2(i-1,1,n-1,k)) / 2)^2, gamma * ((U_LW_2(i,2,n-1,k) + U_LW_2(i-1,2,n-1,k)) / 2) / ((U_LW_2(i,1,n-1,k) + U_LW_2(i-1,1,n-1,k)) / 2)];
                
                A_i_plus_half = [A_i_plus_half_1; A_i_plus_half_2; A_i_plus_half_3];
                A_i_minus_half = [A_i_minus_half_1; A_i_minus_half_2; A_i_minus_half_3];

                %calculate new U
                U_LW_2(i,:,n,k) = (U_LW_2(i,:,n-1,k).' - delta_t_LW_2(k,1) / (2 * delta_x_2) .* ((E_LW_2(i+1,:,n-1,k).' - epsilon .* (U_LW_2(i+1,:,n-1,k).' - U_LW_2(i,:,n-1,k).') .* ((u_plus_c_LW_2(i) + u_plus_c_LW_2(i+1)) / 2) .* (abs(p_LW_2(i+1,n-1,k) - 2*p_LW_2(i,n-1,k) + p_LW_2(i-1,n-1,k)) / (p_LW_2(i+1,n-1,k) + 2*p_LW_2(i,n-1,k) + p_LW_2(i-1,n-1,k)))) - (E_LW_2(i-1,:,n-1,k).' - epsilon .* (U_LW_2(i,:,n-1,k).' - U_LW_2(i-1,:,n-1,k).') .* ((u_plus_c_LW_2(i) + u_plus_c_LW_2(i-1)) / 2) .* (abs(p_LW_2(i+1,n-1,k) - 2*p_LW_2(i,n-1,k) + p_LW_2(i-1,n-1,k)) / (p_LW_2(i+1,n-1,k) + 2*p_LW_2(i,n-1,k) + p_LW_2(i-1,n-1,k))))) + delta_t_LW_2(k,1)^2 / (2 * delta_x_2^2) .* (A_i_plus_half * ( (E_LW_2(i+1,:,n-1,k).' - epsilon .* (U_LW_2(i+1,:,n-1,k).' - U_LW_2(i,:,n-1,k).') .* ((u_plus_c_LW_2(i) + u_plus_c_LW_2(i+1)) / 2) .* (abs(p_LW_2(i+1,n-1,k) - 2*p_LW_2(i,n-1,k) + p_LW_2(i-1,n-1,k)) / (p_LW_2(i+1,n-1,k) + 2*p_LW_2(i,n-1,k) + p_LW_2(i-1,n-1,k)))) - (E_LW_2(i,:,n-1,k).' - epsilon .* (U_LW_2(i,:,n-1,k).' - U_LW_2(i-1,:,n-1,k).') .* ((u_plus_c_LW_2(i) + u_plus_c_LW_2(i-1))/ 2) .* (abs(p_LW_2(i+1,n-1,k) - 2*p_LW_2(i,n-1,k) + p_LW_2(i-1,n-1,k)) / (p_LW_2(i+1,n-1,k) + 2*p_LW_2(i,n-1,k) + p_LW_2(i-1,n-1,k))) ) ) - A_i_minus_half * ( (E_LW_2(i,:,n-1,k).' - epsilon .* (U_LW_2(i+1,:,n-1,k).' - U_LW_2(i,:,n-1,k).') .* ((u_plus_c_LW_2(i) + u_plus_c_LW_2(i+1)) / 2) .* (abs(p_LW_2(i+1,n-1,k) - 2*p_LW_2(i,n-1,k) + p_LW_2(i-1,n-1,k)) / (p_LW_2(i+1,n-1,k) + 2*p_LW_2(i,n-1,k) + p_LW_2(i-1,n-1,k)))) - (E_LW_2(i-1,:,n-1,k).' - epsilon .* (U_LW_2(i,:,n-1,k).' - U_LW_2(i-1,:,n-1,k).') .* ((u_plus_c_LW_2(i)+u_plus_c_LW_2(i-1)) / 2) .* (abs(p_LW_2(i+1,n-1,k) - 2*p_LW_2(i,n-1,k) + p_LW_2(i-1,n-1,k)) / (p_LW_2(i+1,n-1,k) + 2*p_LW_2(i,n-1,k) + p_LW_2(i-1,n-1,k))) ) ) )).';

                %update E based on new U
                E_LW_2(i,1,n,k) = U_LW_2(i,2,n,k);
                E_LW_2(i,2,n,k) = (U_LW_2(i,2,n,k))^2 / U_LW_2(i,1,n,k) + (gamma - 1) * (U_LW_2(i,3,n,k) - (U_LW_2(i,2,n,k))^2 / (2*U_LW_2(i,1,n,k)));
                E_LW_2(i,3,n,k) = (U_LW_2(i,2,n,k)*U_LW_2(i,3,n,k)) / U_LW_2(i,1,n,k) + (U_LW_2(i,2,n,k) / U_LW_2(i,1,n,k)) * ((gamma - 1) * (U_LW_2(i,3,n,k) - (U_LW_2(i,2,n,k))^2 / (2*U_LW_2(i,1,n,k))) );
        
                %calculate new pressure based on U
                p_LW_2(i,n,k) = (gamma - 1) * (U_LW_2(i,3,n,k) - (U_LW_2(i,2,n,k))^2 / (2*U_LW_2(i,1,n,k)));
    
                %Calculate new speed of sound and max(|u|+c)
                c_LW_2(i,1) = real(sqrt(gamma * ((gamma - 1) * (U_LW_2(i,3,n,k) - (U_LW_2(i,2,n,k))^2 / (2*U_LW_2(i,1,n,k)))) / (U_LW_2(i,1,n,k))));
                
                
            end
            u_plus_c_LW_2 = abs(U_LW_2(:, 2, n, k) ./ U_LW_2(:, 1, n, k)) + c_LW_2;
    
            %Calculate new time step and total time elapsed
            delta_t_LW_2(k,1) = nu(k) * delta_x_2 / max(u_plus_c_LW_2);
            t_total_LW_2(k,1) = t_total_LW_2(k,1) + delta_t_LW_2(k,1);
        end
    
end

%analytical for LW grid 2
for k = 1:length(nu)
    %Note: U_anlyt stores rho, u, and p (not the normal U parameters of
    %rho, rho*u, and e)
    for i = 1:IL(2)
        if x_arr_2(i) >= u_S * t_total_LW_2(k)
            U_LW_2_anlyt(i, 1, k) = rho_R;
            U_LW_2_anlyt(i, 2, k) = u_R;
            U_LW_2_anlyt(i, 3, k) = p_R;
            
        elseif (x_arr_2(i) > u_2 * t_total_LW_2(k)) && (x_arr_2(i) < u_S * t_total_LW_2(k)) 
            U_LW_2_anlyt(i, 1, k) = rho_2;
            U_LW_2_anlyt(i, 2, k) = u_2;
            U_LW_2_anlyt(i, 3, k) = p_2;

        elseif (x_arr_2(i) > (u_3 - c_3) * t_total_LW_2(k)) && (x_arr_2(i) < u_3 * t_total_LW_2(k))
            U_LW_2_anlyt(i, 1, k) = rho_3;
            U_LW_2_anlyt(i, 2, k) = u_3;
            U_LW_2_anlyt(i, 3, k) = p_3;

        elseif (x_arr_2(i) > (u_L - c_L) * t_total_LW_2(k)) && (x_arr_2(i) < (u_3 - c_3) * t_total_LW_2(k))
            u_5 = 2 / (gamma + 1) * ( (x_arr_2(i) / t_total_LW_2(k)) + c_L + (gamma - 1) / 2 * u_L);
            c_5 = u_5 - x_arr_2(i) / t_total_LW_2(k);
            p_5 = p_L * (c_5 / c_L)^((2 * gamma) / (gamma - 1));
            rho_5 = (p_5 / p_L)^(1 / gamma) * rho_L;
            U_LW_2_anlyt(i, 1, k) = rho_5;
            U_LW_2_anlyt(i, 2, k) = u_5;
            U_LW_2_anlyt(i, 3, k) = p_5;

        else
            U_LW_2_anlyt(i, 1, k) = rho_L;
            U_LW_2_anlyt(i, 2, k) = u_L;
            U_LW_2_anlyt(i, 3, k) = p_L;
        end
    end
end

%Steger-Warming calculations for Grid 1
for k = 1:length(nu)
    for i = 1:IL(1)
    %Populate initial conditions
        if x_arr_1(i) <= 0
            %Left half
                U_SW_1(i,1,1,k) = rho_L;
                U_SW_1(i,2,1,k) = rho_L * u_L;
                U_SW_1(i,3,1,k) = p_L / (gamma - 1) + (rho_L / 2) * u_L^2;
    
                E_SW_1(i,1,1,k) = U_SW_1(i,2,1,k);
                E_SW_1(i,2,1,k) = (U_SW_1(i,2,1,k))^2 / U_SW_1(i,1,1,k) + (gamma - 1) * (U_SW_1(i,3,1,k) - (U_SW_1(i,2,1,k))^2 / (2*U_SW_1(i,1,1,k)));
                E_SW_1(i,3,1,k) = (U_SW_1(i,2,1,k)*U_SW_1(i,3,1,k)) / U_SW_1(i,1,1,k) + (U_SW_1(i,2,1,k) / U_SW_1(i,1,1,k)) * ((gamma - 1) * (U_SW_1(i,3,1,k) - (U_SW_1(i,2,1,k))^2 / (2*U_SW_1(i,1,1,k))) );

                c_SW_1(i,1) = sqrt(gamma * ((gamma - 1) * (U_SW_1(i,3,1,k) - (U_SW_1(i,2,1,k))^2 / (2*U_SW_1(i,1,1,k)))) / (U_SW_1(i,1,1,k)));

        else
            %Right half
                U_SW_1(i,1,1,k) = rho_R;
                U_SW_1(i,2,1,k) = rho_R * u_R;
                U_SW_1(i,3,1,k) = p_R / (gamma - 1) + (rho_R / 2) * u_R^2;
    
                E_SW_1(i,1,1,k) = U_SW_1(i,2,1,k);
                E_SW_1(i,2,1,k) = (U_SW_1(i,2,1,k))^2 / U_SW_1(i,1,1,k) + (gamma - 1) * (U_SW_1(i,3,1,k) - (U_SW_1(i,2,1,k))^2 / (2*U_SW_1(i,1,1,k)));
                E_SW_1(i,3,1,k) = (U_SW_1(i,2,1,k)*U_SW_1(i,3,1,k)) / U_SW_1(i,1,1,k) + (U_SW_1(i,2,1,k) / U_SW_1(i,1,1,k)) * ((gamma - 1) * (U_SW_1(i,3,1,k) - (U_SW_1(i,2,1,k))^2 / (2*U_SW_1(i,1,1,k))) );
    
                c_SW_1(i,1) = sqrt(gamma * ((gamma - 1) * (U_SW_1(i,3,1,k) - (U_SW_1(i,2,1,k))^2 / (2*U_SW_1(i,1,1,k)))) / (U_SW_1(i,1,1,k)));

        end        
    end
end


%Steger-Warming Grid 1
for i = 1:IL(1)
    %Populate initial conditions
        if x_arr_1(i) <= 0
            %Left half
                U1(1,i,1) = rho_L;
                U2(1,i,1) = rho_L * u_L;
                U3(1,i,1) = p_L / (gamma - 1) + (rho_L / 2) * u_L^2;

                U1_nu2(1,i,1) = rho_L;
                U2_nu2(1,i,1) = rho_L * u_L;
                U3_nu2(1,i,1) = p_L / (gamma - 1) + (rho_L / 2) * u_L^2;
    
        else
            %Right half
                U1(1,i,1) = rho_R;
                U2(1,i,1) = rho_R * u_R;
                U3(1,i,1) = p_R / (gamma - 1) + (rho_R / 2) * u_R^2;

                U1_nu2(1,i,1) = rho_R;
                U2_nu2(1,i,1) = rho_R * u_R;
                U3_nu2(1,i,1) = p_R / (gamma - 1) + (rho_R / 2) * u_R^2;
    
        end        
end

for n = 1:t_steps(1)
   
    %initialize positive and negative eigenvalue matrices
    eigenN = zeros(3,3,IL(1));
    eigenP = zeros(3,3,IL(1));
    
    %calculate timestep size at each timestep
    %need to loop through eigenvalues at every x location to find maximum
    %eigenvalues
    
    %Recalculate dt at every new timestep 
    for j = 1:IL(1) %for each point in the domain
        At(:,:,j) = [U2(:,j,n)./U1(:,j,n) U1(:,j,n)  0;...
                     0 U2(:,j,n)./U1(:,j,n) 1/U1(:,j,n);...
                     0 gamma*(gamma-1)*(U3(:,j,n)-(U2(:,j,n).^2)./(2*U1(:,j,n))) U2(:,j,n)./U1(:,j,n);];

       egn(:,:,j) = eig(At(:,:,j));
       dt_check(j) = nu_1*delta_x_1/max(abs(egn(:,:,j)));

       At_nu2(:,:,j) = [U2_nu2(:,j,n)./U1_nu2(:,j,n) U1_nu2(:,j,n)  0;...
                     0 U2_nu2(:,j,n)./U1_nu2(:,j,n) 1/U1_nu2(:,j,n);...
                     0 gamma*(gamma-1)*(U3_nu2(:,j,n)-(U2_nu2(:,j,n).^2)./(2*U1_nu2(:,j,n))) U2_nu2(:,j,n)./U1_nu2(:,j,n);];

       egn_nu2(:,:,j) = eig(At_nu2(:,:,j));
       dt_check_nu2(j) = nu_2*delta_x_1/max(abs(egn_nu2(:,:,j)));
    end
    %timestep
    dt(n) = min(dt_check); %choose the minimum dt at each timestep
    t_total_SW_1(1,1) = t_total_SW_1(1,1) + dt(n);

    dt_nu2(n) = min(dt_check_nu2); %choose the minimum dt at each timestep
    t_total_SW_1(2,1) = t_total_SW_1(2,1) + dt_nu2(n);
     
     
         %Pressure 
         Pnd(:,:,n) = (gamma-1)*(U3(:,:,n)-(U2(:,:,n).^2)./(2*U1(:,:,n)) );
         
         Pnd_nu2(:,:,n) = (gamma-1)*(U3_nu2(:,:,n)-(U2_nu2(:,:,n).^2)./(2*U1_nu2(:,:,n)) );

         %density is U1
         %get speed of sound c
         cS(:,:,n) = sqrt(gamma*Pnd(:,:,n)./U1(:,:,n));

         cS_nu2(:,:,n) = sqrt(gamma*Pnd_nu2(:,:,n)./U1_nu2(:,:,n));


         %Populate P and P_inv matrices for each x location
    for i = 1:IL(1)
        
        Pvec(:,:,i) = [1 1/(2*cS(:,i,n)^2) 1/(2*cS(:,i,n)^2);...
            U2(:,i,n)/U1(:,i,n) ((U2(:,i,n)/U1(:,i,n))+cS(:,i,n))/(2*cS(:,i,n)^2) ((U2(:,i,n)/U1(:,i,n))-cS(:,i,n))/(2*cS(:,i,n)^2);...
            (1/2)*(U2(:,i,n)/U1(:,i,n))^2 (1/(4*cS(:,i,n)^2))*(U2(:,i,n)/U1(:,i,n))^2+(U2(:,i,n)/U1(:,i,n))/(2*cS(:,i,n))+1/(2*(gamma-1)) (1/(4*cS(:,i,n)^2))*(U2(:,i,n)/U1(:,i,n))^2-(U2(:,i,n)/U1(:,i,n))/(2*cS(:,i,n))+1/(2*(gamma-1));];

        Pvec_inv(:,:,i) = [1-(gamma-1)*(1/(2*cS(:,i,n)^2))*(U2(:,i,n)/U1(:,i,n))^2 (gamma-1)*(U2(:,i,n)/U1(:,i,n))/cS(:,i,n)^2 -(gamma-1)/cS(:,i,n)^2;...
            -(U2(:,i,n)/U1(:,i,n))*cS(:,i,n)+(gamma-1)*(1/2)*(U2(:,i,n)/U1(:,i,n))^2 cS(:,i,n)-(gamma-1)*(U2(:,i,n)/U1(:,i,n)) (gamma-1);...
             (U2(:,i,n)/U1(:,i,n))*cS(:,i,n)+(gamma-1)*(1/2)*(U2(:,i,n)/U1(:,i,n))^2 -cS(:,i,n)-(gamma-1)*(U2(:,i,n)/U1(:,i,n)) (gamma-1);]; 
         
        Pvec_nu2(:,:,i) = [1 1/(2*cS_nu2(:,i,n)^2) 1/(2*cS_nu2(:,i,n)^2);...
            U2_nu2(:,i,n)/U1_nu2(:,i,n) ((U2_nu2(:,i,n)/U1_nu2(:,i,n))+cS_nu2(:,i,n))/(2*cS_nu2(:,i,n)^2) ((U2_nu2(:,i,n)/U1_nu2(:,i,n))-cS_nu2(:,i,n))/(2*cS_nu2(:,i,n)^2);...
            (1/2)*(U2_nu2(:,i,n)/U1_nu2(:,i,n))^2 (1/(4*cS_nu2(:,i,n)^2))*(U2_nu2(:,i,n)/U1_nu2(:,i,n))^2+(U2_nu2(:,i,n)/U1_nu2(:,i,n))/(2*cS_nu2(:,i,n))+1/(2*(gamma-1)) (1/(4*cS_nu2(:,i,n)^2))*(U2_nu2(:,i,n)/U1_nu2(:,i,n))^2-(U2_nu2(:,i,n)/U1_nu2(:,i,n))/(2*cS_nu2(:,i,n))+1/(2*(gamma-1));];

        Pvec_inv_nu2(:,:,i) = [1-(gamma-1)*(1/(2*cS_nu2(:,i,n)^2))*(U2_nu2(:,i,n)/U1_nu2(:,i,n))^2 (gamma-1)*(U2_nu2(:,i,n)/U1_nu2(:,i,n))/cS_nu2(:,i,n)^2 -(gamma-1)/cS_nu2(:,i,n)^2;...
            -(U2_nu2(:,i,n)/U1_nu2(:,i,n))*cS_nu2(:,i,n)+(gamma-1)*(1/2)*(U2_nu2(:,i,n)/U1_nu2(:,i,n))^2 cS_nu2(:,i,n)-(gamma-1)*(U2_nu2(:,i,n)/U1(:,i,n)) (gamma-1);...
             (U2_nu2(:,i,n)/U1_nu2(:,i,n))*cS_nu2(:,i,n)+(gamma-1)*(1/2)*(U2_nu2(:,i,n)/U1_nu2(:,i,n))^2 -cS_nu2(:,i,n)-(gamma-1)*(U2_nu2(:,i,n)/U1(:,i,n)) (gamma-1);]; 
         
    end
    
     %flux vector terms
            E1(:,:,n) = U2(:,:,n);
            E2(:,:,n) = (U2(:,:,n).^2)./U1(:,:,n)+(gamma-1)*(U3(:,:,n)-(U2(:,:,n).^2)./(2*U1(:,:,n)) );
            E3(:,:,n) = U2(:,:,n).*U3(:,:,n)./U1(:,:,n)+(U2(:,:,n)./U1(:,:,n)).*((gamma-1).*(U3(:,:,n)-(U2(:,:,n).^2)./(2*U1(:,:,n)) )); 
    
    
        for i = 1:IL(1)
            
             %conservative vector
             U(:,:,i) = [U1(:,i,n); U2(:,i,n); U3(:,i,n);];
             U_nu2(:,:,i) = [U1_nu2(:,i,n); U2_nu2(:,i,n); U3_nu2(:,i,n);];
            
                
             if i < IL(1)
                %i+1 conservative vector needed for Ei+1 
                U(:,:,i+1) = [U1(:,i+1,n); U2(:,i+1,n); U3(:,i+1,n);];
                U_nu2(:,:,i+1) = [U1_nu2(:,i+1,n); U2_nu2(:,i+1,n); U3_nu2(:,i+1,n);];
                
             end               
            
            %populate positive and negative eigenvalue matrices at each x
            %location
            
            %%E+ and E- at location i
            for d = 1:3
                eigenP(d,d,i) = 0.5*(egn(d,:,i)+abs(egn(d,:,i)));
                eigenN(d,d,i) = 0.5*(egn(d,:,i)-abs(egn(d,:,i)));

                eigenP_nu2(d,d,i) = 0.5*(egn_nu2(d,:,i)+abs(egn_nu2(d,:,i)));
                eigenN_nu2(d,d,i) = 0.5*(egn_nu2(d,:,i)-abs(egn_nu2(d,:,i)));
            end
            
             %calculate E+ and E- at i
                W(:,:,i) = Pvec_inv(:,:,i)*U(:,:,i);
                
                Wp(:,:,i) = eigenP(:,:,i)*W(:,:,i);
                Wn(:,:,i) = eigenN(:,:,i)*W(:,:,i);
                
                Ep(:,:,i) = Pvec(:,:,i)*Wp(:,:,i);
                En(:,:,i) = Pvec(:,:,i)*Wn(:,:,i);

                W_nu2(:,:,i) = Pvec_inv_nu2(:,:,i)*U_nu2(:,:,i);
                
                Wp_nu2(:,:,i) = eigenP_nu2(:,:,i)*W_nu2(:,:,i);
                Wn_nu2(:,:,i) = eigenN_nu2(:,:,i)*W_nu2(:,:,i);
                
                Ep_nu2(:,:,i) = Pvec_nu2(:,:,i)*Wp_nu2(:,:,i);
                En_nu2(:,:,i) = Pvec_nu2(:,:,i)*Wn_nu2(:,:,i);
                
                

            %apply dirichlet boundary conditions
            if i == 1
                
                U1(:,i,n+1) =  rho_L;
                U2(:,i,n+1) = rho_L*u_L;
                U3(:,i,n+1) = p_L/(gamma-1)+rho_L/2*(u_L^2);
                
                U1_nu2(:,i,n+1) =  rho_L;
                U2_nu2(:,i,n+1) = rho_L*u_L;
                U3_nu2(:,i,n+1) = p_L/(gamma-1)+rho_L/2*(u_L^2);
                
            elseif i == IL(1)

                U1(:,i,n+1) =  rho_R;
                U2(:,i,n+1) = rho_R*u_R;
                U3(:,i,n+1) = p_R/(gamma-1)+rho_R/2*(u_R^2);
                
                U1_nu2(:,i,n+1) =  rho_R;
                U2_nu2(:,i,n+1) = rho_R*u_R;
                U3_nu2(:,i,n+1) = p_R/(gamma-1)+rho_R/2*(u_R^2);
                
            else
             
             %% E+ and E- at i+1
            for d = 1:3
                eigenP(d,d,i+1) = 0.5*(egn(d,:,i+1)+abs(egn(d,:,i+1)));
                eigenN(d,d,i+1) = 0.5*(egn(d,:,i+1)-abs(egn(d,:,i+1)));  

                eigenP_nu2(d,d,i+1) = 0.5*(egn_nu2(d,:,i+1)+abs(egn_nu2(d,:,i+1)));
                eigenN_nu2(d,d,i+1) = 0.5*(egn_nu2(d,:,i+1)-abs(egn_nu2(d,:,i+1)));  
            end
                                            
                
                %Calculate E- at i+1
                W(:,:,i+1) = Pvec_inv(:,:,i+1)*U(:,:,i+1); %W = (P^-1)U
                Wn(:,:,i+1) = eigenN(:,:,i+1)*W(:,:,i+1); %Wn = (lambda^-)W
                En(:,:,i+1) = Pvec(:,:,i+1)*Wn(:,:,i+1);% En = E^- = (P)Wn

                W_nu2(:,:,i+1) = Pvec_inv_nu2(:,:,i+1)*U_nu2(:,:,i+1); %W = (P^-1)U
                Wn_nu2(:,:,i+1) = eigenN_nu2(:,:,i+1)*W_nu2(:,:,i+1); %Wn = (lambda^-)W
                En_nu2(:,:,i+1) = Pvec_nu2(:,:,i+1)*Wn_nu2(:,:,i+1);% En = E^- = (P)Wn
                
                %calculation
                U1(:,i,n+1) = U1(:,i,n) - (dt(n)/(delta_x_1))*( (Ep(1,:,i)-Ep(1,:,i-1))+(En(1,:,i+1)-En(1,:,i))  );
                U2(:,i,n+1) = U2(:,i,n) - (dt(n)/(delta_x_1))*( (Ep(2,:,i)-Ep(2,:,i-1))+(En(2,:,i+1)-En(2,:,i))  );
                U3(:,i,n+1) = U3(:,i,n) - (dt(n)/(delta_x_1))*( (Ep(3,:,i)-Ep(3,:,i-1))+(En(3,:,i+1)-En(3,:,i))  );

                U1_nu2(:,i,n+1) = U1_nu2(:,i,n) - (dt_nu2(n)/(delta_x_1))*( (Ep_nu2(1,:,i)-Ep_nu2(1,:,i-1))+(En_nu2(1,:,i+1)-En_nu2(1,:,i))  );
                U2_nu2(:,i,n+1) = U2_nu2(:,i,n) - (dt_nu2(n)/(delta_x_1))*( (Ep_nu2(2,:,i)-Ep_nu2(2,:,i-1))+(En_nu2(2,:,i+1)-En_nu2(2,:,i))  );
                U3_nu2(:,i,n+1) = U3_nu2(:,i,n) - (dt_nu2(n)/(delta_x_1))*( (Ep_nu2(3,:,i)-Ep_nu2(3,:,i-1))+(En_nu2(3,:,i+1)-En_nu2(3,:,i))  );
         
            end
            
        end  
end

%SW Grid 2
for i = 1:IL(2)
    %Populate initial conditions
        if x_arr_2(i) <= 0
            %Left half
                U1_2(1,i,1) = rho_L;
                U2_2(1,i,1) = rho_L * u_L;
                U3_2(1,i,1) = p_L / (gamma - 1) + (rho_L / 2) * u_L^2;
    
                U1_2_nu2(1,i,1) = rho_L;
                U2_2_nu2(1,i,1) = rho_L * u_L;
                U3_2_nu2(1,i,1) = p_L / (gamma - 1) + (rho_L / 2) * u_L^2;
        else
            %Right half
                U1_2(1,i,1) = rho_R;
                U2_2(1,i,1) = rho_R * u_R;
                U3_2(1,i,1) = p_R / (gamma - 1) + (rho_R / 2) * u_R^2;

                U1_2_nu2(1,i,1) = rho_R;
                U2_2_nu2(1,i,1) = rho_R * u_R;
                U3_2_nu2(1,i,1) = p_R / (gamma - 1) + (rho_R / 2) * u_R^2;

        end        
end

for n = 1:t_steps(2)
   
    %initialize positive and negative eigenvalue matrices
    eigenN_2 = zeros(3,3,IL(2));
    eigenP_2 = zeros(3,3,IL(2));
    
    eigenN_2_nu2 = zeros(3,3,IL(2));
    eigenP_2_nu2 = zeros(3,3,IL(2));

    %calculate timestep size at each timestep
    %need to loop through eigenvalues at every x location to find maximum
    %eigenvalues
    
    %Recalculate dt at every new timestep 
    for j = 1:IL(2) %for each point in the domain
        At_2(:,:,j) = [U2_2(:,j,n)./U1_2(:,j,n) U1_2(:,j,n)  0;...
                     0 U2_2(:,j,n)./U1_2(:,j,n) 1/U1_2(:,j,n);...
                     0 gamma*(gamma-1)*(U3_2(:,j,n)-(U2_2(:,j,n).^2)./(2*U1_2(:,j,n))) U2_2(:,j,n)./U1_2(:,j,n);];

       egn_2(:,:,j) = eig(At_2(:,:,j));
       dt_check_2(j) = nu_1*delta_x_2/max(abs(egn_2(:,:,j)));

       At_2_nu2(:,:,j) = [U2_2_nu2(:,j,n)./U1_2_nu2(:,j,n) U1_2_nu2(:,j,n)  0;...
                     0 U2_2_nu2(:,j,n)./U1_2_nu2(:,j,n) 1/U1_2_nu2(:,j,n);...
                     0 gamma*(gamma-1)*(U3_2_nu2(:,j,n)-(U2_2_nu2(:,j,n).^2)./(2*U1_2_nu2(:,j,n))) U2_2_nu2(:,j,n)./U1_2_nu2(:,j,n);];

       egn_2_nu2(:,:,j) = eig(At_2_nu2(:,:,j));
       dt_check_2_nu2(j) = nu_2*delta_x_2/max(abs(egn_2_nu2(:,:,j)));
    end
    %timestep
    dt_2(n) = min(dt_check_2); %choose the minimum dt at each timestep
    t_total_SW_2(1,1) = t_total_SW_2(1,1) + dt_2(n);

    dt_2_nu2(n) = min(dt_check_2_nu2); %choose the minimum dt at each timestep
    t_total_SW_2_nu2(2,1) = t_total_SW_2(2,1) + dt_2_nu2(n);
     
     
         %Pressure 
         Pnd_2(:,:,n) = (gamma-1)*(U3_2(:,:,n)-(U2_2(:,:,n).^2)./(2*U1_2(:,:,n)) );
         Pnd_2_nu2(:,:,n) = (gamma-1)*(U3_2_nu2(:,:,n)-(U2_2_nu2(:,:,n).^2)./(2*U1_2_nu2(:,:,n)) );

         %density is U1
         %get speed of sound c
         cS_2(:,:,n) = sqrt(gamma*Pnd_2(:,:,n)./U1_2(:,:,n));
         cS_2_nu2(:,:,n) = sqrt(gamma*Pnd_2_nu2(:,:,n)./U1_2_nu2(:,:,n));


         %Populate P and P_inv matrices for each x location
    for i = 1:IL(2)
        
        Pvec_2(:,:,i) = [1 1/(2*cS_2(:,i,n)^2) 1/(2*cS_2(:,i,n)^2);...
            U2_2(:,i,n)/U1_2(:,i,n) ((U2_2(:,i,n)/U1_2(:,i,n))+cS_2(:,i,n))/(2*cS_2(:,i,n)^2) ((U2_2(:,i,n)/U1_2(:,i,n))-cS_2(:,i,n))/(2*cS_2(:,i,n)^2);...
            (1/2)*(U2_2(:,i,n)/U1_2(:,i,n))^2 (1/(4*cS_2(:,i,n)^2))*(U2_2(:,i,n)/U1_2(:,i,n))^2+(U2_2(:,i,n)/U1_2(:,i,n))/(2*cS_2(:,i,n))+1/(2*(gamma-1)) (1/(4*cS_2(:,i,n)^2))*(U2_2(:,i,n)/U1_2(:,i,n))^2-(U2_2(:,i,n)/U1_2(:,i,n))/(2*cS_2(:,i,n))+1/(2*(gamma-1));];

        Pvec_inv_2(:,:,i) = [1-(gamma-1)*(1/(2*cS_2(:,i,n)^2))*(U2_2(:,i,n)/U1_2(:,i,n))^2 (gamma-1)*(U2_2(:,i,n)/U1_2(:,i,n))/cS_2(:,i,n)^2 -(gamma-1)/cS_2(:,i,n)^2;...
            -(U2_2(:,i,n)/U1_2(:,i,n))*cS_2(:,i,n)+(gamma-1)*(1/2)*(U2_2(:,i,n)/U1_2(:,i,n))^2 cS_2(:,i,n)-(gamma-1)*(U2_2(:,i,n)/U1_2(:,i,n)) (gamma-1);...
             (U2_2(:,i,n)/U1_2(:,i,n))*cS_2(:,i,n)+(gamma-1)*(1/2)*(U2_2(:,i,n)/U1_2(:,i,n))^2 -cS_2(:,i,n)-(gamma-1)*(U2_2(:,i,n)/U1_2(:,i,n)) (gamma-1);]; 
         
        Pvec_2_nu2(:,:,i) = [1 1/(2*cS_2_nu2(:,i,n)^2) 1/(2*cS_2_nu2(:,i,n)^2);...
            U2_2_nu2(:,i,n)/U1_2_nu2(:,i,n) ((U2_2_nu2(:,i,n)/U1_2_nu2(:,i,n))+cS_2_nu2(:,i,n))/(2*cS_2_nu2(:,i,n)^2) ((U2_2_nu2(:,i,n)/U1_2_nu2(:,i,n))-cS_2_nu2(:,i,n))/(2*cS_2_nu2(:,i,n)^2);...
            (1/2)*(U2_2_nu2(:,i,n)/U1_2_nu2(:,i,n))^2 (1/(4*cS_2_nu2(:,i,n)^2))*(U2_2_nu2(:,i,n)/U1_2_nu2(:,i,n))^2+(U2_2_nu2(:,i,n)/U1_2_nu2(:,i,n))/(2*cS_2_nu2(:,i,n))+1/(2*(gamma-1)) (1/(4*cS_2_nu2(:,i,n)^2))*(U2_2_nu2(:,i,n)/U1_2_nu2(:,i,n))^2-(U2_2_nu2(:,i,n)/U1_2_nu2(:,i,n))/(2*cS_2_nu2(:,i,n))+1/(2*(gamma-1));];

        Pvec_inv_2_nu2(:,:,i) = [1-(gamma-1)*(1/(2*cS_2_nu2(:,i,n)^2))*(U2_2_nu2(:,i,n)/U1_2_nu2(:,i,n))^2 (gamma-1)*(U2_2_nu2(:,i,n)/U1_2_nu2(:,i,n))/cS_2_nu2(:,i,n)^2 -(gamma-1)/cS_2_nu2(:,i,n)^2;...
            -(U2_2_nu2(:,i,n)/U1_2_nu2(:,i,n))*cS_2_nu2(:,i,n)+(gamma-1)*(1/2)*(U2_2_nu2(:,i,n)/U1_2_nu2(:,i,n))^2 cS_2_nu2(:,i,n)-(gamma-1)*(U2_2_nu2(:,i,n)/U1_2_nu2(:,i,n)) (gamma-1);...
             (U2_2_nu2(:,i,n)/U1_2_nu2(:,i,n))*cS_2_nu2(:,i,n)+(gamma-1)*(1/2)*(U2_2_nu2(:,i,n)/U1_2_nu2(:,i,n))^2 -cS_2_nu2(:,i,n)-(gamma-1)*(U2_2_nu2(:,i,n)/U1_2_nu2(:,i,n)) (gamma-1);]; 
         
    end
   
    
        for i = 1:IL(2)
            
             %conservative vector
             U_2(:,:,i) = [U1_2(:,i,n); U2_2(:,i,n); U3_2(:,i,n);];
             U_2_nu2(:,:,i) = [U1_2_nu2(:,i,n); U2_2_nu2(:,i,n); U3_2_nu2(:,i,n);];
            
                
             if i < IL(2)
                %i+1 conservative vector needed for Ei+1 
                U_2(:,:,i+1) = [U1_2(:,i+1,n); U2_2(:,i+1,n); U3_2(:,i+1,n);];
                U_2_nu2(:,:,i+1) = [U1_2_nu2(:,i+1,n); U2_2_nu2(:,i+1,n); U3_2_nu2(:,i+1,n);];
                
             end               
            
            %populate positive and negative eigenvalue matrices at each x
            %location
            
            %%E+ and E- at location i
            for d = 1:3
                eigenP_2(d,d,i) = 0.5*(egn_2(d,:,i)+abs(egn_2(d,:,i)));
                eigenN_2(d,d,i) = 0.5*(egn_2(d,:,i)-abs(egn_2(d,:,i)));

                eigenP_2_nu2(d,d,i) = 0.5*(egn_2_nu2(d,:,i)+abs(egn_2_nu2(d,:,i)));
                eigenN_2_nu2(d,d,i) = 0.5*(egn_2_nu2(d,:,i)-abs(egn_2_nu2(d,:,i)));
            end
            
             %calculate E+ and E- at i
                W_2(:,:,i) = Pvec_inv_2(:,:,i)*U_2(:,:,i);
                
                Wp_2(:,:,i) = eigenP_2(:,:,i)*W_2(:,:,i);
                Wn_2(:,:,i) = eigenN_2(:,:,i)*W_2(:,:,i);
                
                Ep_2(:,:,i) = Pvec_2(:,:,i)*Wp_2(:,:,i);
                En_2(:,:,i) = Pvec_2(:,:,i)*Wn_2(:,:,i);
                
                W_2_nu2(:,:,i) = Pvec_inv_2_nu2(:,:,i)*U_2_nu2(:,:,i);
                
                Wp_2_nu2(:,:,i) = eigenP_2_nu2(:,:,i)*W_2_nu2(:,:,i);
                Wn_2_nu2(:,:,i) = eigenN_2_nu2(:,:,i)*W_2_nu2(:,:,i);
                
                Ep_2_nu2(:,:,i) = Pvec_2_nu2(:,:,i)*Wp_2_nu2(:,:,i);
                En_2_nu2(:,:,i) = Pvec_2_nu2(:,:,i)*Wn_2_nu2(:,:,i);

            %apply dirichlet boundary conditions
            if i == 1
                
                U1_2(:,i,n+1) =  rho_L;
                U2_2(:,i,n+1) = rho_L*u_L;
                U3_2(:,i,n+1) = p_L/(gamma-1)+rho_L/2*(u_L^2);

                U1_2_nu2(:,i,n+1) =  rho_L;
                U2_2_nu2(:,i,n+1) = rho_L*u_L;
                U3_2_nu2(:,i,n+1) = p_L/(gamma-1)+rho_L/2*(u_L^2);
                
            elseif i == IL(2)
                
                U1_2(:,i,n+1) =  rho_R;
                U2_2(:,i,n+1) = rho_R*u_R;
                U3_2(:,i,n+1) = p_R/(gamma-1)+rho_R/2*(u_R^2);

                U1_2_nu2(:,i,n+1) =  rho_R;
                U2_2_nu2(:,i,n+1) = rho_R*u_R;
                U3_2_nu2(:,i,n+1) = p_R/(gamma-1)+rho_R/2*(u_R^2);
                
            else
             
             %% E+ and E- at i+1
            for d = 1:3
                eigenP_2(d,d,i+1) = 0.5*(egn_2(d,:,i+1)+abs(egn_2(d,:,i+1)));
                eigenN_2(d,d,i+1) = 0.5*(egn_2(d,:,i+1)-abs(egn_2(d,:,i+1)));  

                eigenP_2_nu2(d,d,i+1) = 0.5*(egn_2_nu2(d,:,i+1)+abs(egn_2_nu2(d,:,i+1)));
                eigenN_2_nu2(d,d,i+1) = 0.5*(egn_2_nu2(d,:,i+1)-abs(egn_2_nu2(d,:,i+1)));  
            end
                                            
                
                %Calculate E- at i+1
                W_2(:,:,i+1) = Pvec_inv_2(:,:,i+1)*U_2(:,:,i+1); %W = (P^-1)U
                Wn_2(:,:,i+1) = eigenN_2(:,:,i+1)*W_2(:,:,i+1); %Wn = (lambda^-)W
                En_2(:,:,i+1) = Pvec_2(:,:,i+1)*Wn_2(:,:,i+1);% En = E^- = (P)Wn

                W_2_nu2(:,:,i+1) = Pvec_inv_2_nu2(:,:,i+1)*U_2_nu2(:,:,i+1); %W = (P^-1)U
                Wn_2_nu2(:,:,i+1) = eigenN_2_nu2(:,:,i+1)*W_2_nu2(:,:,i+1); %Wn = (lambda^-)W
                En_2_nu2(:,:,i+1) = Pvec_2_nu2(:,:,i+1)*Wn_2_nu2(:,:,i+1);% En = E^- = (P)Wn
                
                %calculation
                U1_2(:,i,n+1) = U1_2(:,i,n) - (dt_2(n)/(delta_x_2))*( (Ep_2(1,:,i)-Ep_2(1,:,i-1))+(En_2(1,:,i+1)-En_2(1,:,i))  );
                U2_2(:,i,n+1) = U2_2(:,i,n) - (dt_2(n)/(delta_x_2))*( (Ep_2(2,:,i)-Ep_2(2,:,i-1))+(En_2(2,:,i+1)-En_2(2,:,i))  );
                U3_2(:,i,n+1) = U3_2(:,i,n) - (dt_2(n)/(delta_x_2))*( (Ep_2(3,:,i)-Ep_2(3,:,i-1))+(En_2(3,:,i+1)-En_2(3,:,i))  );

                U1_2_nu2(:,i,n+1) = U1_2_nu2(:,i,n) - (dt_2_nu2(n)/(delta_x_2))*( (Ep_2_nu2(1,:,i)-Ep_2_nu2(1,:,i-1))+(En_2_nu2(1,:,i+1)-En_2_nu2(1,:,i))  );
                U2_2_nu2(:,i,n+1) = U2_2_nu2(:,i,n) - (dt_2_nu2(n)/(delta_x_2))*( (Ep_2_nu2(2,:,i)-Ep_2_nu2(2,:,i-1))+(En_2_nu2(2,:,i+1)-En_2_nu2(2,:,i))  );
                U3_2_nu2(:,i,n+1) = U3_2_nu2(:,i,n) - (dt_2_nu2(n)/(delta_x_2))*( (Ep_2_nu2(3,:,i)-Ep_2_nu2(3,:,i-1))+(En_2_nu2(3,:,i+1)-En_2_nu2(3,:,i))  );
         
            end
            
        end  
end


%analytical for SW grid 1
for k = 1:length(nu)
    %Note: U_anlyt stores rho, u, and p (not the normal U parameters of
    %rho, rho*u, and e)
    for i = 1:IL(1)
        if x_arr_1(i) >= u_S * t_total_SW_1(k)
            U_SW_1_anlyt(i, 1, k) = rho_R;
            U_SW_1_anlyt(i, 2, k) = u_R;
            U_SW_1_anlyt(i, 3, k) = p_R;
            
        elseif (x_arr_1(i) > u_2 * t_total_SW_1(k)) && (x_arr_1(i) < u_S * t_total_SW_1(k)) 
            U_SW_1_anlyt(i, 1, k) = rho_2;
            U_SW_1_anlyt(i, 2, k) = u_2;
            U_SW_1_anlyt(i, 3, k) = p_2;

        elseif (x_arr_1(i) > (u_3 - c_3) * t_total_SW_1(k)) && (x_arr_1(i) < u_3 * t_total_SW_1(k))
            U_SW_1_anlyt(i, 1, k) = rho_3;
            U_SW_1_anlyt(i, 2, k) = u_3;
            U_SW_1_anlyt(i, 3, k) = p_3;

        elseif (x_arr_1(i) > (u_L - c_L) * t_total_SW_1(k)) && (x_arr_1(i) < (u_3 - c_3) * t_total_SW_1(k))
            u_5 = 2 / (gamma + 1) * ( (x_arr_1(i) / t_total_SW_1(k)) + c_L + (gamma - 1) / 2 * u_L);
            c_5 = u_5 - x_arr_1(i) / t_total_SW_1(k);
            p_5 = p_L * (c_5 / c_L)^((2 * gamma) / (gamma - 1));
            rho_5 = (p_5 / p_L)^(1 / gamma) * rho_L;
            U_SW_1_anlyt(i, 1, k) = rho_5;
            U_SW_1_anlyt(i, 2, k) = u_5;
            U_SW_1_anlyt(i, 3, k) = p_5;

        else
            U_SW_1_anlyt(i, 1, k) = rho_L;
            U_SW_1_anlyt(i, 2, k) = u_L;
            U_SW_1_anlyt(i, 3, k) = p_L;
        end
    end
end

%analytical for SW grid 2
for k = 1:length(nu)
    %Note: U_anlyt stores rho, u, and p (not the normal U parameters of
    %rho, rho*u, and e)
    for i = 1:IL(2)
        if x_arr_2(i) >= u_S * t_total_SW_2(k)
            U_SW_2_anlyt(i, 1, k) = rho_R;
            U_SW_2_anlyt(i, 2, k) = u_R;
            U_SW_2_anlyt(i, 3, k) = p_R;
            
        elseif (x_arr_2(i) > u_2 * t_total_SW_2(k)) && (x_arr_2(i) < u_S * t_total_SW_2(k)) 
            U_SW_2_anlyt(i, 1, k) = rho_2;
            U_SW_2_anlyt(i, 2, k) = u_2;
            U_SW_2_anlyt(i, 3, k) = p_2;

        elseif (x_arr_2(i) > (u_3 - c_3) * t_total_SW_2(k)) && (x_arr_2(i) < u_3 * t_total_SW_2(k))
            U_SW_2_anlyt(i, 1, k) = rho_3;
            U_SW_2_anlyt(i, 2, k) = u_3;
            U_SW_2_anlyt(i, 3, k) = p_3;

        elseif (x_arr_2(i) > (u_L - c_L) * t_total_SW_2(k)) && (x_arr_2(i) < (u_3 - c_3) * t_total_SW_2(k))
            u_5 = 2 / (gamma + 1) * ( (x_arr_2(i) / t_total_SW_2(k)) + c_L + (gamma - 1) / 2 * u_L);
            c_5 = u_5 - x_arr_2(i) / t_total_SW_2(k);
            p_5 = p_L * (c_5 / c_L)^((2 * gamma) / (gamma - 1));
            rho_5 = (p_5 / p_L)^(1 / gamma) * rho_L;
            U_SW_2_anlyt(i, 1, k) = rho_5;
            U_SW_2_anlyt(i, 2, k) = u_5;
            U_SW_2_anlyt(i, 3, k) = p_5;

        else
            U_SW_2_anlyt(i, 1, k) = rho_L;
            U_SW_2_anlyt(i, 2, k) = u_L;
            U_SW_2_anlyt(i, 3, k) = p_L;
        end
    end
end


%% Figure/Plots
figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_1, U_MC_1(:,1, end, 1), 69,'b', 'd');
scatter(x_arr_1, U_MC_1(:,1,end, 2), 69,'g', '*');
plot(x_arr_1, U_MC_1_anlyt(:, 1, 2), '-g', 'LineWidth', 2);
plot(x_arr_1, U_MC_1_anlyt(:, 1, 1), '-b', 'LineWidth', 2);
title('MacCormack Scheme with Artificial Dissipation','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(1)) ' at Time $t = $ ' num2str(t_total_MC_1(1,1)) ' s'],'Interpreter','LaTeX') 
ylim([0 1]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Density $\rho(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$',  '$\rho_{\mathrm{exact}}(x)$ for $\nu=1.3$', '$\rho_{\mathrm{exact}}(x)$ for $\nu=0.8$','Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_1, U_MC_1(:, 2, end, 1) ./ U_MC_1(:, 1, end, 1), 69,'b', 'd');
scatter(x_arr_1, U_MC_1(:, 2, end, 2) ./ U_MC_1(:, 1, end, 2), 69,'g', '*');
plot(x_arr_1, U_MC_1_anlyt(:, 2, 2), '-g', 'LineWidth', 2);
plot(x_arr_1, U_MC_1_anlyt(:, 2, 1), '-b', 'LineWidth', 2);
title('MacCormack Scheme with Artificial Dissipation','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(1)) ' at Time $t = $ ' num2str(t_total_MC_1(1,1)) ' s'],'Interpreter','LaTeX')
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Velocity $u(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$',  '$u_{\mathrm{exact}}(x)$ for $\nu=1.3$','$u_{\mathrm{exact}}(x)$ for $\nu=0.8$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);
ylim([0 500])

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_1, (gamma - 1) .* (U_MC_1(:, 3, end, 1) - (U_MC_1(:, 2, end, 1)).^2 ./ (2 .* U_MC_1(:, 1, end, 1))), 69,'b', 'd');
scatter(x_arr_1, (gamma - 1) .* (U_MC_1(:, 3, end, 2) - (U_MC_1(:, 2, end, 2)).^2 ./ (2 .* U_MC_1(:, 1, end, 2))), 69,'g', '*');
plot(x_arr_1, U_MC_1_anlyt(:, 3, 2), '-g', 'LineWidth', 2);
plot(x_arr_1, U_MC_1_anlyt(:, 3, 1), '-b', 'LineWidth', 2);
title('MacCormack Scheme with Artificial Dissipation','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(1)) ' at Time $t = $ ' num2str(t_total_MC_1(1,1)) ' s'],'Interpreter','LaTeX')
ylim([p_R p_L]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Pressure $p(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$',  '$p_{\mathrm{exact}}(x)$ for $\nu=1.3$','$p_{\mathrm{exact}}(x)$ for $\nu=0.8$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_2, U_MC_2(:,1, end, 1), 69,'b', 'd');
scatter(x_arr_2, U_MC_2(:,1,end, 2), 69,'g', '*');
plot(x_arr_2, U_MC_2_anlyt(:, 1, 2), '-g', 'LineWidth', 2);
plot(x_arr_2, U_MC_2_anlyt(:, 1, 1), '-b', 'LineWidth', 2);
title('MacCormack Scheme with Artificial Dissipation','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(2)) ' at Time $t = $ ' num2str(t_total_MC_2(1,1)) ' s'],'Interpreter','LaTeX') 
ylim([0 1]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Density $\rho(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$',  '$\rho_{\mathrm{exact}}(x)$ for $\nu=1.3$', '$\rho_{\mathrm{exact}}(x)$ for $\nu=0.8$','Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_2, U_MC_2(:, 2, end, 1) ./ U_MC_2(:, 1, end, 1), 69,'b', 'd');
scatter(x_arr_2, U_MC_2(:, 2, end, 2) ./ U_MC_2(:, 1, end, 2), 69,'g', '*');
plot(x_arr_2, U_MC_2_anlyt(:, 2, 2), '-g', 'LineWidth', 2);
plot(x_arr_2, U_MC_2_anlyt(:, 2, 1), '-b', 'LineWidth', 2);
title('MacCormack Scheme with Artificial Dissipation','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(2)) ' at Time $t = $ ' num2str(t_total_MC_2(1,1)) ' s'],'Interpreter','LaTeX')
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Velocity $u(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$',  '$u_{\mathrm{exact}}(x)$ for $\nu=1.3$','$u_{\mathrm{exact}}(x)$ for $\nu=0.8$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);
ylim([0 500])

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_2, (gamma - 1) .* (U_MC_2(:, 3, end, 1) - (U_MC_2(:, 2, end, 1)).^2 ./ (2 .* U_MC_2(:, 1, end, 1))), 69,'b', 'd');
scatter(x_arr_2, (gamma - 1) .* (U_MC_2(:, 3, end, 2) - (U_MC_2(:, 2, end, 2)).^2 ./ (2 .* U_MC_2(:, 1, end, 2))), 69,'g', '*');
plot(x_arr_2, U_MC_2_anlyt(:, 3, 2), '-g', 'LineWidth', 2);
plot(x_arr_2, U_MC_2_anlyt(:, 3, 1), '-b', 'LineWidth', 2);
title('MacCormack Scheme with Artificial Dissipation','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(2)) ' at Time $t = $ ' num2str(t_total_MC_2(1,1)) ' s'],'Interpreter','LaTeX')
ylim([p_R p_L]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Pressure $p(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$',  '$p_{\mathrm{exact}}(x)$ for $\nu=1.3$', '$p_{\mathrm{exact}}(x)$ for $\nu=0.8$','Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_1, U_LW_1(:,1, end, 1), 69,'b', 'd');
scatter(x_arr_1, U_LW_1(:,1, end, 2), 69,'g', '*');
plot(x_arr_1, U_LW_1_anlyt(:, 1, 2), '-g', 'LineWidth', 2);
plot(x_arr_1, U_LW_1_anlyt(:, 1, 1), '-b', 'LineWidth', 2);
title('Lax-Wendroff Scheme with Artificial Dissipation','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(1)) ' at Time $t = $ ' num2str(t_total_LW_1(1,1)) ' s'],'Interpreter','LaTeX') 
ylim([0 1]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Density $\rho(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$',  '$\rho_{\mathrm{exact}}(x)$ for $\nu=1.3$', '$\rho_{\mathrm{exact}}(x)$ for $\nu=0.8$','Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_1, U_LW_1(:, 2, end, 1) ./ U_LW_1(:, 1, end, 1), 69,'b', 'd');
scatter(x_arr_1, U_LW_1(:, 2, end, 2) ./ U_LW_1(:, 1, end, 2), 69,'g', '*');
plot(x_arr_1, U_LW_1_anlyt(:, 2, 2), '-g', 'LineWidth', 2);
plot(x_arr_1, U_LW_1_anlyt(:, 2, 1), '-b', 'LineWidth', 2);
title('Lax-Wendroff Scheme with Artificial Dissipation','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(1)) ' at Time $t = $ ' num2str(t_total_LW_1(1,1)) ' s'],'Interpreter','LaTeX')
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Velocity $u(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$',  '$u_{\mathrm{exact}}(x)$ for $\nu=1.3$', '$u_{\mathrm{exact}}(x)$ for $\nu=0.8$','Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);
ylim([0 500])

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_1, (gamma - 1) .* (U_LW_1(:, 3, end, 1) - (U_LW_1(:, 2, end, 1)).^2 ./ (2 .* U_LW_1(:, 1, end, 1))), 69,'b', 'd');
scatter(x_arr_1, (gamma - 1) .* (U_LW_1(:, 3, end, 2) - (U_LW_1(:, 2, end, 2)).^2 ./ (2 .* U_LW_1(:, 1, end, 2))), 69,'g', '*');
plot(x_arr_1, U_LW_1_anlyt(:, 3, 2), '-g', 'LineWidth', 2);
plot(x_arr_1, U_LW_1_anlyt(:, 3, 1), '-b', 'LineWidth', 2);
title('Lax-Wendroff Scheme with Artificial Dissipation','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(1)) ' at Time $t = $ ' num2str(t_total_LW_1(1,1)) ' s'],'Interpreter','LaTeX')
ylim([p_R p_L]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Pressure $p(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$',  '$p_{\mathrm{exact}}(x)$ for $\nu=1.3$','$p_{\mathrm{exact}}(x)$ for $\nu=0.8$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_2, U_LW_2(:,1, end, 1), 69,'b', 'd');
scatter(x_arr_2, U_LW_2(:,1, end, 2), 69,'g', '*');
plot(x_arr_2, U_LW_2_anlyt(:, 1, 2), '-g', 'LineWidth', 2);
plot(x_arr_2, U_LW_2_anlyt(:, 1, 1), '-b', 'LineWidth', 2);
title('Lax-Wendroff Scheme with Artificial Dissipation','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(2)) ' at Time $t = $ ' num2str(t_total_LW_2(1,1)) ' s'],'Interpreter','LaTeX') 
ylim([0 1]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Density $\rho(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$',  '$\rho_{\mathrm{exact}}(x)$ for $\nu=1.3$', '$\rho_{\mathrm{exact}}(x)$ for $\nu=0.8$','Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_2, U_LW_2(:, 2, end, 1) ./ U_LW_2(:, 1, end, 1), 69,'b', 'd');
scatter(x_arr_2, U_LW_2(:, 2, end, 2) ./ U_LW_2(:, 1, end, 2), 69,'g', '*');
plot(x_arr_2, U_LW_2_anlyt(:, 2, 2), '-g', 'LineWidth', 2);
plot(x_arr_2, U_LW_2_anlyt(:, 2, 1), '-b', 'LineWidth', 2);
title('Lax-Wendroff Scheme with Artificial Dissipation','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(2)) ' at Time $t = $ ' num2str(t_total_LW_2(1,1)) ' s'],'Interpreter','LaTeX')
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Velocity $u(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$',  '$u_{\mathrm{exact}}(x)$ for $\nu=1.3$','$u_{\mathrm{exact}}(x)$ for $\nu=0.8$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);
ylim([0 500])

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_2, (gamma - 1) .* (U_LW_2(:, 3, end, 1) - (U_LW_2(:, 2, end, 1)).^2 ./ (2 .* U_LW_2(:, 1, end, 1))), 69,'b', 'd');
scatter(x_arr_2, (gamma - 1) .* (U_LW_2(:, 3, end, 2) - (U_LW_2(:, 2, end, 2)).^2 ./ (2 .* U_LW_2(:, 1, end, 2))), 69,'g', '*');
plot(x_arr_2, U_LW_2_anlyt(:, 3, 2), '-g', 'LineWidth', 2);
plot(x_arr_2, U_LW_2_anlyt(:, 3, 1), '-b', 'LineWidth', 2);
title('Lax-Wendroff Scheme with Artificial Dissipation','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(2)) ' at Time $t = $ ' num2str(t_total_LW_2(1,1)) ' s'],'Interpreter','LaTeX')
ylim([p_R p_L]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Pressure $p(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$',  '$p_{\mathrm{exact}}(x)$ for $\nu=1.3$','$p_{\mathrm{exact}}(x)$ for $\nu=0.8$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);


figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_1, U1(:,:,end), 69,'b', 'd');
scatter(x_arr_1, U1_nu2(:,:,end), 69,'g', '*');
plot(x_arr_1, U_SW_1_anlyt(:, 1, 2), '-g', 'LineWidth', 2);
plot(x_arr_1, U_SW_1_anlyt(:, 1, 1), '-b', 'LineWidth', 2);
title('Steger-Warming Flux Vector Splitting Scheme','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(1)) ' at Time $t = $ ' num2str(t_total_SW_1(1,1)) ' s'],'Interpreter','LaTeX') 
ylim([0 1]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Density $\rho(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$','$\rho_{\mathrm{exact}}(x)$ for $\nu=1.3$','$\rho_{\mathrm{exact}}(x)$ for $\nu=0.8$','Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_1, U2(:,:,end)./U1(:,:,end), 69,'b', 'd');
scatter(x_arr_1, U2_nu2(:,:,end)./U1_nu2(:,:,end), 69,'g', '*');
plot(x_arr_1, U_SW_1_anlyt(:, 2, 2), '-g', 'LineWidth', 2);
plot(x_arr_1, U_SW_1_anlyt(:, 2, 1), '-b', 'LineWidth', 2);
title('Steger-Warming Flux Vector Splitting Scheme','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(1)) ' at Time $t = $ ' num2str(t_total_SW_1(1,1)) ' s'],'Interpreter','LaTeX')
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Velocity $u(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$','$u_{\mathrm{exact}}(x)$ for $\nu=1.3$','$u_{\mathrm{exact}}(x)$ for $\nu=0.8$','Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);
ylim([0 500])

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_1, (gamma - 1) .* ( U3(:,:,end) - ( U2(:,:,end)).^2 ./ (2 .*  U1(:,:,end))), 69,'b', 'd');
scatter(x_arr_1, (gamma - 1) .* ( U3_nu2(:,:,end) - ( U2_nu2(:,:,end)).^2 ./ (2 .*  U1_nu2(:,:,end))), 69,'g', '*');
plot(x_arr_1, U_SW_1_anlyt(:, 3, 2), '-g', 'LineWidth', 2);
plot(x_arr_1, U_SW_1_anlyt(:, 3, 1), '-b', 'LineWidth', 2);
title('Steger-Warming Flux Vector Splitting Scheme','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(1)) ' at Time $t = $ ' num2str(t_total_SW_1(1,1)) ' s'],'Interpreter','LaTeX')
ylim([p_R p_L]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Pressure $p(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$','$p_{\mathrm{exact}}(x)$ for $\nu=1.3$','$p_{\mathrm{exact}}(x)$ for $\nu=0.8$','Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_2, U1_2(:,:,end), 69,'b', 'd');
scatter(x_arr_2, U1_2_nu2(:,:,end), 69,'g', '*');
plot(x_arr_2, U_SW_2_anlyt(:, 1, 2), '-g', 'LineWidth', 2);
plot(x_arr_2, U_SW_2_anlyt(:, 1, 1), '-b', 'LineWidth', 2);
title('Steger-Warming Flux Vector Splitting Scheme','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(2)) ' at Time $t = $ ' num2str(t_total_SW_2(1,1)) ' s'],'Interpreter','LaTeX') 
ylim([0 1]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Density $\rho(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$','$\rho_{\mathrm{exact}}(x)$ for $\nu=1.3$','$\rho_{\mathrm{exact}}(x)$ for $\nu=0.8$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_2, U2_2(:,:,end)./U1_2(:,:,end), 69,'b', 'd');
scatter(x_arr_2, U2_2_nu2(:,:,end)./U1_2_nu2(:,:,end), 69,'g', '*');
plot(x_arr_2, U_SW_2_anlyt(:, 2, 2), '-g', 'LineWidth', 2);
plot(x_arr_2, U_SW_2_anlyt(:, 2, 1), '-b', 'LineWidth', 2);
title('Steger-Warming Flux Vector Splitting Scheme','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(2)) ' at Time $t = $ ' num2str(t_total_SW_2(1,1)) ' s'],'Interpreter','LaTeX')
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Velocity $u(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$','$u_{\mathrm{exact}}(x)$ for $\nu=1.3$','$u_{\mathrm{exact}}(x)$ for $\nu=0.8$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);
ylim([0 500])

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_2, (gamma - 1) .* ( U3_2(:,:,end) - ( U2_2(:,:,end)).^2 ./ (2 .*  U1_2(:,:,end))), 69,'b', 'd');
scatter(x_arr_2, (gamma - 1) .* ( U3_2_nu2(:,:,end) - ( U2_2_nu2(:,:,end)).^2 ./ (2 .*  U1_2_nu2(:,:,end))), 69,'g', '*');
plot(x_arr_2, U_SW_2_anlyt(:, 3, 2), '-g', 'LineWidth', 2);
plot(x_arr_2, U_SW_2_anlyt(:, 3, 1), '-b', 'LineWidth', 2);
title('Steger-Warming Flux Vector Splitting Scheme','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(2)) ' at Time $t = $ ' num2str(t_total_SW_2(1,1)) ' s'],'Interpreter','LaTeX')
ylim([p_R p_L]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Pressure $p(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$','$p_{\mathrm{exact}}(x)$ for $\nu=1.3$','$p_{\mathrm{exact}}(x)$ for $\nu=0.8$',  'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);
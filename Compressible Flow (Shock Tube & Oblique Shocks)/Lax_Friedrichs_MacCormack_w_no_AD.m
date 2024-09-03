%==========================================================================
% AUTHOR: David L. Tran

% DESCRIPTION: Computes the numerical and analytical solution of Sod's
% shock tube using the Lax-Friedrichs method and MacCormack method for
% specified conditions.
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
P = 3.031;          %pressure ratio p_2/p_R
alpha = (gamma + 1) / (gamma - 1);
                    
%grid properties
delta_x = 6.25e-2;  %grid spacing in [m]
IL = 300;           %number of grid cells
nu_1 = 0.8;         %first Courant number
nu_2 = 1.3;         %second Courant number
nu = [nu_1; nu_2];  %array of Courant numbers
t_steps = 70;       %number of time steps
x_0 = 0;            %center point
t_total_LF = zeros(length(nu), 1);   
                    %time elapsed in Lax-Friedrichs simulation in [s]
t_total_MC = zeros(length(nu), 1);  

%spatial grids
x_right = linspace(0+delta_x/2, (0+delta_x/2 * IL), IL/2);
x_left = linspace(0-delta_x/2, (0-delta_x/2 * IL), IL/2);
x_arr = [flip(x_left), x_right];

%U and E arrays
U_LF = zeros(IL, 3, t_steps+1, length(nu));
E_LF = zeros(IL, 3, t_steps+1, length(nu));
U_anlyt_LF = zeros(IL, 3, length(nu));
E_anlyt_LF = zeros(IL, 3, length(nu));

U_MC = zeros(IL, 3, t_steps+1, length(nu));
E_MC = zeros(IL, 3, t_steps+1, length(nu));
U_MC_p = zeros(IL, 3, t_steps+1, length(nu));
U_anlyt_MC = zeros(IL, 3, length(nu));
E_anlyt_MC = zeros(IL, 3, length(nu));



%% Loops

for k = 1:length(nu)
    for i = 1:IL
        if x_arr(i) < 0
            %left of diaphragm
            U_LF(i, 1, 1, k) = rho_L;
            U_LF(i, 2, 1, k) = rho_L * u_L;
            U_LF(i, 3, 1, k) = p_L / (gamma - 1) + (rho_L) / 2 * u_L^2;
            E_LF(i, 1, 1, k) = U_LF(i, 2, 1, k);
            E_LF(i, 2, 1, k) = (U_LF(i, 2, 1, k)^2 / U_LF(i, 1, 1, k)) + (gamma - 1) * (U_LF(i, 3, 1, k) - (U_LF(i, 2, 1, k)^2) / (2 * U_LF(i,1, 1, k)));
            E_LF(i, 3, 1, k) = (U_LF(i, 2, 1, k) * U_LF(i, 3, 1, k)) / (U_LF(i, 1, 1, k)) + (U_LF(i, 2, 1, k) / U_LF(i, 1, 1, k)) * (gamma - 1) * (U_LF(i, 3, 1, k) - (U_LF(i, 2, 1, k)^2) / (2 * U_LF(i, 1, 1, k)));
            c(i) = sqrt( (gamma * (gamma - 1) * (U_LF(i, 3, 1, k) - (U_LF(i, 1, 1, k) / 2) * (U_LF(i, 2, 1, k)^2 / U_LF(i, 1, 1, k)^2))) / U_LF(i, 1, 1, k));
            u_plus_c(i) = abs(U_LF(i, 2, 1, k) / U_LF(i, 1, 1, k)) + c(i);

            U_MC(i, 1, 1, k) = rho_L;
            U_MC(i, 2, 1, k) = rho_L * u_L;
            U_MC(i, 3, 1, k) = p_L / (gamma - 1) + (rho_L) / 2 * u_L^2;
            E_MC(i, 1, 1, k) = U_MC(i, 2, 1, k);
            E_MC(i, 2, 1, k) = (U_MC(i, 2, 1, k)^2 / U_MC(i, 1, 1, k)) + (gamma - 1) * (U_MC(i, 3, 1, k) - (U_MC(i, 2, 1, k)^2) / (2 * U_MC(i,1, 1, k)));
            E_MC(i, 3, 1, k) = (U_MC(i, 2, 1, k) * U_MC(i, 3, 1, k)) / (U_MC(i, 1, 1, k)) + (U_MC(i, 2, 1, k) / U_MC(i, 1, 1, k)) * (gamma - 1) * (U_MC(i, 3, 1, k) - (U_MC(i, 2, 1, k)^2) / (2 * U_MC(i, 1, 1, k)));
            c_MC(i) = sqrt( (gamma * (gamma - 1) * (U_MC(i, 3, 1, k) - (U_MC(i, 1, 1, k) / 2) * (U_MC(i, 2, 1, k)^2 / U_MC(i, 1, 1, k)^2))) / U_MC(i, 1, 1, k));
            u_plus_c_MC(i) = abs(U_MC(i, 2, 1, k) / U_MC(i, 1, 1, k)) + c_MC(i);
        else
            %right of diaphragm
            U_LF(i, 1, 1, k) = rho_R;
            U_LF(i, 2, 1, k) = rho_R * u_R;
            U_LF(i, 3, 1, k) = p_R / (gamma - 1) + (rho_R) / 2 * u_R^2;
            E_LF(i, 1, 1, k) = U_LF(i, 2, 1, k);
            E_LF(i, 2, 1, k) = (U_LF(i, 2, 1, k)^2 / U_LF(i, 1, 1, k)) + (gamma - 1) * (U_LF(i, 3, 1, k) - (U_LF(i, 2, 1, k)^2) / (2 * U_LF(i, 1, 1, k)));
            E_LF(i, 3, 1, k) = (U_LF(i, 2, 1, k) * U_LF(i, 3, 1, k)) / (U_LF(i, 1, 1, k)) + (U_LF(i, 2, 1, k) / U_LF(i, 1, 1, k)) * (gamma - 1) * (U_LF(i, 3, 1, k) - (U_LF(i, 2, 1, k)^2) / (2 * U_LF(i, 1, 1, k)));
            c(i) = sqrt( (gamma * (gamma - 1) * (U_LF(i, 3, 1, k) - (U_LF(i, 1, 1, k) / 2) * (U_LF(i, 2, 1, k)^2 / U_LF(i, 1, 1, k)^2))) / U_LF(i, 1, 1, k));
            u_plus_c(i) = abs(U_LF(i, 2, 1, k) / U_LF(i, 1, 1, k)) + c(i);

            U_MC(i, 1, 1, k) = rho_R;
            U_MC(i, 2, 1, k) = rho_R * u_R;
            U_MC(i, 3, 1, k) = p_R / (gamma - 1) + (rho_R) / 2 * u_R^2;
            E_MC(i, 1, 1, k) = U_MC(i, 2, 1, k);
            E_MC(i, 2, 1, k) = (U_MC(i, 2, 1, k)^2 / U_MC(i, 1, 1, k)) + (gamma - 1) * (U_MC(i, 3, 1, k) - (U_MC(i, 2, 1, k)^2) / (2 * U_MC(i,1, 1, k)));
            E_MC(i, 3, 1, k) = (U_MC(i, 2, 1, k) * U_MC(i, 3, 1, k)) / (U_MC(i, 1, 1, k)) + (U_MC(i, 2, 1, k) / U_MC(i, 1, 1, k)) * (gamma - 1) * (U_MC(i, 3, 1, k) - (U_MC(i, 2, 1, k)^2) / (2 * U_MC(i, 1, 1, k)));
            c_MC(i) = sqrt( (gamma * (gamma - 1) * (U_MC(i, 3, 1, k) - (U_MC(i, 1, 1, k) / 2) * (U_MC(i, 2, 1, k)^2 / U_MC(i, 1, 1, k)^2))) / U_MC(i, 1, 1, k));
            u_plus_c_MC(i) = abs(U_MC(i, 2, 1, k) / U_MC(i, 1, 1, k)) + c_MC(i);
        end
    end
    delta_t = nu(k) * delta_x /  max(u_plus_c);
    t_total_LF(k) = delta_t;

    delta_t_MC = nu(k) * delta_x / max(u_plus_c_MC);
    t_total_MC(k) = delta_t_MC;

    %Assign BCs for all of time
    U_LF(1, 1, :, k) = rho_L;
    U_LF(1, 2, :, k) = U_LF(1, 2, 1, k);
    U_LF(1, 3, :, k) = U_LF(1, 3, 1, k);
    U_LF(end, 1, :, k) = rho_R;
    U_LF(end, 2, :, k) = U_LF(end, 2, 1, k);
    U_LF(end, 3, :, k) = U_LF(end, 3, 1, k);
    E_LF(1, 1, :, k) = E_LF(1, 1, 1, k);
    E_LF(1, 2, :, k) = E_LF(1, 2, 1, k);
    E_LF(1, 3, :, k) = E_LF(1, 3, 1, k);
    E_LF(end, 1, :, k) = E_LF(end, 1, 1, k);
    E_LF(end, 2, :, k) = E_LF(end, 2, 1, k);
    E_LF(end, 3, :, k) = E_LF(end, 3, 1, k);
    
    U_MC(1, 1, :, k) = rho_L;
    U_MC(1, 2, :, k) = U_MC(1, 2, 1, k);
    U_MC(1, 3, :, k) = U_MC(1, 3, 1, k);
    U_MC(end, 1, :, k) = rho_R;
    U_MC(end, 2, :, k) = U_MC(end, 2, 1, k);
    U_MC(end, 3, :, k) = U_MC(end, 3, 1, k);
    U_MC_p(1, 1, :, k) = rho_L; 
    U_MC_p(1, 2, :, k) = U_MC(1, 2, 1, k); 
    U_MC_p(1, 3, :, k) = U_MC(1, 3, 1, k);
    U_MC_p(end, 1, :, k) = rho_R;
    U_MC_p(end, 2, :, k) = U_MC(end, 2, 1, k); 
    U_MC_p(end, 3, :, k) = U_MC(end, 3, 1, k); 
    E_MC(1, 1, :, k) = E_MC(1, 1, 1, k);
    E_MC(1, 2, :, k) = E_MC(1, 2, 1, k);
    E_MC(1, 3, :, k) = E_MC(1, 3, 1, k);
    E_MC(end, 1, :, k) = E_MC(end, 1, 1, k);
    E_MC(end, 2, :, k) = E_MC(end, 2, 1, k);
    E_MC(end, 3, :, k) = E_MC(end, 3, 1, k);
    E_MC_p(1, 1, :, k) = E_MC(1, 1, :, k);
    E_MC_p(1, 2, :, k) = E_MC(1, 2, :, k);
    E_MC_p(1, 3, :, k) = E_MC(1, 3, :, k);
    E_MC_p(end, 1, :, k) = E_MC(end, 1, 1, k);
    E_MC_p(end, 2, :, k) = E_MC(end, 2, 1, k);
    E_MC_p(end, 3, :, k) = E_MC(end, 3, 1, k);


    for n = 2:t_steps+1
        for j = 2:IL-1
            %Lax-Friedrichs
            U_LF(j, 1, n, k) = U_LF(j, 1, n-1, k) - delta_t / (2 * delta_x) * (E_LF(j+1, 1, n-1, k) - E_LF(j-1, 1, n-1, k)) + 0.5 * (U_LF(j+1, 1, n-1, k) - 2 * U_LF(j, 1, n-1, k) + U_LF(j-1, 1, n-1, k));
            U_LF(j, 2, n, k) = U_LF(j, 2, n-1, k) - delta_t / (2 * delta_x) * (E_LF(j+1, 2, n-1, k) - E_LF(j-1, 2, n-1, k)) + 0.5 * (U_LF(j+1, 2, n-1, k) - 2 * U_LF(j, 2, n-1, k) + U_LF(j-1, 2, n-1, k));
            U_LF(j, 3, n, k) = U_LF(j, 3, n-1, k) - delta_t / (2 * delta_x) * (E_LF(j+1, 3, n-1, k) - E_LF(j-1, 3, n-1, k)) + 0.5 * (U_LF(j+1, 3, n-1, k) - 2 * U_LF(j, 3, n-1, k) + U_LF(j-1, 3, n-1, k));
            E_LF(j, 1, n, k) = U_LF(j, 2, n, k);
            E_LF(j, 2, n, k) = (U_LF(j, 2, n, k)^2 / U_LF(j, 1, n, k)) + (gamma - 1) * (U_LF(j, 3, n, k) - (U_LF(j, 2, n, k)^2 / (2 * U_LF(j, 1, n, k))) );
            E_LF(j, 3, n, k) = (U_LF(j, 2, n, k) * U_LF(j, 3, n, k)) / (U_LF(j, 1, n, k)) + ((U_LF(j, 2, n, k)) / (U_LF(j, 1, n, k))) * ((gamma - 1) * (U_LF(j, 3, n, k) - U_LF(j, 2, n, k)^2 / (2 * U_LF(j, 1, n, k))) );
            c(j) = sqrt( (gamma * (gamma - 1) * (U_LF(j, 3, n, k) - (U_LF(j, 1, n, k) / 2) * (U_LF(j, 2, n, k)^2 / U_LF(j, 1, n, k)^2))) / U_LF(j, 1, n, k));
            u_plus_c(j) = abs(U_LF(j, 2, n, k) / U_LF(j, 1, n, k)) + c(j);

            %MacCormack
            U_MC_p(j, 1, n, k) = U_MC(j, 1, n-1, k) - (delta_t_MC / delta_x) * (E_MC(j+1, 1, n-1, k) - E_MC(j, 1, n-1, k));
            U_MC_p(j, 2, n, k) = U_MC(j, 2, n-1, k) - (delta_t_MC / delta_x) * (E_MC(j+1, 2, n-1, k) - E_MC(j, 2, n-1, k));
            U_MC_p(j, 3, n, k) = U_MC(j, 3, n-1, k) - (delta_t_MC / delta_x) * (E_MC(j+1, 3, n-1, k) - E_MC(j, 3, n-1, k));

            E_MC_p(j, 1, n, k) = U_MC_p(j, 2, n, k);
            E_MC_p(j, 2, n, k) = (U_MC_p(j, 2, n, k)^2 / U_MC_p(j, 1, n, k)) + (gamma - 1) * (U_MC_p(j, 3, n, k) - (U_MC_p(j, 2, n, k)^2 / (2*U_MC_p(j, 1, n, k))));
            E_MC_p(j, 3, n, k) = U_MC_p(j, 2, n, k) * U_MC_p(j, 3, n, k) / U_MC_p(j, 1, n, k) + (U_MC_p(j, 2, n, k) / U_MC_p(j, 1, n, k)) * ((gamma - 1) * (U_MC_p(j, 3, n, k) - (U_MC_p(j, 2, n, k)^2 / (2*U_MC_p(j, 1, n, k)))));

            U_MC(j, 1, n, k) = 0.5 * (U_MC(j, 1, n-1, k) + U_MC_p(j, 1, n, k) - (delta_t_MC / delta_x) * (E_MC_p(j, 1, n, k) - E_MC_p(j-1, 1, n, k)) );
            U_MC(j, 2, n, k) = 0.5 * (U_MC(j, 2, n-1, k) + U_MC_p(j, 2, n, k) - (delta_t_MC / delta_x) * (E_MC_p(j, 2, n, k) - E_MC_p(j-1, 2, n, k)) );
            U_MC(j, 3, n, k) = 0.5 * (U_MC(j, 3, n-1, k) + U_MC_p(j, 3, n, k) - (delta_t_MC / delta_x) * (E_MC_p(j, 3, n, k) - E_MC_p(j-1, 3, n, k)) );

            if U_MC(j, 2, n, k) < 0
                U_MC(j, 1, n, k) = rho_L;
                U_MC(j, 2, n, k) = rho_L * u_L;
                U_MC(j, 3, n, k) = U_MC(j, 3, n-1, k);
            end
            

            E_MC(j, 1, n, k) = U_MC(j, 2, n, k);
            E_MC(j, 2, n, k) = (U_MC(j, 2, n, k)^2 / U_MC(j, 1, n, k)) + (gamma - 1) * (U_MC(j, 3, n, k) - (U_MC(j, 2, n, k)^2 / (2*U_MC(j, 1, n, k))));
            E_MC(j, 3, n, k) = U_MC(j, 2, n, k) * U_MC(j, 3, n, k) / U_MC(j, 1, n, k) + (U_MC(j, 2, n, k) / U_MC(j, 1, n, k)) * ((gamma - 1) * (U_MC(j, 3, n, k) - (U_MC(j, 2, n, k)^2 / (2*U_MC(j, 1, n, k)))));
            
            c_MC(j) = sqrt(gamma * ((gamma - 1) * (U_MC(j, 3, n, k) - (U_MC(j, 2, n, k))^2 / (2 * U_MC(j, 1, n, k)))) / U_MC(j, 1, n, k));
            u_plus_c_MC(j) = abs(U_MC(j, 2, n, k) / U_MC(j, 1, n, k)) + real(c_MC(j));
        end
        delta_t = nu(k) * delta_x / max(u_plus_c);
        t_total_LF(k) = t_total_LF(k) + delta_t;

        delta_t_MC = nu(k) * delta_x / max(u_plus_c_MC);
        t_total_MC(k) = t_total_MC(k) + delta_t_MC;

    end
end

%Analytical calculation
c_L = sqrt(gamma * p_L / rho_L);
x_L = (u_L - c_L) .* t_total_LF;

c_R = sqrt(gamma * p_R / rho_R);

p_2 = p_R * P;
rho_2 = rho_R * (1 + alpha * P) / (alpha + P);
u_S = c_R * sqrt((gamma - 1) / (2 * gamma) * (1 + alpha * P)) + u_R;
u_2 = (P - 1) / (sqrt(1 + alpha * P))*(1 / sqrt(gamma * (gamma - 1) / 2)) * c_R + u_R;

p_3 = p_2;
u_3 = u_2;
rho_3 = (p_3 / p_L)^(1 / gamma) * rho_L;
c_3 = sqrt(gamma * p_3 / rho_3);

%analytical for LF
for k = 1:length(nu)
    %Note: U_anlyt stores rho, u, and p (not the normal U parameters of
    %rho, rho*u, and e)
    for i = 1:IL
        if x_arr(i) >= u_S * t_total_LF(k)
            U_anlyt_LF(i, 1, k) = rho_R;
            U_anlyt_LF(i, 2, k) = u_R;
            U_anlyt_LF(i, 3, k) = p_R;
            
        elseif (x_arr(i) > u_2 * t_total_LF(k)) && (x_arr(i) < u_S * t_total_LF(k)) 
            U_anlyt_LF(i, 1, k) = rho_2;
            U_anlyt_LF(i, 2, k) = u_2;
            U_anlyt_LF(i, 3, k) = p_2;

        elseif (x_arr(i) > (u_3 - c_3) * t_total_LF(k)) && (x_arr(i) < u_3 * t_total_LF(k))
            U_anlyt_LF(i, 1, k) = rho_3;
            U_anlyt_LF(i, 2, k) = u_3;
            U_anlyt_LF(i, 3, k) = p_3;

        elseif (x_arr(i) > (u_L - c_L) * t_total_LF(k)) && (x_arr(i) < (u_3 - c_3) * t_total_LF(k))
            u_5 = 2 / (gamma + 1) * ( (x_arr(i) / t_total_LF(k)) + c_L + (gamma - 1) / 2 * u_L);
            c_5 = u_5 - x_arr(i) / t_total_LF(k);
            p_5 = p_L * (c_5 / c_L)^((2 * gamma) / (gamma - 1));
            rho_5 = (p_5 / p_L)^(1 / gamma) * rho_L;
            U_anlyt_LF(i, 1, k) = rho_5;
            U_anlyt_LF(i, 2, k) = u_5;
            U_anlyt_LF(i, 3, k) = p_5;

        else
            U_anlyt_LF(i, 1, k) = rho_L;
            U_anlyt_LF(i, 2, k) = u_L;
            U_anlyt_LF(i, 3, k) = p_L;
        end
    end
end

x_L_MC = (u_L - c_L) .* t_total_MC;
%analytical for MacCormack
for k = 1:length(nu)
    %Note: U_anlyt stores rho, u, and p (not the normal U parameters of
    %rho, rho*u, and e)
    for i = 1:IL
        if x_arr(i) >= u_S * t_total_MC(k)
            U_anlyt_MC(i, 1, k) = rho_R;
            U_anlyt_MC(i, 2, k) = u_R;
            U_anlyt_MC(i, 3, k) = p_R;
            
        elseif (x_arr(i) > u_2 * t_total_MC(k)) && (x_arr(i) < u_S * t_total_MC(k)) 
            U_anlyt_MC(i, 1, k) = rho_2;
            U_anlyt_MC(i, 2, k) = u_2;
            U_anlyt_MC(i, 3, k) = p_2;

        elseif (x_arr(i) > (u_3 - c_3) * t_total_MC(k)) && (x_arr(i) < u_3 * t_total_MC(k))
            U_anlyt_MC(i, 1, k) = rho_3;
            U_anlyt_MC(i, 2, k) = u_3;
            U_anlyt_MC(i, 3, k) = p_3;

        elseif (x_arr(i) > (u_L - c_L) * t_total_MC(k)) && (x_arr(i) < (u_3 - c_3) * t_total_MC(k))
            u_5 = 2 / (gamma + 1) * ( (x_arr(i) / t_total_MC(k)) + c_L + (gamma - 1) / 2 * u_L);
            c_5 = u_5 - x_arr(i) / t_total_MC(k);
            p_5 = p_L * (c_5 / c_L)^((2 * gamma) / (gamma - 1));
            rho_5 = (p_5 / p_L)^(1 / gamma) * rho_L;
            U_anlyt_MC(i, 1, k) = rho_5;
            U_anlyt_MC(i, 2, k) = u_5;
            U_anlyt_MC(i, 3, k) = p_5;

        else
            U_anlyt_MC(i, 1, k) = rho_L;
            U_anlyt_MC(i, 2, k) = u_L;
            U_anlyt_MC(i, 3, k) = p_L;
        end
    end
end


%% Plots/Figures
figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr, U_LF(:,1,end, 1), 69,'b', 'd');
scatter(x_arr, U_LF(:,1,end, 2), 69,'g', '*');
plot(x_arr, U_anlyt_LF(:, 1, 1), '-b', 'LineWidth', 2);
plot(x_arr, U_anlyt_LF(:, 1, 2), '-g', 'LineWidth', 2);
title("Lax-Friedrichs Scheme",'Interpreter','LaTeX');  
ylim([0 1]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Density $\rho(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$', '$\rho_{\mathrm{exact}}(x)$ for $\nu=0.8$', '$\rho_{\mathrm{exact}}(x)$ for $\nu=1.3$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr, U_LF(:,2,end, 1) ./ U_LF(:,1,end, 1), 69,'b', 'd');
scatter(x_arr, U_LF(:,2,end, 2) ./ U_LF(:,1,end, 2), 69,'g', '*');
plot(x_arr, U_anlyt_LF(:, 2, 1), '-b', 'LineWidth', 2);
plot(x_arr, U_anlyt_LF(:, 2, 2), '-g', 'LineWidth', 2);
title("Lax-Friedrichs Scheme",'Interpreter','LaTeX');               
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Velocity $u(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$', '$u_{\mathrm{exact}}(x)$ for $\nu=0.8$', '$u_{\mathrm{exact}}(x)$ for $\nu=1.3$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr, (gamma - 1) .* (U_LF(:, 3, end, 1) - (U_LF(:, 2, end, 1).^2) ./ (2 .* U_LF(:, 1, end, 1))), 69,'b', 'd');
scatter(x_arr, (gamma - 1) .* (U_LF(:, 3, end, 2) - (U_LF(:, 2, end, 2).^2) ./ (2 .* U_LF(:, 1, end, 2))), 69,'g', '*');
plot(x_arr, U_anlyt_LF(:, 3, 1), '-b', 'LineWidth', 2);
plot(x_arr, U_anlyt_LF(:, 3, 2), '-g', 'LineWidth', 2);
title("Lax-Friedrichs Scheme",'Interpreter','LaTeX');               
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Pressure $p(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$', '$p_{\mathrm{exact}}(x)$ for $\nu=0.8$', '$p_{\mathrm{exact}}(x)$ for $\nu=1.3$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr, U_MC(:,1,end, 1), 69,'b', 'd');
scatter(x_arr, U_MC(:,1,end, 2), 69,'g', '*');
plot(x_arr, U_anlyt_MC(:, 1, 1), '-b', 'LineWidth', 2);
plot(x_arr, U_anlyt_MC(:, 1, 2), '-g', 'LineWidth', 2);
ylim([0 1]);
title("MacCormack Scheme",'Interpreter','LaTeX');               
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Density $\rho(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$', '$\rho_{\mathrm{exact}}(x)$ for $\nu=0.8$', '$\rho_{\mathrm{exact}}(x)$ for $\nu=1.3$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr, U_MC(:,2,end, 1) ./ U_MC(:,1,end,1), 69,'b', 'd');
scatter(x_arr, U_MC(:,2,end, 2) ./ U_MC(:,1,end,2), 69,'g', '*');
plot(x_arr, U_anlyt_MC(:, 2, 1), '-b', 'LineWidth', 2);
plot(x_arr, U_anlyt_MC(:, 2, 2), '-g', 'LineWidth', 2);
ylim([0 500]);
title("MacCormack Scheme",'Interpreter','LaTeX');               
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Velocity $u(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$', '$u_{\mathrm{exact}}(x)$ for $\nu=0.8$', '$u_{\mathrm{exact}}(x)$ for $\nu=1.3$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr, (gamma - 1) .* (U_MC(:,3,end, 1) - (U_MC(:,2,end,1).^2 ./ (2 .* U_MC(:,1,end,1)))), 69,'b', 'd');
scatter(x_arr, (gamma - 1) .* (U_MC(:,3,end, 2) - (U_MC(:,2,end,2).^2 ./ (2 .* U_MC(:,1,end,2)))), 69,'g', '*');
plot(x_arr, U_anlyt_MC(:, 3, 1), '-b', 'LineWidth', 2);
plot(x_arr, U_anlyt_MC(:, 3, 2), '-g', 'LineWidth', 2);
ylim([p_R p_L]);
title("MacCormack Scheme",'Interpreter','LaTeX');               
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Pressure $p(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\nu=1.3$', '$p_{\mathrm{exact}}(x)$ for $\nu=0.8$', '$p_{\mathrm{exact}}(x)$ for $\nu=1.3$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

%==========================================================================
% AUTHOR: David L. Tran
%
% Roe's flux-difference splitting scheme ONLY
%
% DESCRIPTION: Numerically solves the conservative 1D Euler equation in a
% finite volume formulation using the Roe scheme. Initial conditions
% correspond to those of Sod's shock tube. This is a Riemann problem
% wherein sharp discontinuities are present in the exact solution.
%
%==========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                          %
%                       Clear Cache                        %
%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all; clearvars;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                          %
%                        Variables                         %
%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%CFL conditions
nu_1 = 0.8;         %first CFL number
nu_2 = 1.3;         %second CFL number
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

x0 = delta_x_1 / 2;
x_arr_1_left = -x0 - (0:(IL(1)/2)-1)*delta_x_1;
x_arr_1_right = x0 + (0:(IL(1)/2)-1)*delta_x_1;
x_arr_1 = [flip(x_arr_1_left), x_arr_1_right];

U_Roe1 = zeros(IL(1), t_steps(1), 3, length(nu));

t_total_1 = zeros(length(nu), 1);
delta_t_1 = zeros(length(nu), 1);

%GRID 2 ARRAYS (time steps = 140, IL = 600)
delta_x_2 = delta_x_1 / 2;        %grid 2 spacing in [m]

x0_2 = delta_x_2 / 2;
x_arr_2_left = -x0_2 - (0:(IL(2)/2)-1)*delta_x_2;
x_arr_2_right = x0_2 + (0:(IL(2)/2)-1)*delta_x_2;
x_arr_2 = [flip(x_arr_2_left), x_arr_2_right];

U_Roe2 = zeros(IL(2), t_steps(2), 3, length(nu));

t_total_2 = zeros(length(nu), 1);
delta_t_2 = zeros(length(nu), 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                          %
%                     Main Body/Loops                      %
%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Roe's Method for Grid 1
for k = 1:length(nu)
    %Initialize IC
    [U_Roe1(:,1,1,k),U_Roe1(:,1,2,k),U_Roe1(:,1,3,k)] = initcond(IL(1), x_arr_1, gamma, rho_L, p_L, u_L, rho_R, p_R, u_R);
    
    %Calculate speed of sound and |u|+c for each grid point
    c = sos(gamma, pressure(gamma, U_Roe1(:,1,3,k), U_Roe1(:,1,1,k), U_Roe1(:,1,2,k)./U_Roe1(:,1,1,k)), U_Roe1(:,1,1,k) );
    u_plus_c = abs(U_Roe1(:, 1, 2, k) ./ U_Roe1(:, 1, 1, k)) + c;

    %Calculate initial time step size
    delta_t_1(k,1) = nu(k) * delta_x_1 / max(u_plus_c);
    t_total_1(k,1) = t_total_1(k,1) + delta_t_1(k,1);

    for n = 1:t_steps(1)
        %Loop through all of time
        %Apply Dirichlet BCs
        U_Roe1(1,n+1,1,k) = rho_L;
        U_Roe1(1,n+1,2,k) = rho_L * u_L;
        U_Roe1(1,n+1,3,k) = p_L / (gamma - 1) + rho_L / 2 * u_L^2;
    
        U_Roe1(IL(1),n+1,1,k) = rho_R;
        U_Roe1(IL(1),n+1,2,k) = rho_R * u_R;
        U_Roe1(IL(1),n+1,3,k) = p_R / (gamma - 1) + rho_R / 2 * u_R^2;
    
        %Apply Dirichlet BCs for speed of sound and |u|+c
        c(1) = sos(gamma, p_L, rho_L);
        c(IL(1)) = sos(gamma, p_R, rho_R);
        u_plus_c(1) = abs(u_L) + c(1);
        u_plus_c(IL(1)) = abs(u_R) + c(end);

        for i = 2:IL(1)-1
            %Loop through all domain points
            R_i_phalf = sqrt(U_Roe1(i+1,n,1,k) / U_Roe1(i,n,1,k));
            R_i_mhalf = sqrt(U_Roe1(i,n,1,k) / U_Roe1(i-1,n,1,k));

            u_im1 = U_Roe1(i-1,n,2,k)/U_Roe1(i-1,n,1,k);
            u_i = U_Roe1(i,n,2,k)/U_Roe1(i,n,1,k);
            u_ip1 = U_Roe1(i+1,n,2,k)/U_Roe1(i+1,n,1,k);

            H_im1 = gamma ./ (gamma - 1) .* (pressure(gamma, U_Roe1(i-1,n,3,k), U_Roe1(i-1,n,1,k), u_im1) ./ U_Roe1(i-1,n,1,k)) + (U_Roe1(i-1,n,2,k).^2 ./ (2.*U_Roe1(i-1,n,1,k).^2));
            H_i = gamma ./ (gamma - 1) .* (pressure(gamma, U_Roe1(i,n,3,k), U_Roe1(i,n,1,k), u_i) ./ U_Roe1(i,n,1,k)) + (U_Roe1(i,n,2,k).^2 ./ (2.*U_Roe1(i,n,1,k).^2));
            H_ip1 = gamma ./ (gamma - 1) .* (pressure(gamma, U_Roe1(i+1,n,3,k), U_Roe1(i+1,n,1,k), u_ip1) ./ U_Roe1(i+1,n,1,k)) + (U_Roe1(i+1,n,2,k).^2 ./ (2.*U_Roe1(i+1,n,1,k).^2));


            %Compute Roe averages for i+1/2
            rho_avg = R_i_phalf * U_Roe1(i,n,1,k);
            u_avg = ((R_i_phalf*u_ip1 + u_i)) / (R_i_phalf + 1);
            H_avg = (R_i_phalf * H_ip1 + H_i) / (R_i_phalf + 1);
            c_avg = sqrt( (gamma - 1) * (H_avg - u_avg^2 / 2) );

            %Compute Roe averages for i-1/2
            rho_mavg = R_i_mhalf * U_Roe1(i-1,n,1,k);
            u_mavg = ( (R_i_mhalf*u_i + u_im1)) / (R_i_mhalf + 1);
            H_mavg = (R_i_mhalf * H_i + H_im1) / (R_i_mhalf + 1);
            c_mavg = sqrt( (gamma - 1) * (H_mavg - u_mavg^2 / 2) );

            %Compute matrix |Lambda_{i+1/2}|
            Lambda_i_half = [abs(u_avg),    0,               0;...
                              0,    abs(u_avg + c_avg),      0;...
                              0,         0,        abs(u_avg - c_avg)];

            %Compute matrix |Lambda_{i-1/2}|
            Lambda_i_mhalf = [abs(u_mavg),    0,               0;...
                              0,    abs(u_mavg + c_mavg),      0;...
                              0,         0,        abs(u_mavg - c_mavg)];
        
            %Compute matrix P_{i+1/2}
            P_i_half = [1,                1/(2*c_avg^2),                                        1/(2*c_avg^2);...
                        u_avg,      u_avg/(2*c_avg^2) + 1/(2*c_avg),                        u_avg/(2*c_avg^2)-1/(2*c_avg);...
                        u_avg^2/2, u_avg^2/(4*c_avg^2)+u_avg/(2*c_avg)+1/(2*(gamma-1)), u_avg^2/(4*c_avg^2)-u_avg/(2*c_avg)+1/(2*(gamma-1))];

            %Compute matrix P_{i-1/2}
            P_i_mhalf = [1,                1/(2*c_mavg^2),                                        1/(2*c_mavg^2);...
                        u_mavg,      u_mavg/(2*c_mavg^2) + 1/(2*c_mavg),                        u_mavg/(2*c_mavg^2)-1/(2*c_mavg);...
                        u_mavg^2/2, u_mavg^2/(4*c_mavg^2)+u_mavg/(2*c_mavg)+1/(2*(gamma-1)), u_mavg^2/(4*c_mavg^2)-u_mavg/(2*c_mavg)+1/(2*(gamma-1))];

            %Compute matrix P^{-1}_{i+1/2}
            P_inv_i_half = [(2*c_avg^2-gamma*u_avg^2+u_avg^2)/(2*c_avg^2),      u_avg*(gamma-1)/c_avg^2,        (1-gamma)/c_avg^2;...
                        0.5*u_avg*(u_avg*(gamma-1) - 2*c_avg),                 c_avg - gamma*u_avg+u_avg,           gamma-1;...
                        0.5*u_avg*(u_avg*(gamma-1) + 2*c_avg),                -c_avg - gamma*u_avg+u_avg ,          gamma-1];

            %Compute matrix P^{-1}_{i-1/2}
            P_inv_i_mhalf = [(2*c_mavg^2-gamma*u_mavg^2+u_mavg^2)/(2*c_mavg^2),      u_mavg*(gamma-1)/c_mavg^2,        (1-gamma)/c_mavg^2;...
                        0.5*u_mavg*(u_mavg*(gamma-1) - 2*c_mavg),                 c_mavg - gamma*u_mavg+u_mavg,           gamma-1;...
                        0.5*u_mavg*(u_mavg*(gamma-1) + 2*c_mavg),                -c_mavg - gamma*u_mavg+u_mavg ,          gamma-1];

            %Compute |A_{i+1/2}|
            A_i_half = P_i_half * Lambda_i_half;
            A_i_half = A_i_half * P_inv_i_half;

%             A_i_half_test = [0,    1,     0;...
%                             (gamma - 3)/2*u_avg^2,  (3-gamma)*u_avg,  gamma-1;...
%                             0.5*(gamma-1)*u_avg^3-u_avg*H_avg, H_avg-(gamma-1)*u_avg^2, gamma*u_avg];

            %Compute |A_{i-1/2}|
            A_i_mhalf = P_i_mhalf * Lambda_i_mhalf;
            A_i_mhalf = A_i_mhalf * P_inv_i_mhalf;      

%             A_i_mhalf_test = [0,    1,     0;...
%                             (gamma - 3)/2*u_mavg^2,  (3-gamma)*u_mavg,  gamma-1;...
%                             0.5*(gamma-1)*u_mavg^3-u_mavg*H_mavg, H_mavg-(gamma-1)*u_mavg^2, gamma*u_mavg];

            %Compute flux vectors E_{i-1}
            [Em1, Em2, Em3] = flux(U_Roe1(i-1,n,1,k), U_Roe1(i-1,n,2,k), U_Roe1(i-1,n,3,k), gamma);

            %Compute flux vectors E_i
            [E1, E2, E3] = flux(U_Roe1(i,n,1,k), U_Roe1(i,n,2,k), U_Roe1(i,n,3,k), gamma);
            
            %Compute flux vectors E_{i+1}
            [Ep1, Ep2, Ep3] = flux(U_Roe1(i+1,n,1,k), U_Roe1(i+1,n,2,k), U_Roe1(i+1,n,3,k), gamma);

            %Compute Roe flux vectors E*_{i-1/2}
            [Erm1, Erm2, Erm3] = roe_flux(U_Roe1(i-1,n,1,k), U_Roe1(i-1,n,2,k), U_Roe1(i-1,n,3,k), U_Roe1(i,n,1,k), U_Roe1(i,n,2,k), U_Roe1(i,n,3,k), Em1, Em2, Em3, E1, E2, E3, A_i_mhalf);

            %Compute Roe flux vectors E*_{i+1/2}
            [Erp1, Erp2, Erp3] = roe_flux(U_Roe1(i,n,1,k), U_Roe1(i,n,2,k), U_Roe1(i,n,3,k), U_Roe1(i+1,n,1,k), U_Roe1(i+1,n,2,k), U_Roe1(i+1,n,3,k), E1, E2, E3, Ep1, Ep2, Ep3, A_i_half);

            Erp = [Erp1; Erp2; Erp3];
            Erm = [Erm1; Erm2; Erm3];
            
            %Update conservative variables           
            U_Roe1(i,n+1,:,k) = squeeze(U_Roe1(i,n,:,k)) - (delta_t_1(k,1) / delta_x_1) * ( Erp - Erm );

%             if pressure(gamma, U_Roe1(i,n+1,3,k), U_Roe1(i,n+1,1,k), U_Roe1(i,n+1,2,k)./U_Roe1(i,n+1,1,k)) < 0
%                 fprintf('Negative found at %d', i);
%             end
            
            %Calculate speed of sound and |u|+c for each step
            c(i) = sos(gamma, pressure(gamma, U_Roe1(i,n+1,3,k), U_Roe1(i,n+1,1,k), U_Roe1(i,n+1,2,k)./U_Roe1(i,n+1,1,k)), U_Roe1(i,n+1,1,k) );
            u_plus_c(i) = abs(U_Roe1(i,n+1,2,k) / U_Roe1(i,n+1,1,k)) + c(i);
            
        end
        %Calculate variable time step size and update time count
        delta_t_1(k,1) = nu(k) * delta_x_1 / max(u_plus_c);
        t_total_1(k,1) = t_total_1(k,1) + delta_t_1(k,1);
        
    end

end

%analytical for SW grid 1
for k = 1:length(nu)
    %Note: U_anlyt stores rho, u, and p (not the normal U parameters of
    %rho, rho*u, and e)
    for i = 1:IL(1)
        if x_arr_1(i) >= u_S * t_total_1(k)
            U_Roe1_anlyt(i, 1, k) = rho_R;
            U_Roe1_anlyt(i, 2, k) = u_R;
            U_Roe1_anlyt(i, 3, k) = p_R;
            
        elseif (x_arr_1(i) > u_2 * t_total_1(k)) && (x_arr_1(i) < u_S * t_total_1(k)) 
            U_Roe1_anlyt(i, 1, k) = rho_2;
            U_Roe1_anlyt(i, 2, k) = u_2;
            U_Roe1_anlyt(i, 3, k) = p_2;

        elseif (x_arr_1(i) > (u_3 - c_3) * t_total_1(k)) && (x_arr_1(i) < u_3 * t_total_1(k))
            U_Roe1_anlyt(i, 1, k) = rho_3;
            U_Roe1_anlyt(i, 2, k) = u_3;
            U_Roe1_anlyt(i, 3, k) = p_3;

        elseif (x_arr_1(i) > (u_L - c_L) * t_total_1(k)) && (x_arr_1(i) < (u_3 - c_3) * t_total_1(k))
            u_5 = 2 / (gamma + 1) * ( (x_arr_1(i) / t_total_1(k)) + c_L + (gamma - 1) / 2 * u_L);
            c_5 = u_5 - x_arr_1(i) / t_total_1(k);
            p_5 = p_L * (c_5 / c_L)^((2 * gamma) / (gamma - 1));
            rho_5 = (p_5 / p_L)^(1 / gamma) * rho_L;
            U_Roe1_anlyt(i, 1, k) = rho_5;
            U_Roe1_anlyt(i, 2, k) = u_5;
            U_Roe1_anlyt(i, 3, k) = p_5;

        else
            U_Roe1_anlyt(i, 1, k) = rho_L;
            U_Roe1_anlyt(i, 2, k) = u_L;
            U_Roe1_anlyt(i, 3, k) = p_L;
        end
    end
end

%Roe's method for Grid 2
for k = 1:length(nu)
    %Initialize IC
    [U_Roe2(:,1,1,k),U_Roe2(:,1,2,k),U_Roe2(:,1,3,k)] = initcond(IL(2), x_arr_2, gamma, rho_L, p_L, u_L, rho_R, p_R, u_R);
    
    %Calculate speed of sound and |u|+c for each grid point
    c = sos(gamma, pressure(gamma, U_Roe2(:,1,3,k), U_Roe2(:,1,1,k), U_Roe2(:,1,2,k)./U_Roe2(:,1,1,k)), U_Roe2(:,1,1,k) );
    u_plus_c = abs(U_Roe2(:, 1, 2, k) ./ U_Roe2(:, 1, 1, k)) + c;

    %Calculate initial time step size
    delta_t_2(k,1) = nu(k) * delta_x_2 / max(u_plus_c);
    t_total_2(k,1) = t_total_2(k,1) + delta_t_2(k,1);

    for n = 1:t_steps(2)
        %Loop through all of time
        %Apply Dirichlet BCs
        U_Roe2(1,n+1,1,k) = rho_L;
        U_Roe2(1,n+1,2,k) = rho_L * u_L;
        U_Roe2(1,n+1,3,k) = p_L / (gamma - 1) + rho_L / 2 * u_L^2;
    
        U_Roe2(IL(2),n+1,1,k) = rho_R;
        U_Roe2(IL(2),n+1,2,k) = rho_R * u_R;
        U_Roe2(IL(2),n+1,3,k) = p_R / (gamma - 1) + rho_R / 2 * u_R^2;
    
        %Apply Dirichlet BCs for speed of sound and |u|+c
        c(1) = sos(gamma, p_L, rho_L);
        c(IL(2)) = sos(gamma, p_R, rho_R);
        u_plus_c(1) = abs(u_L) + c(1);
        u_plus_c(IL(2)) = abs(u_R) + c(end);

        for i = 2:IL(2)-1
            %Loop through all domain points
            R_i_phalf = sqrt(U_Roe2(i+1,n,1,k) / U_Roe2(i,n,1,k));
            R_i_mhalf = sqrt(U_Roe2(i,n,1,k) / U_Roe2(i-1,n,1,k));

            u_im1 = U_Roe2(i-1,n,2,k)/U_Roe2(i-1,n,1,k);
            u_i = U_Roe2(i,n,2,k)/U_Roe2(i,n,1,k);
            u_ip1 = U_Roe2(i+1,n,2,k)/U_Roe2(i+1,n,1,k);

            H_im1 = gamma ./ (gamma - 1) .* (pressure(gamma, U_Roe2(i-1,n,3,k), U_Roe2(i-1,n,1,k), u_im1) ./ U_Roe2(i-1,n,1,k)) + (U_Roe2(i-1,n,2,k).^2 ./ (2.*U_Roe2(i-1,n,1,k).^2));
            H_i = gamma ./ (gamma - 1) .* (pressure(gamma, U_Roe2(i,n,3,k), U_Roe2(i,n,1,k), u_i) ./ U_Roe2(i,n,1,k)) + (U_Roe2(i,n,2,k).^2 ./ (2.*U_Roe2(i,n,1,k).^2));
            H_ip1 = gamma ./ (gamma - 1) .* (pressure(gamma, U_Roe2(i+1,n,3,k), U_Roe2(i+1,n,1,k), u_ip1) ./ U_Roe2(i+1,n,1,k)) + (U_Roe2(i+1,n,2,k).^2 ./ (2.*U_Roe2(i+1,n,1,k).^2));


            %Compute Roe averages for i+1/2
            rho_avg = R_i_phalf * U_Roe2(i,n,1,k);
            u_avg = ((R_i_phalf*u_ip1 + u_i)) / (R_i_phalf + 1);
            H_avg = (R_i_phalf * H_ip1 + H_i) / (R_i_phalf + 1);
            c_avg = sqrt( (gamma - 1) * (H_avg - u_avg^2 / 2) );

            %Compute Roe averages for i-1/2
            rho_mavg = R_i_mhalf * U_Roe2(i-1,n,1,k);
            u_mavg = ( (R_i_mhalf*u_i + u_im1)) / (R_i_mhalf + 1);
            H_mavg = (R_i_mhalf * H_i + H_im1) / (R_i_mhalf + 1);
            c_mavg = sqrt( (gamma - 1) * (H_mavg - u_mavg^2 / 2) );

            %Compute matrix |Lambda_{i+1/2}|
            Lambda_i_half = [abs(u_avg),    0,               0;...
                              0,    abs(u_avg + c_avg),      0;...
                              0,         0,        abs(u_avg - c_avg)];

            %Compute matrix |Lambda_{i-1/2}|
            Lambda_i_mhalf = [abs(u_mavg),    0,               0;...
                              0,    abs(u_mavg + c_mavg),      0;...
                              0,         0,        abs(u_mavg - c_mavg)];
        
            %Compute matrix P_{i+1/2}
            P_i_half = [1,                1/(2*c_avg^2),                                        1/(2*c_avg^2);...
                        u_avg,      u_avg/(2*c_avg^2) + 1/(2*c_avg),                        u_avg/(2*c_avg^2)-1/(2*c_avg);...
                        u_avg^2/2, u_avg^2/(4*c_avg^2)+u_avg/(2*c_avg)+1/(2*(gamma-1)), u_avg^2/(4*c_avg^2)-u_avg/(2*c_avg)+1/(2*(gamma-1))];

            %Compute matrix P_{i-1/2}
            P_i_mhalf = [1,                1/(2*c_mavg^2),                                        1/(2*c_mavg^2);...
                        u_mavg,      u_mavg/(2*c_mavg^2) + 1/(2*c_mavg),                        u_mavg/(2*c_mavg^2)-1/(2*c_mavg);...
                        u_mavg^2/2, u_mavg^2/(4*c_mavg^2)+u_mavg/(2*c_mavg)+1/(2*(gamma-1)), u_mavg^2/(4*c_mavg^2)-u_mavg/(2*c_mavg)+1/(2*(gamma-1))];

            %Compute matrix P^{-1}_{i+1/2}
            P_inv_i_half = [(2*c_avg^2-gamma*u_avg^2+u_avg^2)/(2*c_avg^2),      u_avg*(gamma-1)/c_avg^2,        (1-gamma)/c_avg^2;...
                        0.5*u_avg*(u_avg*(gamma-1) - 2*c_avg),                 c_avg - gamma*u_avg+u_avg,           gamma-1;...
                        0.5*u_avg*(u_avg*(gamma-1) + 2*c_avg),                -c_avg - gamma*u_avg+u_avg ,          gamma-1];

            %Compute matrix P^{-1}_{i-1/2}
            P_inv_i_mhalf = [(2*c_mavg^2-gamma*u_mavg^2+u_mavg^2)/(2*c_mavg^2),      u_mavg*(gamma-1)/c_mavg^2,        (1-gamma)/c_mavg^2;...
                        0.5*u_mavg*(u_mavg*(gamma-1) - 2*c_mavg),                 c_mavg - gamma*u_mavg+u_mavg,           gamma-1;...
                        0.5*u_mavg*(u_mavg*(gamma-1) + 2*c_mavg),                -c_mavg - gamma*u_mavg+u_mavg ,          gamma-1];

            %Compute |A_{i+1/2}|
            A_i_half = P_i_half * Lambda_i_half;
            A_i_half = A_i_half * P_inv_i_half;

%             A_i_half_test = [0,    1,     0;...
%                             (gamma - 3)/2*u_avg^2,  (3-gamma)*u_avg,  gamma-1;...
%                             0.5*(gamma-1)*u_avg^3-u_avg*H_avg, H_avg-(gamma-1)*u_avg^2, gamma*u_avg];

            %Compute |A_{i-1/2}|
            A_i_mhalf = P_i_mhalf * Lambda_i_mhalf;
            A_i_mhalf = A_i_mhalf * P_inv_i_mhalf;      

%             A_i_mhalf_test = [0,    1,     0;...
%                             (gamma - 3)/2*u_mavg^2,  (3-gamma)*u_mavg,  gamma-1;...
%                             0.5*(gamma-1)*u_mavg^3-u_mavg*H_mavg, H_mavg-(gamma-1)*u_mavg^2, gamma*u_mavg];

            %Compute flux vectors E_{i-1}
            [Em1, Em2, Em3] = flux(U_Roe2(i-1,n,1,k), U_Roe2(i-1,n,2,k), U_Roe2(i-1,n,3,k), gamma);

            %Compute flux vectors E_i
            [E1, E2, E3] = flux(U_Roe2(i,n,1,k), U_Roe2(i,n,2,k), U_Roe2(i,n,3,k), gamma);
            
            %Compute flux vectors E_{i+1}
            [Ep1, Ep2, Ep3] = flux(U_Roe2(i+1,n,1,k), U_Roe2(i+1,n,2,k), U_Roe2(i+1,n,3,k), gamma);

            %Compute Roe flux vectors E*_{i-1/2}
            [Erm1, Erm2, Erm3] = roe_flux(U_Roe2(i-1,n,1,k), U_Roe2(i-1,n,2,k), U_Roe2(i-1,n,3,k), U_Roe2(i,n,1,k), U_Roe2(i,n,2,k), U_Roe2(i,n,3,k), Em1, Em2, Em3, E1, E2, E3, A_i_mhalf);

            %Compute Roe flux vectors E*_{i+1/2}
            [Erp1, Erp2, Erp3] = roe_flux(U_Roe2(i,n,1,k), U_Roe2(i,n,2,k), U_Roe2(i,n,3,k), U_Roe2(i+1,n,1,k), U_Roe2(i+1,n,2,k), U_Roe2(i+1,n,3,k), E1, E2, E3, Ep1, Ep2, Ep3, A_i_half);

            Erp = [Erp1; Erp2; Erp3];
            Erm = [Erm1; Erm2; Erm3];
            
            %Update conservative variables           
            U_Roe2(i,n+1,:,k) = squeeze(U_Roe2(i,n,:,k)) - (delta_t_2(k,1) / delta_x_2) * ( Erp - Erm );

%             if pressure(gamma, U_Roe1(i,n+1,3,k), U_Roe1(i,n+1,1,k), U_Roe1(i,n+1,2,k)./U_Roe1(i,n+1,1,k)) < 0
%                 fprintf('Negative found at %d', i);
%             end
            
            %Calculate speed of sound and |u|+c for each step
            c(i) = sos(gamma, pressure(gamma, U_Roe2(i,n+1,3,k), U_Roe2(i,n+1,1,k), U_Roe2(i,n+1,2,k)./U_Roe2(i,n+1,1,k)), U_Roe2(i,n+1,1,k) );
            u_plus_c(i) = abs(U_Roe2(i,n+1,2,k) / U_Roe2(i,n+1,1,k)) + c(i);
            
        end
        %Calculate variable time step size and update time count
        delta_t_2(k,1) = nu(k) * delta_x_2 / max(u_plus_c);
        t_total_2(k,1) = t_total_2(k,1) + delta_t_2(k,1);
        
    end

end

%analytical for SW grid 2
for k = 1:length(nu)
    %Note: U_anlyt stores rho, u, and p (not the normal U parameters of
    %rho, rho*u, and e)
    for i = 1:IL(2)
        if x_arr_2(i) >= u_S * t_total_2(k)
            U_Roe2_anlyt(i, 1, k) = rho_R;
            U_Roe2_anlyt(i, 2, k) = u_R;
            U_Roe2_anlyt(i, 3, k) = p_R;
            
        elseif (x_arr_2(i) > u_2 * t_total_2(k)) && (x_arr_2(i) < u_S * t_total_2(k)) 
            U_Roe2_anlyt(i, 1, k) = rho_2;
            U_Roe2_anlyt(i, 2, k) = u_2;
            U_Roe2_anlyt(i, 3, k) = p_2;

        elseif (x_arr_2(i) > (u_3 - c_3) * t_total_2(k)) && (x_arr_2(i) < u_3 * t_total_2(k))
            U_Roe2_anlyt(i, 1, k) = rho_3;
            U_Roe2_anlyt(i, 2, k) = u_3;
            U_Roe2_anlyt(i, 3, k) = p_3;

        elseif (x_arr_2(i) > (u_L - c_L) * t_total_2(k)) && (x_arr_2(i) < (u_3 - c_3) * t_total_2(k))
            u_5 = 2 / (gamma + 1) * ( (x_arr_2(i) / t_total_2(k)) + c_L + (gamma - 1) / 2 * u_L);
            c_5 = u_5 - x_arr_2(i) / t_total_2(k);
            p_5 = p_L * (c_5 / c_L)^((2 * gamma) / (gamma - 1));
            rho_5 = (p_5 / p_L)^(1 / gamma) * rho_L;
            U_Roe2_anlyt(i, 1, k) = rho_5;
            U_Roe2_anlyt(i, 2, k) = u_5;
            U_Roe2_anlyt(i, 3, k) = p_5;

        else
            U_Roe2_anlyt(i, 1, k) = rho_L;
            U_Roe2_anlyt(i, 2, k) = u_L;
            U_Roe2_anlyt(i, 3, k) = p_L;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                          %
%                      Plots/Figures                       %
%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Roe Grid 1
figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_1, U_Roe1(:,end, 1, 1), 69,'b', 'd');
plot(x_arr_1, U_Roe1_anlyt(:, 1, 1), '-b', 'LineWidth', 2);
title('Roe Flux-Difference Splitting Scheme','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(1)) ' at Time $t = $ ' num2str(t_total_1(1,1)) ' s'],'Interpreter','LaTeX') 
ylim([0 1]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Density $\rho(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\rho_{\mathrm{exact}}(x)$ for $\nu=0.8$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_1, U_Roe1(:,end, 1, 2), 69,'g', '*');
plot(x_arr_1, U_Roe1_anlyt(:, 1, 2), '-g', 'LineWidth', 2);
title('Roe Flux-Difference Splitting Scheme','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(1)) ' at Time $t = $ ' num2str(t_total_1(2,1)) ' s'],'Interpreter','LaTeX') 
ylim([0 1]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Density $\rho(x)$','Interpreter','LaTeX');  
legend('$\nu=1.3$','$\rho_{\mathrm{exact}}(x)$ for $\nu=1.3$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_1, U_Roe1(:,end,2,1)./U_Roe1(:,end,1,1), 69,'b', 'd');
plot(x_arr_1, U_Roe1_anlyt(:, 2, 1), '-b', 'LineWidth', 2);
title('Roe Flux-Difference Splitting Scheme','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(1)) ' at Time $t = $ ' num2str(t_total_1(1,1)) ' s'],'Interpreter','LaTeX')
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Velocity $u(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$u_{\mathrm{exact}}(x)$ for $\nu=0.8$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);
ylim([-50 500]);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_1, U_Roe1(:,end,2,2)./U_Roe1(:,end,1,2), 69,'g', '*');
plot(x_arr_1, U_Roe1_anlyt(:, 2, 2), '-g', 'LineWidth', 2);
title('Roe Flux-Difference Splitting Scheme','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(1)) ' at Time $t = $ ' num2str(t_total_1(2,1)) ' s'],'Interpreter','LaTeX')
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Velocity $u(x)$','Interpreter','LaTeX');  
legend('$\nu=1.3$','$u_{\mathrm{exact}}(x)$ for $\nu=1.3$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);
ylim([-50 500]);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_1, pressure(gamma, U_Roe1(:,end,3,1), U_Roe1(:,end,1,1), U_Roe1(:,end,2,1)./U_Roe1(:,end,1,1)), 69,'b', 'd');
plot(x_arr_1, U_Roe1_anlyt(:, 3, 1), '-b', 'LineWidth', 2);
title('Roe Flux-Difference Splitting Scheme','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(1)) ' at Time $t = $ ' num2str(t_total_1(1,1)) ' s'],'Interpreter','LaTeX')
ylim([0 p_L]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Pressure $p(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$p_{\mathrm{exact}}(x)$ for $\nu=0.8$',  'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_1, pressure(gamma, U_Roe1(:,end,3,2), U_Roe1(:,end,1,2), U_Roe1(:,end,2,2)./U_Roe1(:,end,1,2)), 69,'g', '*');
plot(x_arr_1, U_Roe1_anlyt(:, 3, 2), '-g', 'LineWidth', 2);
title('Roe Flux-Difference Splitting Scheme','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(1)) ' at Time $t = $ ' num2str(t_total_1(2,1)) ' s'],'Interpreter','LaTeX')
ylim([0 p_L]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Pressure $p(x)$','Interpreter','LaTeX');  
legend('$\nu=1.3$','$p_{\mathrm{exact}}(x)$ for $\nu=1.3$',  'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

%Roe Grid 2
figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_2, U_Roe2(:,end, 1, 1), 69,'b', 'd');
plot(x_arr_2, U_Roe2_anlyt(:, 1, 1), '-b', 'LineWidth', 2);
title('Roe Flux-Difference Splitting Scheme','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(2)) ' at Time $t = $ ' num2str(t_total_2(1,1)) ' s'],'Interpreter','LaTeX') 
ylim([0 1]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Density $\rho(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$\rho_{\mathrm{exact}}(x)$ for $\nu=0.8$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_2, U_Roe2(:,end, 1, 2), 69,'g', '*');
plot(x_arr_2, U_Roe2_anlyt(:, 1, 2), '-g', 'LineWidth', 2);
title('Roe Flux-Difference Splitting Scheme','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(2)) ' at Time $t = $ ' num2str(t_total_2(2,1)) ' s'],'Interpreter','LaTeX') 
ylim([0 1]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Density $\rho(x)$','Interpreter','LaTeX');  
legend('$\nu=1.3$','$\rho_{\mathrm{exact}}(x)$ for $\nu=1.3$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_2, U_Roe2(:,end,2,1)./U_Roe2(:,end,1,1), 69,'b', 'd');
plot(x_arr_2, U_Roe2_anlyt(:, 2, 1), '-b', 'LineWidth', 2);
title('Roe Flux-Difference Splitting Scheme','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(2)) ' at Time $t = $ ' num2str(t_total_2(1,1)) ' s'],'Interpreter','LaTeX')
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Velocity $u(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$u_{\mathrm{exact}}(x)$ for $\nu=0.8$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);
ylim([-50 500]);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_2, U_Roe2(:,end,2,2)./U_Roe2(:,end,1,2), 69,'g', '*');
plot(x_arr_2, U_Roe2_anlyt(:, 2, 2), '-g', 'LineWidth', 2);
title('Roe Flux-Difference Splitting Scheme','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(2)) ' at Time $t = $ ' num2str(t_total_2(2,1)) ' s'],'Interpreter','LaTeX')
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Velocity $u(x)$','Interpreter','LaTeX');  
legend('$\nu=1.3$','$u_{\mathrm{exact}}(x)$ for $\nu=1.3$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);
ylim([-100 500]);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_2, pressure(gamma, U_Roe2(:,end,3,1), U_Roe2(:,end,1,1), U_Roe2(:,end,2,1)./U_Roe2(:,end,1,1)), 69,'b', 'd');
plot(x_arr_2, U_Roe2_anlyt(:, 3, 1), '-b', 'LineWidth', 2);
title('Roe Flux-Difference Splitting Scheme','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(2)) ' at Time $t = $ ' num2str(t_total_2(1,1)) ' s'],'Interpreter','LaTeX')
ylim([0 p_L]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Pressure $p(x)$','Interpreter','LaTeX');  
legend('$\nu=0.8$','$p_{\mathrm{exact}}(x)$ for $\nu=0.8$',  'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_2, pressure(gamma, U_Roe2(:,end,3,2), U_Roe2(:,end,1,2), U_Roe2(:,end,2,2)./U_Roe2(:,end,1,2)), 69,'g', '*');
plot(x_arr_2, U_Roe2_anlyt(:, 3, 2), '-g', 'LineWidth', 2);
title('Roe Flux-Difference Splitting Scheme','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(2)) ' at Time $t = $ ' num2str(t_total_2(2,1)) ' s'],'Interpreter','LaTeX')
ylim([0 p_L]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Pressure $p(x)$','Interpreter','LaTeX');  
legend('$\nu=1.3$','$p_{\mathrm{exact}}(x)$ for $\nu=1.3$',  'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                          %
%                        Functions                         %
%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U1, U2, U3] = initcond(IL, x_arr, gamma, rho_L, p_L, u_L, rho_R, p_R, u_R)
%This function assigns the initial conditions.

for i = 1:IL
    %Populate initial conditions
        if x_arr(i) <= 0
            %Left half
                U1(1,i) = rho_L;
                U2(1,i) = rho_L * u_L;
                U3(1,i) = p_L / (gamma - 1) + (rho_L / 2) * u_L^2;
    
        else
            %Right half
                U1(1,i) = rho_R;
                U2(1,i) = rho_R * u_R;
                U3(1,i) = p_R / (gamma - 1) + (rho_R / 2) * u_R^2;
    
        end        
end

end

function [c] = sos(gamma, p, rho)
%This functions calculates the speed of sound, c, for a calorically perfect
%gas using INPUT PARAMS: (1) heat capacity ratio, gamma; (2) pressure, p;
%(3) density, rho.

c = sqrt(gamma .* p ./ rho);

end

function [p] = pressure(gamma, e, rho, u)
%This function calculates the pressure

p = (gamma - 1) .* (e - (rho ./ 2) .* u.^2);

end

function [E1, E2, E3] = flux(U1, U2, U3, gamma)
%This function calculates the conservative flux vector, E.

E1 = U2;
E2 = U2.^2./U1 + pressure(gamma, U3, U1, U2./U1);
E3 = (U3 + pressure(gamma, U3, U1, U2./U1)) .* (U2./U1);

end

function [E1r, E2r, E3r] = roe_flux(U1, U2, U3, U1p1, U2p1, U3p1, E1, E2, E3, E1p1, E2p1, E3p1, A)
%This function calculates the conservative flux vector, E.

U = [U1; U2; U3];
Up1 = [U1p1; U2p1; U3p1];
E = [E1; E2; E3];
Ep1 = [E1p1; E2p1; E3p1];

Er = 0.5 * (E + Ep1) - 0.5 * A * (Up1 - U);

E1r = Er(1);
E2r = Er(2);
E3r = Er(3);

end
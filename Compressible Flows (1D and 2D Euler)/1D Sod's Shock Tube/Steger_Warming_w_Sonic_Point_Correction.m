%==========================================================================
% AUTHOR: David L. Tran
%
% Steger-Warming scheme with sonic correction
%
% DESCRIPTION: Numerically solves the conservative 1D Euler equation in a
% finite volume formulation using the Steger-Warming scheme with a 
% sonic-point correction. Initial conditions correspond to those of Sod's 
% shock tube. This is a Riemann problem wherein sharp discontinuities are 
% present in the exact solution.
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

%Sonic-point correction
epsilon_prime = 0.06;

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
IL_1 = 600;         %number of grid points for sim 1
IL_2 = 600;         %number of grid points for sim 2
IL = [IL_1; IL_2];

t_steps_1 = 140;     %number of time steps for sim 1
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

U_SW1 = zeros(3, IL(1));
E_SW1 = zeros(3, IL(1));

t_total = 0;
delta_t = 0;

%GRID 2 ARRAYS (time steps = 140, IL = 600)
delta_x_2 = delta_x_1 / 2;        %grid 2 spacing in [m]

x0_2 = delta_x_2 / 2;
x_arr_2_left = -x0_2 - (0:(IL(2)/2)-1)*delta_x_2;
x_arr_2_right = x0_2 + (0:(IL(2)/2)-1)*delta_x_2;
x_arr_2 = [flip(x_arr_2_left), x_arr_2_right];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                          %
%                     Main Body/Loops                      %
%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Steger-Warming with Sonic-Point Correction

%Initialize ICs
[U(1,:), U(2,:), U(3,:)] = initcond(IL(1), x_arr_1, gamma, rho_L, p_L, u_L, rho_R, p_R, u_R);
u = U(2,:)./U(1,:);
p = pressure(gamma, U(3,:), U(1,:), u);
[E(1,:), E(2,:), E(3,:)] = flux(U(1,:), U(2,:), U(3,:), gamma);
c = sos(gamma, p, U(1,:));

%Set initial time step size
delta_t = nu(1) * delta_x_1 / max(abs(u) + c);
t_total = t_total + delta_t;

for n = 1:t_steps(1)
    %loop through all of time
    
    Uold = U; %assign dummy variable

    for i = 2:IL(1)-1
        [A_plus_i, A_minus_i] = geteig1D(Uold(:,i), gamma, epsilon_prime);
        [A_plus_ip1, A_minus_ip1] = geteig1D(Uold(:,i+1), gamma, epsilon_prime);
        [A_plus_im1, A_minus_im1] = geteig1D(Uold(:,i-1), gamma, epsilon_prime);

        %Update U
        U(:,i) = Uold(:,i) - (delta_t / delta_x_1) * (A_plus_i*Uold(:,i) ...
            + A_minus_ip1*Uold(:,i+1) - A_minus_i*Uold(:,i) - A_plus_im1*Uold(:,i-1));
    end
    u = U(2,:)./U(1,:);
    p = pressure(gamma, U(3,:), U(1,:), u);
    [E(1,:), E(2,:), E(3,:)] = flux(U(1,:), U(2,:), U(3,:), gamma);

    c = sos(gamma, p, U(1,:));

    %Calculate new time step size and accumulated time
    delta_t = nu(1) * delta_x_1 / max(abs(u) + c);
    t_total = t_total + delta_t;

end


%analytical for SW grid 1
    %Note: U_anlyt stores rho, u, and p (not the normal U parameters of
    %rho, rho*u, and e)
    for i = 1:IL(1)
        if x_arr_1(i) >= u_S * t_total
            U_anlyt(i, 1) = rho_R;
            U_anlyt(i, 2) = u_R;
            U_anlyt(i, 3) = p_R;
            
        elseif (x_arr_1(i) > u_2 * t_total) && (x_arr_1(i) < u_S * t_total) 
            U_anlyt(i, 1) = rho_2;
            U_anlyt(i, 2) = u_2;
            U_anlyt(i, 3) = p_2;

        elseif (x_arr_1(i) > (u_3 - c_3) * t_total) && (x_arr_1(i) < u_3 * t_total)
            U_anlyt(i, 1) = rho_3;
            U_anlyt(i, 2) = u_3;
            U_anlyt(i, 3) = p_3;

        elseif (x_arr_1(i) > (u_L - c_L) * t_total) && (x_arr_1(i) < (u_3 - c_3) * t_total)
            u_5 = 2 / (gamma + 1) * ( (x_arr_1(i) / t_total) + c_L + (gamma - 1) / 2 * u_L);
            c_5 = u_5 - x_arr_1(i) / t_total;
            p_5 = p_L * (c_5 / c_L)^((2 * gamma) / (gamma - 1));
            rho_5 = (p_5 / p_L)^(1 / gamma) * rho_L;
            U_anlyt(i, 1) = rho_5;
            U_anlyt(i, 2) = u_5;
            U_anlyt(i, 3) = p_5;

        else
            U_anlyt(i, 1) = rho_L;
            U_anlyt(i, 2) = u_L;
            U_anlyt(i, 3) = p_L;
        end
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                          %
%                      Plots/Figures                       %
%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Steger-Warming with Sonic-Point Correction
figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_1, U(1, :), 69,'b', 'd');
plot(x_arr_1, U_anlyt(:, 1), '-b', 'LineWidth', 2);
title('Steger-Warming Flux-Vector Splitting Scheme with Sonic-Point Correction','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(2)) ' at Time $t = $ ' num2str(t_total) ' s'],'Interpreter','LaTeX') 
ylim([0 1]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Density $\rho(x)$','Interpreter','LaTeX');  
% legend('$\nu=0.8$','$\rho_{\mathrm{exact}}(x)$ for $\nu=0.8$', 'Interpreter','latex','location','best');
legend('$\nu=1.3$','$\rho_{\mathrm{exact}}(x)$ for $\nu=1.3$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_1, U(2,:)./U(1,:), 69,'b', 'd');
plot(x_arr_1, U_anlyt(:, 2), '-b', 'LineWidth', 2);
title('Steger-Warming Flux-Vector Splitting Scheme with Sonic-Point Correction','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(2)) ' at Time $t = $ ' num2str(t_total) ' s'],'Interpreter','LaTeX')
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Velocity $u(x)$','Interpreter','LaTeX');  
% legend('$\nu=0.8$','$u_{\mathrm{exact}}(x)$ for $\nu=0.8$', 'Interpreter','latex','location','best');
legend('$\nu=1.3$','$u_{\mathrm{exact}}(x)$ for $\nu=1.3$', 'Interpreter','latex','location','best');
set(gca,'LineWidth',2,'FontSize',18);
ylim([-50 500]);

figure;
hold on;
set(gcf,'Position',[0 25 1200 600]);
scatter(x_arr_1, pressure(gamma, U(3,:), U(1,:), U(2,:)./U(1,:)), 69,'b', 'd');
plot(x_arr_1, U_anlyt(:, 3), '-b', 'LineWidth', 2);
title('Steger-Warming Flux-Vector Splitting Scheme with Sonic-Point Correction','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL(2)) ' at Time $t = $ ' num2str(t_total) ' s'],'Interpreter','LaTeX')
ylim([0 p_L]);
xlabel('Position $x$','Interpreter','LaTeX');
ylabel('Pressure $p(x)$','Interpreter','LaTeX');  
% legend('$\nu=0.8$','$p_{\mathrm{exact}}(x)$ for $\nu=0.8$',  'Interpreter','latex','location','best');
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

function [Aplus, Aminus] = geteig1D(U, gamma, epsilon)
%This function obtains the A+ and A- arrays using the Steger-Warming
%flux-vector splitting scheme with a sonic point-correction.
u = U(2) ./ U(1);
p = pressure(gamma, U(3), U(1), u);
c = sos(gamma, p, U(1));

%Compute auxiliary variables
alpha = (u^2) / 2;
beta = gamma - 1;

%Initialize arrays
lambda = zeros(3,3);
lambda_plus = zeros(3,3);
lambda_minus = zeros(3,3);
P = zeros(3,3);
P_inv = zeros(3,3);

lambda = [u,    0,     0;...
         0,    u+c,    0;...
         0,     0,    u-c];

%Row 1 of P
P(1,1) = 1;
P(1,2) = 1/(2*c^2);
P(1,3) = 1/(2*c^2);

%Row 2 of P
P(2,1) = u;
P(2,2) = u/(2*c^2) + 1/(2*c);
P(2,3) = u/(2*c^2) - 1/(2*c);

%Row 3 of P
P(3,1) = alpha;
P(3,2) = alpha/(2*c^2) + u/(2*c) + 1/(2*beta);
P(3,3) = alpha/(2*c^2) - u/(2*c) + 1/(2*beta);

%Row 1 of P^{-1}
P_inv(1,1) = 1 - (alpha*beta)/(c^2);
P_inv(1,2) = u*beta/(c^2);
P_inv(1,3) = -beta/(c^2);

%Row 2 of P^{-1}
P_inv(2,1) = alpha*beta - u*c;
P_inv(2,2) = -u*beta + c;
P_inv(2,3) = beta;

%Row 3 of P^{-1}
P_inv(3,1) = alpha*beta + u*c;
P_inv(3,2) = -u*beta - c;
P_inv(3,3) = beta;

%Perform sonic correction
for i = 1:3
    lambda_plus(i,i) = (lambda(i,i) + sqrt(lambda(i,i)^2 + epsilon^2) )/ 2;
    lambda_minus(i,i) = (lambda(i,i) - sqrt(lambda(i,i)^2 + epsilon^2) )/ 2;
end

Aplus = P * lambda_plus * P_inv;
Aminus = P * lambda_minus * P_inv;


end

%==========================================================================
%
% Copyright (c) 2023, User assumes full liability from results. Repository
% must be cited and no commercial uses permitted.
%
% AUTHOR: David L. Tran
%
% DESCRIPTION: Solve the 1D Diffusion Heat Equation to obtain the
% temperature distribution across a bar of specified user inputs.
%
%==========================================================================

%% Clear Cache
clc; close all; clear all;


%% Constants
%Physical Constants
k = 100;                   %Thermal conductivity in [W/(m*K)]
T_A = 373;                  %Left boundary temperature in [K]
T_B = 473;                  %Right boundary temperature in [K]
q_gen = 1000;                 %Uniform Volumetric Heat Generation Source in [W/m^3]

%CFD Grid Domain Geometric Constants
N = 62;                     %number of grids
L = 5;                     %bar length in [m]
A = 0.1;                   %bar cross-sectional area in [m^2]


%% Call Functions

[x, T_Analytical] = calculateAnalytical(k, T_A, T_B, q_gen, L, A);
[x_CFD, T_CFD] = calculateCFD(k, T_A, T_B, q_gen, N, L, A);
plotResults(k, N, T_CFD, T_Analytical, q_gen, x, x_CFD);


%% Functions
function [x_CFD, T_CFD] = calculateCFD(k, T_A, T_B, q_gen, N, L, A)
%{
% Calculate the CFD temperature distribution of the 1D bar
% INPUTS:
    % k:            Thermal conductivity in [W/(m*K)]
    % T_A:          Left boundary temperature in [K]
    % T_B:          Right boundary temperature in [K]
    % q_gen:        Uniform Volumetric Heat Generation Source in [W/m^3]
    % L:            length of bar in [m]
    % A:            cross-sectional area of bar in [m^2]

% OUTPUTS:
    % x_CFD:        Array of x coordinates of CFD temp dist in [m] 
    % T_CFD:        Array of CFD temperature distribution in [K]
%}



% Mesh constants
L_cell = L / N;         %length of each cell in [m]
d_LP = L_cell;          %distance bt left and point cell in consideration of centroids in [m]
d_PR = L_cell;          %distance bt point and right cell in consideration of centroids in [m]
d = d_LP;

x_CFD = L_cell/2:L_cell:L;
x_CFD = [0; x_CFD.'; L];

% Calculate material properties
DA = k * A / d;         %material parameter in [W/K]
Dl_Al = DA;
Dr_Ar = DA;
SV = q_gen * A * L_cell;%heat source per unit volume per cell in [W]

A_mat = zeros(N,N);
S_mat = zeros(N,1);

% Assemble matrix
for i = 1:N
    for j = 1:N
        if i == 1
            %Left boundary cell
            a_L = 0;
            a_R = Dr_Ar;

            S_p = -2 * Dl_Al;
            S_u = (T_A - 273) * (2 * Dl_Al) + q_gen * A * L_cell;

            a_p = a_L + a_R - S_p;

            A_mat(i, 1) = a_p;
            A_mat(i, 2) = -a_R;
            if j > 2
                A_mat(i, j) = 0;
            end
            S_mat(i, 1) = S_u;
        elseif i == N
            %Right boundary cell
            a_L = Dl_Al;
            a_R = 0;

            S_p = -2 * Dr_Ar;
            S_u = (T_B - 273)*(2 * Dr_Ar) + q_gen * A * L_cell;

            a_p = a_L + a_R - S_p;
            
            A_mat(i, N-1) = -a_L;
            A_mat(i, N) = a_p;
            if j < N - 1
                A_mat(i, j) = 0;
            end
            S_mat(i, 1) = S_u;
        else
            %Interior cells
            a_L = Dl_Al;
            a_R = Dr_Ar;

            S_p = 0;
            S_u = q_gen * A * L_cell;

            a_p = a_L + a_R - S_p;

            A_mat(i, i-1) = -a_L;
            A_mat(i, i) = a_p;
            A_mat(i, i+1) = -a_R;
            if j > (i + 1)
                A_mat(i, j) = 0;
            end

            S_mat(i, 1) = S_u;
        end
    end
end

%Solve system
T_CFD = linsolve(A_mat, S_mat);
T_CFD = [T_A-273; T_CFD; T_B-273];

end


function [x, T_Analytical] = calculateAnalytical(k, T_A, T_B, q_gen, L, A)
%{
% Calculate and plot the analytical solution of the temperature
distribution 
% INPUTS:
    % k:            Thermal conductivity in [W/(m*K)]
    % T_A:          Left boundary temperature in [K]
    % T_B:          Right boundary temperature in [K]
    % q_gen:        Uniform Volumetric Heat Generation Source in [W/m^3]
    % L:            length of bar in [m]
    % A:            cross-sectional area of bar in [m^2]

% OUTPUTS:
    % x:            Array of x coordinates of analytical temp dist in [m] 
    % T_Analytical: Array of analytical temperature distribution in [K]
%}

%Initialize x array
x = linspace(0, L, 5000);

%Calcualte analytical temperature distribution
T_Analytical = T_A + x ./ L .* (T_B - T_A) + q_gen ./ (2 * k) .* x .* (L - x) - 273;
end

function plotResults(k, N, T_CFD, T_Analytical, q_gen, x, x_CFD)
%{
% This functions plots the analytical temperature distribution across a
% plane wall with uniform heat generation q_gen.
% INPUTS:
    % k:            Thermal conductivity in [W/(m*K)]
    % N:            number of cells
    % T_CFD:        CFD temperature distribution in [K]
    % T_Analytical: Analytical temperature distribution in [K]
    % x:            array of x points for analytical soln in [m]
    % x_CFD:        array of x points for CFD soln in [m]
    
%}


figure;
scatter(x_CFD, T_CFD, 'filled'); hold on;
plot(x, T_Analytical, '-b', 'LineWidth', 2);
hold off;
xlabel('$x$ [m]', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$T$ [$^\circ$C]', 'Interpreter', 'latex', 'FontSize', 16);
title(['1D Bar Conduction Simulation with $N = $ ',num2str(N), ' cells'], 'Interpreter', 'latex', 'FontSize', 18)
legend('CFD','Analytical','Location', 'best', 'Interpreter', 'latex');
end


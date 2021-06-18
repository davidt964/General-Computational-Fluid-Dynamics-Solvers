%==========================================================================
%David Tran
%UCLA ID: 205-492-874
%{
Description: This script implemnts a finite-difference algorithm to
simulate potential flow over the upper half of the NACA 0012 airfoil. Plots
of the pressure coefficient c_p will be generated along the airfoil surface
and compared to my previous panel method, experimental, and academically
published panel method results. Further, 2D contours of c_p in the entire
flow field will be created along with the corresponding streamlines.
%}

%=============================Clear Cache==================================
clc; close all; clear all;


%=============================Constants====================================
V_inf = 1;                                  %Freestream speed in [m/s]
c = 1;                                      %Chord length in [m]
alpha = 0;                                  %Angle of attack in [degrees]
alpha_rad = alpha * pi/180;                 %Angle of attack in [radians]
L = 11 * c;                                 %Total length of the domain in [m]
H = 4 * c;                                  %Height of the domain in [m]
N = 220;                                    %Number of points in the x direction
M = 40;                                     %Number of points in the y direction
delta_x = L / N;                            %Grid size in x direction in [m]
delta_y = H / M;                            %Grid size in y direction in [m]
epsilon = 1e-6;                             %Iteration convergence criterion
MAX_IT = 20000;                             %Maximum number of iterations

psi = zeros(N+1, M+1);                      %Stores values of streamfunction
tempPsi = zeros(N+1, M+1);                  %Stores temporary values of psi that gets updated after each iteration

[xgrid, ygrid] = gridGen(L, H, N, M);       %Generate the grid
[xA, yA] = plotAirfoil(c, N, delta_x);      %Plot the airfoil

tempPsi = applyBCs(tempPsi, H, N, M, yA, V_inf);     
                                            %Obtain initial psi (k) for only the boundary conditions

error = 1;                                  %Temporary error value to enter loop
count = 1;                                  %Count of the number of iterations

h = waitbar(0,'Data loading...');
while (error > epsilon && count < MAX_IT) 
    %Only do this for interior points
    for i=2:N
        for j=2:M
            %Equation 26 (RHS is old psi and LHS is new psi)
            psi(i, j) = delta_y^2 / (2*(delta_x^2 + delta_y^2)) * (tempPsi(i+1,j) + tempPsi(i-1,j)) ...
                + delta_x^2 / (2*(delta_x^2 + delta_y^2)) * (tempPsi(i,j+1) + tempPsi(i,j-1));
        end
    end
    
    %Apply BCs for updated Psi (k + 1)
    psi = applyBCs(psi, H, N, M, yA, V_inf);
    
    %Calculate the error
    error = max(max(abs(psi - tempPsi)));

    %Update tempPsi with the newly calculated iteration
    tempPsi = psi;
    waitbar(count/MAX_IT);
    
    count = count + 1;
end
delete(h);

Cp_surf = calcSurf(psi, N, delta_x, delta_y, V_inf);
Cp = calcPress(psi, N, M, delta_x, delta_y, V_inf, Cp_surf); %Obtain Cp

%==================================Plots===================================
%Experimental and Computational Data Files
Cp_comp = importdata("Cpsurf_Mason.dat");
xread = Cp_comp(:,1);
Cpread = Cp_comp(:,2);

Cp_comp2 = importdata("Cpsurf_GregoryOReilly.dat");
xread2 = Cp_comp2(:,1);
Cpread2 = Cp_comp2(:,2);

%My source panel method for N = 400 pts
Cp_comp3 = importdata("Cp_400.dat");
xread3 = Cp_comp3(:,1);
Cpread3 = Cp_comp3(:,2);

%Plot the grid
figure(1);
hold on;
xlabel('$x$ (m)', 'Interpreter', 'latex');
ylabel('$y$ (m)', 'Interpreter', 'latex');
title('Grid of the Computational Domain');
xlim([0 L]);
ylim([0 H]);

%Plot the vertical lines
for i = 1:N+1
    plot(xgrid(i, :), ygrid(i, :), 'k');
end

%Plot the horizontal lines
for j = 1:M+1
    plot(xgrid(:, j), ygrid(:, j), 'k');
end
hold on;
axis equal;

%Plot the airfoil
plot(xA, yA, 'm', 'LineWidth', 2);
hold off;

%Plot the streamlines
figure(2);
contour(xgrid, ygrid, psi, M, 'LineWidth', 1.5, 'LineColor', 'b');
hold on;
plot(xA, yA, 'k', 'LineWidth', 2); 
xlabel('$x$ (m)', 'Interpreter', 'latex');
ylabel('$y$ (m)', 'Interpreter', 'latex');
title('Streamlines of NACA 0012 Airfoil');
xlim([0 L]);
ylim([0 H]);
axis equal;
hold off;

%Plot the contours of Cp
figure(3);
contourf(xgrid, ygrid, Cp);
hold on;
plot(xA, yA, 'k', 'LineWidth', 2); 
h = colorbar;
ylabel(h, '$C_p$', 'Interpreter', 'latex');
xlabel('$x$ (m)', 'Interpreter', 'latex');
ylabel('$y$ (m)', 'Interpreter', 'latex');
title('$C_p$ Contours', 'Interpreter', 'latex');
axis equal;
hold off;

%Plot my panel method results
x_end = c/delta_x+1;
x_af_start = 5*c/delta_x+1;
x_af_end = 6*c/delta_x+1;
figure(4);
plot(xread3(1:200), Cpread3(1:200), 'r', 'LineWidth', 2); %200 is N/2, where N = 400 since we only want the top surface
hold on;
plot(xA(1:x_end),Cp_surf(x_af_start:x_af_end), '-.om'); %Finite difference
xlabel('$x$ (m)', 'Interpreter', 'latex');
ylabel('$C_p$', 'Interpreter', 'latex');
title('Finite Difference Method vs. Source Panel Method');
xlim([0 c]);
legend('Source Panel Method', 'Finite Difference');
grid on;
hold off;

%Plot the experimental results
figure(5);
plot(xread2, Cpread2, 'r', 'LineWidth', 2);
hold on;
plot(xA(1:x_end),Cp_surf(x_af_start:x_af_end), '-.om'); %Finite difference
xlabel('$x$ (m)', 'Interpreter', 'latex');
ylabel('$C_p$', 'Interpreter', 'latex');
title('Finite Difference Method vs. Experimental Results');
xlim([0 c])
legend('Experimental', 'Finite Difference');
grid on;
hold off;

%Plot published panel method results
figure(6);
plot(xread, Cpread, 'r', 'LineWidth', 2);
hold on;
plot(xA(1:x_end),Cp_surf(x_af_start:x_af_end), '-.om'); %Finite difference
xlabel('$x$ (m)', 'Interpreter', 'latex');
ylabel('$C_p$', 'Interpreter', 'latex');
title('Finite Difference Method vs. Published Results');
xlim([0 c]);
legend('Mason', 'Finite Difference');
grid on;
hold off;


%===============================Functions==================================
function [xgrid, ygrid] = gridGen(L, H, N, M)
%{
This function generates the grid taking in the length L and height H of the
computational domain, spacing out the grid lines based on the number of
points in the x direction N and number of points in the y direction M.
%}
delta_x = L / N;                            %Grid size in x direction in [m]
delta_y = H / M;                            %Grid size in y direction in [m]

for i=1:N+1
    for j = 1:M+1
        xgrid(i, j) = (i - 1) * delta_x;
        ygrid(i, j) = (j - 1) * delta_y;
    end
end
end

function [xA, yA] = plotAirfoil(c, N, delta_x)
%{
Plots the airfoil with input parameters: chord length c, number of points
in the x direction N, and spacing between x grid points delta_x. Returns
the x and y coordinates of the airfoil [xA, yA].
%}
t = 0.12;                               %NACA 00XX, t is XX/100
xA = zeros(1,N+1);                      %array to store x values
yA = zeros(1,N+1);                      %array to store y values

for k=1:N+1
    if k < (5 * c / delta_x)
        xA(k) = (k - 1) * delta_x;
        yA(k) = 0;
    elseif (k >= 5 * c / delta_x + 1) && (k <= 6 * c / delta_x + 1)
        xA(k) = delta_x * (k - 1);
        yA(k) = t/(0.2) * c * (0.2969 * sqrt((xA(k) - 5*c) / c) - 0.1260 * ((xA(k) - 5*c)/c) - 0.3516 * ((xA(k) - 5*c)/c)^2 + 0.2843 * ((xA(k) - 5*c)/c)^3 - 0.1015 * ((xA(k) - 5*c)/c)^4 );
    else
        xA(k) = delta_x * (k - 1);
        yA(k) = 0;
    end
end

%For some reason, MATLAB kept converting my yA array to a complex double
%array with an imaginary element of +0.0000i. In other words, there was no
%imaginary element for any of the values of yA but MATLAB made it imaginary?
%Hence, the inclusion of this line.
yA = real(yA);
end

function psi = applyBCs(psi, H, N, M, yA, V_inf)
%{
Applies the boundary conditions at the bottom wall (including airfoil),
left, right, and top walls. Returns the streamfunction psi at these points
%}

%Top Wall BC
for i = 1:N+1
    psi(i, M+1) = V_inf * H;
end

%Left Wall BC
for j = 1:M+1
    psi(1, j) = V_inf * (H / M) * (j - 1);
end

%Right Wall
for j = 1:M+1
    psi(N+1, j) = psi(N, j);
end

%Lower Boundary BC
for i = 1:N+1
    psi(i, 1) = -(yA(i) / ((H/M) - yA(i))) * psi(i, 2);
end
end

function Cp_surf = calcSurf(psi, N, delta_x, delta_y, V_inf)
%{
Calculates the pressure coefficient on the surface of the airfoil
%}
for i = 1:N+1
    if i > 1 && i < N +1
        u(i) = (psi(i, 2) - psi(i, 1)) / delta_y;
        v(i) = -(psi(i+1, 1) - psi(i-1, 1)) / (2 * delta_x);
    else
        u(i) = V_inf;  %only a horizontal component
        v(i) = 0;      %no normal (y) component
    end
end

Cp_surf = 1 - (u.^2 + v.^2) / V_inf^2;

end

function Cp = calcPress(psi, N, M, delta_x, delta_y, V_inf, Cp_surf)
%{
Calculates the pressure coefficient for the interior points in the
computational grid domain with inputs streamfunction psi, number of points 
in x direction N, number of points in y direction M, change in x delta_x, 
change in y delta_y, and freestream velocity V_inf. Returns the pressure 
coefficient for the interior points only.
%}

%Calculate velocity at each interior grid point
for i = 2:N
    for j = 2:M
        %Only for interior points
        u(i, j) = (psi(i,j+1) - psi(i,j-1)) / (2 * delta_y);
        v(i, j) = -(psi(i+1, j) - psi(i-1, j)) / (2 * delta_x);
    end
end

%Apply velocity boundary conditions
%Left wall
for j = 1:M+1
    u(1, j) = V_inf;
    v(1, j) = 0;
end

%Top wall
for i = 1:N+1
    u(i, M+1) = V_inf;
    v(i, M+1) = 0;
end

%Right wall
for j = 1:M+1
    u(N+1, j) = V_inf;
    v(N+1, j) = 0;
end

%Calculate Cp
for i = 1:N+1
    for j = 1:M+1
        if j == 1
            Cp(i, j) = Cp_surf(i);
        else
            Cp(i, j) = 1 - (u(i,j)^2 + v(i,j)^2)/(V_inf^2);
        end
    end
end

end
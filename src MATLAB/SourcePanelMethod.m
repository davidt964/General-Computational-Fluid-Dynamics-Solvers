%==========================================================================
%Source Panel Method
%David Tran
%UCLA ID: 205-492-874
%{
Description:
This script seeks to implement the source panel method to generate a plot
of the pressure coefficient Cp as a function of x/c, where c is the chord
length. Results from this method will be compared to provided experimental
and computational results to assess its accuracy. Further, data will be
created for panel counts of 50, 100, 200, and 400 to demonstrate the effect
of increasing the number of panels.
%}

%=============================Clear Cache==================================
clc; close all; clear all;

%=============================Constants====================================
alpha = 0;                              %Angle of attack in degrees
N_1 = 50;                               %Number of panels for set 1
N_2 = 100;                              %Number of panels for set 2
N_3 = 200;                              %Number of panels for set 3
N_4 = 400;                              %Number of panels for set 4
V_infty = 1;                            %speed of freestream velocity in [m/s]
c = 1;                                  %chord length in [m]

%Plot the airfoil
[X, Y] = plotAirfoil();

%%First set of data points for N = 50 panels
A1 = zeros(N_1, N_1);
A1_t = zeros(N_1, N_1);
B1 = zeros(N_1, 1);
B1_t = zeros(N_1, 1);
beta_i1 = zeros(1, N_1);
phi_i1 = zeros(1, N_1);
Sj_1 = zeros(1, N_1);
[X1, Y1] = panels(N_1);                  %Obtain panel endpoints
[xc, yc] = controlPoint(N_1, X1, Y1);    %Obtain panel control points

%Populate the A and B matrices
for i=1:N_1
    phi_i1(i) = phi_calc(X1(i), Y1(i), X1(i+1), Y1(i+1));
    beta_i1(i) = phi_i1(i) + pi/2;
    B1(i) = -V_infty * cos(beta_i1(i)); 
    Sj_1(i) = sqrt((X1(i+1) - X1(i))^2 + (Y1(i+1) - Y1(i))^2);
    B1_t(i) = V_infty * sin(beta_i1(i));
    for j=1:N_1
        if i ~= j
            phi_i1(j) = phi_calc(X1(j), Y1(j), X1(j+1), Y1(j+1));
            A1(i,j) = 1 / (2*pi) * integral(xc(i), yc(i), X1(j), Y1(j), X1(j+1), Y1(j+1), phi_i1(i), phi_i1(j));
            A1_t(i, j) = 1 / (2*pi) * integral_j(xc(i), yc(i), X1(j), Y1(j), X1(j+1), Y1(j+1), phi_i1(i), phi_i1(j));
        elseif i == j
            A1(i, j) = 0.5;
            A1_t(i, j) = 0;
        end
    end
end

%Solve the Ax = B system and verify that Eq. 3.157 holds (sum of i-th
%source/sink multiplied by i-th panel length for all panels is 0)
lambda_1 = A1\B1;
sourceSum(lambda_1, Sj_1, N_1);

%Obtain the tangential velocity
V_i1 = tangentialVelocity(A1_t, lambda_1, B1_t, N_1);

%Calcuclate the pressure coefficient
Cp_1 = Cp(V_i1, V_infty, N_1);

%%Second set of data points for N = 100 panels
A2 = zeros(N_2, N_2);
A2_t = zeros(N_2, N_2);
B2 = zeros(N_2, 1);
B2_t = zeros(N_2, 1);
beta_i2 = zeros(1, N_2);
phi_i2 = zeros(1, N_2);
Sj_2 = zeros(1, N_2);
[X2, Y2] = panels(N_2);                    %Obtain panel endpoints
[xc2, yc2] = controlPoint(N_2, X2, Y2);    %Obtain panel control points

%Populate the A and B matrices
for i=1:N_2
    phi_i2(i) = phi_calc(X2(i), Y2(i), X2(i+1), Y2(i+1));
    beta_i2(i) = phi_i2(i) + pi/2;
    B2(i) = -V_infty * cos(beta_i2(i)); 
    Sj_2(i) = sqrt((X2(i+1) - X2(i))^2 + (Y2(i+1) - Y2(i))^2);
    B2_t(i) = V_infty * sin(beta_i2(i));
    for j=1:N_2
        if i ~= j
            phi_i2(j) = phi_calc(X2(j), Y2(j), X2(j+1), Y2(j+1));
            A2(i,j) = 1 / (2*pi) * integral(xc2(i), yc2(i), X2(j), Y2(j), X2(j+1), Y2(j+1), phi_i2(i), phi_i2(j));
            A2_t(i, j) = 1 / (2*pi) * integral_j(xc2(i), yc2(i), X2(j), Y2(j), X2(j+1), Y2(j+1), phi_i2(i), phi_i2(j));
        elseif i == j
            A2(i, j) = 0.5;
            A2_t(i, j) = 0;
        end
    end
end

%Solve the Ax = B system and verify that Eq. 3.157 holds (sum of i-th
%source/sink multiplied by i-th panel length for all panels is 0)
lambda_2 = A2\B2;
sourceSum(lambda_2, Sj_2, N_2);

%Obtain the tangential velocity
V_i2 = tangentialVelocity(A2_t, lambda_2, B2_t, N_2);

%Calcuclate the pressure coefficient
Cp_2 = Cp(V_i2, V_infty, N_2);

%%Third set of data points for N = 200
A3 = zeros(N_3, N_3);
A3_t = zeros(N_3, N_3);
B3 = zeros(N_3, 1);
B3_t = zeros(N_3, 1);
beta_i3 = zeros(1, N_3);
phi_i3 = zeros(1, N_3);
Sj_3 = zeros(1, N_3);
[X3, Y3] = panels(N_3);                    %Obtain panel endpoints
[xc3, yc3] = controlPoint(N_3, X3, Y3);    %Obtain panel control points

%Populate the A and B matrices
for i=1:N_3
    phi_i3(i) = phi_calc(X3(i), Y3(i), X3(i+1), Y3(i+1));
    beta_i3(i) = phi_i3(i) + pi/2;
    B3(i) = -V_infty * cos(beta_i3(i)); 
    Sj_3(i) = sqrt((X3(i+1) - X3(i))^2 + (Y3(i+1) - Y3(i))^2);
    B3_t(i) = V_infty * sin(beta_i3(i));
    for j=1:N_3
        if i ~= j
            phi_i3(j) = phi_calc(X3(j), Y3(j), X3(j+1), Y3(j+1));
            A3(i,j) = 1 / (2*pi) * integral(xc3(i), yc3(i), X3(j), Y3(j), X3(j+1), Y3(j+1), phi_i3(i), phi_i3(j));
            A3_t(i, j) = 1 / (2*pi) * integral_j(xc3(i), yc3(i), X3(j), Y3(j), X3(j+1), Y3(j+1), phi_i3(i), phi_i3(j));
        elseif i == j
            A3(i, j) = 0.5;
            A3_t(i, j) = 0;
        end
    end
end

%Solve the Ax = B system and verify that Eq. 3.157 holds (sum of i-th
%source/sink multiplied by i-th panel length for all panels is 0)
lambda_3 = A3\B3;
sourceSum(lambda_3, Sj_3, N_3);

%Obtain the tangential velocity
V_i3 = tangentialVelocity(A3_t, lambda_3, B3_t, N_3);

%Calcuclate the pressure coefficient
Cp_3 = Cp(V_i3, V_infty, N_3);

%%Fourth set of data points for N = 400 panels
A4 = zeros(N_4, N_4);
A4_t = zeros(N_4, N_4);
B4 = zeros(N_4, 1);
B4_t = zeros(N_4, 1);
beta_i4 = zeros(1, N_4);
phi_i4 = zeros(1, N_4);
Sj_4 = zeros(1, N_4);
[X4, Y4] = panels(N_4);                    %Obtain panel endpoints
[xc4, yc4] = controlPoint(N_4, X4, Y4);    %Obtain panel control points

%Populate the A and B matrices
for i=1:N_4
    phi_i4(i) = phi_calc(X4(i), Y4(i), X4(i+1), Y4(i+1));
    beta_i4(i) = phi_i4(i) + pi/2;
    B4(i) = -V_infty * cos(beta_i4(i)); 
    Sj_4(i) = sqrt((X4(i+1) - X4(i))^2 + (Y4(i+1) - Y4(i))^2);
    B4_t(i) = V_infty * sin(beta_i4(i));
    for j=1:N_4
        if i ~= j
            phi_i4(j) = phi_calc(X4(j), Y4(j), X4(j+1), Y4(j+1));
            A4(i,j) = 1 / (2*pi) * integral(xc4(i), yc4(i), X4(j), Y4(j), X4(j+1), Y4(j+1), phi_i4(i), phi_i4(j));
            A4_t(i, j) = 1 / (2*pi) * integral_j(xc4(i), yc4(i), X4(j), Y4(j), X4(j+1), Y4(j+1), phi_i4(i), phi_i4(j));
        elseif i == j
            A4(i, j) = 0.5;
            A4_t(i, j) = 0;
        end
    end
end

%Solve the Ax = B system and verify that Eq. 3.157 holds (sum of i-th
%source/sink multiplied by i-th panel length for all panels is 0)
lambda_4 = A4\B4;
sourceSum(lambda_4, Sj_4, N_4);

%Obtain the tangential velocity
V_i4 = tangentialVelocity(A4_t, lambda_4, B4_t, N_4);

%Calcuclate the pressure coefficient
Cp_4 = Cp(V_i4, V_infty, N_4);

%================================Plots=====================================
%Experimental and Computational Data Files
Cp_comp = importdata("Cpsurf_Mason.dat");
xread = Cp_comp(:,1);
Cpread = Cp_comp(:,2);

Cp_comp2 = importdata("Cpsurf_GregoryOReilly.dat");
xread2 = Cp_comp2(:,1);
Cpread2 = Cp_comp2(:,2);

%Plots for N = 50
set(gcf, 'Position',  [400, 100, 800, 800]);
figure(1);
subplot(3,1,1);
plot(X, Y, 'k', 'LineWidth', 2);
hold on;
plot(X1, Y1, '-.or', 'LineWidth', 1.25);
title('NACA 0012 Geometry', 'Interpreter', 'latex');
xlabel('$x/c$', 'Interpreter', 'latex');
ylabel('$y/c$', 'Interpreter', 'latex');
N1_legend = (sprintf('%d Source Panels', N_1));
legend('Airfoil', N1_legend, 'Interpreter', 'latex');
grid on;
axis equal;

subplot(3,1,2);
Cp1_title = (sprintf('NACA 0012 $C_p$ vs. Published Data for %d panels', N_1));
hold on;
title(Cp1_title, 'Interpreter', 'latex');
xlabel('$x/c$', 'Interpreter', 'latex');
ylabel('$C_p$', 'Interpreter', 'latex');
plot(xc(1:N_1/2), Cp_1(1:N_1/2), '-.om', 'LineWidth', 2);
plot(xread, Cpread, 'r', 'LineWidth', 2);
xlim([0 c]);
ylim([-2 1]);
legend('Source Panel Method', 'Mason Panel Method', 'Interpreter', 'latex');
grid on;

subplot(3,1,3);
Cp1_title1 = (sprintf('NACA 0012 $C_p$ vs. Experimental Data for %d panels', N_1));
hold on;
title(Cp1_title1, 'Interpreter', 'latex');
xlabel('$x/c$', 'Interpreter', 'latex');
ylabel('$C_p$', 'Interpreter', 'latex');
plot(xc(1:N_1/2), Cp_1(1:N_1/2), '-.om', 'LineWidth', 2);
plot(xread2, Cpread2, 'b', 'LineWidth', 2);
xlim([0 c]);
ylim([-2 1]);
legend('Source Panel Method', 'Experimental Data', 'Interpreter', 'latex');
grid on;

%Plots for N = 100
figure(2);
set(gcf, 'Position',  [400, 100, 800, 800]);
subplot(3,1,1);
plot(X, Y, 'k', 'LineWidth', 2);
hold on;
plot(X2, Y2, '-.or', 'LineWidth', 1.25);
title('NACA 0012 Geometry', 'Interpreter', 'latex');
xlabel('$x/c$', 'Interpreter', 'latex');
ylabel('$y/c$', 'Interpreter', 'latex');
N2_legend = (sprintf('%d Source Panels', N_2));
legend('Airfoil', N2_legend, 'Interpreter', 'latex');
grid on;
axis equal;

subplot(3,1,2);
Cp2_title = (sprintf('NACA 0012 $C_p$ vs. Published Data for %d panels', N_2));
hold on;
title(Cp2_title, 'Interpreter', 'latex');
xlabel('$x/c$', 'Interpreter', 'latex');
ylabel('$C_p$', 'Interpreter', 'latex');
plot(xc2(1:N_2/2), Cp_2(1:N_2/2), '-.om', 'LineWidth', 2);
plot(xread, Cpread, 'r', 'LineWidth', 2);
xlim([0 c]);
ylim([-2 1]);
legend('Source Panel Method', 'Mason Panel Method', 'Interpreter', 'latex');
grid on;

subplot(3,1,3);
Cp2_title1 = (sprintf('NACA 0012 $C_p$ vs. Experimental Data for %d panels', N_2));
hold on;
title(Cp2_title1, 'Interpreter', 'latex');
xlabel('$x/c$', 'Interpreter', 'latex');
ylabel('$C_p$', 'Interpreter', 'latex');
plot(xc2(1:N_2/2), Cp_2(1:N_2/2), '-.om', 'LineWidth', 2);
plot(xread2, Cpread2, 'b', 'LineWidth', 2);
xlim([0 c]);
ylim([-2 1]);
legend('Source Panel Method', 'Experimental Data', 'Interpreter', 'latex');
grid on;

%Plots for N = 200
figure(3);
set(gcf, 'Position',  [400, 100, 800, 800]);
subplot(3,1,1);
plot(X, Y, 'k', 'LineWidth', 2);
hold on;
plot(X3, Y3, '-.or', 'LineWidth', 1.25);
title('NACA 0012 Geometry', 'Interpreter', 'latex');
xlabel('$x/c$', 'Interpreter', 'latex');
ylabel('$y/c$', 'Interpreter', 'latex');
N3_legend = (sprintf('%d Source Panels', N_3));
legend('Airfoil', N3_legend, 'Interpreter', 'latex');
grid on;
axis equal;

subplot(3,1,2);
Cp3_title = (sprintf('NACA 0012 $C_p$ vs. Published Data for %d panels', N_3));
hold on;
title(Cp3_title, 'Interpreter', 'latex');
xlabel('$x/c$', 'Interpreter', 'latex');
ylabel('$C_p$', 'Interpreter', 'latex');
plot(xc3(1:N_3/2), Cp_3(1:N_3/2), '-.om', 'LineWidth', 2);
plot(xread, Cpread, 'r', 'LineWidth', 2);
xlim([0 c]);
ylim([-2 1]);
legend('Source Panel Method', 'Mason Panel Method', 'Interpreter', 'latex');
grid on;

subplot(3,1,3);
Cp3_title1 = (sprintf('NACA 0012 $C_p$ vs. Experimental Data for %d panels', N_3));
hold on;
title(Cp3_title1, 'Interpreter', 'latex');
xlabel('$x/c$', 'Interpreter', 'latex');
ylabel('$C_p$', 'Interpreter', 'latex');
plot(xc3(1:N_3/2), Cp_3(1:N_3/2), '-.om', 'LineWidth', 2);
plot(xread2, Cpread2, 'b', 'LineWidth', 2);
xlim([0 c]);
ylim([-2 1]);
legend('Source Panel Method', 'Experimental Data', 'Interpreter', 'latex');
grid on;

%Plots for N = 400
figure(4);
set(gcf, 'Position',  [400, 100, 800, 800]);
subplot(3,1,1);
plot(X, Y, 'k', 'LineWidth', 2);
hold on;
plot(X4, Y4, '-.or', 'LineWidth', 1.25);
title('NACA 0012 Geometry', 'Interpreter', 'latex');
xlabel('$x/c$', 'Interpreter', 'latex');
ylabel('$y/c$', 'Interpreter', 'latex');
N4_legend = (sprintf('%d Source Panels', N_4));
legend('Airfoil', N4_legend, 'Interpreter', 'latex');
grid on;
axis equal;

subplot(3,1,2);
Cp4_title = (sprintf('NACA 0012 $C_p$ vs. Published Data for %d panels', N_4));
hold on;
title(Cp4_title, 'Interpreter', 'latex');
xlabel('$x/c$', 'Interpreter', 'latex');
ylabel('$C_p$', 'Interpreter', 'latex');
plot(xc4(1:N_4/2), Cp_4(1:N_4/2), '-.om', 'LineWidth', 2);
plot(xread, Cpread, 'r', 'LineWidth', 2);
xlim([0 c]);
ylim([-2 1]);
legend('Source Panel Method', 'Mason Panel Method', 'Interpreter', 'latex');
grid on;

subplot(3,1,3);
Cp4_title1 = (sprintf('NACA 0012 $C_p$ vs. Experimental Data for %d panels', N_4));
hold on;
title(Cp4_title1, 'Interpreter', 'latex');
xlabel('$x/c$', 'Interpreter', 'latex');
ylabel('$C_p$', 'Interpreter', 'latex');
plot(xc4(1:N_4/2), Cp_4(1:N_4/2), '-.om', 'LineWidth', 2);
plot(xread2, Cpread2, 'b', 'LineWidth', 2);
xlim([0 c]);
ylim([-2 1]);
legend('Source Panel Method', 'Experimental Data', 'Interpreter', 'latex');
grid on;

%%NOTE: Also, although I am aware that it says to declare functions in
%%separate MATLAB files, I find that to be too clunky. When submitting the
%%final product, it becomes annoying submitting all of the separate files
%%and for graders to download them all.
%=============================Functions====================================
function [xA, yA] = plotAirfoil()
%This function computes the coordinates of a NACA XXXX airfoil with input
%parameters N, which is the number of panels, and topStatus, which tells
%the function if the points are on the top or bottom surface (1 for top and
%2 for bottom). Function returns an array of the x and y coordinates.
t = 0.12;                               %NACA 00XX, t is XX/100
c = 1;                                  %chord length
N = 1000;                               %Number of points
xA_upper = linspace(0, c, (N/2)+1);     %precalculated evenly spaced x values from [0, c]
xA_lower = linspace(c, 0, N/2);         %precalculated evenly spaced x values from [c, 0]
xA = zeros(1,N+1);                      %array to store x values
yA = zeros(1,N+1);                      %array to store y values

for k=1:N+1
    if k <= N/2 + 1
        xA(k) = xA_upper(k);
        yA(k) = t/(0.2) * c * (0.2969 * sqrt(xA(k) / c) - 0.1260 * (xA(k)/c) - 0.3516 * (xA(k)/c)^2 + 0.2843 * (xA(k)/c)^3 - 0.1015 * (xA(k)/c)^4);
    else
        xA(k) = xA_lower(k-(N/2) - 1);
        yA(k) = -t/(0.2) * c * (0.2969 * sqrt(xA(k) / c) - 0.1260 * (xA(k)/c) - 0.3516 * (xA(k)/c)^2 + 0.2843 * (xA(k)/c)^3 - 0.1015 * (xA(k)/c)^4);
    end
end

end

function [X, Y] = panels(N)
%This function computes the coordinates of a NACA XXXX airfoil with input
%parameters N, which is the number of panels, and topStatus, which tells
%the function if the points are on the top or bottom surface (1 for top and
%2 for bottom). Function returns an array of the x and y coordinates.
t = 0.12;                               %NACA 00XX, t is XX/100
c = 1;                                  %chord length
xA_upper = linspace(0, c, (N/2)+1);     %precalculated evenly spaced x values from [0, c]
xA_lower = linspace(c, 0, N/2);         %precalculated evenly spaced x values from [c, 0]
X = zeros(1,N+1);                       %array to store x values
Y = zeros(1,N+1);                       %array to store y values

for k=1:N+1
    if k <= N/2 + 1
        X(k) = xA_upper(k);
        Y(k) = t/(0.2) * c * (0.2969 * sqrt(X(k) / c) - 0.1260 * (X(k)/c) - 0.3516 * (X(k)/c)^2 + 0.2843 * (X(k)/c)^3 - 0.1015 * (X(k)/c)^4);
    else
        X(k) = xA_lower(k-(N/2) - 1);
        Y(k) = -t/(0.2) * c * (0.2969 * sqrt(X(k) / c) - 0.1260 * (X(k)/c) - 0.3516 * (X(k)/c)^2 + 0.2843 * (X(k)/c)^3 - 0.1015 * (X(k)/c)^4);
    end
end
end

function [xc, yc] = controlPoint(N, x, y)
%Calculates the control points of each panel from the x and y arrays.
%Takes in inputs N (number of panels), x and y (coordinates of the
%endpoints of the panels)

for k = 1:N
    xc(k) = (x(k) + x(k+1)) / 2;
    yc(k) = (y(k) + y(k+1)) / 2;
end

end

function Iij = integral(xi, yi, Xj, Yj, Xjp1, Yjp1, phi_i, phi_j)
%Calculate the ij integral component using the integral given by Eq. 3.162
%and 3.163.
A = -(xi - Xj) * cos(phi_j) - (yi - Yj) * sin(phi_j);
B = (xi - Xj)^2 + (yi - Yj)^2;
C = sin(phi_i - phi_j);
D = (yi - Yj) * cos(phi_i) - (xi - Xj) * sin(phi_i);
E = sqrt(B - A^2);
Sj = sqrt((Xjp1 - Xj)^2 + (Yjp1 - Yj)^2);

Iij = (C / 2) * log((Sj^2 + 2*A*Sj + B) / B) + (D - A*C)/E * (atan2(Sj + A, E) - atan2(A, E));

end

function Jij = integral_j(xi, yi, Xj, Yj, Xjp1, Yjp1, phi_i, phi_j)
%Calculate the integral over the jth panel given by equation 3.165
A = -(xi - Xj) * cos(phi_j) - (yi - Yj) * sin(phi_j);
B = (xi - Xj)^2 + (yi - Yj)^2;
C = sin(phi_i - phi_j);
D = (yi - Yj) * cos(phi_i) - (xi - Xj) * sin(phi_i);
E = sqrt(B - A^2);
Sj = sqrt((Xjp1 - Xj)^2 + (Yjp1 - Yj)^2);

Jij = (D - A*C)/(2*E) * log((Sj^2 + 2*A*Sj + B) / B) - C * (atan2(Sj + A, E) - atan2(A, E));

end

function phi = phi_calc(X, Y, Xp1, Yp1)
%This function calculates phi (the angle of the panel with respect to the
%horizontal). Has four inputs, the start points of the panel and the end
%points of the panel.
delta_X = Xp1 - X;
delta_Y = Yp1 - Y;

phi = atan2(delta_Y, delta_X);

end

function Cpi = Cp(Vi, V_infty, N)
%This function calculate the pressure coefficient Cp with two inputs: the
%tangential velocity component, Vi, and the freestream velocity, V_infty.
for i=1:N
    Cpi(i) = 1 - (Vi(i) / V_infty)^2;
end

end

function sourceSum(lambda, sj, N)
%lambda and sj are arrays that store the source strength of each panel and
%length of each panel. It has no outputs and prints the sum of the source
%strength multiplied by the panel length.
product = zeros(1, N);

for i=1:N
    product(i) = lambda(i) * sj(i);
end

zeroCheck = sum(product, 'all');

fprintf('Sum of all sources/sinks for %d panels: %f\n', N, zeroCheck);

end

function Vi = tangentialVelocity(A_t, lambda, B_t, N)
%Calculates the tangential velocity of each panel, where B_t is defined by
%B_t = V_infty * sin(beta_i). The first three parameters are arrays.
Vi = zeros(N, 1);

for i = 1:N
    Vi(i) = Vi(i) + B_t(i);
    for j = 1:N
        Vi(i) = Vi(i) + A_t(i, j) * lambda(j);
    end
end
end

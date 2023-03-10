%==========================================================================
%Author: David L. Tran
%Description: Creates a rectangular grid and calculates the streamfunction
%and various flow properties (e.g., Cp) at each point. The domain is
%defined as a having a trigonometric curve on the top and bottom walls.
%}

%=============================Clear Cache==================================
clc; close all; clear all;


%=============================Constants====================================
V_inf = 1;                                  %flow speed in [m/s]
L = 10;                                     %domain length in [m]
H = 1;                                      %domain height in [m]
N = 50;                                     %number of horizontal lines
M = 50;                                     %number of vertical lines
xE = L / 3;                                 %start of curved segments  in [m]
xH = xE;                                    
xF = 2 * L / 3;                             %end of curved segments in [m]
xG = xF;
delta_x = L / N;                            %Grid size in x direction in [m]
delta_y = H / M;                            %Grid size in y direction in [m]
epsilon = 1e-6;                             %Iteration convergence criterion
MAX_IT = 20000;                             %Maximum number of iterations

psi = zeros(N+1, M+1);                      %Stores values of streamfunction
tempPsi = zeros(N+1, M+1);                  %Stores temporary values of psi that gets updated after each iteration

[xgrid, ygrid] = gridGen(L, H, N, M);       %Generate the grid
[xTop, yTop] = plotTop(L, H, N, xH);
[xBot, yBot] = plotBottom(L, H, N, xE);

error = 1;                                  %Temporary error value to enter loop
count = 1;                                  %Count of the number of iterations

h = waitbar(0,'Data loading...');

tempPsi = applyBCs(psi, H, N, M, yTop, yBot, V_inf);    
                                            %Obtain initial psi (k) for only the boundary conditions
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
    psi = applyBCs(psi, H, N, M, yTop, yBot, V_inf);
    
    %Calculate the error
    error = max(max(abs(psi - tempPsi)));

    %Update tempPsi with the newly calculated iteration
    tempPsi = psi;
    waitbar(count/MAX_IT);
    
    count = count + 1;
end
waitbar(MAX_IT/MAX_IT);
delete(h);

cp_topsurf = calcTopSurf(psi, N, M, delta_x, delta_y, V_inf);
cp_botsurf = calcBotSurf(psi, N, delta_x, delta_y, V_inf);

Cp = calcPress(psi, cp_botsurf, cp_topsurf, V_inf, L, H, N, M);


%==================================Plots===================================
figure(1);
set(gcf, 'Position',  [400, 100, 800, 800]);
%THIS IS JUST A SANITY CHECK TO SEE IF I AM DOING EVERYTHING RIGHT; IT
%LOOKS GOOD TO ME I THINK?
hold on;
contour(xgrid, ygrid, psi, M/3, 'LineWidth', 1.5, 'LineColor', 'b');
plot(xTop, yTop, 'k', 'LineWidth', 2);
plot(xBot, yBot, 'k', 'LineWidth', 2);
xlabel('$x$ (m)', 'Interpreter', 'latex');
ylabel('$y$ (m)', 'Interpreter', 'latex');
title('Streamlines');
axis equal;
hold off;

%Plot the contours of Cp
figure(2);
set(gcf, 'Position',  [400, 100, 800, 800]);
hold on;
contourf(xgrid, ygrid, Cp);
h = colorbar;
ylabel(h, '$C_p$', 'Interpreter', 'latex');
xlabel('$x$ (m)', 'Interpreter', 'latex');
ylabel('$y$ (m)', 'Interpreter', 'latex');
title('$C_p$ Contours', 'Interpreter', 'latex');
axis equal;
hold off;


%Plot the bottom and top cp vs. x
%Bottom surface
figure(3);
set(gcf, 'Position',  [400, 100, 800, 800]);
hold on;
plot(xBot, cp_botsurf, 'k', 'LineWidth', 2);
xlabel('$x$ (m)', 'Interpreter', 'latex');
ylabel('$C_p$', 'Interpreter', 'latex');
title('$C_p$ vs. $x$ for Bottom Surface', 'Interpreter', 'latex');
grid on;
hold off;

disp('%========Bottom Surface========');
disp('i         xi            Cpi');
for i = 1:N+1
    fprintf('%d       %f           %f \n', i, xBot(i), cp_botsurf(i));
end

%Top surface
figure(4);
set(gcf, 'Position',  [400, 100, 800, 800]);
hold on;
plot(xTop, cp_topsurf, 'k', 'LineWidth', 2);
xlabel('$x$ (m)', 'Interpreter', 'latex');
ylabel('$C_p$', 'Interpreter', 'latex');
title('$C_p$ vs. $x$ for Top Surface', 'Interpreter', 'latex');
grid on;
hold off;

disp('%========Top Surface========');
disp('i         xi            Cpi');
for i = 1:N+1
    fprintf('%d       %f           %f \n', i, xTop(i), cp_topsurf(i));
end

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
        ygrid(i, j) = (j - 1) * delta_y - H/2;
    end
end
end

function [xTop, yTop] = plotTop(L, H, N, xH)
delta_x = L / N;                            %Grid size in x direction in [m]

for i = 1:N+1
    if i <= floor((L/3)/delta_x) + 1
        xTop(i) = delta_x * (i - 1);
        yTop(i) = H/2;
    elseif i > floor((L/3)/delta_x) + 1 && i <= floor((2*L/3)/delta_x) + 1
        xTop(i) = delta_x * (i - 1);
        yTop(i) = (H/2) + (H/100) * sin((6 * pi *(xTop(i) - xH)) / L);
    else
        xTop(i) = delta_x * (i - 1);
        yTop(i) = H/2;
    end
end

end

function [xBot, yBot] = plotBottom(L, H, N, xE)
delta_x = L / N;                            %Grid size in x direction in [m]

for i = 1:N+1
    if i <= floor((L/3)/delta_x) + 1
        xBot(i) = delta_x * (i - 1);
        yBot(i) = -H/2;
    elseif i > floor((L/3)/delta_x) + 1 && i <= floor((2*L/3)/delta_x) + 1
        xBot(i) = delta_x * (i - 1);
        yBot(i) = -(H/2) + (H/120) * sin((6 * pi *(xBot(i) - xE)) / L);
    else
        xBot(i) = delta_x * (i - 1);
        yBot(i) = -H/2;
    end
end
end

function psi = applyBCs(psi, H, N, M, yTop, yBot, V_inf)
delta_y = H/M;

%Top Wall BC
for i = 1:N+1
    psi(i, M+1) = (V_inf * H * delta_y + 2 * psi(i, M) * yTop(i) - psi(i, M)*H) / (2 * delta_y + 2*yTop(i) - H);
end

%Left Wall BC
for j = 1:M+1
    psi(1, j) = -V_inf * (H / 2) + V_inf * (j - 1) * delta_y;
end

%Right Wall
for j = 1:M+1
    psi(N+1, j) = psi(N, j);
end

%Lower Boundary BC
for i = 1:N+1
    psi(i, 1) = (-V_inf * H * delta_y - 2 * psi(i, 2) * yBot(i) - H*psi(i,2)) / (2 * delta_y - 2*yBot(i) - H);
end
end

function Cp_topsurf = calcTopSurf(psi, N, M, delta_x, delta_y, V_inf)
for i = 1:N+1
    if i > 1 && i < N + 1
        u(i) = (psi(i, M+1) - psi(i, M)) / delta_y;
        v(i) = -(psi(i, M+1) - psi(i-1, M+1)) / (delta_x);
%         v(i) = (psi(i+1, M+1) - psi(i - 1, M+1)) / (2 * delta_x); %test to see
%         if central difference yields same results, which it does
    else
        u(i) = V_inf;  %only a horizontal component
        v(i) = 0;      %no normal (y) component
    end
end

Cp_topsurf = 1 - (u.^2 + v.^2) / V_inf^2;

end

function Cp_botsurf = calcBotSurf(psi, N, delta_x, delta_y, V_inf)
for i = 1:N+1
    if i > 1 && i < N + 1
        u(i) = (psi(i, 2) - psi(i, 1)) / delta_y;
        v(i) = -(psi(i, 1) - psi(i-1, 1)) / (delta_x);
%         v(i) = -(psi(i+1, 1) - psi(i-1, 1)) / (2 * delta_x); %test to see
%         if central difference yields same results, which it does
    else
        u(i) = V_inf;  %only a horizontal component
        v(i) = 0;      %no normal (y) component
    end
end

Cp_botsurf = 1 - (u.^2 + v.^2) / V_inf^2;

end

function Cp = calcPress(psi, cp_botsurf, cp_topsurf, V_inf, L, H, N, M)
delta_x = L / N;                            %Grid size in x direction in [m]
delta_y = H / M;                            %Grid size in y direction in [m]

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

%Right wall
for j = 1:M+1
    u(N+1, j) = V_inf;
    v(N+1, j) = 0;
end

for i = 1:N+1
    for j = 1:M+1
        if j == 1
            Cp(i, j) = cp_botsurf(i);
        elseif j == M+1
            Cp(i, j) = cp_topsurf(i);
        else
            Cp(i, j) = 1 - (u(i,j)^2 + v(i,j)^2)/(V_inf^2);
        end
    end
end
end


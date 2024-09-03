%==========================================================================
% AUTHOR: David L. Tran
%
% DESCRIPTION: An elliptic grid is generated for the NACA 0012 airfoil
% with P(zeta, eta) = Q(zeta, eta) = 0. A uniformly generated grid is used
% as an initial guess and two two types of boundaries are generated: (1) flat
% boundaries with upper half of NACA 0012 airfoil and (2) asymmetrical wavy
% walls. To solve the nonlinear coupled elliptic grid generation equations,
% the Gauss-Seidel method is employed: just calculated values of x_ij^(n+1)
% and y_ij^(n+1) replace the old values x_ij^(n) and y_ij^(n), where (n) is
% the iteration counter. Then, we implement an exponential grid method for
% the same two aforementioned boundaries to compare the two grid generation
% techniques. Newton's method is used to solve for alpha, a parameter which
% controls grid spacing.
%
%==========================================================================

%% Clear Cache
clc; close all; clearvars;

%% Variables

epsilon = 1e-5;         %convergence criteria
MAX_IT = 2000;          %maximum number of iterations in loop before termination

%NACA 0012 airfoil characteristics
t = 0.12;               %maximum airfoil thickness
IE = 21;
IT = 41;
IL = 61;
JL = 31;
JL_wavy = 41;

x_start = 0;            %start coordinate for x for NACA 0012
x_end = 1;              %end coordinate for x for NACA 0012
H = 2;                  %height of physical domain
y_0 = 0;                %start coordinate for y

                        
Delta_ymin_H = 0.005;   %constant for exponential grid stretching function in J direction (Delta y_min/H)

y_naca0012_u = @(x) (t / 0.2) .* (0.2969 .* sqrt(x) - 0.1260 .* x - 0.3537 .* x.^2 + 0.2843 .* x.^3 - 0.1015 .* x.^4);
                        %equation to plot upper half of NACA 0012 airfoil
y_u = @(x) (7 + cos(4 .* pi .* x)) ./ 4;
                        %equation to plot the upper wavy wall
y_l = @(x) (1 - cos(2 .* pi .* x)) ./ 4;
                        %equation to plot the lower wavy wall

%Flat wall
x_naca0012_pos = linspace(x_start, x_end, IT-IE+1);
                        %array for x coordinates of NACA 0012 airfoil
y_naca0012 = y_naca0012_u(x_naca0012_pos);
                        %array for y coordinates of upperhalf of NACA 0012 airfoil
delta_x = x_naca0012_pos(2) - x_naca0012_pos(1);
                        %grid spacing in x direction
x_left = x_start - (IE - 1) * delta_x;
                        %end position to the left of the airfoil
x_right = x_end + (IL - IT) * delta_x;
                        %end position to the right of the airfoil
x_left_arr = linspace(x_left, x_start, IE);
                        %array for the domain to the left of the airfoil
x_right_arr = linspace(x_end, x_right, IL-IT+1);
                        %array for the domain to the right of the airfoil 
x_domain = [x_left_arr(1:end-1), x_naca0012_pos, x_right_arr(2:end)];
                        %array for x coordinates of physical domain
y_1 = [zeros(IE-1,1); y_naca0012'; zeros(IL-IT,1)];                        
y_domain = linspace(0, H, JL);
                        %array for y coordinates of physical domain
y_domain_wavy = linspace(0, H, JL_wavy);

%Wavy Wall
x_start_wavy = 0;       %start x-coordinate of wavy wall
L = 1;                  %length of wavy domain
x_domain_wavy = linspace(x_start_wavy, L, IL);
y_domain_lower = y_l(x_domain_wavy);
y_domain_upper = y_u(x_domain_wavy);

F_alpha = @(alpha) Delta_ymin_H - (exp(alpha / JL) - 1) / (exp(alpha) - 1);
F_prime_alpha = @(alpha) -(JL * exp(alpha) - (JL - 1) * exp((JL + 1) * alpha / (JL)) - exp(alpha / JL)) / (JL * (exp(alpha) -1)^2);

F_alpha_wavy = @(alpha) Delta_ymin_H - (exp(alpha / JL_wavy) - 1) / (exp(alpha) - 1);
F_prime_alpha_wavy = @(alpha) -(JL_wavy * exp(alpha) - (JL_wavy - 1) * exp((JL_wavy + 1) * alpha / (JL_wavy)) - exp(alpha / JL_wavy)) / (JL_wavy * (exp(alpha) - 1)^2);

alpha_guess = 0.9;
alpha_1 = [alpha_guess];
alpha_wavy = [alpha_guess];

delta_y = zeros(JL,1);
delta_y(1) = Delta_ymin_H * H;
delta_y_wavy = zeros(JL_wavy,1);
delta_y_wavy(1) = min(y_domain_upper) * Delta_ymin_H;

y_flat_exp = zeros(IL, JL);
y_flat_exp(:,1) = y_1;
y_flat_exp(:,end) = H;
y_wavy_exp = zeros(IL, JL_wavy);
y_wavy_exp(:,1) = y_domain_lower;
y_wavy_exp(:,end) = y_domain_upper;

%Elliptic Grid
x_elliptic = zeros(IL,JL);
y_elliptic = zeros(IL,JL);
delta_eta = 1;
delta_zeta = 1;

%BVs of x coordinates
x_elliptic(1,:) = x_left;       %Left Boundary
x_elliptic(IL,:) = x_right;     %Right Boundary
x_elliptic(:,1) = x_domain;     %Bottom Boundary
x_elliptic(:,JL) = x_domain;    %Top Boundary    

%BVs of y coordinates
y_elliptic(1,:) = y_domain;     %Left Boundary
y_elliptic(IL,:) = y_domain;    %Right Boundary
y_elliptic(:,1) = y_1';         %Bottom Boundary
y_elliptic(:,JL) = H;           %Top Boundary

%Initial Guess (Uniform Grid)
for i = 2:IL-1
    for j = 2:JL-1
        x_elliptic(i,j) = x_domain(i);
        if x_domain(i) >= x_start && x_domain(i) <= x_end
            y_elliptic(i,j) = y_domain(j) + y_naca0012(i - IE + 1);
        else
            y_elliptic(i,j) = y_domain(j);
        end
    end
end

%Wavy Elliptic Grid
x_elliptic_wavy = zeros(IL,JL_wavy);
y_elliptic_wavy = zeros(IL,JL_wavy);

%BVs of x coordinates (Wavy)
x_elliptic_wavy(1,:) = x_left;       %Left Boundary
x_elliptic_wavy(IL,:) = x_right;     %Right Boundary
x_elliptic_wavy(:,1) = x_domain;     %Bottom Boundary
x_elliptic_wavy(:,JL_wavy) = x_domain;    %Top Boundary


%BVs of y coordinates (Wavy)
y_elliptic_wavy(1,:) = y_domain_wavy;     %Left Boundary
y_elliptic_wavy(IL,:) = y_domain_wavy;    %Right Boundary
y_elliptic_wavy(:,1) = y_domain_lower;    %Bottom Boundary
y_elliptic_wavy(:,JL_wavy) = y_domain_upper;           
                                          %Top Boundary                           

%% Loops

%Elliptic Grid Generation
tol_x = 1;
tol_y = 1;
ellip_count = 1;
while ((tol_x > epsilon) || (tol_y > epsilon)) && (ellip_count < MAX_IT)
    %set temporary storing variable to calculate residual
    x_temp = x_elliptic;
    y_temp = y_elliptic;
    for j = 2:JL-1
        for i = 2:IL-1
            alpha_n = ( (x_elliptic(i, j+1) - x_elliptic(i,j-1)) / (2) )^2 + ( (y_elliptic(i, j+1) - y_elliptic(i,j-1)) / (2) )^2;
            beta_n = ((x_elliptic(i+1,j) - x_elliptic(i-1,j)) / (2)) * ((x_elliptic(i, j+1) - x_elliptic(i, j-1)) / (2)) + ...
                ((y_elliptic(i+1,j) - y_elliptic(i-1,j)) / (2)) * ((y_elliptic(i, j+1) - y_elliptic(i, j-1)) / (2));
            gamma_n = ((x_elliptic(i+1,j) - x_elliptic(i-1,j)) / (2))^2 + ((y_elliptic(i+1, j) - y_elliptic(i-1, j)) / (2))^2;
            x_zetazeta = (x_elliptic(i+1,j) - 2 * x_elliptic(i,j) + x_elliptic(i-1,j)) / (delta_zeta^2);
            x_zetaeta = (x_elliptic(i+1,j+1) - x_elliptic(i-1,j+1) - x_elliptic(i+1, j-1) + x_elliptic(i-1,j-1)) / (4 * delta_zeta * delta_eta);
            x_etaeta = (x_elliptic(i,j+1) - 2 * x_elliptic(i,j) + x_elliptic(i,j-1)) / (delta_eta^2);
            y_zetazeta = (y_elliptic(i+1,j) - 2 * y_elliptic(i,j) + y_elliptic(i-1,j)) / (delta_zeta^2);
            y_zetaeta = (y_elliptic(i+1,j+1) - y_elliptic(i-1,j+1) - y_elliptic(i+1, j-1) + y_elliptic(i-1,j-1)) / (4 * delta_zeta * delta_eta);
            y_etaeta = (y_elliptic(i,j+1) - 2 * y_elliptic(i,j) + y_elliptic(i,j-1)) / (delta_eta^2);

            %immediately store new calculated values for x_ij and y_ij into
            %loop and overwrite old values by adding resid_x and resid_y to
            %the old values
            x_elliptic(i,j) = (1 /(-2*(alpha_n + gamma_n))) * (-alpha_n * (x_elliptic(i+1,j) + x_elliptic(i-1,j)) + 2 * beta_n * ((x_elliptic(i+1,j+1) - x_elliptic(i-1,j+1) - x_elliptic(i+1,j-1) + x_elliptic(i-1,j-1)) / 4) - gamma_n * (x_elliptic(i,j+1) + x_elliptic(i,j-1)));
            y_elliptic(i,j) = (1 /(-2*(alpha_n + gamma_n))) * (-alpha_n * (y_elliptic(i+1,j) + y_elliptic(i-1,j)) + 2 * beta_n * ((y_elliptic(i+1,j+1) - y_elliptic(i-1,j+1) - y_elliptic(i+1,j-1) + y_elliptic(i-1,j-1)) / 4) - gamma_n * (y_elliptic(i,j+1) + y_elliptic(i,j-1)));
            resid_x(i-1,j-1) = abs(x_elliptic(i,j) - x_temp(i,j));
            resid_y(i-1,j-1) = abs(y_elliptic(i,j) - y_temp(i,j));
        end
               
    end
    %calculate residual
    tol_x = max(abs(resid_x), [], "all");
    tol_y = max(abs(resid_y), [], "all");
    ellip_count = ellip_count + 1;
end

tol_x_wavy = 1;
tol_y_wavy = 1;
ellip_count_wavy = 1;
while ((tol_x_wavy > epsilon) || (tol_y_wavy > epsilon)) && (ellip_count_wavy < MAX_IT)
    %set temporary storing variable to calculate residual
    x_temp = x_elliptic_wavy;
    y_temp = y_elliptic_wavy;
    for j = 2:JL_wavy-1
        for i = 2:IL-1
            alpha_n = ( (x_elliptic_wavy(i, j+1) - x_elliptic_wavy(i,j-1)) / (2) )^2 + ( (y_elliptic_wavy(i, j+1) - y_elliptic_wavy(i,j-1)) / (2) )^2;
            beta_n = ((x_elliptic_wavy(i+1,j) - x_elliptic_wavy(i-1,j)) / (2)) * ((x_elliptic_wavy(i, j+1) - x_elliptic_wavy(i, j-1)) / (2)) + ...
                ((y_elliptic_wavy(i+1,j) - y_elliptic_wavy(i-1,j)) / (2)) * ((y_elliptic_wavy(i, j+1) - y_elliptic_wavy(i, j-1)) / (2));
            gamma_n = ((x_elliptic_wavy(i+1,j) - x_elliptic_wavy(i-1,j)) / (2))^2 + ((y_elliptic_wavy(i+1, j) - y_elliptic_wavy(i-1, j)) / (2))^2;
            x_zetazeta = (x_elliptic_wavy(i+1,j) - 2 * x_elliptic_wavy(i,j) + x_elliptic_wavy(i-1,j)) / (delta_zeta^2);
            x_zetaeta = (x_elliptic_wavy(i+1,j+1) - x_elliptic_wavy(i-1,j+1) - x_elliptic_wavy(i+1, j-1) + x_elliptic_wavy(i-1,j-1)) / (4 * delta_zeta * delta_eta);
            x_etaeta = (x_elliptic_wavy(i,j+1) - 2 * x_elliptic_wavy(i,j) + x_elliptic_wavy(i,j-1)) / (delta_eta^2);
            y_zetazeta = (y_elliptic_wavy(i+1,j) - 2 * y_elliptic_wavy(i,j) + y_elliptic_wavy(i-1,j)) / (delta_zeta^2);
            y_zetaeta = (y_elliptic_wavy(i+1,j+1) - y_elliptic_wavy(i-1,j+1) - y_elliptic_wavy(i+1, j-1) + y_elliptic_wavy(i-1,j-1)) / (4 * delta_zeta * delta_eta);
            y_etaeta = (y_elliptic_wavy(i,j+1) - 2 * y_elliptic_wavy(i,j) + y_elliptic_wavy(i,j-1)) / (delta_eta^2);

            %immediately store new calculated values for x_ij and y_ij into
            %loop and overwrite old values by adding resid_x and resid_y to
            %the old values
            x_elliptic_wavy(i,j) = (1 /(-2*(alpha_n + gamma_n))) * (-alpha_n * (x_elliptic_wavy(i+1,j) + x_elliptic_wavy(i-1,j)) + 2 * beta_n * ((x_elliptic_wavy(i+1,j+1) - x_elliptic_wavy(i-1,j+1) - x_elliptic_wavy(i+1,j-1) + x_elliptic_wavy(i-1,j-1)) / 4) - gamma_n * (x_elliptic_wavy(i,j+1) + x_elliptic_wavy(i,j-1)));
            y_elliptic_wavy(i,j) = (1 /(-2*(alpha_n + gamma_n))) * (-alpha_n * (y_elliptic_wavy(i+1,j) + y_elliptic_wavy(i-1,j)) + 2 * beta_n * ((y_elliptic_wavy(i+1,j+1) - y_elliptic_wavy(i-1,j+1) - y_elliptic_wavy(i+1,j-1) + y_elliptic_wavy(i-1,j-1)) / 4) - gamma_n * (y_elliptic_wavy(i,j+1) + y_elliptic_wavy(i,j-1)));
            resid_x_wavy(i-1,j-1) = abs(x_elliptic_wavy(i,j) - x_temp(i,j));
            resid_y_wavy(i-1,j-1) = abs(y_elliptic_wavy(i,j) - y_temp(i,j));
        end
               
    end
    %calculate residual
    tol_x_wavy = max(abs(resid_x_wavy), [], "all");
    tol_y_wavy = max(abs(resid_y_wavy), [], "all");
    ellip_count_wavy = ellip_count_wavy + 1;
end


%Newton's method to solve for alpha (flat and wavy walls)
resid_exp_1 = abs(F_alpha(alpha_guess));
resid_exp_wavy = abs(F_alpha_wavy(alpha_guess));
count = 1;
count_wavy = 1;

while resid_exp_1 > epsilon && count_wavy < MAX_IT
    count = count + 1;

    alpha_1(count) = alpha_1(count-1) - ((F_alpha(alpha_1(count-1))) / (F_prime_alpha(alpha_1(count-1))));

    resid_exp_1 = abs(alpha_1(count) - alpha_1(count-1));
end

fprintf('Alpha for Flat Wall: %f\n', alpha_1(end));

while resid_exp_wavy > epsilon && count_wavy < MAX_IT
    count_wavy = count_wavy + 1;

    alpha_wavy(count_wavy) = alpha_wavy(count_wavy - 1) - (F_alpha_wavy(alpha_wavy(count_wavy-1)) ) / (F_prime_alpha_wavy(alpha_wavy(count_wavy-1)) );

    resid_exp_wavy = abs(alpha_wavy(count_wavy) - alpha_wavy(count_wavy-1));
end

fprintf('Alpha for Wavy Wall: %f\n', alpha_wavy(end));

% Loop for Delta_y of each point for exp grid
for j = 2:JL-1
    delta_y(j) = H * (exp(j * alpha_1(end) / JL) - 1) / (exp(alpha_1(end))-1);
    for i = 1:IL
        
        if x_domain(i) >= x_start && x_domain(i) <= x_end

            y_flat_exp(i,j) = y_naca0012(i - IE + 1) + delta_y(j-1);
        else
            y_flat_exp(i,j) =  delta_y(j-1);
        end
    end
end

%replace 

for j = 2:JL_wavy-1
    delta_y_wavy(j) = min(y_domain_upper) * (exp(j * alpha_wavy(end) / JL_wavy) - 1) / (exp(alpha_wavy(end)));    
    for i = 1:IL
        y_wavy_exp(i,j) = y_domain_lower(i) + delta_y_wavy(j-1);  
    end
end


%% Plots/Figures

%Elliptic grid (flat wall)
figure;
hold on;
set(gcf,'Position',[0 25 1200 800]);
for i = 1:IL
%     plot(x_elliptic(i,:),y_flat_exp(:,i), 'b');

    for j = 1:JL
            x_plot_e(i,j) = x_elliptic(i,j);
            y_plot_e(i,j) = y_elliptic(i,j);
%             plot(x_plot, y_flat_exp(i,:), 'k');
    
    end
    plot(x_plot_e(i,:), y_plot_e(i,:),'k');
end
for i = 1:JL
    plot(x_plot_e(:,i), y_plot_e(:,i),'b');
end
title("Flat Wall (Elliptic Grid)",'Interpreter','LaTeX');               
xlabel('$x$','Interpreter','LaTeX');
ylabel('$y$','Interpreter','LaTeX');  
set(gca,'FontSize',18);

%Elliptic grid (wavy wall)
figure;
hold on;
set(gcf,'Position',[0 25 1200 800]);
for i = 1:IL
    for j = 1:JL_wavy
            x_plot_e(i,j) = x_elliptic_wavy(i,j);
            y_plot_e(i,j) = y_elliptic_wavy(i,j);
%             plot(x_plot, y_flat_exp(i,:), 'k');
    
    end
    plot(x_plot_e(i,:), y_plot_e(i,:),'k');
end
for i = 1:JL_wavy
    plot(x_plot_e(:,i), y_plot_e(:,i),'b');
end
title("Wavy Wall (Elliptic Grid)",'Interpreter','LaTeX');               
xlabel('$x$','Interpreter','LaTeX');
ylabel('$y$','Interpreter','LaTeX');  
set(gca,'FontSize',18);


%Exponential grid (flat wall)
figure;
hold on;
set(gcf,'Position',[0 25 1200 800]);
for i = 1:JL
    plot(x_domain,y_flat_exp(:,i), 'b');
end
for i = 1:IL
    if x_domain(i) > x_start && x_domain(i) < x_end
        x_plot = x_domain(i) .* ones(JL,1);
        plot(x_plot, y_flat_exp(i,:), 'k');
    else    
        xline(x_domain(i));
    end
end
title("Flat Wall (Exponential Grid)",'Interpreter','LaTeX');               
xlabel('$x$','Interpreter','LaTeX');
ylabel('$y$','Interpreter','LaTeX');  
set(gca,'FontSize',18);

%Exponential grid (wavy wall)
figure;
hold on;
set(gcf,'Position',[0 25 1200 800]);
for i = 1:JL_wavy
    plot(x_domain_wavy,y_wavy_exp(:,i), 'b');
end
for i = 1:IL
        x_plot = x_domain_wavy(i) .* ones(JL_wavy,1);
        plot(x_plot, y_wavy_exp(i,:), 'k');
end
title("Wavy Wall (Exponential Grid)",'Interpreter','LaTeX');               
xlabel('$x$','Interpreter','LaTeX');
ylabel('$y$','Interpreter','LaTeX');  
set(gca,'FontSize',18);
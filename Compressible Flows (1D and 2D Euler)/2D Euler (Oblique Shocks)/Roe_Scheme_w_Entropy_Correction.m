%==========================================================================
% AUTHOR: David L. Tran
%
% Roe's method with sonic-point entropy correction ONLY for this script
%
% DESCRIPTION: Solves Roe's scheme with a sonic-pt entropy correction for 
% a 2D supersonic engine inlet with prescribed supersonic inlet conditions. 
% A grid is generated using the differential equation method which uses an
% algebraically generated grid as an initial guess. The numerical
% approximation is repeated until convergence within a set tolerance is
% desired.
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

%Gas parameters
gamma = 1.4;    %heat capacity ratio
R = 287.1;      %specific gas constant of air in [J/(kg*K)]

%Engine Geometry/Dimensions
H = 1;                  %Height of the inlet in [m]
L = 3.2;                %Length of the section in which top wall converges in [m]

%Inlet conditions
p_1 = 1e5;                                %pressure at inlet in [Pa]
rho_1 = 1;                                %density in [kg/m^3]
M_1 = 2.9;                                %inlet Mach number

%CFL number
nu = 0.75;       %Courant number for simulation

% Grid 1 (uncomment when in use)

% IL = 42;      %number of control volumes in x-dir.
% JL = 22;      %number of control volumes in y-dir.
% IS = 6;       %control volume in which engine inlet begins to converge

%Grid 2 (uncomment when in use)

% IL = 82;      %number of control volumes in x-dir.
% JL = 42;      %number of control volumes in y-dir.
% IS = 12;       %control volume in which engine inlet begins to converge

%Grid 3 (uncomment when in use)

IL = 162;     %number of control volumes in x-dir.
JL = 82;      %number of control volumes in y-dir.
IS = 24;      %control volume in which engine inlet begins to converge

%Solver Criterion

tol = 1e-5;     %convergence criteria
IT_MAX = 3000;  %maximum number of iterations before program termination
count = 1;      %counter for number of iterations run in loop

%Roe arrays
nxN=zeros(IL,JL);
nxS=zeros(IL,JL);
nxW=zeros(IL,JL);
nxE=zeros(IL,JL);
nyN=zeros(IL,JL);
nyS=zeros(IL,JL);
nyW=zeros(IL,JL);
nyE=zeros(IL,JL);
V=zeros(IL,JL);
dx=zeros(IL,JL);
dy=zeros(IL,JL);
cx=zeros(IL,JL);
cy=zeros(IL,JL);
dt1=zeros(IT_MAX,1);
sN=zeros(IL,JL);
sS=zeros(IL,JL);
sW=zeros(IL,JL);
sE=zeros(IL,JL);
U=zeros(IL,JL);
unew=zeros(IL,JL);
unew1=zeros(IL,JL);
unew2=zeros(IL,JL);
unew3=zeros(IL,JL);
unew4=zeros(IL,JL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                          %
%                     Main Body/Loops                      %
%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Euler equations                    %%
%%           U_t+E_x+F_y=0                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Conservative variables
%u1:density (rho)
%u2:x-momentum
%u3:y-momentum
%u4:energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Flux variables in terms of conservative variables   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%E vector
%e1=rho*u
%e2=rho*u^2+P
%e3=rho*u*v
%e4=(e+P)*u
%F vector
%f1=rho*v
%f2=rho*u*v
%f3=rho*v^2+P
%f4=(e+P)*v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize initial conditions to U arrays

[u1,u2,u3,u4] = initcond(IL, JL, gamma, p_1, rho_1, M_1);
p = pressure(gamma,u1,u2,u3,u4);  %Pressure
u = u_vel(u1,u2);                 %x-velocity
v = v_vel(u1,u3);                 %y-velocity
h = (u4 + p) ./ u1;               %enthalpy

%Generate grid
[x,y,t1] = grid_gen4(IL, JL, IS, L, H);

%Calculate unit normals and volumes of the cells
for j=1:JL
   for i=1:IL
      
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%   Unit normals at each face of control volume       %%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       [nxE(i,j),nyE(i,j),sE(i,j)]=norm_comp_x(y(i+1,j+1),y(i+1,j)...
           ,x(i+1,j+1),x(i+1,j));
       [nxW(i,j),nyW(i,j),sW(i,j)]=norm_comp_x(y(i,j+1),y(i,j)...
           ,x(i,j+1),x(i,j));
       [nxN(i,j),nyN(i,j),sN(i,j)]=norm_comp_y(y(i,j+1),y(i+1,j+1)...
           ,x(i,j+1),x(i+1,j+1));
       [nxS(i,j),nyS(i,j),sS(i,j)]=norm_comp_y(y(i,j),y(i+1,j)...
           ,x(i,j),x(i+1,j));
      
      
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%             Area of control volume  V         %%
       %%          change in x (dx),change in y (dy)    %%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       [V(i,j),dx(i,j),dy(i,j),cx(i,j),cy(i,j)]=volume(x(i,j),x(i+1,j),x(i+1,j+1),x(i,j+1)...
           ,y(i,j),y(i+1,j),y(i+1,j+1),y(i,j+1));
  
  
   end
end

dx1=max(max(dx));
dy1=max(max(dy));

%Time marching loop until steady-state
err = 1; %to enter while loop

while (count < IT_MAX) && (err > tol)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%      Boundary Condition     %%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [u1,u2,u3,u4,p,u,v]=...
       boundarycondition(IL,JL,IS,t1,...
       gamma,u1,u,v,p);
    c=sos(gamma,p,u1);              %Speed of sound
    h = (u4 + p) ./ u1;               %enthalpy
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %      Time step of each C.V.  %%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   dt=nu./((abs(u)./dx)+(abs(v)./dy)...
       +c.*sqrt((1./dx.^2)+(1./dy.^2))); %Time step determined by CFL cond.
  
   dt1(count)=min(min(dt));
  
  
   [e1,e2,e3,e4]=...
       Eflux(gamma,u1,u2,u3,u4);   %Flux vectors in conservative variables
   [f1,f2,f3,f4]=...
       Fflux(gamma,u1,u2,u3,u4);   %Flux vectors in conservative variables
  
   for j=2:JL-1
       for i=2:IL-1
           %Calculate Roe's averages for each of the four faces
           [rho_bar_iph, u_bar_iph, v_bar_iph, h_bar_iph, c_bar_iph] = RoeAvg(u1(i,j), u1(i+1,j), u(i,j), u(i+1,j), v(i,j), v(i+1,j), h(i,j), h(i+1,j), gamma);
           [rho_bar_imh, u_bar_imh, v_bar_imh, h_bar_imh, c_bar_imh] = RoeAvg(u1(i-1,j), u1(i,j), u(i-1,j), u(i,j), v(i-1,j), v(i,j), h(i-1,j), h(i,j), gamma);
           [rho_bar_jph, u_bar_jph, v_bar_jph, h_bar_jph, c_bar_jph] = RoeAvg(u1(i,j), u1(i,j+1), u(i,j), u(i,j+1), v(i,j), v(i,j+1), h(i,j), h(i,j+1), gamma);
           [rho_bar_jmh, u_bar_jmh, v_bar_jmh, h_bar_jmh, c_bar_jmh] = RoeAvg(u1(i,j-1), u1(i,j), u(i,j-1), u(i,j), v(i,j-1), v(i,j), h(i,j-1), h(i,j), gamma);
           
           %Calculate |A| for each of the four faces using Roe's averages
           A_iph = calcA(u_bar_iph, v_bar_iph, rho_bar_iph, c_bar_iph, nxE(i,j), nyE(i,j), gamma); %i+1/2
           A_imh = calcA(u_bar_imh, v_bar_imh, rho_bar_imh, c_bar_imh, nxW(i,j), nyW(i,j), gamma); %i-1/2
           A_jph = calcA(u_bar_jph, v_bar_jph, rho_bar_jph, c_bar_jph, nxN(i,j), nyN(i,j), gamma); %j+1/2
           A_jmh = calcA(u_bar_jmh, v_bar_jmh, rho_bar_jmh, c_bar_jmh, nxS(i,j), nyS(i,j), gamma); %j-1/2

           %Calculate prime fluxes of each point in stencil     
           Eprime_ip1 = Fluxprime(e1(i+1,j),e2(i+1,j),e3(i+1,j),e4(i+1,j)...
               ,f1(i+1,j),f2(i+1,j),f3(i+1,j),f4(i+1,j),nxE(i,j),nyE(i,j));
           Eprime_w_ip1 = Fluxprime(e1(i,j),e2(i,j),e3(i,j),e4(i,j)...
               ,f1(i,j),f2(i,j),f3(i,j),f4(i,j),nxE(i,j),nyE(i,j)); %E'' with i+1
           Eprime_w_im1 = Fluxprime(e1(i,j),e2(i,j),e3(i,j),e4(i,j)...
               ,f1(i,j),f2(i,j),f3(i,j),f4(i,j),nxW(i,j),nyW(i,j)); %E'' with i-1
           Eprime_im1 = Fluxprime(e1(i-1,j),e2(i-1,j),e3(i-1,j),e4(i-1,j)...
               ,f1(i-1,j),f2(i-1,j),f3(i-1,j),f4(i-1,j),nxW(i,j),nyW(i,j));
           Fprime_jp1=Fluxprime(e1(i,j+1),e2(i,j+1),e3(i,j+1),e4(i,j+1)...
               ,f1(i,j+1),f2(i,j+1),f3(i,j+1),f4(i,j+1),nxN(i,j),nyN(i,j));
           Fprime_w_jp1=Fluxprime(e1(i,j),e2(i,j),e3(i,j),e4(i,j)...
               ,f1(i,j),f2(i,j),f3(i,j),f4(i,j),nxN(i,j),nyN(i,j)); %F'' with j+1
           Fprime_w_jm1=Fluxprime(e1(i,j),e2(i,j),e3(i,j),e4(i,j)...
               ,f1(i,j),f2(i,j),f3(i,j),f4(i,j),nxS(i,j),nyS(i,j)); %F'' with j-1
           Fprime_jm1=Fluxprime(e1(i,j-1),e2(i,j-1),e3(i,j-1),e4(i,j-1)...
               ,f1(i,j-1),f2(i,j-1),f3(i,j-1),f4(i,j-1),nxS(i,j),nyS(i,j));
          
           %Store U and store the i-1, i+1, j-1, and j+1 values
           U = [u1(i,j);u2(i,j);u3(i,j);u4(i,j)];
           UIp=[u1(i+1,j);u2(i+1,j);u3(i+1,j);u4(i+1,j)];
           UIm=[u1(i-1,j);u2(i-1,j);u3(i-1,j);u4(i-1,j)];
           UJp=[u1(i,j+1);u2(i,j+1);u3(i,j+1);u4(i,j+1)];
           UJm=[u1(i,j-1);u2(i,j-1);u3(i,j-1);u4(i,j-1)];

           %Store the S terms for each face
           S=[sE(i,j);sW(i,j);sN(i,j);sS(i,j)];

           %Store the u and v velocity components
           uVel=[u(i,j);u(i+1,j);u(i,j+1)];
           vVel=[v(i,j);v(i+1,j);v(i,j+1)];
           C=[c(i,j);c(i+1,j);c(i,j+1)];

           %normal vectors of each face
           Nor=[nxE(i,j);nyE(i,j);nxW(i,j);nyW(i,j);nxN(i,j);nyN(i,j);...
               nxS(i,j);nyS(i,j)];
              
           %Calculate Roe's fluxes
           E_iph_Roe = RoeFlux(Eprime_w_ip1, Eprime_ip1, U, UIp, A_iph);
           E_imh_Roe = RoeFlux(Eprime_im1, Eprime_w_im1, UIm, U, A_imh);
           F_jph_Roe = RoeFlux(Fprime_w_jp1, Fprime_jp1, U, UJp, A_jph);
           F_jmh_Roe = RoeFlux(Fprime_jm1, Fprime_w_jm1, UJm, U, A_jmh);

           %Update the step
           [U] = Roe(dt1(count), V(i,j), S, U, E_imh_Roe, E_iph_Roe, F_jmh_Roe, F_jph_Roe);
           
           %Store just calculated values in new array
           unew1(i,j)=U(1);
           unew2(i,j)=U(2);
           unew3(i,j)=U(3);
           unew4(i,j)=U(4);
       end
   end   
 
   u1=unew1;
   u2=unew2;
   u3=unew3;
   u4=unew4;
  
   pold = p;                       %store old pressure
   p=pressure(gamma,u1,u2,u3,u4);  %Pressure
   u=u_vel(u1,u2);                 %x-velocity
   v=v_vel(u1,u3);                 %y-velocity
  
   c=sos(gamma,p,u1);              %Speed of sound

   %Update counter
   count = count + 1;

   %Calculate new error based on pressures
   err = max(max(abs((p-pold)./p)));
   err_arr(count-1) = err;
  
end

count_arr = 1:1:count-1;
M = sqrt(u.^2 + v.^2) ./ c;
T=sum(dt1);


%Exact shock calculation

shock_count = 1;

%Mach numbers
M_2 = 2.378;
M_3 = 1.942;
M_4 = 1.551;

M_exact = zeros(IL,1);

%Angles
beta_1 = 29;
beta_2 = 34.22;
beta_3 = 41.61;

beta_1_rad = beta_1 * pi/180;
beta_2_rad = beta_2 * pi/180;
beta_3_rad = beta_3 * pi/180;

%collision point of oblique shock 1 with lower wall
x_shock_1 = H / tan(beta_1_rad);

for i = 2:IL-1

    if shock_count == 1 && y(i,floor(JL/4)) >= H - x(i,floor(JL/4)) * tan(beta_1_rad)
        shock_count = shock_count + 1;
        M_exact(1:i,1) = M_1;
        i_shock_1 = i;
        fprintf('%d\n', i);

    elseif shock_count == 2 && (y(i,floor(JL/4)) <= (x(i,floor(JL/4)) - x_shock_1) * tan(beta_2_rad)) && (x(i,floor(JL/4)) >= x_shock_1)   
        shock_count = shock_count + 1;
        M_exact(i_shock_1+1:i,1) = M_2;
        i_shock_2 = i;
        x_shock_2 = x(i,floor(JL/4));
        fprintf('%d', i);
% 
%     elseif shock_count == 3 && (y(i,JL/2) <= H - x(i,JL/2) * tan(beta_1)) && (x(i,JL/2) >= x_shock_2)   
%         shock_count = shock_count + 1;
%         M_exact(i_shock_1+1:i,1) = M_2;    

    end
    

end

M_exact(i_shock_2:end,1) = M_3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                          %
%                      Plots/Figures                       %
%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on;
contourf(cx,cy,p);
c = colorbar;
c.Label.String = 'Pressure (Pa)';
c.Label.Interpreter = 'latex';
axis equal;
set(gcf,'Position',[0 25 1200 600]);
title('Roe Flux-Difference Splitting Scheme w/ Sonic-Point Entropy Correction','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL) ', $JL = $ ' num2str(JL) ', and $IS = $ ' num2str(IS)],'Interpreter','LaTeX')
xlabel('Horizontal Position $x$ (m)','Interpreter','LaTeX');
ylabel('Vertical Position $y$ (m)','Interpreter','LaTeX'); 
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
contourf(cx,cy,M);
c = colorbar;
c.Label.String = 'Mach Number';
c.Label.Interpreter = 'latex';
axis equal;
set(gcf,'Position',[0 25 1200 600]);
title('Roe Flux-Difference Splitting Scheme w/ Sonic-Point Entropy Correction','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL) ', $JL = $ ' num2str(JL) ', and $IS = $ ' num2str(IS)],'Interpreter','LaTeX')
xlabel('Horizontal Position $x$ (m)','Interpreter','LaTeX');
ylabel('Vertical Position $y$ (m)','Interpreter','LaTeX'); 
set(gca,'LineWidth',2,'FontSize',18);

figure;
semilogy(count_arr, err_arr, '-b', 'LineWidth', 2);
hold;
set(gcf,'Position',[0 25 1200 600]);
title('Roe Flux-Difference Splitting Scheme w/ Sonic-Point Entropy Correction','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL) ', $JL = $ ' num2str(JL) ', and $IS = $ ' num2str(IS)],'Interpreter','LaTeX');
xlabel('Iteration Count','Interpreter','LaTeX');
ylabel('Residual','Interpreter','LaTeX'); 
set(gca,'LineWidth',2,'FontSize',18);

figure;
hold on;
plot(x(2:IL-1,floor(JL/4)), M(2:IL-1,floor(JL/4)), '-b', 'LineWidth', 2);
plot(x(2:end,floor(JL/4)), M_exact, '-g', 'LineWidth', 2);
axis equal;
set(gcf,'Position',[0 25 1200 600]);
title('Roe Flux-Difference Splitting Scheme w/ Sonic-Point Entropy Correction','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL) ', $JL = $ ' num2str(JL) ', and $IS = $ ' num2str(IS)],'Interpreter','LaTeX')
xlabel('Horizontal Position $x$ (m)','Interpreter','LaTeX');
ylabel('Mach Number','Interpreter','LaTeX'); 
legend('Approximation', 'Exact', 'Interpreter', 'latex');
set(gca,'LineWidth',2,'FontSize',18);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                          %
%                        Functions                         %
%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Eprime=Fluxprime(e1,e2,e3,e4,f1,f2,f3,f4,nx,ny)
E=[e1;e2;e3;e4];
F=[f1;f2;f3;f4];
Eprime=E*nx+F*ny;
end

function [newFlux] = RoeFlux(Flux_i, Flux_ip1, U_i, U_ip1, A)
%This function calculates the Roe flux.
newFlux = 0.5 * (Flux_i + Flux_ip1) - 0.5 * A * (U_ip1 - U_i);
end

function [rho_bar, u_bar, v_bar, h_bar, c_bar] = RoeAvg(rho, rho_p1, u, u_p1, v, v_p1, h, h_p1, gamma)
%This function calculates the various Roe's averages for an arbitrary state
%and the plus 1 state (can be i-1 and i, i and i+1, j-1 and j, & j and j+1)
R_avg = sqrt(rho_p1 / rho);
rho_bar = R_avg * rho;
u_bar = (u + R_avg * u_p1) / (1 + R_avg);
v_bar = (v + R_avg * v_p1) / (1 + R_avg);
h_bar = (h + R_avg * h_p1) / (1 + R_avg);
c_bar = sqrt( (gamma - 1) * (h_bar - (u_bar^2 + v_bar^2) / 2 ) );
end

function [A] = calcA(u, v, rho, c, n_x, n_y, gamma)
%Calculates |A| for the Roe scheme

eps_corr = 0.06; %sonic-point entropy correction factor

%Auxiliary variables
alpha = (u^2 + v^2) / 2;
beta = gamma - 1;
S_x = n_x;
S_y = n_y;
u_prime = S_x * u + S_y * v;
v_prime = -S_y * u + S_x * v;

%populate row 1 of P
P(1,1) = 1;
P(1,2) = 1 / (2*c^2);
P(1,3) = 0;
P(1,4) = 1 / (2*c^2);
%populate row 2 of P
P(2,1) = u;
P(2,2) = u / (2*c^2) + S_x / (2*c);
P(2,3) = -S_y * rho;
P(2,4) = u / (2*c^2) - S_x / (2*c);
%populate row 3 of P
P(3,1) = v;
P(3,2) = v / (2*c^2) + S_y / (2*c);
P(3,3) = S_x * rho;
P(3,4) = v / (2*c^2) - S_y / (2*c);
%populate row 4 of P
P(4,1) = alpha;
P(4,2) = ( alpha / (2*c^2) ) + ( u_prime / (2*c) ) + ( 1 / (2*beta) );
P(4,3) = rho * v_prime;
P(4,4) = ( alpha / (2*c^2) ) - ( u_prime / (2*c) ) + ( 1 / (2*beta) );

%populate row 1 of P^{-1}
P_inv(1,1) = 1 - (alpha * beta) / (c^2);
P_inv(1,2) = u * beta / (c^2);
P_inv(1,3) = v * beta / (c^2);
P_inv(1,4) = -beta / (c^2);
%populate row 2 of P^{-1}
P_inv(2,1) = alpha * beta - u_prime * c;
P_inv(2,2) = -u * beta + S_x * c;
P_inv(2,3) = -v * beta + S_y * c;
P_inv(2,4) = beta;
%populate row 3 of P^{-1}
P_inv(3,1) = -v_prime / rho;
P_inv(3,2) = -S_y / rho;
P_inv(3,3) = S_x / rho;
P_inv(3,4) = 0;
%populate row 4 of P^{-1}
P_inv(4,1) = alpha * beta + u_prime * c;
P_inv(4,2) = -u * beta - S_x * c;
P_inv(4,3) = -v * beta - S_y * c;
P_inv(4,4) = beta;

%populate |Lambda|

Lambda(1,1) = abs(u_prime);
Lambda(2,2) = abs(u_prime + c);
Lambda(3,3) = abs(u_prime);
Lambda(4,4) = abs(u_prime - c);

if abs(Lambda(1,1)) < 2 * eps_corr
    %Perform correction for 1st and 3rd eigenvalues (identical)
    Lambda(1,1) = (Lambda(1,1)^2 + 4 * eps_corr^2) / (4 * eps_corr);
    Lambda(3,3) = (Lambda(3,3)^2 + 4 * eps_corr^2) / (4 * eps_corr);

elseif abs(Lambda(2,2)) < 2 * eps_corr
    %Perform correction for 2nd eigenvalue
    Lambda(2,2) = (Lambda(2,2)^2 + 4 * eps_corr^2) / (4 * eps_corr);    
       
elseif abs(Lambda(4,4)) < 2 * eps_corr
    %Perform correction for 4th eigenvalue
    Lambda(4,4) = (Lambda(4,4)^2 + 4 * eps_corr^2) / (4 * eps_corr);

end 

%Calculate |A|
A = P * Lambda * P_inv;

end

function [U_new] = Roe(dt, V, S, U, E_imh, E_iph, F_jmh, F_jph)
%This function updates the Roe solver
s1=S(1); %S_i+1/2
s2=S(2); %S_i-1/2
s3=S(3); %S_j+1/2
s4=S(4); %S_j-1/2
%Calculate updated step
U_new = U - dt / V .* ( (E_iph .* s1 - E_imh .* s2) + (F_jph .* s3 - F_jmh .* s4) );
end

function [u1, u2, u3, u4] = initcond(IL, JL, gamma, p_1, rho_1, M_1)
%This function assigns the initial conditions to the four components of the
%2D conservative vector, U = [rho, rho*u, rho*v, e].
%Inlet conditions
c = sos(gamma,p_1,rho_1);                 %speed of sound of inlet in [m/s]
u_1 = M_1*c;                              %initial horizontal velocity component in [m/s]
v_1 = 0;                                  %initial vertical velocity component in [m]
[e_1] = energy(gamma,rho_1,u_1,v_1,p_1);   %initial total energy per unit volume in [J/m^3]
%%Use Inlet Condition as Initial Condition
u1(1:IL,1:JL) = rho_1;
u2(1:IL,1:JL) = rho_1*u_1;
u3(1:IL,1:JL) = rho_1*v_1;
u4(1:IL,1:JL) = e_1;
end

function [u1,u2,u3,u4,p,u,v]=boundarycondition(IL,JL,IS,t1,gamma,u1,u,v,p)
%This function assigns the boundary conditions for the grid
%Inlet conditions
p1=10^5;
rho1=1;
V1=0;
M1=2.9;
c1 =sos(gamma,p1,rho1);
U1=M1*c1;
tht=-t1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Inlet           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
u1(1,2:JL-1)=rho1;
p(1,2:JL-1)=p1;
u(1,2:JL-1)=U1;
v(1,2:JL-1)=V1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Outlet             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Zero-order extrapolation where variables at the cell center of each cell
%in the last (JL) column are set to the adjacent cell (JL - 1) values.
u1(IL,2:JL-1)=u1(IL-1,2:JL-1);
p(IL,2:JL-1)=p(IL-1,2:JL-1);
u(IL,2:JL-1)=u(IL-1,2:JL-1);
v(IL,2:JL-1)=v(IL-1,2:JL-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Lower Wall         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set value at j = 1 to be equal to j = 2
u1(1:IL,1)=u1(1:IL,2);
p(1:IL,1)=p(1:IL,2);
u(1:IL,1)=u(1:IL,2);
v(1:IL,1)=-v(1:IL,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Top Wall           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Before initial shock
u1(1:IS-1,JL)=u1(1:IS-1,JL-1);
p(1:IS-1,JL)=p(1:IS-1,JL-1);
u(1:IS-1,JL)=u(1:IS-1,JL-1);
v(1:IS-1,JL)=-v(1:IS-1,JL-1);
%After initial shock
u1(IS:IL,JL)=u1(IS:IL,JL-1);
p(IS:IL,JL)=p(IS:IL,JL-1);
u(IS:IL,JL)=u(IS:IL,JL-1)*cos(2*tht)+v(IS:IL,JL-1)*sin(2*tht);
v(IS:IL,JL)=u(IS:IL,JL-1)*sin(2*tht)-v(IS:IL,JL-1)*cos(2*tht);
u2=u1.*u;
u3=u1.*v;
u4=energy(gamma,u1,u,v,p);
end

function [c] = sos(gamma, p, rho)
%This functions calculates the speed of sound, c, for a calorically perfect
%gas using INPUT PARAMS: (1) heat capacity ratio, gamma; (2) pressure, p;
%(3) density, rho.
c = sqrt(gamma .* p ./ rho);
end

function [p] = pressure(gamma, U_1, U_2, U_3, U_4)
%This functions calculates the pressure using the components of the 2D
%conservative Euler vector, U, which are U_1, U_2, U_3, and U_4
p = (gamma-1) .* ( U_4 - (U_2.^2 + U_3.^2) ./ (2.*U_1) );
end

function [e]= energy(gamma,rho,u,v,p)
%This function calculates the total energy per unit volume, e, using the
%various components of U which correspond to flow properties.
e = (p./(gamma-1))+((rho).*(u.^2+v.^2)/2);
end

function [x4,y4,t1]=grid_gen4(IL, JL, IS, L, H)
% IS:         %Control volume in which upper wall starts converging
% IL:         %Number of control volumes in x-direction
% JL:         %Number of control volumes in y-direction
dx = L/(IL-IS);         %x-direction grid spacing between control volumes in [m]
dy = H/(JL-2);          %y-drection grid spacing between the control volumes of inlet in [m]
tht = 10.94;            %Angle at which wall deflects in [deg]
t1 = tht*pi/180;        %Angle at which wall deflects in [rad]
H2 = H-(L)*tan(t1);     %Height of the outlet in [m]
dy2 = H2/(JL-2);        %y-drection grid spacing between the control volumes of outlet
dr = zeros(IL-IS+1,1);  %array for CVs in coverging section
%Initialize
x = zeros(IL-1,JL-1);
y = zeros(IL-1,JL-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Before Shock                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inlet
x1(1,1:JL-1)=-(IS-2)*dx;
y1(1,1:JL-1)=0:dy:H;
%Right bdy of the control volumes before converging wall
x1(IS-1,1:JL-1)=0;
y1(IS-1,1:JL-1)=0:dy:H;
%Top Wall
x1(1:IS-1,JL-1)=-(IS-2)*dx:dx:0;
y1(1:IS-1,JL-1)=H;
%Lower wall
x1(1:IS-1,1)=-(IS-2)*dx:dx:0;
y1(1:IS-1,1)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial guess to the interior points using algebric method %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=2:JL-2
   for i=2:IS-2
       x1(i,j)=x1(1,j)+(i-1)*dx;
   end
end
for i=2:IS-2
   for j=2:JL-2
       y1(i,j)=y1(i,1)+(j-1)*(H-y1(i,1))/(JL-2);
   end
end
%
% hold on
% for j=1:JL+1
%     plot(x1(1:IS+1,j),y1(1:IS+1,j),'bx-')
% end
%  
% for i=1:IS+1
%     plot(x1(i,1:JL+1),y1(i,1:JL+1),'bx-')
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate interior grid points using differential method (Laplace eq.) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x1,y1]=differential_method(x1,y1,IS-1,JL-1);
hold on
for j=1:JL-1
   plot(x1(1:IS-1,j),y1(1:IS-1,j),'bx-')
end
 for i=1:IS-1
   plot(x1(i,1:JL-1),y1(i,1:JL-1),'bx-')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              After Shock                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Left boundary
x2(1,1:JL-1)=0;
y2(1,1:JL-1)=0:dy:H;
%Top Wall
x2(2:IL-IS+1,JL-1)=dx:dx:L;
y2(2:IL-IS+1,JL-1)=H-(x2(2:IL-IS+1,JL-1))*tan(t1);
%Outlet
x2(IL-IS+1,1:JL-1)=L;
y2(IL-IS+1,1:JL-1)=0:dy2:H2;
%Lower wall
x2(1:IL-IS+1,1)=0:dx:L;
y2(1:IL-IS+1,1)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial guess to the interior points using algebric method %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=2:JL-2
   for i=2:IL-IS
       x2(i,j)=x2(1,j)+(i-1)*(L/(IL-IS));
   end
end
for i=2:IL-IS
   for j=2:JL-2
       y2(i,j)=y2(i,1)+(j-1)*(y2(i,JL-1)-y2(i,1))/(JL-2);
   end
end
hold on
for j=1:JL-1
   plot(x2(1:IL-IS,j),y2(1:IL-IS,j),'bx-')
end
 for i=1:IL-IS
   plot(x2(i,1:JL-1),y2(i,1:JL-1),'bx-')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate interior grid points using differential method (Laplace eq.)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x2,y2]=differential_method(x2,y2,IL-IS+1,JL-1);
x(1:IS-1,:)=x1;
y(1:IS-1,:)=y1;
x(IS-1:IL-1,:)=x2;
y(IS-1:IL-1,:)=y2;
%%%%%%%%%%%%%%%%%%
%  Ghost Cell    %
%%%%%%%%%%%%%%%%%%
%Inlet
x3(1,1:JL+1)=-(IS-1)*dx;
y3(1,1:JL+1)=-dy:dy:H+dy;
%Top Wall (Before shock)
x3(1:IS,JL+1)=-(IS-1)*dx:dx:0;
y3(1:IS,JL+1)=H+dy;
%Lower Wall (Before shock)
x3(1:IS,1)=-(IS-1)*dx:dx:0;
y3(1:IS,1)=-dy;
%Outlet
H3=(H-(L+dx)*tan(t1));
dy3=H3/(JL-2);
x3(IL+1,1:JL+1)=L+dx;
y3(IL+1,1:JL+1)=-dy3:dy3:H3+dy3;
%Lower wall and Top wall (After shock)
x3(IS:IL,1)=0:dx:L;
x3(IS:IL,JL+1)=0:dx:L;
for i=IS-1:IL-1
  dr(i,1)=y(i,JL-1)/(JL-2);
  y3(i+1,1)=-dr(i,1);
  y3(i+1,JL+1)=y(i,JL-1)+dr(i,1);
end
x4=x3;
y4=y3;
x4(2:IL,2:JL)=x;
y4(2:IL,2:JL)=y;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                 Plot Grid                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
for j=1:JL+1
   plot(x4(1:IL+1,j),y4(1:IL+1,j),'bx-')
end
 for i=1:IL+1
   plot(x4(i,1:JL+1),y4(i,1:JL+1),'bx-')
end
xlabel('$x$ (m)','Interpreter','LaTeX');
ylabel('$y$ (m)','Interpreter','LaTeX');
title('Numerical grid of Supersonic Engine Inlet','Interpreter','LaTeX');
subtitle(['for $IL = $ ' num2str(IL) ', $JL = $ ' num2str(JL) ', and $IS = $ ' num2str(IS)],'Interpreter','LaTeX');
set(gcf,'Position',[0 25 1200 600]);
set(gca,'LineWidth',2,'FontSize',18);
end

function [x,y] = differential_method(x,y,IL,JL)
%IL-number of grid points in x-direction
%JL-number of grid points in y-direction
norm1 = 1;
norm2 = 1;
epsilon = 1e-5; %Convergence criteria
kmax = 1e5;  %Maximum number of iteration
k = 0;
xold = x;
yold = y;
%Gauss Seidel to compute the interior coordinates of the computational
%domain
while (k < kmax) && ((norm1 > epsilon) || (norm2 > epsilon))
  
   k=k+1;   
  
   for j=2:JL-1
       for i=2:IL-1   
           alpha=((x(i,j+1)-x(i,j-1))/2)^2+((y(i,j+1)-y(i,j-1))/2)^2;
          
           beta =((x(i+1,j)-x(i-1,j))/2)*((x(i,j+1)-x(i,j-1))/2)...
               +((y(i+1,j)-y(i-1,j))/2)*((y(i,j+1)-y(i,j-1))/2);
          
           gamma=((x(i+1,j)-x(i-1,j))/2)^2+((y(i+1,j)-y(i-1,j))/2)^2;
           RHS1=-alpha*(x(i+1,j)+x(i-1,j))...
               +0.5*beta*(x(i+1,j+1)-x(i-1,j+1)-x(i+1,j-1)+x(i-1,j-1))...
               -gamma*(x(i,j+1)+x(i,j-1));
          
           RHS2=-alpha*(y(i+1,j)+y(i-1,j))...
               +0.5*beta*(y(i+1,j+1)-y(i-1,j+1)-y(i+1,j-1)+y(i-1,j-1))...
               -gamma*(y(i,j+1)+y(i,j-1));
          
           x(i,j)=-RHS1/(2*(alpha+gamma));
           y(i,j)=-RHS2/(2*(alpha+gamma));             
       end
   end
  
   norm1=max(max(abs(x-xold)));
   norm2=max(max(abs(y-yold)));
  
   yold=y;
   xold=x;
  
end
end

function [e1,e2,e3,e4]=Eflux(gamma,u1,u2,u3,u4)
%This function calculates the E flux.
% e1=U(:,2);
% e2=(U(:,2).^2./U(:,1))+(gamma-1).*(U(:,4)-(U(:,2).^2+U(:,3).^2)./(2.*U(:,1)));
% e3=U(:,2).*U(:,3)./U(:,1);
% e4=(gamma.*U(:,4)-(gamma-1).*(U(:,2).^2+U(:,3).^2)./(2.*U(:,1))).*(U(:,2)./U(:,1));
e1=u2;
e2=(u2.^2./u1)+(gamma-1).*(u4-(u2.^2+u3.^2)./(2.*u1));
e3=u2.*u3./u1;
e4=(gamma.*u4-(gamma-1).*(u2.^2+u3.^2)./(2.*u1)).*(u2./u1);
end

function [f1,f2,f3,f4]=Fflux(gamma,u1,u2,u3,u4)
%This function calculates the F flux.
% f1=U(:,3);
% f2=U(:,2).*U(:,3)./U(:,1);
% f3=(U(:,3).^2./U(:,1))+(gamma-1).*(U(:,4)-(U(:,2).^2+U(:,3).^2)./(2.*U(:,1)));
% f4=(gamma.*U(:,4)-(gamma-1).*(U(:,2).^2+U(:,3).^2)./(2.*U(:,1))).*(U(:,3)./U(:,1));
%
f1=u3;
f2=u2.*u3./u1;
f3=(u3.^2./u1)+(gamma-1).*(u4-(u2.^2+u3.^2)./(2.*u1));
f4=(gamma.*u4-(gamma-1).*(u2.^2+u3.^2)./(2.*u1)).*(u3./u1);
end

function [V,dx,dy,cx,cy]=volume(x1,x2,x3,x4,y1,y2,y3,y4)
%This function calculates the volume of a cell.
S1=0.5*abs((x3-x1)*(y4-y1)-(x4-x1)*(y3-y1));
S2=0.5*abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
V=S1+S2;
x=[x1,x2,x3,x4];
y=[y1,y2,y3,y4];
dx=max(x)-min(x);
dy=max(y)-min(y);
x(5)=x(1);
y(5)=y(1);
A=0;
for i=1:4
   A=A+0.5*(x(i)*y(i+1)-x(i+1)*y(i));
end
cx=0;
for i=1:4
   cx=cx+(1/(6*A))*(x(i)+x(i+1))*(x(i)*y(i+1)-x(i+1)*y(i));
end
cy=0;
for i=1:4
   cy=cy+(1/(6*A))*(y(i)+y(i+1))*(x(i)*y(i+1)-x(i+1)*y(i));
end
end

function [nx1,ny1,s1]=norm_comp_x(ypp,ypm,xpp,xpm)
dy1=ypp-ypm;
dx1=xpp-xpm;
s1=sqrt(dx1^2+dy1^2);
nx1=dy1/s1;
ny1=-dx1/s1;
end

function [nx2,ny2,s2]=norm_comp_y(ymp,ypp,xmp,xpp)
dy2=ymp-ypp;
dx2=xmp-xpp;
s2=sqrt(dx2^2+dy2^2);
nx2=dy2/s2;
ny2=-dx2/s2;
end

function u=u_vel(u1,u2)
%This function calculates the u velocity component.
u=u2./u1;
end

function v=v_vel(u1,u3)
%This function calculates the v velocity component.
v=u3./u1;
end




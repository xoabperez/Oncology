% perez_x_hw_2.m
% Created on: 2/17/18
% Created by: Xoab Perez

%% Problem 3a AND 3b
% Implement an avascular tumor growth model that considers a mechanically
% coupled RD model. This is done by reducing the diffusion coefficient
% as stresses increase. Stresses are calculated according to the gradient 
% of the number of cells.

clear; close all; 

% Load Initial Condition
load('hw_initial_condition','N');

% Initialize and assign constants
[sy,sx] = size(N);          % Number of iterations in each direction
dx = .1; dy = .1; dt = .01; % Assigned grid spacing
carcap = 100;               % Carrying capacity in cells
D0 = .1;                    % Diffusion coefficient in mm^2/day
k = 5*ones(sy,sx);          % Proliferation constant in /days
E = 2.5;                    % Young's modulus in N/mm
nu = .45;                   % Poisson value
lam1 = 1;                   % Parameter relating stress gradient to N
lam2 = [0 1e-3 1e-2 5e-2 1e-1 5e-1 1];  % Range of parameters for diffusion

% Build boundary condition matrix: M(div stress) = lambda1(grad N)
[M] = mech_matrix_build(sy,sx,dy,dx,E, nu); 
Um = zeros(sy,sx,2);                    % Divergence of stress tensor
tumor_area = zeros(100,length(lam2));   % Keep track os tumor areas
tumor_area(1,:) = sum(sum(N>25));       % Initial area

% Simulate tumor growth for different values of lambda2
for l = 1:length(lam2)
    load('hw_initial_condition','N');
    for t = 2:100
        %(2) Calculate Tumor Cell Gradient
        gradN = grad(N,dx,dy);

        %(3) Calculate Displacement
        U = M\(lam1*gradN);
        Um(:,:,1) = reshape(U(1:sx*sy),[sy,sx])';      % x-disp matrix
        Um(:,:,2) = reshape(U(sx*sy+1:end),[sy,sx])';  % y-disp matrix

        %(4) Calculate Stress & Strain
        [epsxx,epsyy,epsxy] = strains(Um,dx,dy);
        [sigmaxx,sigmayy,sigmaxy,sigmavm]=stresses(epsxx,epsyy,epsxy,E,nu);

        %(5) Update D
        D = D0*exp(-sigmavm*lam2(l));

        %(1) Calculate Tumor Cell Number
        N = rd_fdm_center_v1(N,D,k,carcap,[dx dy],dt);
        tumor_area(t,l) = sum(sum(N>25));
    end
    
    figure(3)
    plot(N(50,:)); hold on              % Plot cross-section
    legendinfo{l} = [num2str(lam2(l))];  
end
legend(legendinfo);
title('Cell Distribution at y = 5mm, Iteration 100');
title(legend, '\lambda_2')
xlabel('x position');
ylabel('Number of cells');

figure(1)
plot(tumor_area);
legend(legendinfo);
title('Tumor Area with Time');
title(legend, '\lambda_2')
xlabel('Iteration');
ylabel('# of cells with N>25');

figure(2)
plot(log(lam2),tumor_area(end,:));
title('Relationship between Tumor Area and \lambda_2');
xlabel('log(\lambda_2)');
ylabel('Tumor Area at Iteration 100');
 
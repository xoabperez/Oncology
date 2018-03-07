 % Homework 3
%% Problem 2a
clear; clc;close all
dx = 0.1; 
dy =0.1;
sx = 10/dx; 
sy = 10/dy;
dt = 0.02;
carcap = 1; % carrying capacity
D = zeros(sy,sx); % mm^2/day

N_o = zeros(sy,sx); % aerobic cells
N_g = zeros(sy,sx); % glycolysis cells
N_o(:,:) = 0.75*carcap; % initialize N_o
L = zeros(sy,sx);  % initialize lactate
G = 2*ones(sy,sx);  % initialize glucose
chi_L = 0*G;
chi_G = 0*G;
Dno = D; Dng = D;
Dg = 5*D; Dl = D; 

% Initialize model parameters
kg  = 0*5;  % N_g growth rate
ko  = 0*1;  % N_o growth rate
kog = 15; % conversion from N_o to N_g 
kgo = 3;  % conversion from N_g to N_o
gamma = 75; % sensitivity switch parameter
Gmin  = 0.2;
Lstar = .25; % Lactate threshold
Gstar = 0.5; 
Mstar = 0.1;
Nstar = 1;
alpha_G = 1;
alpha_L = 10;
beta_G =1.15;
beta_L = 3.2;
beta_o = 1.6;

%Plotting parameters
No = zeros(1,5/dt); Ng = No; Gtot = No; Ltot = No;
No(1) = sum(sum(N_o)); Ng(1) = sum(sum(N_g));
G_init = sum(sum(G));
Gtot(1) = G_init-sum(sum(G)); Ltot(1) = sum(sum(L));

for t = 2:(5/dt)
    chi_G(G>Gmin) = 1; chi_G(G<Gmin) = 0; chi_G(G==Gmin) = 0;
    chi_L = .5*(1+tanh(gamma*(L-Lstar)));
    No_new = rd_fdm(N_o,Dno,[dx,dy],dt)+dt*(ko.*N_o.*(1-(N_o+N_g)/carcap)+...
        kgo.*N_g.*chi_L-kog.*N_o.*(1-chi_L).*chi_G);
    Ng_new = rd_fdm(N_g,Dng,[dx,dy],dt)+dt*(kg.*N_g.*(1-(N_o+N_g)/carcap)-...
        kgo.*N_g.*chi_L+kog.*N_o.*(1-chi_L).*chi_G);
    G_new = rd_fdm(G,Dg,[dx,dy],dt)+dt*(-beta_o*alpha_G.*G.*N_o./(alpha_G*G+...
        alpha_L*L+Nstar)-beta_G.*G.*N_g./(G+Gstar));
    L_new = rd_fdm(L,Dl,[dx,dy],dt)+dt*(-beta_L*alpha_L.*L.*N_o./(alpha_G*G+...
        alpha_L*L+Mstar)+2*beta_G.*G.*N_g./(G+Gstar));    

    N_o = No_new; N_g = Ng_new; G = G_new; L  = L_new;

    No(t) = sum(sum(N_o)); Ng(t) = sum(sum(N_g));
    Gtot(t) = G_init-sum(sum(G)); Ltot(t) = sum(sum(L));
end
figure(1)
plot(1:5/dt, No, 1:5/dt, Ng);
xlabel('Iteration');
ylabel('Number of Cells');
title('Oxidative vs Glycolitic Cells');
legend('Oxidative','Glycolitic');
saveas(gcf,'2a1','jpg');

figure(2)
plot(1:5/dt, Gtot, 1:5/dt, Ltot);
xlabel('Iteration');
ylabel('Amount');
title('Glucose Consumption and Lactate Level');
legend('Glucose Consumed','Lactate');
saveas(gcf,'2a2','jpg');

%% Problem 2b
clear; clc;close all
dx = 0.1; 
dy =0.1;
sx = 10/dx; 
sy = 10/dy;
dt = 0.02;
carcap = 1; % carrying capacity
D = zeros(sy,sx); % mm^2/day

N_o = zeros(sy,sx); % aerobic cells
N_g = zeros(sy,sx); % glycolysis cells
N_o(:,:) = 0.75*carcap; % initialize N_o
L = zeros(sy,sx);  % initialize lactate
G = 2*ones(sy,sx);  % initialize glucose
chi_L = 0*G;
chi_G = 0*G;
Dno = D; Dng = D;
Dg = 5*D; Dl = D; 

% Initialize model parameters
kg  = 0*5;  % N_g growth rate
ko  = 0*1;  % N_o growth rate
kog = 15; % conversion from N_o to N_g 
kgo = 3;  % conversion from N_g to N_o
gamma = 75; % sensitivity switch parameter
Gmin  = 0.2;
Lstar = .25; % Lactate threshold
Gstar = 0.5; 
Mstar = 0.1;
Nstar = 1;
alpha_G = 1;
alpha_L = 10;
beta_G =1.15;
beta_L = 3.2;
beta_o = 1.6;

%Plotting parameters
No = zeros(1,5/dt+1); Ng = No; Gtot = No; Ltot = No;
No(1) = sum(sum(N_o)); Ng(1) = sum(sum(N_g));
G_init = sum(sum(G));
Gtot(1) = G_init-sum(sum(G)); Ltot(1) = sum(sum(L));

for t1 = 1:5
    for t = ((t1-1)/dt+2):(t1/dt+1)
        chi_G(G>Gmin) = 1; chi_G(G<Gmin) = 0; chi_G(G==Gmin) = 0;
        chi_L = .5*(1+tanh(gamma*(L-Lstar)));
        No_new = rd_fdm(N_o,Dno,[dx,dy],dt)+dt*(ko.*N_o.*(1-(N_o+N_g)/carcap)+...
            kgo.*N_g.*chi_L-kog.*N_o.*(1-chi_L).*chi_G);
        Ng_new = rd_fdm(N_g,Dng,[dx,dy],dt)+dt*(kg.*N_g.*(1-(N_o+N_g)/carcap)-...
            kgo.*N_g.*chi_L+kog.*N_o.*(1-chi_L).*chi_G);
        G_new = rd_fdm(G,Dg,[dx,dy],dt)+dt*(-beta_o*alpha_G.*G.*N_o./(alpha_G*G+...
            alpha_L*L+Nstar)-beta_G.*G.*N_g./(G+Gstar));
        L_new = rd_fdm(L,Dl,[dx,dy],dt)+dt*(-beta_L*alpha_L.*L.*N_o./(alpha_G*G+...
            alpha_L*L+Mstar)+2*beta_G.*G.*N_g./(G+Gstar));    

        N_o = No_new; N_g = Ng_new; G = G_new; L  = L_new;

        No(t) = sum(sum(N_o)); Ng(t) = sum(sum(N_g));
        Gtot(t) = G_init-sum(sum(G)); Ltot(t) = sum(sum(L));
    end
    L = zeros(sy,sx);  % reinitialize lactate
    G = 1*ones(sy,sx);  % reinitialize glucose
end
for t1 = 1:4
    Gtot(t1/dt+2:(t1+1)/dt+1) = Gtot(t1/dt+2:(t1+1)/dt+1)+Gtot(t1/dt+1)-...
        sum(sum(ones(sy,sx)));
end
figure(1)
plot(1:5/dt+1, No, 1:5/dt+1, Ng);
xlabel('Iteration');
ylabel('Number of Cells');
title('Oxidative vs Glycolitic Cells');
legend('Oxidative','Glycolitic');
xlim([0 250]); ylim([0 9000]);
saveas(gcf,'2b1','jpg');

figure(2)
plot(1:5/dt+1, Gtot, 1:5/dt+1, Ltot);
xlabel('Iteration');
ylabel('Amount');
title('Glucose Consumption and Lactate Level');
legend('Glucose Consumed','Lactate');
xlim([0 250]);% ylim([0 18000]);
saveas(gcf,'2b2','jpg');
    
%% Problem 2c
% clear; clc; close all
% dx = 0.1; 
% dy =0.1;
% sx = 10/dx; 
% sy = 10/dy;
% dt = 0.005;
% carcap = 1; % carrying capacity
% D = 0.005*ones(sy,sx); % mm^2/day
% 
% N_o = zeros(sy,sx); % aerobic cells
% N_g = zeros(sy,sx); % glycolysis cells
% N_o(40:60,40:60) = 0.75*carcap; % initialize N_o
% L = zeros(sy,sx);  % initialize lactate
% G = 2*ones(sy,sx);  % initialize glucose
% chi_L = 0*G;
% chi_G = 0*G;
% Dno = D; Dng = D;
% Dg = 5*D; Dl = D; 
% 
% % Initialize model parameters
% kg  = 5;  % N_g growth rate
% ko  = 1;  % N_o growth rate
% kog = 15; % conversion from N_o to N_g 
% kgo = 3;  % conversion from N_g to N_o
% gamma = 75; % sensitivity switch parameter
% Gmin  = 0.2;
% Lstar = .25; % Lactate threshold
% Gstar = 0.5; 
% Mstar = 0.1;
% Nstar = 1;
% alpha_G = 1;
% alpha_L = 10;
% beta_G =10;%1.15;
% beta_L = 3.2;
% beta_o = 1.6;
% 
% %Plotting parameters
% No = zeros(1,5/dt+1); Ng = No; Gtot = No; Ltot = No;
% No(1) = sum(sum(N_o)); Ng(1) = sum(sum(N_g));
% G_init = sum(sum(G));
% Gtot(1) = G_init-sum(sum(G)); Ltot(1) = sum(sum(L));
fig = figure('position',[100 100 1600 900]);
set(fig,'Visible','off');
writerObj = VideoWriter('newvile10.avi','Motion JPEG AVI');
writerObj.Quality = 100;
writerObj.FrameRate = 10;
open(writerObj); 

for t = 901:(5/dt)
    chi_G(G>Gmin) = 1; chi_G(G<Gmin) = 0; chi_G(G==Gmin) = 0;
    chi_L = .5*(1+tanh(gamma*(L-Lstar)));
    No_new = rd_fdm(N_o,Dno,[dx,dy],dt)+dt*(ko.*N_o.*(1-(N_o+N_g)/carcap)+...
        kgo.*N_g.*chi_L-kog.*N_o.*(1-chi_L).*chi_G);
    Ng_new = rd_fdm(N_g,Dng,[dx,dy],dt)+dt*(kg.*N_g.*(1-(N_o+N_g)/carcap)-...
        kgo.*N_g.*chi_L+kog.*N_o.*(1-chi_L).*chi_G);
    G_new = rd_fdm(G,Dg,[dx,dy],dt)+dt*(-beta_o*alpha_G.*G.*N_o./(alpha_G*G+...
        alpha_L*L+Nstar)-beta_G.*G.*N_g./(G+Gstar));
    L_new = rd_fdm(L,Dl,[dx,dy],dt)+dt*(-beta_L*alpha_L.*L.*N_o./(alpha_G*G+...
        alpha_L*L+Mstar)+2*beta_G.*G.*N_g./(G+Gstar));    

    N_o = No_new; N_g = Ng_new; G = G_new; L  = L_new;
    
    No(t) = sum(sum(N_o)); Ng(t) = sum(sum(N_g));
    Gtot(t) = G_init-sum(sum(G)); Ltot(t) = sum(sum(L));
    
%     subplot(2,3,1)
%     imagesc(N_o(:,:));
%     title('Concentration of Oxidative Cells');
%     xlabel('x position');
%     ylabel('y position');
%     colorbar
% 
%     subplot(2,3,2)
%     imagesc(N_g(:,:));
%     title('Concentration of Glycolitic Cells');
%     xlabel('x position');
%     ylabel('y position');
%     colorbar

    subplot(2,3,1)
    mesh(N_o(:,:));
    zlim([0 1]);
    title('Concentration of Oxidative Cells');
    xlabel('x position');
    ylabel('y position');
    zlabel('Normalized amount of cells');
    
    subplot(2,3,2)
    mesh(G(:,:));
    zlim([0 2]);
    title('Glucose Concentration');
    xlabel('x position');
    ylabel('y position');
    zlabel('Normalized amount');
    
    subplot(2,3,4)
    mesh(N_g(:,:));
    zlim([0 1]);
    title('Concentration of Glycolitic Cells');
    xlabel('x position');
    ylabel('y position');
    zlabel('Normalized amount of cells');
    
    subplot(2,3,5)
    mesh(L(:,:));
    zlim([0 1]);
    title('Lactate Concentration');
    xlabel('x position');
    ylabel('y position');
    zlabel('Normalized amount');

    subplot(2,3,3)
    plot(1:t, No(1:t), 1:t, Ng(1:t));
    xlabel('Iteration');
    ylabel('Number of Cells');
    title('Oxidative vs Glycolitic Cells');
    legend('Oxidative','Glycolitic');
    xlim([0 1000]); ylim([0 5000]);

    subplot(2,3,6)
    plot(1:t, Gtot(1:t), 1:t, Ltot(1:t));
    xlabel('Iteration');
    ylabel('Amount');
    title('Glucose Consumption and Lactate Level');
    legend('Glucose Consumed','Lactate');
    xlim([0 1000]); ylim([0 5000]);
    
    if mod(t,5)==0
        F = getframe(fig);
        writeVideo(writerObj,F);
    end
end
close(fig)
close(writerObj)

%% Problem 3a (1D optimization)
clear;
load N_1d
ktrue = k;
Dtrue = D0;


        dx = 0.1; dt = 0.02;
    sx = 100;
    
% Noiseless
N_s_inf = N_s;

% SNR 16
N_s_snr_16 = N_s.*random('Normal',1,1/16,[100 10]);

% SNR 4
N_s_snr_4 = N_s.*random('Normal',1,1/4,[100 10]);

% (1) Optimize Parameters

% (2) Predict future tumor growth

% (3) Calculate parameter error


%% Problem 3b (Parameter Optimization Noiseless)
clear all; clc
load N_2d
warning off;
ktrue = k;
Dtrue = D0;
dx = 0.1; % in mm
dy = 0.1; % in mm
sx = 10/dx;  %no. grid points
sy = 10/dy; %no. grid points
dt = 0.01; %in days
    
LB = zeros(26,1);  % Lower Bounds
UB = 10*ones(26,1); % Upper Bounds
UB(end) = (1/4)/(1/dx^2 + 1/dy^2)/dt;
params = ones(26,1); params(end) = 0.5*UB(end); % Initial Guess...

 options = optimset('TolFun',1e-12,'Tolx',1e-12,'MaxIter',1000,'Display','off');



    
 
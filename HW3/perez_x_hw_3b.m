 % Homework 3
%% Problem 3a (1D optimization)
clear all; clc; close all;
load N_1d
ktrue = k;
Dtrue = .05;
dx = 0.1; dt = 0.02;

% Noiselsss and signal-to-noise ratios of 16, 4
N = {N_s; N_s.*random('Normal',1,1/16,[100 10]);...
    N_s.*random('Normal',1,1/4,[100 10])};

for x = 1:3
N_true = cell2mat(N(x));

% (1) Optimize Parameters
kg = ones(100,1);   % Initial guess for k
Dg = .0625;         % Initial guess for D
tol = 1e-5;         % Error tolerance
maxiter = 1000;     % Max # of iterations
lambda = 2;         
carcap = 1;         % Carrying Capacity
kbounds = [0 10];   % Bounds on k
Dbounds = [0 .125]; % Bounds on D
gamma = 1e-6;       % Difference for derivative

N_2 = rd_fdm_1d_v1(N_true(:,1),Dg,kg,carcap,dx,dt,100); % N_2 estimate
N_3 = rd_fdm_1d_v1(N_2(:,1),Dg,kg,carcap,dx,dt,100);    % N_3 estimate

err = sum((N_2-N_true(:,2)).^2)+sum((N_3-N_true(:,3)).^2); % Initial error
iter = 0;
erriter = 0;

while err > tol && iter < maxiter
    iter = iter+1;
    % Calculate Jacobian
    for n = 1:100
        kt = kg;
        kt(n) = kt(n)+gamma;            % Change slightly
        Nt2 = rd_fdm_1d_v1(N_true(:,1),Dg,kt,carcap,dx,dt,100);
        Nt3 = rd_fdm_1d_v1(Nt2,Dg,kt,carcap,dx,dt,100);
        J(1:100,n) = (Nt2-N_2)/gamma;   % "Derivative" of N_2 estimate
        J(101:200,n) = (Nt3-N_3)/gamma; % "Derivative" of N_3 estimate       
    end
    Dt = Dg+gamma;                     
    Nt2 = rd_fdm_1d_v1(N_true(:,1),Dt,kg,carcap,dx,dt,100);
    Nt3 = rd_fdm_1d_v1(Nt2,Dt,kg,carcap,dx,dt,100);
    J(1:100,n+1) = (Nt2-N_2)/gamma;
    J(101:200,n+1) = (Nt3-N_3)/gamma;
    
    % Calculate change in parameters
    err_v = [N_true(:,2)-N_2; N_true(:,3)-N_3];       
    del = (J'*J+lambda*diag(diag(J'*J)))\(J'*err_v); % Increment
    
    % Update parameters
    kt = kg+del(1:end-1);
    kt(kt<kbounds(1)) = kbounds(1); % Enforce bounds
    kt(kt>kbounds(2)) = kbounds(2);
    Dt = Dt+del(n+1);
    if Dt < Dbounds(1)
        Dt = Dbounds(1);
    elseif Dt > Dbounds(2)
        Dt = Dbounds(2);
    end
    
    % Evaluate the model
    Nt2 = rd_fdm_1d_v1(N_true(:,1),Dt,kt,carcap,dx,dt,100); % New guess
    Nt3 = rd_fdm_1d_v1(Nt2,Dt,kt,carcap,dx,dt,100);
    errt = sum((Nt2-N_true(:,2)).^2)+sum((Nt3-N_true(:,3)).^2); % New err
    
    if errt < err   % Keep new parameters if error decreases
        erriter(iter) = errt;
        err = errt;
        kg = kt;
        Dg = Dt;
        lambda = lambda/2;
        N_2 = Nt2;
        N_3 = Nt3;
    else            % Otherwise, change lambda
        erriter(iter) = errt;
        lambda = lambda*4;
    end
    
end
figure(1)
plot(1:1:100,kt,'k'); hold on
plot(1:1:100,k,'r');
title('Proliferation Rate Comparison')
xlabel('x-position')
legend('Estimated k','Actual k');
saveas(gcf,strcat('3b_k_',num2str(x)),'png');
    
% (2) Predict future tumor growth
N_guess(:,1) = N_true(:,1);
for n = 2:10
    N_guess(:,n) = rd_fdm_1d_v1(N_true(:,1),Dg,kg,carcap,dx,dt,(n-1)*100);
end

figure(2)
set(gcf,'pos',[680 554 915 428])
subplot(1,2,1)
plot(N_true)
title('True data')
subplot(1,2,2)
plot(N_guess)
title('Estimated data')
saveas(gcf,strcat('3b_timepoints_N_x_',num2str(x)),'png');

figure(3)
z = 1:1:100;
plot(z,N_s(:,3),'ko',z,N_s(:,6),'bo'); hold on
plot(z,N_true(:,3),'rx',z,N_true(:,6),'mx');
title('True vs. Estimated Data')
xlabel('x-position')
ylabel('Number of Cells')
legend('True, 72 hrs','True, 144 hrs','Est., 72 hrs','Est., 144 hrs',...
    'Location','South')
saveas(gcf,strcat('3b_timepoints_36_N_x_',num2str(x)),'png');


% (3) Calculate parameter error
k_error(x) = sum((ktrue-kg).^2);
D_error(x) = (Dtrue-Dg)^2;

end
disp(k_error)
disp(D_error)

%% Problem 3b (Parameter Optimization Noiseless)
clear all; clc
load N_2d
warning off;
ktrue = k;
Dtrue = .05;
dx = 0.1; % in mm
dy = 0.1; % in mm
sx = 10/dx;  %no. grid points
sy = 10/dy; %no. grid points
dt = 0.01; %in days
carcap = 100;
    
LB = zeros(26,1);  % Lower Bounds
UB = 2*ones(26,1); % Upper Bounds
UB(end) = .125;
params = ones(26,1); params(end) = 0.5*UB(end); % Initial Guess...

options = optimset('TolFun',1e-12,'Tolx',1e-12,...
    'MaxIter',1000,'Display','off');

N = {N_snr_Inf; N_snr_16; N_snr_4};

for x = 1:3

% (1) Optimize Parameters using lsqnonlin
Nx = cell2mat(N(x));
fun = @(p) [reshape(rd_fdm_center_steps_v1(Nx(:,:,1),p(end),...
    prolif_func(p(1:25)),carcap,[dx dy],dt,100)-Nx(:,:,2),10000,1);...
    reshape(rd_fdm_center_steps_v1(Nx(:,:,1),p(end),...
    prolif_func(p(1:25)),carcap,[dx dy],dt,200)-Nx(:,:,3),10000,1)];

p_est = lsqnonlin(fun,params,LB,UB,options);

% (2) Predict future tumor growth
sserr(1) = 0;
for n = 2:10
    N_guess = rd_fdm_center_steps_v1(Nx(:,:,1),p_est(end),...
    prolif_func(p_est(1:25)),carcap,[dx dy],dt,(n-1)*100);
    sserr(n) = sum((sum(N_guess-N_snr_Inf(:,:,n)).^2));
end
figure(1)
semilogy(sserr(1:10))
title('Forward solve error')
xlabel('Time step')
saveas(gcf,strcat('4_l2_err_',num2str(x)),'png');

% (3) Calculate parameter error
k_error(x) = sum(sum((ktrue-prolif_func(p_est(1:25))).^2));
D_error(x) = (Dtrue-p_est(end))^2;

figure(2)
set(gcf,'pos',[680 554 915 428])
subplot(1,2,1)
mesh(ktrue)
title('True k')
subplot(1,2,2)
mesh(prolif_func(p_est(1:25)));
title('Estimated k')
saveas(gcf,strcat('4_k_',num2str(x)),'png');

figure(3)
set(gcf,'pos',[680 554 915 428])
subplot(1,2,1)
mesh(N_snr_Inf(:,:,10))
title('True N')
subplot(1,2,2)
mesh(N_guess);
title('Estimated N')
saveas(gcf,strcat('4_N_',num2str(x)),'png');

end
    
 
 % Homework 3
%% Problem 2a
clear all; clc;
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
    
    for t = 2:(5/dt)

        % Put code here!
        
        % or Here!
        
        
    end
%% Problem 2b
clear all; clc;
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
    count = 0;
    for t = 2:(5/dt)
    
        % you know what to do
        
    end
%% Problem 2c
clear all; clc;
    dx = 0.1; 
    dy =0.1;
    sx = 10/dx; 
    sy = 10/dy;
    dt = 0.005;
    carcap = 1; % carrying capacity
    D = 0.005*ones(sy,sx); % mm^2/day
    
    N_o = zeros(sy,sx); % aerobic cells
    N_g = zeros(sy,sx); % glycolysis cells
    N_o(40:60,40:60) = 0.75*carcap; % initialize N_o
    L = zeros(sy,sx);  % initialize lactate
    G = 2*ones(sy,sx);  % initialize glucose

    Dno = D; Dng = D;
    Dg = 5*D; Dl = D; 
    
    % Initialize model parameters
    kg  = 5;  % N_g growth rate
    ko  = 1;  % N_o growth rate
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
    beta_G =10;%1.15;
    beta_L = 3.2;
    beta_o = 1.6;
    
    for t = 2:(5/dt)
     
         % empty space for creativity
         
    end

%% Problem 3a (1D optimization)
clear all; clc;
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
    


    
 
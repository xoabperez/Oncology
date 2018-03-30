% perez_x_hw_1.m
% Created on: Feb. 4, 2018
% Created by: Xoab Perez

%% Problem 3a
% k = 0 (no proliferation)
% Initialize your parameters and simulation domain
clc; clear; close all;
X =10; Y =10; dx = .1; dy = .1; % Total and step lengths
dt = 0.01;                      % time step
D0 = .05;                       % diffusion coefficient, mm^2/day
k = 0;                          % proliferation rate /day
carcap = 100;                   % carrying capacity, theta
days = 5;                       % Total time
N=zeros(X/dx+1,Y/dy+1,days/dt); % Initialize tumor cell distribution
N(45:55,45:55,1) = 0.75*carcap; % Initial condition
cells = zeros(days/dt,1);       % Number of cells each day
cells(1) = sum(sum(N(:,:,1)));  % First day

for t = 2:(days/dt) 
    N(:,:,t) = forward(N(:,:,t-1),dx,dy,dt,D0,k,carcap);
    cells(t) = sum(sum(N(:,:,t))); % Total number of cells at time step
end

%Display change in total number of cells with time
plot(cells-cells(1),'Linewidth',2);
xlabel('Iteration')
ylabel('(Total Cell Number) - (Initial)')
title('No Proliferation Coeff. (k = 0)')
set(gca,'LineWidth',1.5,'FontSize',10)

%% Problem 3b
% k = 0, various time step sizes
% Initialize your parameters and simulation domain
clc; clear; close all;
X = 10; Y = 10; dx = .1; dy = .1;   % grid spacing
dt = [0.01 .03 .05 .07 .09 .11];    % time steps
D0 = .05;                           % diffusion coefficient, mm^2/day
k = 0;                              % proliferation rate /day
carcap = 100;                       % carrying capacity, theta
N = zeros(X/dx+1,Y/dy+1,400);       % Initialize tumor cell distribution
N(45:55,45:55,1) = 0.75*carcap;     % Initial condition
cells = zeros(400,length(dt));      % Number of cells each day
cells(1,:) = sum(sum(N(:,:,1)));    % First day

% Forward Evaulation
for n = 1:length(dt)
    for t = 2:400 
        N(:,:,t) = forward(N(:,:,t-1),dx,dy,dt(n),D0,k,carcap);
        cells(t,n) = sum(sum(N(:,:,t))); % Total number of cells at time step
    end
end

%Display change in total number of cells with time
figure(1)
for n = 1:length(dt)
    subplot(3,2,n)
    plot(cells(:,n))
    xlabel('Iteration')
    ylabel('Number of Tumor Cells')
    title(strcat('k = 0, dt = ',num2str(dt(n))))
    set(gca,'LineWidth',1.5,'FontSize',10)
    if (cells(end,n)-cells(1,n)<100)
        ytickformat('%4.2f')
    else
        ytickformat('auto')
    end
end

%% Problem 3c
% k = .003, varying diffusion coefficients
% Initialize your parameters and simulation domain
clc; clear; close all;
X = 10; Y = 10; dx = .1; dy = .1;   % grid spacing
dt = 0.01;                          % time step
k = .003;                           % proliferation rate /day
D0 = k*[.5 1 2 4 8 16 32];          % diffusion coefficient, mm^2/day
carcap = 100;                       % carrying capacity, theta
N = zeros(X/dx+1,Y/dy+1,400);   % Initialize tumor cell distribution
N(45:55,45:55,1) = 0.75*carcap;     % Initial condition
profile = zeros(X/dx+1,length(D0)); % Matrix for profile data
legendInfo = cell(length(D0),1);

% Forward Evaulation
for n = 1:length(D0)
    for t = 2:400 
        N(:,:,t) = forward(N(:,:,t-1),dx,dy,dt,D0(n),k,carcap);
    end
    profile(:,n) = N(:,50,400);     % Profile data at y=5mm
    legendInfo{n} = num2str(D0(n)/k);
end
figure(1)
plot(profile);
legend(legendInfo)
title(legend, 'D_0/k (mm^2)')
xlim([0 100])
xlabel('x-coordinate')
ylabel('Number of Tumor Cells')
title({'Varying Diffusion Coefficients, k = .003 day^{-1}',...
    '(Cross-section at y = 5mm, t = 4 days)'})
set(gca,'LineWidth',1.5,'FontSize',10)

%% Problem 3d
% k = 2.5, original diffusion coeff. examine growth with time
% Initialize your parameters and simulation domain
clc; clear; %close all;
X = 10; Y = 10; dx = .1; dy = .1; % grid spacing
dt = 0.01;                      % time step
D0 = .05;                       % diffusion coefficient, mm^2/day
k = 2.5;                        % proliferation rate /day
carcap = 100;                   % carrying capacity, theta
days = 5;                       % Total time
N=zeros(X/dx+1,Y/dy+1,days/dt); % Initialize tumor cell distribution
N(45:55,45:55,1) = 0.75*carcap; % Initial condition

% Forward Evaulation
for t = 2:(days/dt)
    N(:,:,t) = forward(N(:,:,t-1),dx,dy,dt,D0,k,carcap);
end

figure(1)
for n = 1:5
    plot(N(:,50,n*100)); hold on 
    legendInfo{n} = ['Day ' num2str(n)];
end
legend(legendInfo)
xlim([0 100])
set(gca,'LineWidth',1.5,'FontSize',10);
xlabel('x-coordinate')
ylabel('Number of Tumor Cells')
title({'Tumor Growth','(Cross-section at y = 5mm)'})
%sum(sum(N(:,:,t)));

%% 1D Reaction Diffusion Equation
%{
Garg paper
D = 0.001 cm^2/day
E = 2100 Pa;
poissons ratio = 0.45;

Predicted Distribution at days 200, 300, 400
sx = 0 to 10 cm
dx = 0.01;


%}
%
clear all; clc;
dx = 0.1; % grid spacing
dt = 0.02;  % time step
sx = 100;   % no. of grid points
D0 = 0.1;      % diffusion coefficient
k = 5;   % proliferation rate
carcap = 1; % carrying capacity

% Explicit time differentiation
    N = zeros(100,1);
    D = D0*ones(100,1); 
    
    %initialize tumor cell distribution
    N(45:55,1) = 1*carcap;
    Nm = N;

    % Build mechanical model coefficient matrix
        E = 1*ones(100,1); %youngs modulus
        E(1:49) = 1; %;
        lam_1 = 1; % coupling constant 1
        lam_2 = 5; % coupling constant 2
        M = zeros(100,100); % coefficient matrix
        U = zeros(100,1);    % Displacement vector
        gradN = zeros(100,1); % gradN
        for x = 1:sx
           if x == 1
               M(1,1) = 1;
               gradN(1) = 0;
           elseif x == sx
               M(sx,sx) = 1;
               gradN(sx) = 0;
           else
               if x == 2
                   M(x,x-1) = 0;
               else
                   M(x,x-1) = E(x)*(1/dx^2);
               end
               M(x,x) = -2*E(x)/dx^2;
               if x == sx-1
                   M(x,x+1) = 0;
               else
                   M(x,x+1) = E(x)/dx^2;
               end
           end
           
        end
        
    for t = 2:(30/dt)
        % (1) 1D reaction-diffusion equation
            % Non-mechanically coupled RD 
            for x = 1:sx
                if x ==1 
                    d2Ndx2 = (2*N(x+1,t-1)-2*N(x,t-1))/(dx^2); %boundary condition
                    d1Ndx1 = 0;
                    d1Ddx1 = 0;
                elseif x == sx
                    d2Ndx2 = (2*N(x-1,t-1)-2*N(x,t-1))/(dx^2); % boundary condition
                    d1Ndx1 = 0;
                    d1Ddx1 = 0;
                else
                    d2Ndx2 = (N(x+1,t-1)-2*N(x,t-1)+N(x-1,t-1))/(dx^2);
                    d1Ndx1 =  (N(x+1,t-1)-N(x-1,t-1))/(2*dx);
                    d1Ddx1 = 0; %(D(x+1,t-1)-D(x-1,t-1))/(2*dx);
                end
                N(x,t) = N(x,t-1) + dt*(D0*d2Ndx2 +k*N(x,t-1)*(1-N(x,t-1)/carcap));
            end
            
            % mechanically coupled RD model
            for x = 1:sx
            if x ==1 
                d2Ndx2 =(2*Nm(x+1,t-1)-2*Nm(x,t-1))/(dx^2);
                d1Ndx1 = 0;
                d1Ddx1 = 0;
            elseif x == sx
                d2Ndx2 = (2*Nm(x-1,t-1)-2*Nm(x,t-1))/(dx^2);
                d1Ndx1 =  0;
                d1Ddx1 =0;
            else
                d2Ndx2 = (Nm(x+1,t-1)-2*Nm(x,t-1)+Nm(x-1,t-1))/(dx^2);
                d1Ndx1 =  (Nm(x+1,t-1)-Nm(x-1,t-1))/(2*dx);
                d1Ddx1 = (D(x+1)-D(x-1))/(2*dx);
            end
                Nm(x,t) = Nm(x,t-1)+dt*(D(x)*d2Ndx2 + d1Ddx1*d1Ndx1+k*Nm(x,t-1)*(1-Nm(x,t-1)/carcap));
            end
        
        % (2) Calculate Displacement
            gradN(1) = 0;
            gradN(100) = 0;
            gradN(2:99) = (Nm(3:100,t)-Nm(1:98,t))/(2*dx);
        
            U = M\(lam_1*gradN);
        
        % (3) Calculate Strain
           strain(1) = (U(2)-U(1))/dx;
           strain(100) = (U(100)-U(99))/dx;
           strain(2:99) = (U(3:100)-U(1:98))/(2*dx);
           
        % (4) Calculate Stress
            stress = E(x)*strain;
        
        % (5) Update D(x,t)
            D = D0*exp(-abs(stress)*lam_2);
            
        subplot(2,2,1)
        plot(1:100,Nm(:,t),'g',1:100,N(:,t),'b',[20 20],[0 200],'--k',[40 40],[0 200],'--k',[50 50],[0 200],'r'...
            ,[60 60],[0 200],'--k',[80 80],[0 200],'--k','LineWidth',1.5); axis square
        axis([0 100 0 1*1.05])
        set(gca,'LineWidth',1.5,'FontSize',20);
        title('Cells','FontSize',20)
        
        subplot(2,2,2)
        plot(1:100,U,'LineWidth',1.5)
        axis([0 100 -10 20])
        set(gca,'LineWidth',1.5,'FontSize',20);
        title('Displacement','FontSize',20)
        
        subplot(2,2,3)
        plot(1:100,strain,'LineWidth',1.5)
        axis([0 100 -1 1])
        set(gca,'LineWidth',1.5,'FontSize',20);
        title('Strain','FontSize',20)
       
        subplot(2,2,4)
        plot(1:100,stress,'LineWidth',1.5)
        axis([0 100 -1 1])
        set(gca,'LineWidth',1.5,'FontSize',20);
        title('stress','FontSize',20)
      
        drawnow
    end

%% 1D Reaction Diffusion Equation
%{
Garg paper
D = 0.001 cm^2/day
E = 2100 Pa;
poissons ratio = 0.45;

Predicted Distribution at days 200, 300, 400
sx = 0 to 10 cm
dx = 0.01;


%}
%
clear all; clc;
dx = 0.1; % grid spacing
dt = 0.02;  % time step
sx = 100;   % no. of grid points
D0 = 0.1;      % diffusion coefficient
k = 5;   % proliferation rate
carcap = 1; % carrying capacity

% Explicit time differentiation
    N = zeros(100,1);
    D = D0*ones(100,1); 
    
    %initialize tumor cell distribution
    N(45:55,1) = 1*carcap;
    Nm = N;

    % Build mechanical model coefficient matrix
        E = 1*ones(100,1); %youngs modulus
        E(1:49) =1; %;
        lam_1 = 1; % coupling constant 1
        lam_2 = 5; % coupling constant 2
        M = zeros(100,100); % coefficient matrix
        U = zeros(100,1);    % Displacement vector
        gradN = zeros(100,1); % gradN
        for x = 1:sx
           if x == 1
               M(1,1) = 1;
               gradN(1) = 0;
           elseif x == sx
               M(sx,sx) = 1;
               gradN(sx) = 0;
           else
               if x == 2
                   M(x,x-1) = 0;
               else
                   M(x,x-1) = E(x)*(1/dx^2);
               end
               M(x,x) = -2*E(x)/dx^2;
               if x == sx-1
                   M(x,x+1) = 0;
               else
                   M(x,x+1) = E(x)/dx^2;
               end
           end
           
        end
        
    for t = 2:(30/dt)
        % (1) 1D reaction-diffusion equation
            for x = 1:sx
                if x ==1 
                d2Ndx2 = (2*Nm(x+1,t-1)-2*Nm(x,t-1))/dx^2;  %boundary condition
                d1Ndx1 = 0; d1Ddx1 = 0;
            elseif x == sx
                d2Ndx2 = (2*Nm(x-1,t-1)-2*Nm(x,t-1))/dx^2; % boundary condition
                d1Ndx1 = 0; d1Ddx1 = 0;
            else
                d2Ndx2 = (Nm(x+1,t-1)-2*Nm(x,t-1)+Nm(x-1,t-1))/dx^2;
                d1Ndx1 = (Nm(x+1,t-1)-Nm(x-1,t-1))/(2*dx); 
                d1Ddx1 = (D(x+1)-D(x-1))/(2*dx); 
                end
                Nm(x,t) = Nm(x,t-1) + dt*(D(x)*d2Ndx2+d1Ndx1*d1Ddx1+k*Nm(x,t-1)*(1-Nm(x,t-1)/carcap));
            end
            
                        for x = 1:sx
            if x ==1 
                d2Ndx2 = (2*N(x+1,t-1)-2*N(x,t-1))/dx^2;  %boundary condition
                d1Ndx1 = 0; d1Ddx1 = 0;
            elseif x == sx
                d2Ndx2 = (2*N(x-1,t-1)-2*N(x,t-1))/dx^2; % boundary condition
                d1Ndx1 = 0; d1Ddx1 = 0;
            else
                d2Ndx2 = (N(x+1,t-1)-2*N(x,t-1)+N(x-1,t-1))/dx^2;
            end
                N(x,t) = N(x,t-1) + dt*(D0*d2Ndx2+k*N(x,t-1)*(1-N(x,t-1)/carcap));
            end
        
        % (2) Calculate Displacement
            gradN(1) = 0; gradN(100) = 0;
            gradN(2:99) = (Nm(3:100,t)-Nm(1:98,t))/(2*dx);
        
            U = M\(lam_1*gradN);
        
        % (3) Calculate Strain
           strain(1) = (U(2)-U(1))/dx; strain(100) = (U(100)-U(99))/dx;
           strain(2:99) = (U(3:100)-U(1:98))/(2*dx);
           
        % (4) Calculate Stress
            stress = E(x)*strain;
        
        % (5) Update D(x,t)
            D = D0*exp(-abs(stress)*lam_2);
            %% 
        subplot(2,2,1)
        plot(1:100,Nm(:,t),'g',1:100,N(:,t),'b',[20 20],[0 200],'--k',[40 40],[0 200],'--k',[50 50],[0 200],'r'...
            ,[60 60],[0 200],'--k',[80 80],[0 200],'--k','LineWidth',1.5); axis square
        axis([0 100 0 1*1.05])
             set(gca,'LineWidth',1.5,'FontSize',20);
        title('Cells: MC-RD = green, RD = blue','FontSize',20)
        
        subplot(2,2,2)
        plot(1:100,U,'LineWidth',1.5)
        axis([0 100 -10 20])
        set(gca,'LineWidth',1.5,'FontSize',20);
        title('Displacement','FontSize',20)
        
        subplot(2,2,3)
        plot(1:100,strain,'LineWidth',1.5)
        axis([0 100 -1 1])
        set(gca,'LineWidth',1.5,'FontSize',20);
        title('Strain','FontSize',20)
       
        subplot(2,2,4)
        plot(1:100,stress,'LineWidth',1.5)
        axis([0 100 -1 1])
        set(gca,'LineWidth',1.5,'FontSize',20);
        title('stress','FontSize',20)
      
        drawnow
    end



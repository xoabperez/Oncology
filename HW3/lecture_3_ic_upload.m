%lecture_2_ic.m 
% in-class 1D example
clear all; close all; clc;

x = 0:0.1:1;
exampleF = @(x,p) (x-p(1)).*p(2)+(x).^2*p(3);
ptrue = [10 2 3];
N_true = exampleF(x,ptrue)+random('Normal',0,.1,[1 11]);

pg = ones(3,1); 
lambda = 10;
tol = 1e-5;
p_bounds = [0 20];

N_g = exampleF(x,pg);
N_best = N_g;
err = sum((N_g-N_true).^2);
pe = pg;  % best current values of P

it_count = 0;
err_c = 0;
while err > tol & it_count < 2000
    it_count = it_count + 1;
    
     
     % Calculate Jacobian
     
       for n = 1:3
           pt = pe;
           pt(n) = pt(n)+1e-8;
           Nt = exampleF(x,pt);
           J(:,n) = (Nt-N_best)/(1e-8);
       end
     
     % Calculate change in parameters
        err_v = N_true-N_best;
        del=(J'*J+lambda*diag(diag(J'*J)))\((J'*err_v'));
        
    % Update parameters
        pt = pe;
        for n = 1:3
            pt(n) = pe(n) + del(n);
            if pt(n) < p_bounds(1)
                pt(n) = p_bounds(1);
            elseif pt(n) > p_bounds(2);
                pt(n) = p_bounds(2);
            end
        end
        
    % Evaluate the model
        Nt = exampleF(x,pt);
        err_t = sum((Nt-N_true).^2);
    
        if err_t < err
            err_c(it_count) = err_t;
            err = err_t;
            pe = pt;
            lambda = lambda/2;
            N_best = Nt;
            disp(pe)
        else
            err_c(it_count) = err_t;
            lambda = lambda*1.5;
        end
            
            
    
        subplot(1,2,1)
    plot(x,N_true,'.r',x,N_best,'r','LineWidth',1.5,'MarkerSize',20)
    xlabel('Position')
    ylabel('Measure')
    subplot(1,2,2)
    plot(log10(err_c))
    drawnow
    pause(0.2)
    
end


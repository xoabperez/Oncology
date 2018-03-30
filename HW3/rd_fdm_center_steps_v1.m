function [Nout] = rd_fdm_center_steps_v1(N_in,D,k,carcap,dims,dt,steps)
    [sy, sx] = size(N_in);
    Nout = zeros(sy,sx);
    dx = dims(1);
    dy = dims(2);
    % steps = the number of simulation time steps
    
    
    for tt = 1:steps
    
    for x = 1:sx
        for y = 1:sy
            if x == 1
                diff_x = (-2*N_in(y,x)+2*N_in(y,x+1))/(dx^2);
            elseif x == sx
                diff_x = (-2*N_in(y,x)+2*N_in(y,x-1))/(dx^2);
            else
                diff_x = (N_in(y,x+1)-2*N_in(y,x)+N_in(y,x-1))/(dx^2);
            end
            if y == 1
                diff_y = (2*N_in(y+1,x)-2*N_in(y,x))/(dy^2);
            elseif y == sy
                diff_y = (2*N_in(y-1,x)-2*N_in(y,x))/(dy^2);
            else
                diff_y = (N_in(y+1,x)-2*N_in(y,x)+N_in(y-1,x))/(dy^2);
            end
            
            diffusion = D*(diff_y+diff_x);
            Nout(y,x) = N_in(y,x) +dt*(diffusion+k(y,x)*N_in(y,x)*(1-N_in(y,x)/carcap)) ;  
        end
    end


    N_in = Nout;
    end




end
function [Nout] = rd_fdm_center_v1(N_in,D,k,carcap,dims,dt,steps)
    [sy sx] = size(N_in);
    Nout = zeros(sy,sx);
    dx = dims(1);
    dy = dims(2);
    % steps = the number of simulation time steps
    
    
    for tt = 1:steps
    
    for x = 1:sx
        for y = 1:sy
            if x == 1
                diff_x = (-2*N_in(y,x)+2*N_in(y,x+1))/(dx^2);
                diff_1x = 0; diff_Dx = 0;
            elseif x == sx
                diff_x = (-2*N_in(y,x)+2*N_in(y,x-1))/(dx^2);
                diff_1x = 0; diff_Dx = 0;
            else
                diff_x = (N_in(y,x+1)-2*N_in(y,x)+N_in(y,x-1))/(dx^2);
                diff_1x = (1/(2*dx))*(N_in(y,x+1)-N_in(y,x-1));
                diff_Dx = (1/(2*dx))*(D(y,x+1)-D(y,x-1));
            end
            if y == 1
                diff_y = (2*N_in(y+1,x)-2*N_in(y,x))/(dy^2);
                diff_1y = 0; diff_Dy = 0;
            elseif y == sy
                diff_y = (2*N_in(y-1,x)-2*N_in(y,x))/(dy^2);
                diff_1y = 0; diff_Dy = 0;
            else
                diff_y = (N_in(y+1,x)-2*N_in(y,x)+N_in(y-1,x))/(dy^2);
                diff_1y = (1/(2*dy))*(N_in(y+1,x)-N_in(y-1,x));
                diff_Dy = (1/(2*dy))*(D(y+1,x)-D(y-1,x));
            end
            
            
            diffusion = D(y,x)*(diff_y+diff_x)+diff_Dx*diff_1x+diff_Dy*diff_1y;
            Nout(y,x) = N_in(y,x) +dt*(diffusion+k(y,x)*N_in(y,x)*(1-N_in(y,x)/carcap)) ;  
        end
    end


    N_in = Nout;
    end




end
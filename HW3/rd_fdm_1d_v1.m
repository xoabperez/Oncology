function [N_o] = rd_fdm_1d_v1(N,D,k,carcap,dx,dt,steps)
[sx, ~] = size(N);
N_o = N;
for in_count = 1:steps
    for x = 1:sx
                if x ==1 
                    d2Ndx2 = (2*N(x+1)-2*N(x))/dx^2;  %boundary condition
                    d1Ndx1 = 0; d1Ddx1 = 0;
                elseif x == sx
                    d2Ndx2 = (2*N(x-1)-2*N(x))/dx^2; % boundary condition
                    d1Ndx1 = 0; d1Ddx1 = 0;
                else
                    d2Ndx2 = (N(x+1)-2*N(x)+N(x-1))/dx^2;
                    d1Ndx1 = (N(x+1)-N(x-1))/(2*dx); 
                    d1Ddx1 = 0; %(D(x+1)-D(x-1))/(2*dx); 
                end
                    N_o(x) = N(x) + dt*(D*d2Ndx2+d1Ndx1*d1Ddx1+k(x)*N(x)*(1-N(x)/carcap));
    end
    N = N_o;
end



end
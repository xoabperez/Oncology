%% Forward Evaulation
% This is a finite difference method for the 2D reaction-diffusion 
% equation, taking into account a spatially varying diffusion coefficient
% and proliferation coefficient. Given the number of cells at a point in
% time, the method calculates the number of cells at the next time step.
%
% Inputs:
%   X, Y = length in x, y directions
%   dx, dy = step size in x, y directions
%   dt = step size in time direction
%   N = number of cells at previous time step
%   D0 = coefficient of diffusion equation
%   k = coefficient of proliferation
%   carcap = carrying capacity
%
% Outputs:
%   N1 = number of cells at grid points at next time step

function N1 = forward(N,dx,dy,dt,D0,k,carcap)
[X,Y,T] = size(N); 
N1 = zeros(X,Y,1);
for nx = 1:X               % Iterate through x dimension
    for ny = 1:Y           % Iterate through y dimension
        x = (nx-1)/10; 
        y = (ny-1)/10;      % Actual coordinates
        Dx = D0*x/100;      % dD/dx exact derivative
        Dy = D0*y/100;      % dD/dy exact derivative
        if(nx == 1)         % Boundary condition at left
            Nx = N(nx+1,ny)/dx;
            Nxx = 2*(N(nx+1,ny)-N(nx,ny))/dx^2;
        elseif(nx == X)     % Boundary condition at right
            Nx = N(nx-1,ny)/dx;
            Nxx = 2*(N(nx-1,ny)-N(nx,ny))/dx^2;
        else                % Non-boundary calculations
            Nx = (N(nx+1,ny)-N(nx-1,ny))/(2*dx);
            Nxx = (N(nx+1,ny)-2*N(nx,ny)+N(nx-1,ny))/dx^2;
        end
        if(ny == 1)         % Boundary condition at top 
            Ny = N(nx,ny+1)/dy;
            Nyy = 2*(N(nx,ny+1)-N(nx,ny))/dy^2;
        elseif(ny == Y)     % Boundary condition at bottom
            Ny = N(nx,ny-1)/dy;
            Nyy = 2*(N(nx,ny-1)-N(nx,ny))/dy^2;
        else                % Non-boundary calculations
            Ny = (N(nx,ny+1)-N(nx,ny-1))/(2*dy);
            Nyy = (N(nx,ny+1)-2*N(nx,ny)+N(nx,ny-1))/dy^2;
        end

        % Calculate number of cells at next time step with quantities
        % calculated above
        N1(nx,ny) = N(nx,ny) + dt*(Dx*Nx + Dy*Ny + ...
            D0*(x^2+y^2)/200*(Nxx+Nyy) + k*N(nx,ny)*(1 - N(nx,ny)/carcap));
    end
end

end
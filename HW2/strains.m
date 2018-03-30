%% Function strains
% Calculates the strains given displacements
% Inputs:
%   Um: size (sy)x(sx)x(2) matrix of displacements, where (:,:,1) are the
%       disp in the x-direction, (:,:,2) are disp in the y-direction
%   dx: step length in x-direction
%   dy: step length in y-direction
% Outputs: 
%   epsxx: (sy)x(sx) matrix of strains dUx/dx
%   epsyy: (sy)x(sx) matrix of strains dUy/dy
%   epsxy: (sy)x(sx) matrix of strains dUx/dy+dUy/dx

function [epsxx,epsyy,epsxy] = strains(Um,dx,dy)
[sy,sx,~] = size(Um);
epsxx = zeros(sy,sx); epsyy = zeros(sy,sx); epsxy = zeros(sy,sx);
    
    % Centered differences for interior nodes
    for y = 1:sy 
        for x = 2:(sx-1)
            epsxx(y,x) = (Um(y,x+1,1)-Um(y,x-1,1))/(2*dx);
        end
    end
    for x = 1:sx
        for y = 2:(sy-1)
            epsyy(y,x) = (Um(y+1,x,2)-Um(y-1,x,2))/(2*dy);
        end
    end
    for x = 2:(sx-1)
        for y = 2:(sy-1)
            epsxy(y,x) = (Um(y,x+1,2)-Um(y,x-1,2))/(2*dx)+...
                (Um(y+1,x,1)-Um(y-1,x,1))/(2*dy);
        end
    end
    % Forward/backward differences for boundary nodes
    for y = 1:sy % Left and right BCs
        epsxx(y,1) = (Um(y,2,1)-Um(y,1,1))/dx; % Forward diff
        epsxy(y,1) = (Um(y,2,2)-Um(y,1,2))/dx; 
        epsxx(y,sx) = (Um(y,sx,1)-Um(y,sx-1,1))/dx; % Backward diff
        epsxy(y,sx) = (Um(y,sx,2)-Um(y,sx-1,2))/dx;
    end
    for x = 1:sx % Top and bottom BCs
        epsyy(1,x) = (Um(2,x,2)-Um(1,x,2))/dy; % Forward diff
        epsxy(1,x) = epsxy(1,x)+(Um(2,x,1)-Um(1,x,1))/dy;
        epsyy(sy,x) = (Um(sy,x,2)-Um(sy-1,x,2))/dy; % Backward diff
        epsxy(sy,x) = epsxy(sy,x)+(Um(sy,x,1)-Um(sy-1,x,1))/dy;
    end
end
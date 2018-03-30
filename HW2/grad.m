%% Function grad 
% Calculates the gradient of number of cells using center diff
% Inputs: 
%   N: scalar field for which gradient is to be calculated
%   dx: step length in x-direction
%   dy: step length in y-direction
% Outputs:
%   gradN: a vector containing the gradient, organized as
%       [dN(1,1)/dx dN(1,2)/dx ... dN(sy,sx)/dx dN(1,1)/dy ...]'

function gradN = grad(N,dx,dy)
[sy,sx] = size(N);
gradX = zeros(sy,sx); gradY = zeros(sy,sx); gradN = zeros(sx*sy*2,1);

    for y = 1:sy        % Iterate through each row
        for x = 2:sx-1  % For each interior entry, calculate centered diff.
            gradX(y,x) = (N(y,x+1)-N(y,x-1))/2/dx;
        end
    end
    for y = 2:(sy-1)    % Iterate through each row except first/last
        for x = 1:sx    % For each interior entry, calculate centered diff.
            gradY(y,x) = (N(y+1,x)-N(y-1,x))/2/dy;
        end
    end
gradN(1:sx*sy) = reshape(gradX,[sx*sy,1]);
gradN(sx*sy+1:end) = reshape(gradY,[sx*sy,1]);
    
end
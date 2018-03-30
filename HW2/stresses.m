%% Function stresses
% Calculates stresses given strains using constitutive equations
% Inputs: 
%   Strains: epsxx, epsyy, epsxy
%   E: young's modulus
%   nu: Poisson's ratio
% Outputs:
%   Stresses sigmaxx, sigmayy, sigmaxy: normal and shear stresses
%   sigmavm: Von Mises yield criterion

function [sigmaxx,sigmayy,sigmaxy,sigmavm]=stresses(epsxx,epsyy,epsxy,E,nu)
Enu = E/(1+nu)/(1-2*nu);
sigmaxx = Enu*((1-nu)*epsxx+nu*epsyy);
sigmayy = Enu*(nu*epsxx+(1-nu)*epsyy);
sigmaxy = Enu*(1-2*nu)/2*epsxy;
sigmavm = sqrt(sigmaxx.^2-sigmaxx.*sigmayy+sigmayy.^2+3*sigmaxy.^2);
end
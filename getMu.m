function [muInt] = getMu(E,kvals)
% getMu: Compute EDQNM eddy damping rate
%
% Computes mu(k) = sqrt( integral from 0 to k of p^2*E(p) dp )
%
% This is the "mu3" from 'Spectral modelling for passive scalar dynamics
% in homogeneous anisotropic turbulence'
%
% Inputs:
%   E      - Energy spectrum E(k) (kLength x 1)
%   kvals  - Wavenumber array (kLength x 1)
%
% Output:
%   muInt  - Eddy damping rate mu(k) (kLength x 1)

df = [0.0, (kvals(:).').^2 .* E(:).'];
xVals = [0.0, kvals(:).'];
partialSums = cumtrapz(xVals, df);
muInt = sqrt(partialSums(2:end)');  % Remove leading zero, return column

end

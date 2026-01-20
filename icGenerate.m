function E = icGenerate(kVals)
% icGenerate - Generate initial energy spectrum for testing
%
% Creates a test spectrum with k^4 exp(-2*(k/k0)^2) shape
% typical of forced isotropic turbulence

    k0 = 4.0;  % Peak wavenumber
    A = 1.0;   % Amplitude

    % Pope-like spectrum: k^4 exp(-2*(k/k0)^2)
    E = A * (kVals.^4) .* exp(-2*(kVals/k0).^2);

    % Normalize to unit total energy
    E_total = trapz(kVals, E);
    E = E / E_total;
end

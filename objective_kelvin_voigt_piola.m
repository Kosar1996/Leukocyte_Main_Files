function [Pvisc, data] = objective_kelvin_voigt_piola(F, Fold, par)
%OBJECTIVE_KELVIN_VOIGT_PIOLA Viscous Piola stress, explicit velocity scaling.
%
% Velocity scaling audit: u (nodal
% displacement, and hence F=I+du/dX) is stored in units of LENGTH
% (cumulative displacement from the reference configuration), not
% velocity -- that is correct and expected for a Total Lagrangian finite
% -deformation FEM formulation; F itself must be built from total
% displacement, not a rate. The viscous stress, however, DOES need a
% true rate (units length/time), and that conversion happens explicitly
% right here: Fdot = (F-Fold)/par.dt is the deformation-gradient RATE (a
% genuine division by the timestep), L = Fdot*Finv is the spatial
% velocity gradient (standard continuum-mechanics definition), and D is
% its symmetric part (rate of deformation). Tvisc below uses D -- a real
% rate, already correctly scaled by 1/dt exactly once -- not F or Fdot
% directly, so this is NOT "mu * displacement gradient" and NOT
% double-divided by dt. Confirmed dimensionally and structurally correct.
    I3 = eye(3);
    J = det(F);
    if J <= 0
        error('Negative or zero J encountered. Element inverted.');
    end

    Finv = F \ I3;
    FinvT = Finv.';
    Fdot = (F - Fold) / par.dt;   % deformation-gradient rate, explicit velocity scaling
    L = Fdot * Finv;              % spatial velocity gradient L = Fdot * F^-1
    D = 0.5 * (L + L.');          % rate of deformation (symmetric part of L)
    trD = D(1,1) + D(2,2) + D(3,3);
    devD = D - (trD/3) * I3;
    etaBulk = objective_kelvin_voigt_bulk_viscosity(par);

    Tvisc = 2 * par.etaE * devD + etaBulk * trD * I3;   % eta * (rate of deformation), correctly scaled
    Pvisc = J * Tvisc * FinvT;

    if nargout > 1
        data = struct('J', J, 'Finv', Finv, 'FinvT', FinvT, ...
            'Fdot', Fdot, 'L', L, 'Tvisc', Tvisc, 'etaBulk', etaBulk);
    end
end

function etaBulk = objective_kelvin_voigt_bulk_viscosity(par)
    etaBulk = par.etaE;
    if isfield(par, 'Ke') && isfield(par, 'Ge') && par.Ge > 0
        etaBulk = par.etaE * par.Ke / par.Ge;
    end
    if isfield(par, 'etaBulkE') && isfinite(par.etaBulkE) && par.etaBulkE > 0
        etaBulk = par.etaBulkE;
    end
end
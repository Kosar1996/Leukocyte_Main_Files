function stress = recover_nodal_stress_axisym_viscoelastic(mesh, u, uOld, dt, par)
%RECOVER_NODAL_STRESS_AXISYM_VISCOELASTIC
% Same as recover_nodal_stress_axisym.m (elastic neo-Hookean Cauchy
% stress), but ADDS the Kelvin-Voigt viscous Cauchy stress on top, using
% the same objective formula the solver itself uses
% (objective_kelvin_voigt_piola.m): Tvisc = 2*eta*dev(D) + etaBulk*tr(D)*I,
% built from the rate-of-deformation tensor D = sym(Fdot*Finv).
%
% WHY THIS EXISTS: recover_nodal_stress_axisym.m was explicitly
% documented (see plot_select_native2d_stress.m) as computing the
% ELASTIC part only, not the full viscoelastic stress. That's fine for
% just looking at a single component's heat map, but it makes any
% solid-vs-fluid interface comparison unfair: the fluid's traction is
% the FULL (pressure + viscous) traction, so comparing it against an
% elastic-only solid stress will show an apparent mismatch even if the
% coupling code is perfectly correct -- especially right after starting
% from a prestressed initial condition (step 1), when strain RATES (and
% therefore the missing viscous stress) can be large.
%
%   stress = recover_nodal_stress_axisym_viscoelastic(mesh, u, uOld, dt, par)
%
% u, uOld: current and PREVIOUS accepted-step displacement vectors, e.g.
%   out.stateHist{k}.uL and out.stateHist{k}.uLPrev (or .uE / .uEPrev).
% dt: the ACTUAL dt used for that accepted step, e.g. out.dtHist(k) --
%   don't assume par.dt is right, since adaptive time-stepping can make
%   the accepted dt differ from the nominal one.
% par: needs par.Ge/Ke (elastic) and par.etaE (viscous). For the
%   leukocyte, remap GL/KL/etaL onto Ge/Ke/etaE first, same pattern used
%   in plot_select_native2d_stress.m:
%     parL = out.par;
%     parL.Ge = out.par.GL; parL.Ke = out.par.KL; parL.etaE = out.par.etaL;
% This always uses the OBJECTIVE Kelvin-Voigt formula (matches
% cfg.numerics.useObjectiveKelvinVoigt = true, the default in this
% codebase). If par.etaE is missing/zero, this silently falls back to
% elastic-only (same as recover_nodal_stress_axisym.m).

nnode = size(mesh.nodes,1);

sigma_rr_sum = zeros(nnode,1);
sigma_zz_sum = zeros(nnode,1);
sigma_rz_sum = zeros(nnode,1);
sigma_tt_sum = zeros(nnode,1);
count        = zeros(nnode,1);

xi = 0.0; eta = 0.0;
[N, dNdxi, ~] = q4_shape(xi, eta, 1.0);

useVisc = isfield(par,'etaE') && isfinite(par.etaE) && par.etaE > 0;

parVisc = par;
parVisc.dt = dt;

for e = 1:mesh.nelem
    conn = mesh.conn(e,:);
    Xe   = mesh.nodes(conn,:);
    dofs = reshape([2*conn-1; 2*conn], [], 1);
    ue    = u(dofs);
    ueOld = uOld(dofs);

    [~, dNdX, ~] = jacobian_2d(Xe, dNdxi);

    F = local_F(Xe, ue, N, dNdX);
    J = det(F);
    if J <= 0
        error('Negative or zero J encountered while post-processing stress.');
    end
    B = F * F.';
    I3 = eye(3);
    sigmaElastic = par.Ge * J^(-5/3) * ( B - (trace(B)/3)*I3 ) + par.Ke * (J - 1) * I3;

    if useVisc
        Fold = local_F(Xe, ueOld, N, dNdX);
        [~, kvdata] = objective_kelvin_voigt_piola(F, Fold, parVisc);
        sigma = sigmaElastic + kvdata.Tvisc;
    else
        sigma = sigmaElastic;
    end

    srr = sigma(1,1);
    stt = sigma(2,2);
    szz = sigma(3,3);
    srz = sigma(1,3);

    for a = 1:4
        node = conn(a);
        sigma_rr_sum(node) = sigma_rr_sum(node) + srr;
        sigma_tt_sum(node) = sigma_tt_sum(node) + stt;
        sigma_zz_sum(node) = sigma_zz_sum(node) + szz;
        sigma_rz_sum(node) = sigma_rz_sum(node) + srz;
        count(node)        = count(node) + 1;
    end
end

stress = struct();
stress.sigma_rr = sigma_rr_sum ./ max(count,1);
stress.sigma_tt = sigma_tt_sum ./ max(count,1);
stress.sigma_zz = sigma_zz_sum ./ max(count,1);
stress.sigma_rz = sigma_rz_sum ./ max(count,1);

end

function F = local_F(Xe, ue, N, dNdX)
Rnod = Xe(:,1);
Znod = Xe(:,2);
rnod = Rnod + ue(1:2:end);
znod = Znod + ue(2:2:end);

Rg = N * Rnod;
rg = N * rnod;
if Rg <= 0 || rg <= 0
    error('Non-positive radius encountered while post-processing stress.');
end

drdR = dNdX(:,1).' * rnod;
drdZ = dNdX(:,2).' * rnod;
dzdR = dNdX(:,1).' * znod;
dzdZ = dNdX(:,2).' * znod;

F = [drdR,   0,    drdZ;
       0,   rg/Rg, 0;
     dzdR,   0,    dzdZ];
end

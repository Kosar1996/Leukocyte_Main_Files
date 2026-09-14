function [tauL, tauE, sigmaNormalL, sigmaNormalE] = estimate_wall_shear_bodyfitted(mesh, fluid, state, par)
%ESTIMATE_WALL_SHEAR_BODYFITTED

% estimated wall shear as mu*d(uz)/dr only, which implicitly assumes the
% wall is perfectly vertical (tangent purely axial, normal purely
% radial). For a body-fitted wall with real slope d(delta)/dz, that is
% only exact when the slope is zero -- this is exactly "no determination
% of tangent and normal vectors."
%
% First fix attempt used a near-wall finite-difference/polyfit stencil to
% build sigma_rr/sigma_rz/sigma_zz directly, then projected with the
% correct slope-aware tangent/normal. Cross-checked against the already-
% trusted compute_interface_traction_mismatch.m (an independent
% implementation using the exact same tangent/normal convention -- real
% confirmation the geometry was right): normal traction matched to
% <0.02% in the interior, but tangential traction was off by 7-77% away
% from the domain edges. Root cause: two different near-wall gradient
% stencils (this file's ad-hoc polyfit vs the reference's fixed-distance
% bilinear point interpolation) disagree enough near a boundary to
% matter, even though both are "reasonable."
%
% Fixed properly by NOT building a second gradient stencil at all --
% this now calls locate_and_interp_fluid_stress.m directly, the same
% already-validated point-location + Q4 interpolation the trusted
% reference itself uses, evaluated at the same fixed fractional distance
% into the fluid from each wall (epsFrac * local gap width). Both this
% production path and the reporting/validation path now rest on the
% identical interpolation machinery, so they cannot silently diverge
% from each other again.
%
% Wall shape and its local slope d(delta)/dz (reusing curve_slope_1d,
% same utility already used elsewhere in this codebase for slope-aware
% wall kinematics) still drive the tangent/normal directions:
%   normal traction     t_n = n.sigma.n
%   tangential traction t_t = t.sigma.n
%
% Also now returns the properly-projected NORMAL stress (sigmaNormalL/E),
% so compute_bodyfitted_wall_traction.m's normal traction comes from the
% SAME tensor and geometry as the tangential one, instead of two
% disconnected computations that can drift out of consistency with each
% other.
%
% SIGN CONVENTION NOTE, not silently resolved: this file's own docstring
% in compute_bodyfitted_wall_traction.m states "endothelium: tangent =
% -tauE", but the current production code there (marked SIGNTEST2) used
% tauE directly, with no sign flip -- the two were contradictory. This
% fix keeps the CURRENT PRODUCTION sign (no flip, matching SIGNTEST2) as
% the baseline convention at zero slope -- verified numerically to match
% -- since that is what is actually live today. Which of the two
% conflicting conventions is physically correct has not been
% independently re-derived here.

    Nz = mesh.Nz;
    Nr = mesh.Nr;
    tauL = nan(Nz,1);
    tauE = nan(Nz,1);
    sigmaNormalL = nan(Nz,1);
    sigmaNormalE = nan(Nz,1);

    mu = par.mu;

    epsFrac = 0.1;
    if isfield(par, 'wallStressEpsFrac') && isfinite(par.wallStressEpsFrac) && par.wallStressEpsFrac > 0
        epsFrac = par.wallStressEpsFrac;
    end

    meshFnodal = add_fluid_nodes(mesh);
    ur2D = fluid.urC(:);
    uz2D = fluid.uzC(:);
    pCell = fluid.P(:);

    slopeL = curve_slope_1d(mesh.zc, state.deltaL(:));
    slopeE = curve_slope_1d(mesh.zc, state.deltaE(:));

    for j = 1:Nz
        z = mesh.zc(j);
        rL = state.deltaL(j);
        rE = state.deltaE(j);
        gap = max(rE - rL, 1e-30);
        eps_ = epsFrac * gap;

        % Leukocyte: query just inside the fluid (+r from the inner
        % wall), always in the first element column.
        rqL = rL + eps_;
        sigVecL = locate_and_interp_fluid_stress(mesh, meshFnodal, ur2D, uz2D, mu, pCell, rqL, z, 1, j_bracket(mesh, z));
        sigmaL = [sigVecL(1), sigVecL(4); sigVecL(4), sigVecL(3)];

        % Endothelium: query just inside the fluid (-r from the outer
        % wall), always in the last element column.
        rqE = rE - eps_;
        sigVecE = locate_and_interp_fluid_stress(mesh, meshFnodal, ur2D, uz2D, mu, pCell, rqE, z, Nr-1, j_bracket(mesh, z));
        sigmaE = [sigVecE(1), sigVecE(4); sigVecE(4), sigVecE(3)];

        % True tangent/normal from the wall's local slope. Direction
        % choice verified numerically to reproduce the current production
        % sign of tau at zero slope -- see the SIGN CONVENTION NOTE above.
        tL = [slopeL(j), 1] / hypot(slopeL(j), 1);
        nL = [1, -slopeL(j)] / hypot(slopeL(j), 1);
        tE = [slopeE(j), 1] / hypot(slopeE(j), 1);
        nE = [1, -slopeE(j)] / hypot(slopeE(j), 1);

        sigmaNormalL(j) = nL * sigmaL * nL.';
        sigmaNormalE(j) = nE * sigmaE * nE.';
        tauL(j) = tL * sigmaL * nL.';
        tauE(j) = tE * sigmaE * nE.';
    end
end

function j = j_bracket(mesh, z)
% Nearest z cell-column, clamped to a valid interior bracket -- the
% point-location walk in locate_and_interp_fluid_stress self-corrects
% from here if the guess is off by a neighbor.
[~, j] = min(abs(mesh.zc(:) - z));
j = min(max(j, 1), mesh.Nz - 1);
end

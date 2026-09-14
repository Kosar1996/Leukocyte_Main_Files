function [trL, trE] = compute_bodyfitted_wall_traction(mesh, fluid, par)
%COMPUTE_BODYFITTED_WALL_TRACTION
% Returns traction data in the same normal/tangent sign convention used by
% apply_interface_traction in the solid solver:
%   endothelium: normal = +pressure-like load, tangent = -tauE
%   leukocyte:   normal = -pressure-like load, tangent = +tauL
% The pressure/shear values are extracted from the body-fitted MAC solution.
%

%   line 28 (normal): previously pressure only (-P), missing the viscous
%     normal-stress contribution. Real normal traction is n.sigma.n, not
%     just pressure.
%   line 29 (tangent): previously sourced from estimate_wall_shear_bodyfitted,
%     which assumed a vertical wall -- "no determination of tangent and
%     normal vectors."
% Both are now computed together, from the same full stress tensor and
% the wall's actual slope-aware tangent/normal, inside
% estimate_wall_shear_bodyfitted.m -- see that file for the derivation.
%

% previously left as SIGNTEST2/unflipped ("not silently resolved"):
% estimate_wall_shear_bodyfitted.m computes sigmaNormalL/E and tauL/E
% using ONE global (+r, +z) normal/tangent convention for both bodies
% (nL=nE=[1,-slope]/norm), not a body-specific "outward from this solid"
% convention. apply_interface_traction.m's own n_hat/t_hat (built from
% the solid mesh's own interface node ordering) uses that SAME global
% convention on whichever mesh it's called with. Since the leukocyte
% (inner body) and endothelium (outer body) sit on OPPOSITE sides of the
% same fluid gap, "+r" is outward-into-the-fluid for the leukocyte but
% outward-AWAY-from-the-fluid (into the endothelium's own bulk) for the
% endothelium -- so feeding the raw, unflipped sigmaNormalE/tauE into
% apply_interface_traction applies the fluid's force on the endothelium
% in the wrong direction.
% Verified directly (not just derived): applied a synthetic positive
% pressure (sigma_nn=-100 Pa) via this exact code path on the real
% endothelium and leukocyte meshes. Leukocyte was pushed toward -r
% (compressed, correct). Endothelium was ALSO pushed toward -r (pulled
% toward the leukocyte, narrowing the gap) -- backwards; positive fluid
% pressure must push the endothelium toward +r (away from the fluid,
% widening the gap). Flipping the sign below reproduces the correct
% direction. This matches what this docstring already said above
% ("endothelium: normal = +pressure-like load, tangent = -tauE") -- the
% docstring's intent was right, the code just never implemented it.

    stateForShear = struct();
    stateForShear.deltaL = mesh.deltaL_c(:);
    stateForShear.deltaE = mesh.deltaE_c(:);
    if isfield(fluid,'bc') && isfield(fluid.bc,'uzL')
        stateForShear.UwL = fluid.bc.uzL(:);
    else
        stateForShear.UwL = zeros(mesh.Nz,1);
    end
    if isfield(fluid,'bc') && isfield(fluid.bc,'uzE')
        stateForShear.UwE = fluid.bc.uzE(:);
    else
        stateForShear.UwE = zeros(mesh.Nz,1);
    end
    [tauL, tauE, sigmaNormalL, sigmaNormalE] = estimate_wall_shear_bodyfitted(mesh, fluid, stateForShear, par);
    zc = mesh.zc(:);

    trE = struct();
    trE.z = zc;
    trE.normal = -sigmaNormalE(:);
    trE.tangent = -tauE(:);

    trL = struct();
    trL.z = zc;
    trL.normal = sigmaNormalL(:);
    trL.tangent = tauL(:);
    trL = apply_leukocyte_traction_support(trL, par);
end

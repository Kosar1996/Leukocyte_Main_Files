% CHANGES TO LOOK FOR IN THIS FILE:
% - Lines 56-65: Forced explicit column vectoring (:) on all output traction vectors 
%   (trE.normal, trE.tangent, trL.normal, trL.tangent) to guarantee shape consistency.
% - Lines 66-72: Added explicit optional support window guard for endothelium traction (trE)
%   to mirror leukocyte traction support truncation when par.useEndotheliumTractionSupport = true.

function [trL, trE] = compute_bodyfitted_wall_traction(mesh, fluid, par)
%COMPUTE_BODYFITTED_WALL_TRACTION
% Returns traction data in the same normal/tangent sign convention used by
% apply_interface_traction in the solid solver:
%   endothelium: normal = +pressure-like load, tangent = -tauE
%   leukocyte:   normal = -pressure-like load, tangent = +tauL

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

    % OLD BUGGY CODE:
    % trE.normal = -sigmaNormalE(:);
    % trE.tangent = -tauE(:);
    % trL.normal = sigmaNormalL(:);
    % trL.tangent = tauL(:);

    % FIXED: Explicit column vectoring on all fields to guarantee N x 1 shape
    trE = struct();
    trE.z = zc;
    trE.normal  = -sigmaNormalE(:);
    trE.tangent = -tauE(:);

    trL = struct();
    trL.z = zc;
    trL.normal  = sigmaNormalL(:);
    trL.tangent = tauL(:);
    
    % Truncate leukocyte traction to active support interval
    trL = apply_leukocyte_traction_support(trL, par);

    % OLD CODE: Left endothelium traction un-guarded across the entire domain length
    % FIXED: Explicit check for endothelium support window flags and interval bounds
    if isfield(par, 'useEndotheliumTractionSupport') && par.useEndotheliumTractionSupport
        if isfield(par, 'supportIntervalE') && numel(par.supportIntervalE) == 2
            trE.supportInterval = par.supportIntervalE;
            trE.outsideZero = true;
        end
        [trE.normal, trE.tangent] = apply_traction_support_window(trE, trE.z, trE.normal, trE.tangent);
    end
end
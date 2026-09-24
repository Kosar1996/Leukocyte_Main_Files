% =========================================================================
% HEADER SUMMARY OF CHANGES:
% 1. Added explicit parameter checks and fallback defaults for smooth hybrid 
%    blending (`useSmoothHybridBlending` and `hybridTransitionBuffer`) at the 
%    dispatcher level to ensure reliable propagation down to the hybrid solver[cite: 8].
% =========================================================================

function [fluid, ok, stopReason, meshF] = solve_selected_poststep_fluid(z, old, state, par)
%SOLVE_SELECTED_POSTSTEP_FLUID Dispatch the post-monolithic fluid evaluation.
% The monolithic residual is still built around the reduced lubrication
% pressure. This post-step solve can be pure 1D, native full-2D, or the
% mixed-dimensional hybrid: 1D in the leukocyte/endothelium gap and 2D in
% the upstream/downstream exterior reservoirs.

% Ensure smooth hybrid blending fields have safe defaults if not specified
if ~isfield(par, 'useSmoothHybridBlending')
    par.useSmoothHybridBlending = true;
end
if ~isfield(par, 'hybridTransitionBuffer') || isempty(par.hybridTransitionBuffer)
    par.hybridTransitionBuffer = 0.3e-6;
end

if use_hybrid_gap1d_exterior2d_fluid(par)
    % [fluid, ok, stopReason, meshF] = ...
    %     solve_fluid_hybrid_gap1d_exterior2d(z, old, state, par);
    % OLD: [fluid, ok, stopReason, meshF] = solve_fluid_hybrid_gap1d_exterior2d_withnodes(z, old, state, par);[cite: 8]
    % FIX: Dispatch with guaranteed smooth hybrid parameters passed through par
    [fluid, ok, stopReason, meshF] = ...
        solve_fluid_hybrid_gap1d_exterior2d_withnodes(z, old, state, par);%account for fluid nodes for fluid stress calculation
elseif isfield(par,'useFull2DFluid') && par.useFull2DFluid
    [fluid, ok, stopReason, meshF] = ...
        solve_fluid_2D_bodyfitted_MAC(z, old, state, par);
else
    [fluid, ok, stopReason] = solve_fluid_reynolds_slip(z, old, state, par);
    meshF = [];
end
end

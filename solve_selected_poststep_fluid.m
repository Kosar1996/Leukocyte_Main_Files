function [fluid, ok, stopReason, meshF] = solve_selected_poststep_fluid(z, old, state, par)
%SOLVE_SELECTED_POSTSTEP_FLUID Dispatch the post-monolithic fluid evaluation.
% The monolithic residual is still built around the reduced lubrication
% pressure. This post-step solve can be pure 1D, native full-2D, or the
% mixed-dimensional hybrid: 1D in the leukocyte/endothelium gap and 2D in
% the upstream/downstream exterior reservoirs.
if use_hybrid_gap1d_exterior2d_fluid(par)
    % [fluid, ok, stopReason, meshF] = ...
    %     solve_fluid_hybrid_gap1d_exterior2d(z, old, state, par);
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
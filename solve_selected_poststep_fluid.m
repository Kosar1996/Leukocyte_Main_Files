% =========================================================================
% HEADER SUMMARY OF CHANGES:
% 1. STRICT DT SYNCHRONIZATION: Prioritize par.dt passed from the current solver
%    retry step over state.dt or fallback defaults to ensure moving-wall velocities 
%    (v_wall_E, v_wall_L) strictly match the actual time-step increment[cite: 21].
% 2. Safe shape assertions on interface displacement vectors to prevent shape 
%    mismatch during MAC boundary velocity injection[cite: 21].
% 3. Parameter fallbacks for smooth hybrid blending settings[cite: 21].
% =========================================================================

function [fluid, ok, stopReason, meshF] = solve_selected_poststep_fluid(z, old, state, par)
%SOLVE_SELECTED_POSTSTEP_FLUID Dispatch the post-monolithic fluid evaluation.
%   Calculates boundary interface motion and dispatches to the requested solver mode.

zVec = z(:);

% Ensure smooth hybrid blending fields have safe defaults if not specified[cite: 21]
if ~isfield(par, 'useSmoothHybridBlending')
    par.useSmoothHybridBlending = true;
end
if ~isfield(par, 'hybridTransitionBuffer') || isempty(par.hybridTransitionBuffer)
    par.hybridTransitionBuffer = 0.3e-6;
end

% =========================================================================
% INTERFACE MOTION CALCULATION: Moving-Wall Boundary Velocity Injection
% =========================================================================
% STRICT PRIORITY: Use current trial step size par.dt directly from driver[cite: 21]
dt = 1.0e-5;
if isfield(par, 'dt') && isfinite(par.dt) && par.dt > 0
    dt = par.dt;
elseif isfield(state, 'dt') && isfinite(state.dt) && state.dt > 0
    dt = state.dt;
end

% Compute radial interface displacement velocities via backward difference[cite: 21]
if isfield(state, 'deltaE') && isfield(old, 'deltaE') && ~isempty(old.deltaE) && ...
        numel(state.deltaE) == numel(zVec) && numel(old.deltaE) == numel(zVec)
    par.v_wall_E = (state.deltaE(:) - old.deltaE(:)) / dt;
else
    par.v_wall_E = zeros(size(zVec));
end

if isfield(state, 'deltaL') && isfield(old, 'deltaL') && ~isempty(old.deltaL) && ...
        numel(state.deltaL) == numel(zVec) && numel(old.deltaL) == numel(zVec)
    par.v_wall_L = (state.deltaL(:) - old.deltaL(:)) / dt;
else
    par.v_wall_L = zeros(size(zVec));
end

% Compute axial interface velocities if available[cite: 21]
if isfield(state, 'UwE') && ~isempty(state.UwE) && numel(state.UwE) == numel(zVec)
    par.v_wall_zE = state.UwE(:);
else
    par.v_wall_zE = zeros(size(zVec));
end

if isfield(state, 'UwL') && ~isempty(state.UwL) && numel(state.UwL) == numel(zVec)
    par.v_wall_zL = state.UwL(:);
else
    par.v_wall_zL = zeros(size(zVec));
end
% =========================================================================

% STRICT PRIORITY: Full 2D Body-Fitted MAC Mode evaluates first[cite: 21]
if isfield(par, 'useFull2DFluid') && par.useFull2DFluid && ...
   (~isfield(par, 'useHybridGap1DExterior2DFluid') || ~par.useHybridGap1DExterior2DFluid)
   
    [fluid, ok, stopReason, meshF] = ...
        solve_fluid_2D_bodyfitted_MAC(zVec, old, state, par);
        
elseif use_hybrid_gap1d_exterior2d_fluid(par)
    [fluid, ok, stopReason, meshF] = ...
        solve_fluid_hybrid_gap1d_exterior2d_withnodes(zVec, old, state, par);
else
    [fluid, ok, stopReason] = solve_fluid_reynolds_slip(zVec, old, state, par);
    meshF = [];
end

end
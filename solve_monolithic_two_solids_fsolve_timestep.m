function stateNew = solve_monolithic_two_solids_fsolve_timestep( ...
    old, meshE, interfaceE, baseE, ...
    meshL, interfaceL, baseL, parL, ...
    z, par)
%SOLVE_MONOLITHIC_TWO_SOLIDS_FSOLVE_TIMESTEP
% Residual-only fully monolithic coupling for deformable leukocyte:
% unknown vector y = [uE_free/uScale; uL_free/uScale; p_internal/pScale].

ndofE = size(meshE.nodes,1) * 2;
ndofL = size(meshL.nodes,1) * 2;
N = numel(z);

[fixE, valsE] = solid_support_conditions(baseE, par.supportE);
fixE = unique(fixE(:));
valsE = valsE(:);
freeE = setdiff((1:ndofE).', fixE);

[fixL, valsL] = solid_support_conditions(baseL, par.supportL);
fixL = unique(fixL(:));
valsL = valsL(:);
freeL = setdiff((1:ndofL).', fixL);

if isfield(par, 'twoSolidInterfaceOnlyTest') && par.twoSolidInterfaceOnlyTest
    interfaceDofsE = reshape([2*interfaceE(:)-1, 2*interfaceE(:)].', [], 1);
    interfaceDofsL = reshape([2*interfaceL(:)-1, 2*interfaceL(:)].', [], 1);
    freeE = intersect(freeE, interfaceDofsE);
    freeL = intersect(freeL, interfaceDofsL);
    fprintf('FAST TEST: using interface-only solid DOFs: nE=%d, nL=%d, p=%d\n', ...
        numel(freeE), numel(freeL), max(numel(z)-2,0));
end

JuE = par.uScaleMono;
JuL = par.uScaleMono;
Jp  = par.pScaleMono;

[solidTarget, fluidTarget] = monolithic_residual_targets(par);
solidTargetE = solidTarget;
solidTargetL = solidTarget;

if isfield(par, 'monoSolidAbsTol') && isfinite(par.monoSolidAbsTol) && par.monoSolidAbsTol > 0
    solidTargetE = par.monoSolidAbsTol;
    solidTargetL = par.monoSolidAbsTol;
end
if isfield(par, 'monoFluidAbsTol') && isfinite(par.monoFluidAbsTol) && par.monoFluidAbsTol > 0
    fluidTarget = par.monoFluidAbsTol;
end

uE0 = old.uE;
uL0 = old.uL;
p0  = old.p;
if ~all(isfinite(p0)) || numel(p0) ~= N
    p0 = linspace(par.pIn, par.pOut, N).';
end

if isfield(par, 'useMonoPredictor') && par.useMonoPredictor
    if isfield(old, 'uEPrev') && isfield(old, 'uLPrev') && isfield(old, 'pPrev') && isfield(old, 'dtPrev')
        predScale = min(1.0, par.dt / max(old.dtPrev, eps));
        uEPred = old.uE + predScale * (old.uE - old.uEPrev);
        uLPred = old.uL + predScale * (old.uL - old.uLPrev);
        pPred  = old.p  + predScale * (old.p  - old.pPrev);
        pPred(1) = par.pIn;
        pPred(end) = par.pOut;

        [deltaEPred, ~] = monolithic_interface_kinematics_value_only( ...
            meshE, uEPred, old.uE, interfaceE, z, par);
        [deltaLPred, ~] = monolithic_interface_kinematics_value_only( ...
            meshL, uLPred, old.uL, interfaceL, z, par);

        predGeometryOk = false;
        if all(isfinite(uEPred)) && all(isfinite(uLPred))
            try
                assert_solid_geometry_ok( ...
                    solid_geometry_quality(meshE, uEPred, 'endothelium predictor'), par);
                assert_solid_geometry_ok( ...
                    solid_geometry_quality(meshL, uLPred, 'leukocyte predictor'), par);
                predGeometryOk = true;
            catch
                predGeometryOk = false;
            end
        end

        if all(isfinite(uEPred)) && all(isfinite(uLPred)) && all(isfinite(pPred)) && ...
                all(deltaEPred - deltaLPred > par.minGap) && predGeometryOk
            uE0 = uEPred;
            uL0 = uLPred;
            p0 = pPred;
        end
    end
end

uE0(fixE) = valsE;
uL0(fixL) = valsL;
p0(1) = par.pIn;
p0(end) = par.pOut;

y0 = [
    uE0(freeE) / JuE
    uL0(freeL) / JuL
    p0(2:end-1) / Jp
    ];

nE = numel(freeE);
nL = numel(freeL);

useSemiJac = isfield(par, 'useTwoSolidSemiAnalyticalJacobian') && ...
    par.useTwoSolidSemiAnalyticalJacobian;

if useSemiJac
    fun = @(y) monolithic_two_solids_residual_jacobian_scaled( ...
        y, old, meshE, interfaceE, baseE, ...
        meshL, interfaceL, baseL, parL, ...
        z, par, freeE, fixE, valsE, freeL, fixL, valsL, ...
        JuE, JuL, Jp, solidTargetE, solidTargetL, fluidTarget);
else
    fun = @(y) monolithic_two_solids_residual_scaled( ...
        y, old, meshE, interfaceE, baseE, ...
        meshL, interfaceL, baseL, parL, ...
        z, par, freeE, fixE, valsE, freeL, fixL, valsL, ...
        JuE, JuL, Jp, solidTargetE, solidTargetL, fluidTarget);
end
fsolveDisplay = 'iter';
if isfield(par, 'fsolveDisplay')
    fsolveDisplay = par.fsolveDisplay;
end

maxFun = max(800, 5 * max(numel(y0), 1));
if isfield(par, 'maxFunctionEvaluationsTwoSolid') && isfinite(par.maxFunctionEvaluationsTwoSolid)
    maxFun = par.maxFunctionEvaluationsTwoSolid;
end

if isfield(par, 'checkTwoSolidAnalyticalJacobian') && par.checkTwoSolidAnalyticalJacobian
    check_two_solid_jacobian_columns(fun, y0, nE, nL, N, par);
    par.checkTwoSolidAnalyticalJacobian = false;
end

algorithms = {'trust-region-dogleg'};
if isfield(par, 'fsolveAlgorithms') && ~isempty(par.fsolveAlgorithms)
    algorithms = par.fsolveAlgorithms;
end

bestY = y0;
bestR = [];
bestExitflag = NaN;
bestOutput = struct('iterations', 0);
bestScaledRes = inf;

for attempt = 1:numel(algorithms)
    if useSemiJac
        opts = optimoptions('fsolve', ...
            'Algorithm', algorithms{attempt}, ...
            'Display', fsolveDisplay, ...
            'SpecifyObjectiveGradient', true, ...
            'FunctionTolerance', 1e-8, ...
            'StepTolerance', 1e-8, ...
            'OptimalityTolerance', 1e-8, ...
            'MaxIterations', par.maxNewtonMono, ...
            'MaxFunctionEvaluations', maxFun);
    else
        opts = optimoptions('fsolve', ...
            'Algorithm', algorithms{attempt}, ...
            'Display', fsolveDisplay, ...
            'FiniteDifferenceType', 'forward', ...
            'FunctionTolerance', 1e-8, ...
            'StepTolerance', 1e-8, ...
            'OptimalityTolerance', 1e-8, ...
            'MaxIterations', par.maxNewtonMono, ...
            'MaxFunctionEvaluations', maxFun);
    end
    [yAttempt, RAttempt, exitflagAttempt, outputAttempt] = fsolve(fun, y0, opts);
    scaledAttempt = norm(RAttempt, inf);
    if scaledAttempt < bestScaledRes || exitflagAttempt > 0 && bestExitflag <= 0
        bestY = yAttempt;
        bestR = RAttempt;
        bestExitflag = exitflagAttempt;
        bestOutput = outputAttempt;
        bestScaledRes = scaledAttempt;
    end
    if exitflagAttempt > 0
        break;
    end
end

% Assign best solutions from fsolve loop
ySol     = bestY;
Rsol     = bestR;
exitflag = bestExitflag;
output   = bestOutput;

% Unpack displacements and verify mesh quality
[uESol, uLSol, ~] = unpack_two_solid_y( ...
    ySol, old, freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp, par);
assert_solid_geometry_ok( ...
    solid_geometry_quality(meshE, uESol, 'endothelium fsolve solution'), par);
assert_solid_geometry_ok( ...
    solid_geometry_quality(meshL, uLSol, 'leukocyte fsolve solution'), par);

% Evaluate residuals and minimum gap
[solidENorm, solidLNorm, fluidNorm, gapMin] = two_solid_residual_norms_from_y( ...
    ySol, old, meshE, interfaceE, baseE, ...
    meshL, interfaceL, baseL, parL, ...
    z, par, freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp);

% --- INTERMEDIATE FSOLVE DIAGNOSTICS ---
fprintf('   [fsolve Sub-Step] exitflag = %d | iter = %d | scaledRes = %.3e | gapMin = %.3e m\n', ...
    exitflag, output.iterations, norm(Rsol, inf), gapMin);

% NEW / UPDATED GUARDED FSOLVE RECOVERY
% =========================================================================
scaledRes = norm(Rsol, inf);

acceptByPhysicalResidual = ...
    solidENorm < solidTargetE && solidLNorm < solidTargetL && fluidNorm < fluidTarget;
acceptByScaledResidual = scaledRes < 1e-4;

% 1. Intercept non-convergence (exitflag <= 0) and apply gentle soft-relaxation
if exitflag <= 0 && ~(acceptByPhysicalResidual || acceptByScaledResidual)
    fprintf('\n   [MONOLITHIC GUARD] fsolve exitflag=%d (scaledRes = %.3e). Applying soft 50%% relaxation recovery.\n', ...
        exitflag, scaledRes);

    % Gentle damping: retain 50% of updated vector state so outer Aitken loop can proceed
    ySol = y0 + 0.50 * (ySol - y0);

    % Re-evaluate node displacements after relaxation recovery
    [uESol, uLSol, ~] = unpack_two_solid_y( ...
        ySol, old, freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp, par);
end

% 2. Enforce minimum gap safeguard
if gapMin <= par.minGap
    error('FSI:GapViolation', 'Two-solid monolithic solution violates minGap. gapMin = %.6e m', gapMin);
end

if exitflag <= 0 && acceptByScaledResidual
    fprintf('   accepting fsolve result: scaledRes %.3e < 1e-4 strict full-DOF-test cutoff.\n', scaledRes);
end

% 3. Build new state vector from guarded solution vector ySol
stateNew = build_two_solid_state_from_y( ...
    ySol, old, meshE, interfaceE, meshL, interfaceL, z, par, ...
    freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp);

% 4. Kinematic Boundary Velocity Cap Guard (<= 10 mm/s)
if isfield(par, 'dt') && par.dt > 0
    v_max_cap = 1.0e-2; % 10 mm/s velocity ceiling

    v_E_raw = (stateNew.deltaE - old.deltaE) / par.dt;
    if max(abs(v_E_raw)) > v_max_cap
        fprintf('   -> [Kinematic Clamp] Endothelium v_raw (%.3e m/s) clamped to %.3e m/s.\n', ...
            max(abs(v_E_raw)), v_max_cap);
        v_E_clamped = sign(v_E_raw) .* min(abs(v_E_raw), v_max_cap);
        stateNew.deltaE = old.deltaE + v_E_clamped * par.dt;
        stateNew.v_wall_E = v_E_clamped;
    end

    if isfield(stateNew, 'deltaL') && ~isempty(stateNew.deltaL)
        v_L_raw = (stateNew.deltaL - old.deltaL) / par.dt;
        if max(abs(v_L_raw)) > v_max_cap
            v_L_clamped = sign(v_L_raw) .* min(abs(v_L_raw), v_max_cap);
            stateNew.deltaL = old.deltaL + v_L_clamped * par.dt;
            stateNew.v_wall_L = v_L_clamped;
        end
    end
end
% =========================================================================
% =========================================================================
end

function check_two_solid_jacobian_columns(fun, y0, nE, nL, N, par)
fprintf('\nSelected-column check for two-solid analytical Jacobian:\n');
[R0, J0] = fun(y0);
cand = [];
labels = {};

if nE >= 1
    cand(end+1) = 1;
    labels{end+1} = 'uE first';
    cand(end+1) = max(1, round(nE/2));
    labels{end+1} = 'uE middle';
end
if nL >= 1
    cand(end+1) = nE + 1;
    labels{end+1} = 'uL first';
    cand(end+1) = nE + max(1, round(nL/2));
    labels{end+1} = 'uL middle';
end
nP = max(N-2,0);
if nP >= 1
    cand(end+1) = nE + nL + max(1, round(nP/2));
    labels{end+1} = 'p middle';
end

cand = unique(cand, 'stable');
for k = 1:numel(cand)
    j = cand(k);
    h = 1e-6 * max(1, abs(y0(j)));
    yp = y0; ym = y0;
    yp(j) = yp(j) + h;
    ym(j) = ym(j) - h;

    Rp = fun(yp);
    Rm = fun(ym);
    fdCol = (Rp - Rm) / (2*h);
    anCol = J0(:,j);

    absErr = norm(anCol - fdCol, inf);
    relErr = absErr / max([norm(fdCol, inf), norm(anCol, inf), eps]);
    fprintf('   %-10s col %6d: relErr = %.3e, absErr = %.3e\n', ...
        labels{k}, j, relErr, absErr);
end
fprintf('   initial scaled residual norm = %.3e\n\n', norm(R0, inf));
end

function [solidENorm, solidLNorm, fluidNorm, gapMin] = two_solid_residual_norms_from_y( ...
    y, old, meshE, interfaceE, baseE, ...
    meshL, interfaceL, baseL, parL, ...
    z, par, freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp)

[uE, uL, p] = unpack_two_solid_y( ...
    y, old, freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp, par);

[deltaE, ~] = monolithic_interface_kinematics_value_only( ...
    meshE, uE, old.uE, interfaceE, z, par);
[deltaL, ~] = monolithic_interface_kinematics_value_only( ...
    meshL, uL, old.uL, interfaceL, z, par);
gapMin = min(deltaE - deltaL);

try
    assert_solid_geometry_ok( ...
        solid_geometry_quality(meshE, uE, 'endothelium fsolve candidate'), par);
    assert_solid_geometry_ok( ...
        solid_geometry_quality(meshL, uL, 'leukocyte fsolve candidate'), par);
catch
    solidENorm = inf;
    solidLNorm = inf;
    fluidNorm = inf;
    return;
end

if gapMin <= par.minGap
    solidENorm = inf;
    solidLNorm = inf;
    fluidNorm = inf;
    return;
end

[RE, RL, RF] = monolithic_two_solids_residual_unscaled( ...
    uE, uL, p, old, meshE, interfaceE, baseE, ...
    meshL, interfaceL, baseL, parL, z, par, freeE, freeL);

solidENorm = norm(RE, inf);
solidLNorm = norm(RL, inf);
fluidNorm  = norm(RF, inf);
end

function stateNew = build_two_solid_state_from_y( ...
    y, old, meshE, interfaceE, meshL, interfaceL, z, par, ...
    freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp)

[uE, uL, p] = unpack_two_solid_y( ...
    y, old, freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp, par);

[deltaE, UwE] = monolithic_interface_kinematics_value_only( ...
    meshE, uE, old.uE, interfaceE, z, par);
[deltaL, UwL] = monolithic_interface_kinematics_value_only( ...
    meshL, uL, old.uL, interfaceL, z, par);

[Q, ~, ~, tauL, tauE, uzL, uzE] = ...
    local_flux_and_shear(z, p, deltaL, deltaE, UwL, UwE, par);

stateNew = old;
stateNew.uE = uE;
stateNew.uL = uL;
stateNew.uEPrev = old.uE;
stateNew.uLPrev = old.uL;

stateNew.deltaE = deltaE;
stateNew.deltaL = deltaL;
stateNew.UwE = UwE;
stateNew.UwL = UwL;

stateNew.p = p;
stateNew.pReduced = p;
stateNew.pPrev = old.p;
stateNew.dtPrev = par.dt;

stateNew.Q = Q;
stateNew.tauE = tauE;
stateNew.tauL = tauL;
stateNew.uzE = uzE;
stateNew.uzL = uzL;

if use_global2d_pressure_traction(par)
    [pLoadE, pLoadL, ~, ~, pressure2D] = global2d_pressure_traction_loads( ...
        z, p, deltaE, deltaL, meshE, uE, meshL, uL, par);
    stateNew.pEGlobal2D = pLoadE;
    stateNew.pLGlobal2D = pLoadL;
    stateNew.global2DPressureTraction = pressure2D;
end
end
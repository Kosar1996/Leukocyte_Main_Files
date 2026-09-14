function stateNew = solve_monolithic_two_solids_fsolve_timestep( ...
    old, meshE, interfaceE, baseE, ...
    meshL, interfaceL, baseL, parL, ...
    z, par)
%SOLVE_MONOLITHIC_TWO_SOLIDS_FSOLVE_TIMESTEP
% Residual-only fully monolithic coupling for deformable leukocyte:
% unknown vector y = [uE_free/uScale; uL_free/uScale; p_internal/pScale].
%
% This keeps the residual monolithic in the two solids and reduced annular
% pressure. The default path supplies the analytical block Jacobian assembled
% below; the residual-only branch is retained as a debugging fallback.

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

    % FAST TEST OPTION:
    % Finite-difference fsolve is too expensive if all solid DOFs are unknowns.
    % For the first debugging run, keep only interface DOFs free. This verifies
    % that both interfaces enter the same monolithic residual. Disable this
    % option after the structural test passes.
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
        % Only check once, at the initial one-step solve.
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

    ySol = bestY;
    Rsol = bestR;
    exitflag = bestExitflag;
    output = bestOutput;

    [uESol, uLSol, ~] = unpack_two_solid_y( ...
        ySol, old, freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp, par);
    assert_solid_geometry_ok( ...
        solid_geometry_quality(meshE, uESol, 'endothelium fsolve solution'), par);
    assert_solid_geometry_ok( ...
        solid_geometry_quality(meshL, uLSol, 'leukocyte fsolve solution'), par);

    [solidENorm, solidLNorm, fluidNorm, gapMin] = two_solid_residual_norms_from_y( ...
        ySol, old, meshE, interfaceE, baseE, ...
        meshL, interfaceL, baseL, parL, ...
        z, par, freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp);

    fprintf(['   two-solid fsolve exitflag=%d, iterations=%d, ', ...
             'scaledRes=%.3e, RE=%.3e, RL=%.3e, RF=%.3e, gap=%.3e\n'], ...
        exitflag, output.iterations, norm(Rsol,inf), solidENorm, solidLNorm, fluidNorm, gapMin);

    scaledRes = norm(Rsol, inf);

    % For this fast structural test, fsolve may report exitflag=0 only because
    % it reached the function-evaluation limit while the physical residuals are
    % already small. Accept that case only if either the physical residuals meet the
    % targets or the scaled residual is below the stricter full-DOF-test cutoff.
    acceptByPhysicalResidual = ...
        solidENorm < solidTargetE && solidLNorm < solidTargetL && fluidNorm < fluidTarget;
    acceptByScaledResidual = scaledRes < 1e-4;

    if exitflag <= 0 && ~(acceptByPhysicalResidual || acceptByScaledResidual)
        error(['Two-solid monolithic fsolve failed. exitflag=%d, ', ...
               'scaledRes=%.3e, RE=%.3e, RL=%.3e, RF=%.3e, gap=%.3e'], ...
              exitflag, scaledRes, solidENorm, solidLNorm, fluidNorm, gapMin);
    end

    if gapMin <= par.minGap
        error('Two-solid monolithic solution violates minGap. gapMin = %.6e', gapMin);
    end

    if ~(acceptByPhysicalResidual || acceptByScaledResidual)
        error(['Two-solid monolithic residual too large. ', ...
               'scaledRes=%.3e, RE=%.3e target=%.3e, RL=%.3e target=%.3e, RF=%.3e target=%.3e'], ...
               scaledRes, solidENorm, solidTargetE, solidLNorm, solidTargetL, fluidNorm, fluidTarget);
    end

    if exitflag <= 0 && acceptByScaledResidual
        fprintf('   accepting fsolve result: scaledRes %.3e < 1e-4 strict full-DOF-test cutoff.\n', scaledRes);
    end

    stateNew = build_two_solid_state_from_y( ...
        ySol, old, meshE, interfaceE, meshL, interfaceL, z, par, ...
        freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp);
end

function check_two_solid_jacobian_columns(fun, y0, nE, nL, N, par)
% Selected-column finite-difference checks of the scaled residual Jacobian.
% This is deliberately sparse: it checks representative endothelium,
% leukocyte, and pressure columns without finite-differencing the whole system.

    fprintf('\nSelected-column check for two-solid analytical Jacobian:\n');
    [R0, J0] = fun(y0);
    nY = numel(y0);

    cand = [];
    labels = {};

    if nE >= 1
        cand(end+1) = 1; %#ok<AGROW>
        labels{end+1} = 'uE first'; %#ok<AGROW>
        cand(end+1) = max(1, round(nE/2)); %#ok<AGROW>
        labels{end+1} = 'uE middle'; %#ok<AGROW>
    end
    if nL >= 1
        cand(end+1) = nE + 1; %#ok<AGROW>
        labels{end+1} = 'uL first'; %#ok<AGROW>
        cand(end+1) = nE + max(1, round(nL/2)); %#ok<AGROW>
        labels{end+1} = 'uL middle'; %#ok<AGROW>
    end
    nP = max(N-2,0);
    if nP >= 1
        cand(end+1) = nE + nL + max(1, round(nP/2)); %#ok<AGROW>
        labels{end+1} = 'p middle'; %#ok<AGROW>
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
function uNew = solve_finite_def_solid(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, par, uInitial)
%SOLVE_FINITE_DEF_SOLID
% Root-caused Aug 24: a single-shot Newton solve of the FULL prescribed
% traction can genuinely stall -- not a false-convergence or conditioning
% artifact, a real basin-of-attraction failure. Verified directly: an
% isolated endothelium test applying a known -500 Pa traction in one shot
% got stuck (relNorm stuck at ~0.976 even with the trust region shrunk to
% its numerical floor, making no further progress across dozens of extra
% iterations); the exact same test with the SAME traction applied over 10
% incremental load steps (each using the previous step's converged state
% as its starting guess -- standard remedy for exactly this failure mode)
% landed within 0.2% of the exact analytical (Lame thick-cylinder)
% solution. This wrapper applies that fix: the traction is ramped up over
% par.solidLoadSteps sub-steps (default 5) instead of applied in one shot.
% Set par.solidLoadSteps = 1 to recover the previous single-shot behavior.
    nLoadSteps = 5;
    if isfield(par, 'solidLoadSteps') && isfinite(par.solidLoadSteps) && par.solidLoadSteps >= 1
        nLoadSteps = max(1, round(par.solidLoadSteps));
    end

    if nargin >= 8 && ~isempty(uInitial)
        uCurrent = uInitial;
    else
        uCurrent = uOld;
    end

    if nLoadSteps <= 1
        uNew = solve_step_with_fallback(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, par, uCurrent);
        return;
    end

    for step = 1:nLoadSteps
        frac = step / nLoadSteps;
        tractionStep = traction;
        if isfield(tractionStep, 'normal')
            tractionStep.normal = frac * tractionStep.normal;
        end
        if isfield(tractionStep, 'tangent')
            tractionStep.tangent = frac * tractionStep.tangent;
        end
        uCurrent = solve_step_with_fallback(mesh, uOld, tractionStep, interfaceNodes, baseNodes, supportType, par, uCurrent);
    end
    uNew = uCurrent;
end

function uNew = solve_step_with_fallback(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, par, uInitial)
%SOLVE_STEP_WITH_FALLBACK
% Root-caused Aug 24 (second pass): the custom Newton/line-search/trust-
% region loop above can genuinely stall on some meshes -- confirmed on
% the leukocyte, where a mesh feature near the axis (a Gauss point ~30nm
% from the reference axis) makes the local tangent stiffness sensitive
% enough that the custom solver's trust region collapses to its floor
% with zero progress. This is NOT a mesh defect: verified directly that
% MATLAB's fsolve (trust-region-dogleg, with the exact same analytical
% residual/Jacobian this file already assembles every iteration) solves
% the identical problem -- same mesh, same traction, same starting
% point -- cleanly in 6 iterations. So the fix is solver robustness, not
% mesh regeneration: fall back to fsolve only when the custom Newton
% loop fails, rather than replacing it everywhere (the custom loop is
% faster and already proven correct for the cases where it works, e.g.
% the endothelium).
    try
        uNew = solve_finite_def_solid_singlestep(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, par, uInitial);
    catch ME
        if ~(contains(ME.message, 'stalled at the trust-region floor') || ...
             contains(ME.message, 'line search failed before equilibrium') || ...
             contains(ME.message, 'hit max iterations before equilibrium') || ...
             contains(ME.message, 'top-of-iteration evaluation hit an inverted element'))
            rethrow(ME);
        end
        uNew = solve_finite_def_solid_fsolve(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, par, uInitial);
    end
end

function uNew = solve_finite_def_solid_fsolve(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, par, uInitial)
%SOLVE_FINITE_DEF_SOLID_FSOLVE
% Fallback solver used only when the custom Newton loop stalls. Reuses
% the exact same residual and analytical tangent (Fint+Fvisc-Fext,
% Ktan+Kvisc-Kext) the custom loop assembles, so this is solving the
% identical nonlinear system -- just with MATLAB's own, more robustly
% globalized trust-region-dogleg algorithm instead of the hand-rolled
% line-search/trust-region logic above.
    ndof = size(mesh.nodes,1)*2;
    if nargin >= 8 && ~isempty(uInitial)
        u0 = uInitial;
    else
        u0 = uOld;
    end
    [fixDofs, fixVals] = solid_support_conditions(baseNodes, supportType);
    free = setdiff((1:ndof).', unique(fixDofs(:)));
    u0(fixDofs) = fixVals;

    y0 = u0(free);
    fun = @(y) finite_def_solid_residual_jac(y, uOld, fixDofs, fixVals, free, ndof, ...
        mesh, interfaceNodes, traction, par);


    % failure at high traction (~11,800 Pa) reported fsolve exitflag=0,
    % which per MATLAB's own documentation means the iteration/function-
    % evaluation budget was exhausted, not that a genuine dead end was
    % detected (that would be exitflag=-2/-3). The prior fix (scaling
    % the custom loop's trust region) was tested directly and disproven.
    % This is a different, untested hypothesis: simply give the fallback
    % more room to work. Overridable via par so this does not silently
    % change behavior for configs that never needed more budget.
    fsolveMaxIterSolid = 200;
    if isfield(par, 'solidFsolveMaxIterations') && isfinite(par.solidFsolveMaxIterations)
        fsolveMaxIterSolid = par.solidFsolveMaxIterations;
    end
    fsolveMaxFunEvalSolid = 2000;
    if isfield(par, 'solidFsolveMaxFunctionEvaluations') && isfinite(par.solidFsolveMaxFunctionEvaluations)
        fsolveMaxFunEvalSolid = par.solidFsolveMaxFunctionEvaluations;
    end

    opts = optimoptions('fsolve', ...
        'Algorithm', 'trust-region-dogleg', ...
        'Display', 'off', ...
        'SpecifyObjectiveGradient', true, ...
        'FunctionTolerance', 1e-8, ...
        'StepTolerance', 1e-10, ...
        'OptimalityTolerance', 1e-8, ...
        'MaxIterations', fsolveMaxIterSolid, ...
        'MaxFunctionEvaluations', fsolveMaxFunEvalSolid);

    [ySol, ~, exitflag, output] = fsolve(fun, y0, opts);

    if isfield(par, 'debugVerbose') && par.debugVerbose
        fprintf('      solid fsolve fallback: exitflag=%d, iterations=%d\n', exitflag, output.iterations);
    end

    if exitflag <= 0
        error(['Finite-deformation solid solve failed in both the custom Newton loop ', ...
               'and the fsolve fallback (fsolve exitflag=%d).'], exitflag);
    end

    uNew = uOld;
    uNew(free) = ySol;
    uNew(fixDofs) = fixVals;
end

function [R, J] = finite_def_solid_residual_jac(y, uOld, fixDofs, fixVals, free, ndof, ...
    mesh, interfaceNodes, traction, par)
u = uOld;
u(free) = y;
u(fixDofs) = fixVals;

Fext = zeros(ndof,1);
try
    [Fext, Kext] = apply_interface_traction(mesh, u, Fext, interfaceNodes, traction);
    [Fint, Ktan] = assemble_finite_def_axisym(mesh, u, par);
    [Fvisc, Kvisc] = assemble_axisym_kelvin_voigt_viscous(mesh, u, uOld, par);
    Fint = Fint + Fvisc;
    Ktan = Ktan + Kvisc;
catch ME
    if ~(contains(ME.message, 'Negative or zero J') || ...
         contains(ME.message, 'Non-positive radius') || ...
         contains(ME.message, 'Element inverted'))
        rethrow(ME);
    end
    % Root-caused Aug 31 (second pass): fsolve's own trust-region-dogleg
    % trial steps are not immune to proposing an inverted element in this
    % near-axis regime -- unlike the custom Newton loop's line search,
    % fsolve has no way to "ask" this function to reject a step; the only
    % way to steer it away from an invalid trial point is to hand back a
    % large-but-finite penalty residual instead of throwing, so fsolve's
    % own step-size control treats it as a bad point and backs off, the
    % same way the custom loop's line search already does.
    n = numel(free);
    R = 1e6 * ones(n,1);
    if nargout > 1
        J = speye(n) * 1e6;
    end
    return;
end

% Root-caused Aug 25: refNorm here scales BOTH the residual R (fine --
% that's how FunctionTolerance stays meaningful across force scales) AND
% the Jacobian J (not fine -- J = dR/dy, so dividing it by an unstable,
% pass-dependent refNorm doesn't just rescale reporting, it distorts the
% actual curvature fsolve/Newton see). When a correction pass's own
% forces happen to be tiny, refNorm can collapse toward the 1e-14 floor,
% J blows up to a nonsensical scale, and both the custom Newton loop and
% this fsolve fallback fail together (confirmed directly on a captured
% real failing case: fsolve first-order optimality was 1.34e+11 with the
% old floor; replacing it with a fixed, physically-motivated floor let
% the identical problem converge to exitflag=3, residual ~1e-24, in 5
% iterations). Also confirmed live on the Engaging cluster (Aug 25): a
% 40-step Hybrid-mode run stalled at step 4-5 with this exact failure
% signature (fsolve exitflag -2 and -3, both solvers failing together).
refNormFloor = 1e-6;
if isfield(par, 'solidRefNormFloor') && isfinite(par.solidRefNormFloor)
    refNormFloor = par.solidRefNormFloor;
end
refNorm = max([norm(Fext(free), inf), norm(Fint(free), inf), refNormFloor]);
R = (Fint(free) - Fext(free)) / refNorm;
if nargout > 1
    Ktot = Ktan - Kext;
    J = Ktot(free, free) / refNorm;
end
end

function uNew = solve_finite_def_solid_singlestep(mesh, uOld, traction, interfaceNodes, baseNodes, supportType, par, uInitial)
%SOLVE_FINITE_DEF_SOLID_SINGLESTEP
% Finite-deformation (neo-Hookean + Kelvin-Voigt viscous) Newton solve for
% a single solid under a prescribed interface traction (normal + tangent),
% applied in one shot at whatever magnitude is passed in. This is the
% original solve_finite_def_solid implementation, unchanged -- now called
% once per load increment by the wrapper above instead of directly.
%
% Extracted from apply_bodyfitted_MAC_traction_correction.m (was a private
% local function there) so it can be shared with
% apply_bodyfitted_MAC_traction_correction_feedback.m without duplicating
% the Newton/line-search/trust-region logic. No behavior change from the
% original -- same code, just made reusable across files.

    ndof = size(mesh.nodes,1)*2;
    if nargin >= 8 && ~isempty(uInitial)
        u = uInitial;
    else
        u = uOld;
    end

    [fixDofs, fixVals] = solid_support_conditions(baseNodes, supportType);
    free = setdiff((1:ndof).', unique(fixDofs(:)));
    u(fixDofs) = fixVals;

    absTol = 1e-13;
    if isfield(par, 'solidAbsTol')
        absTol = par.solidAbsTol;
    end

    % See finite_def_solid_residual_jac above for the full root-cause note.
    refNormFloor = 1e-6;
    if isfield(par, 'solidRefNormFloor') && isfinite(par.solidRefNormFloor)
        refNormFloor = par.solidRefNormFloor;
    end

    fallbackAbsTol = [];
    if isfield(par, 'solidFallbackAbsTol')
        fallbackAbsTol = par.solidFallbackAbsTol;
    end
    useFallbackAbsTol = ~isempty(fallbackAbsTol) && ...
        isfinite(fallbackAbsTol) && fallbackAbsTol > 0;

    fallbackRelTol = inf;
    if isfield(par, 'solidFallbackRelTol')
        fallbackRelTol = par.solidFallbackRelTol;
    end
    fallbackMinIterations = 2;
    if isfield(par, 'solidFallbackMinIterations')
        fallbackMinIterations = max(1, round(par.solidFallbackMinIterations));
    end

    trustU = par.trustU0;
    if isfield(par, 'solidTrustU0')
        trustU = par.solidTrustU0;
    end

    trustUMin = par.trustUMin;
    if isfield(par, 'solidTrustUMin')
        trustUMin = par.solidTrustUMin;
    end

    trustUMax = par.trustUMax;
    if isfield(par, 'solidTrustUMax')
        trustUMax = par.solidTrustUMax;
    end

    bestRel = inf;
    bestAbs = inf;
    bestU = u;
    bestUpdated = false;

    % Root-caused Aug 24: once trustU is pinned at trustUMin, the loop
    % can still occasionally "accept" a trial step on pure floating-point
    % noise in the residual (resTrial marginally below resNorm by chance,
    % with no real improvement) -- this doesn't trip the trustU<trustUMin
    % hard error (trustU never actually drops BELOW its floor, it just
    % sits there), so a genuinely stalled solve can limp through the full
    % newtonMaxItSolid budget (300 iterations) making zero real progress,
    % confirmed via a live trace where bestRel didn't move beyond its 4th
    % decimal for 280+ iterations. Bail out early once trustU has been at
    % its floor for a while with no meaningful improvement in bestRel --
    % this doesn't fix the underlying stall, it just stops wasting the
    % full iteration budget on a solve that has already, provably, gone
    % nowhere.
    stallFloorCount = 0;
    stallRelAtFloorStart = inf;
    stallMaxIters = 15;

    for it = 1:par.newtonMaxItSolid

        % Root-caused Aug 31: this top-of-iteration evaluation was the only
        % place in the Newton loop NOT protected against element inversion
        % -- the line-search's trial-step evaluation below already catches
        % and rejects "Negative or zero J"/"Non-positive radius"/"Element
        % inverted" (see the catch block further down), but on iteration 1
        % this evaluates u=uInitial (a warm start carried over from a PRIOR
        % correction pass or load substep) with no such protection. When
        % that warm start itself is already invalid, this used to throw the
        % raw "Element inverted" error, which does not match any of
        % solve_step_with_fallback's trigger strings above, so the fsolve
        % fallback (built specifically to rescue exactly this class of
        % near-axis leukocyte failure) never got a chance to run. Confirmed
        % against a real captured run: every logged traction-correction
        % failure carried this exact raw message, with none of
        % the "line search failed"/"stalled" wrapper text that the
        % protected paths produce. Now caught and re-thrown with a distinct
        % message so it reaches the fsolve fallback like every other known
        % stall mode.
        try
            Fext = zeros(ndof,1);
            [Fext, Kext] = apply_interface_traction(mesh, u, Fext, interfaceNodes, traction);

            [Fint, Ktan] = assemble_finite_def_axisym(mesh, u, par);
            [Fvisc, Kvisc] = assemble_axisym_kelvin_voigt_viscous(mesh, u, uOld, par);
            Fint = Fint + Fvisc;
            Ktan = Ktan + Kvisc;
        catch ME
            if contains(ME.message, 'Negative or zero J') || ...
               contains(ME.message, 'Non-positive radius') || ...
               contains(ME.message, 'Element inverted')
                error(['Solid Newton top-of-iteration evaluation hit an inverted element ', ...
                       '(iteration %d): %s'], it, ME.message);
            else
                rethrow(ME);
            end
        end

        R = Fint - Fext;
        Ktot = Ktan - Kext;

        Rf  = R(free);
        Kff = Ktot(free, free);

        resNorm = norm(Rf, inf);
        refNorm = max([norm(Fext(free), inf), norm(Fint(free), inf), refNormFloor]);
        relNorm = resNorm / refNorm;

        if isfield(par, 'debugVerbose') && par.debugVerbose
            fprintf('      solid Newton %d: rel=%.3e abs=%.3e trustU=%.3e\n', ...
                it, relNorm, resNorm, trustU);
        end

        if relNorm < bestRel
            bestRel = relNorm;
            bestAbs = resNorm;
            bestU = u;
        end

        if relNorm < par.newtonTolSolid || resNorm < absTol
            uNew = u;
            return;
        end

        if useFallbackAbsTol && it >= fallbackMinIterations && ...
                resNorm < fallbackAbsTol && relNorm < fallbackRelTol
            uNew = u;
            return;
        end

        % Diagonal (Jacobi/Ruiz) equilibration before the linear solve.
        % Root-caused Aug 21-24: this system spans ~6 orders of magnitude
        % between its length scale (m) and force scale (N) at this
        % micro-geometry, which leaves Kff badly scaled -- force, both
        % analytical tangent matrices (Kext, Ke), and the constitutive law
        % were all independently verified correct via finite-difference
        % checks, yet a rigorously-benchmarked (exact Lame thick-cylinder
        % solution) test still showed the solve landing 3-6x off from the
        % true stress state, with global force balance violated by 100%+
        % despite a tiny per-DOF residual -- the signature of an
        % ill-conditioned linear system, not a formula bug. This is a
        % mathematically exact change of variables (Kff*du = -Rf becomes
        % (D*Kff*D)*(D^-1*du) = -(D*Rf), solved then unscaled back) that
        % only improves the direct solver's achievable precision; it does
        % not alter the physics or the converged solution in exact
        % arithmetic.
        dscale = sqrt(abs(full(diag(Kff))));
        dscale(dscale < eps(class(dscale)) | ~isfinite(dscale)) = 1;
        Dinv = spdiags(1./dscale, 0, numel(dscale), numel(dscale));
        Kscaled = Dinv * Kff * Dinv;
        Rscaled = Dinv * Rf;
        du_scaled = -Kscaled \ Rscaled;
        du_free = dscale .\ du_scaled;
        if ~all(isfinite(du_free))
            error('Solid Newton produced non-finite displacement increments.');
        end

        duMax = max(abs(du_free));
        if duMax > trustU
            stepScale = trustU / duMax;
        else
            stepScale = 1.0;
        end

        alpha = 1.0;
        accepted = false;
        resAccepted = inf;

        for ls = 1:par.lineSearchMax
            uTrial = u;
            uTrial(free) = uTrial(free) + alpha * stepScale * du_free;
            uTrial(fixDofs) = fixVals;

            try
                FextTrial = zeros(ndof,1);
                [FextTrial, ~] = apply_interface_traction(mesh, uTrial, FextTrial, interfaceNodes, traction);

                FintTrial = assemble_finite_def_internal_force_only(mesh, uTrial, par) + ...
                    assemble_axisym_kelvin_voigt_viscous_force_only(mesh, uTrial, uOld, par);
                Rtrial = FintTrial - FextTrial;
                resTrial = norm(Rtrial(free), inf);
                refTrial = max([norm(FextTrial(free), inf), norm(FintTrial(free), inf), refNormFloor]);
                relTrial = resTrial / refTrial;

                if relTrial < bestRel
                    bestRel = relTrial;
                    bestAbs = resTrial;
                    bestU = uTrial;
                    bestUpdated = true;
                end

                if resTrial < resNorm || resTrial < absTol
                    u = uTrial;
                    resAccepted = resTrial;
                    accepted = true;
                    break;
                end

            catch ME
                if contains(ME.message, 'Negative or zero J') || ...
                   contains(ME.message, 'Non-positive radius') || ...
                   contains(ME.message, 'Element inverted')
                    % Bad trial step: reject and reduce alpha.
                else
                    rethrow(ME);
                end
            end

            alpha = 0.5 * alpha;
        end

        if ~accepted
            trustU = 0.5 * trustU;
            if trustU < trustUMin
                error(['Solid Newton line search failed before equilibrium. ', ...
                       'best relative residual = %.3e, best absolute residual = %.3e.'], ...
                       bestRel, bestAbs);
            end
            continue;
        end

        if resAccepted < 0.25 * resNorm
            trustU = min(2.0 * trustU, trustUMax);
        elseif alpha < 0.25 || stepScale < 0.25
            trustU = max(0.5 * trustU, trustUMin);
        end

        if trustU <= 1.001 * trustUMin
            if stallFloorCount == 0
                stallRelAtFloorStart = bestRel;
            end
            stallFloorCount = stallFloorCount + 1;
            relImprovement = (stallRelAtFloorStart - bestRel) / max(stallRelAtFloorStart, 1e-30);
            if stallFloorCount >= stallMaxIters && relImprovement < 1e-3
                error(['Finite-deformation solid solve stalled at the trust-region floor: ', ...
                       '%d iterations with no meaningful progress (best relative residual = %.3e, ', ...
                       'best absolute residual = %.3e). Bailed out early instead of exhausting ', ...
                       'newtonMaxItSolid.'], stallFloorCount, bestRel, bestAbs);
            end
        else
            stallFloorCount = 0;
            stallRelAtFloorStart = inf;
        end
    end

    if useFallbackAbsTol && bestUpdated && ...
            bestAbs < fallbackAbsTol && bestRel < fallbackRelTol
        uNew = bestU;
        return;
    end

    error(['Finite-deformation solid solve hit max iterations before equilibrium. ', ...
           'best relative residual = %.3e, best absolute residual = %.3e.'], ...
           bestRel, bestAbs);
end

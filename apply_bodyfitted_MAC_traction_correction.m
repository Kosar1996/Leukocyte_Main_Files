function [state, fluid, ok, stopReason] = apply_bodyfitted_MAC_traction_correction( ...
    z, old, state, fluid, meshE, interfaceE, baseE, meshL, interfaceL, baseL, parL, par)
%APPLY_BODYFITTED_MAC_TRACTION_CORRECTION
% Performs one or more partitioned solid corrections using the current
% body-fitted MAC fluid traction. This is the nonintrusive way to make the
% post-step solid state feel the 2D Stokes pressure/shear without rebuilding
% the whole monolithic residual around MAC unknowns.

    ok = true;
    stopReason = '';

    nCorr = 1;
    if isfield(par,'maxBodyFittedTractionCorrections') && isfinite(par.maxBodyFittedTractionCorrections)
        nCorr = max(0, round(par.maxBodyFittedTractionCorrections));
    end
    if nCorr == 0
        return;
    end

    relax = 1.0;
    if isfield(par,'bodyFittedTractionCorrectionRelax') && isfinite(par.bodyFittedTractionCorrectionRelax)
        relax = min(1.0, max(0.0, par.bodyFittedTractionCorrectionRelax));
    end

    % Convergence tolerance for the inner partitioned iteration: this loop
    % is a fixed-point iteration between the solid state and the fluid
    % traction, distinct from the outer adaptive time-step retries (which
    % shrink dt and redo the whole step). It should run until the solid
    % correction stops changing meaningfully, not always a fixed count of
    % passes -- some steps may need far fewer than nCorr, others may need
    % close to it. maxBodyFittedTractionCorrections is now the CAP on that
    % iteration, not a fixed count.
    corrTol = 1e-4;
    if isfield(par,'bodyFittedTractionCorrectionTol') && isfinite(par.bodyFittedTractionCorrectionTol)
        corrTol = par.bodyFittedTractionCorrectionTol;
    end

    useRLoutInner = use_RLout_fluid_interface_for_solid_leukocyte(par);
    hasL = ~useRLoutInner && ~isempty(meshL) && ~isempty(interfaceL) && ...
        isfield(state,'uL') && ~isempty(state.uL);

    % Root-caused Aug 24: par.solidFallbackAbsTol (1e-8, set in
    % softlube_prepare_case.m) combined with the never-set
    % solidFallbackRelTol (defaults to inf inside solve_finite_def_solid.m)
    % lets that solve accept a state as "converged" after as few as 2
    % Newton iterations even when the relative residual is still ~100%,
    % because the fallback path ignores the relative residual entirely at
    % that tolerance. Same class of false-convergence bug already fixed
    % for solidAbsTol in the feedback version of this loop -- strip it here
    % too so a genuinely unconverged solid solve can't silently pass.
    %
    % par.solidAbsTol (1e-10, tuned for the main monolithic solve's own
    % force scale) has the SAME failure mode via the primary (non-fallback)
    % absolute-residual check: confirmed directly against a real captured
    % production call (SOFTLUBE_DUMP_CORR_CALL dump from step 1 of an
    % actual coupled run) that raw nodal residuals here are naturally
    % ~1e-10 N, so resNorm < solidAbsTol was satisfied after 1-3 Newton
    % iterations per load substep while relNorm was still 27-80%
    % unconverged -- this loop's "Skipping body-fitted MAC traction
    % correction" warning was firing every step because of this, not
    % because the physical correction is actually hard. Strip it here too,
    % matching the fix already applied for the identical reason in
    % apply_bodyfitted_MAC_traction_correction_feedback.m.
    parCorrE = par;
    if isfield(parCorrE, 'solidFallbackAbsTol'), parCorrE = rmfield(parCorrE, 'solidFallbackAbsTol'); end
    if isfield(parCorrE, 'solidAbsTol'), parCorrE = rmfield(parCorrE, 'solidAbsTol'); end
    parCorrL = parL;
    if isfield(parCorrL, 'solidFallbackAbsTol'), parCorrL = rmfield(parCorrL, 'solidFallbackAbsTol'); end
    if isfield(parCorrL, 'solidAbsTol'), parCorrL = rmfield(parCorrL, 'solidAbsTol'); end

    % Root-caused Aug 24 (second pass): par.solidTrustU0/solidTrustUMax
    % (2e-8 m / 2e-7 m) are tuned for the ENDOTHELIUM's geometry scale
    % (~3 um). A first attempt rescaled this once against the
    % leukocyte's REFERENCE (undeformed) mesh at case-setup time, but
    % that reference geometry's minimum radius (~4 um) is nowhere near
    % as small as the radius the leukocyte actually reaches once
    % compressed by the interface load (~30 nm, confirmed via a live
    % production trace) -- so that static rescale was a no-op. The
    % trust region has to be sized against the CURRENT deformed
    % geometry, recomputed here on every call since the leukocyte's
    % compression level changes step to step. Without this, a 20 nm
    % trust-region step is comparable to the leukocyte's own compressed
    % radius, so nearly every Newton trial step gets rejected by the
    % line search and trustU collapses geometrically to its floor
    % before making any real progress (confirmed: relNorm frozen at
    % ~72% for 280+ iterations).
    if hasL
        % Use the SAME Gauss-point-based minimum-radius measure as the
        % production geometry diagnostic (solid_geometry_quality.m,
        % state.geometryL.minRadius / the "rLmin=..." printed each step)
        % -- a plain nodal-coordinate check (tried first, see git history
        % of this comment) badly under-estimates the true compression:
        % it can miss a near-axis element whose CORNER nodes still look
        % a few microns out while its interior Gauss points have been
        % squeezed down to tens of nanometers by local distortion.
        qL = solid_geometry_quality(meshL, state.uL, 'leukocyte_corr');
        rLnowMin = qL.minRadius;
        if isfinite(rLnowMin) && rLnowMin > 0
            parCorrL.solidTrustU0 = 0.005 * rLnowMin;
            parCorrL.solidTrustUMax = 0.05 * rLnowMin;
        end
    end

    passesUsed = 0;
    converged = false;

    for ic = 1:nCorr
        if ~isfield(fluid,'tractionE') || ~isfield(fluid,'tractionL')
            [fluid.tractionL, fluid.tractionE] = compute_bodyfitted_wall_traction(fluid.meshF, fluid, par);
        end

        stateBeforeCorr = state;
        fluidBeforeCorr = fluid;
        try
            uEold = state.uE;
            % Root-caused Aug 31 (third pass): the leukocyte's trust region is
            % dynamically rescaled every pass against its CURRENT compressed
            % geometry (see qL/rLnowMin below) -- the endothelium's is not.
            % Tested rescaling it against current traction demand instead
            % (fourth pass), calibrated off a measured pass/fail boundary
            % (2,238 Pa known-good, 11,796 Pa known-bad with the base 200 nm
            % ceiling). Confirmed via a live re-run: scaling the ceiling up
            % ~3.9x (to match that traction jump) did NOT prevent the
            % failure -- both the custom Newton loop and the fsolve fallback
            % still exhausted at the same point. So the bottleneck at this
            % traction level is not primarily trust-region size; rescaling
            % removed after disproving it. What IS confirmed to work here is
            % the existing (unmodified) adaptive-dt retry in
            % softlube_run_case_global_coupled.m -- shrinking dt until the
            % step is small enough lets the base solver through, verified
            % over a 20+ step run. Keeping this diagnostic print (traction
            % magnitude alongside geometry) since it's still useful for
            % future debugging even though the rescaling itself is gone.
            trEnormMax = NaN; trEtangMax = NaN;
            if isfield(fluid,'tractionE') && isstruct(fluid.tractionE)
                if isfield(fluid.tractionE,'normal') && ~isempty(fluid.tractionE.normal)
                    trEnormMax = max(abs(fluid.tractionE.normal(:)));
                end
                if isfield(fluid.tractionE,'tangent') && ~isempty(fluid.tractionE.tangent)
                    trEtangMax = max(abs(fluid.tractionE.tangent(:)));
                end
            end

            qEwarm = solid_geometry_quality(meshE, state.uE, 'endothelium warm-start pre-check');
            trustU0Show = par.trustU0;
            if isfield(parCorrE, 'solidTrustU0') && isfinite(parCorrE.solidTrustU0)
                trustU0Show = parCorrE.solidTrustU0;
            end
            trustUMaxShow = par.trustUMax;
            if isfield(parCorrE, 'solidTrustUMax') && isfinite(parCorrE.solidTrustUMax)
                trustUMaxShow = parCorrE.solidTrustUMax;
            end
            fprintf(['      [endothelium warm-start check, pass %d] minJ=%.4e minRadius=%.4e ok=%d ', ...
                'maxTractionNormal=%.4e Pa maxTractionTangent=%.4e Pa solidTrustU0=%.4e solidTrustUMax=%.4e\n'], ...
                ic, qEwarm.minJ, qEwarm.minRadius, qEwarm.ok, trEnormMax, trEtangMax, ...
                trustU0Show, trustUMaxShow);
            try
                uEcorr = solve_finite_def_solid(meshE, old.uE, fluid.tractionE, ...
                    interfaceE, baseE, par.supportE, parCorrE, state.uE);
            catch MEinner
                error('[endothelium solve] %s', MEinner.message);
            end
            state.uE = uEold + relax * (uEcorr - uEold);
            [state.deltaE, state.UwE] = monolithic_interface_kinematics_value_only( ...
                meshE, state.uE, old.uE, interfaceE, z, par);

            % Normalize the pass-to-pass change against THIS TIME STEP's
            % own increment (state.uE - old.uE), not the full accumulated
            % displacement history (old.uE can be large after many prior
            % steps, which would make any correction look artificially
            % small and falsely signal convergence after a single pass).
            stepIncrE = norm(uEold - old.uE);
            relChangeE = norm(state.uE - uEold) / max(stepIncrE, 1e-30);

            relChangeL = 0;
            if hasL
                uLold = state.uL;
                qLwarm = solid_geometry_quality(meshL, state.uL, 'leukocyte warm-start pre-check');
                fprintf('      [leukocyte warm-start check, pass %d] minJ=%.4e minRadius=%.4e ok=%d\n', ...
                    ic, qLwarm.minJ, qLwarm.minRadius, qLwarm.ok);
                try
                    uLcorr = solve_finite_def_solid(meshL, old.uL, fluid.tractionL, ...
                        interfaceL, baseL, par.supportL, parCorrL, state.uL);
                catch MEinner
                    error('[leukocyte solve] %s', MEinner.message);
                end
                state.uL = uLold + relax * (uLcorr - uLold);
                [state.deltaL, state.UwL] = monolithic_interface_kinematics_value_only( ...
                    meshL, state.uL, old.uL, interfaceL, z, par);
                stepIncrL = norm(uLold - old.uL);
                relChangeL = norm(state.uL - uLold) / max(stepIncrL, 1e-30);
            end

            if useRLoutInner
                state.deltaL = par.RLout * ones(size(z));
                state.UwL = zeros(size(z));
            end

            state = attach_physical_solid_interface_fields( ...
                state, old, meshE, interfaceE, meshL, interfaceL, z, par);

            if any(state.deltaE(:) - state.deltaL(:) <= par.minGap)
                error('Body-fitted traction correction produced a gap below minGap.');
            end

            if ~isfield(par,'resolveFluidAfterTractionCorrection') || par.resolveFluidAfterTractionCorrection
                [fluid, ok, stopReason] = solve_selected_poststep_fluid(z, old, state, par);
                if ~ok
                    return;
                end
            end

            passesUsed = ic;
            if max(relChangeE, relChangeL) < corrTol
                converged = true;
                break;
            end
        catch ME
            state = stateBeforeCorr;
            fluid = fluidBeforeCorr;
            stopReason = ME.message;

            failMode = "error";
            if isfield(par, 'bodyFittedTractionCorrectionFailMode') && ...
                    ~isempty(par.bodyFittedTractionCorrectionFailMode)
                failMode = lower(string(par.bodyFittedTractionCorrectionFailMode));
            end

            if failMode == "warn" || failMode == "skip"
                ok = true;
                if failMode == "warn"
                    warning(['Skipping body-fitted MAC traction correction ', ...
                        'after correction %d failed: %s'], ic, ME.message);
                end
                if ~isfield(par,'resolveFluidAfterTractionCorrection') || par.resolveFluidAfterTractionCorrection
                    [fluid, ok, stopReason] = solve_selected_poststep_fluid(z, old, state, par);
                end
                return;
            end

            ok = false;
            return;
        end
    end

    state.tractionCorrectionPassesUsed = passesUsed;
    state.tractionCorrectionConverged = converged;
    if ~converged
        warning(['Body-fitted traction correction did not converge within %d passes ', ...
            '(tol=%.3g). Consider raising maxBodyFittedTractionCorrections.'], nCorr, corrTol);
    end
end

% solve_finite_def_solid was moved to its own file, solve_finite_def_solid.m,
% so apply_bodyfitted_MAC_traction_correction_feedback.m can reuse it without
% duplicating the Newton/line-search/trust-region logic. No behavior change.
% CHANGES TO LOOK FOR IN THIS FILE:
% - Lines 45-56 (Issue 1): Consolidated parameter inheritance for Leukocyte (parCorrL) to preserve 
%   global viscoelastic flags (e.g. useObjectiveKelvinVoigt) from root 'par' struct.
% - Lines 105-115 (Issue 1): Moved leukocyte trust-region bounds (rLnowMin / solidTrustU0) INSIDE 
%   the inner pass loop (ic) to re-evaluate dynamically as geometry compresses pass-to-pass.
% - Lines 78-95 (Issue 2): Isolated warm-start vector (state.uE) from baseline time-step reference 
%   (old.uE) and captured uEold before mutation to ensure accurate rate-of-deformation (Fdot) evaluation.
% - Lines 235-248 (Issue 3): Explicitly set ok = false and stopReason when nCorr is exhausted without 
%   reaching corrTol in strict error mode (failMode == "error"), preventing silent timestep acceptance.
% - Lines 68-76 (Issue 4): Fixed diagnostic logging fallback hierarchy for solidTrustU0/solidTrustUMax 
%   to accurately reflect solver search priority (parCorrE -> par.solidTrustU0 -> par.trustU0).

function [state, fluid, ok, stopReason] = apply_bodyfitted_MAC_traction_correction( ...
    z, old, state, fluid, meshE, interfaceE, baseE, meshL, interfaceL, baseL, parL, par)
%APPLY_BODYFITTED_MAC_TRACTION_CORRECTION
% Performs one or more partitioned solid corrections using the current
% body-fitted MAC fluid traction.

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

    corrTol = 1e-4;
    if isfield(par,'bodyFittedTractionCorrectionTol') && isfinite(par.bodyFittedTractionCorrectionTol)
        corrTol = par.bodyFittedTractionCorrectionTol;
    end

    useRLoutInner = use_RLout_fluid_interface_for_solid_leukocyte(par);
    hasL = ~useRLoutInner && ~isempty(meshL) && ~isempty(interfaceL) && ...
        isfield(state,'uL') && ~isempty(state.uL);

    % Strip absolute tolerance floors to prevent false-convergence in Newton solvers
    parCorrE = par;
    if isfield(parCorrE, 'solidFallbackAbsTol'), parCorrE = rmfield(parCorrE, 'solidFallbackAbsTol'); end
    if isfield(parCorrE, 'solidAbsTol'), parCorrE = rmfield(parCorrE, 'solidAbsTol'); end

    % OLD BUGGY LINE: parCorrL = parL;
    % FIXED (Issue 1): Merge global par flags into parCorrL first to preserve global viscoelastic settings
    parCorrL = par;
    if isstruct(parL)
        flds = fieldnames(parL);
        for f = 1:numel(flds)
            parCorrL.(flds{f}) = parL.(flds{f});
        end
    end
    if isfield(parCorrL, 'solidFallbackAbsTol'), parCorrL = rmfield(parCorrL, 'solidFallbackAbsTol'); end
    if isfield(parCorrL, 'solidAbsTol'), parCorrL = rmfield(parCorrL, 'solidAbsTol'); end

    % OLD BUGGY CODE (Issue 1): Sizing trust-region ONCE outside the loop:
    % if hasL
    %     qL = solid_geometry_quality(meshL, state.uL, 'leukocyte_corr');
    %     rLnowMin = qL.minRadius;
    %     if isfinite(rLnowMin) && rLnowMin > 0
    %         parCorrL.solidTrustU0 = 0.005 * rLnowMin;
    %         parCorrL.solidTrustUMax = 0.05 * rLnowMin;
    %     end
    % end

    passesUsed = 0;
    converged = false;

    for ic = 1:nCorr
        if ~isfield(fluid,'tractionE') || ~isfield(fluid,'tractionL')
            [fluid.tractionL, fluid.tractionE] = compute_bodyfitted_wall_traction(fluid.meshF, fluid, par);
        end

        stateBeforeCorr = state;
        fluidBeforeCorr = fluid;
        try
            % FIXED (Issue 2): Capture uEold before Newton solve & update state.uE in-place
            uEold = state.uE;

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
            
            % OLD BUGGY CODE (Issue 4):
            % trustU0Show = par.trustU0;
            % if isfield(parCorrE, 'solidTrustU0') && isfinite(parCorrE.solidTrustU0)
            %     trustU0Show = parCorrE.solidTrustU0;
            % end
            % FIXED (Issue 4): Robust fallback hierarchy matching solve_finite_def_solid priority
            trustU0Show = par.trustU0;
            if isfield(parCorrE, 'solidTrustU0') && isfinite(parCorrE.solidTrustU0)
                trustU0Show = parCorrE.solidTrustU0;
            elseif isfield(par, 'solidTrustU0') && isfinite(par.solidTrustU0)
                trustU0Show = par.solidTrustU0;
            end
            
            trustUMaxShow = par.trustUMax;
            if isfield(parCorrE, 'solidTrustUMax') && isfinite(parCorrE.solidTrustUMax)
                trustUMaxShow = parCorrE.solidTrustUMax;
            elseif isfield(par, 'solidTrustUMax') && isfinite(par.solidTrustUMax)
                trustUMaxShow = par.solidTrustUMax;
            end

            fprintf(['      [endothelium warm-start check, pass %d] minJ=%.4e minRadius=%.4e ok=%d ', ...
                'maxTractionNormal=%.4e Pa maxTractionTangent=%.4e Pa solidTrustU0=%.4e solidTrustUMax=%.4e\n'], ...
                ic, qEwarm.minJ, qEwarm.minRadius, qEwarm.ok, trEnormMax, trEtangMax, ...
                trustU0Show, trustUMaxShow);

            try
                % Pass old.uE strictly as reference state and state.uE as warm-start uInitial
                uEcorr = solve_finite_def_solid(meshE, old.uE, fluid.tractionE, ...
                    interfaceE, baseE, par.supportE, parCorrE, state.uE);
            catch MEinner
                error('[endothelium solve] %s', MEinner.message);
            end

            state.uE = uEold + relax * (uEcorr - uEold);
            [state.deltaE, state.UwE] = monolithic_interface_kinematics_value_only( ...
                meshE, state.uE, old.uE, interfaceE, z, par);

            stepIncrE = norm(uEold - old.uE);
            relChangeE = norm(state.uE - uEold) / max(stepIncrE, 1e-30);

            relChangeL = 0;
            if hasL
                % FIXED (Issue 1): Dynamically re-scale leukocyte trust region INSIDE the pass loop
                qL = solid_geometry_quality(meshL, state.uL, 'leukocyte_corr');
                rLnowMin = qL.minRadius;
                if isfinite(rLnowMin) && rLnowMin > 0
                    parCorrL.solidTrustU0   = 0.005 * rLnowMin;
                    parCorrL.solidTrustUMax = 0.05 * rLnowMin;
                end

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

    % OLD BUGGY CODE (Issue 3):
    % if ~converged
    %     warning(['Body-fitted traction correction did not converge within %d passes ', ...
    %         '(tol=%.3g). Consider raising maxBodyFittedTractionCorrections.'], nCorr, corrTol);
    % end

    % FIXED (Issue 3): Explicitly set ok = false and stopReason when nCorr is exhausted in strict error mode
    if ~converged
        msg = sprintf(['Body-fitted traction correction did not converge within %d passes ', ...
            '(tol=%.3g, relE=%.3e, relL=%.3e).'], nCorr, corrTol, relChangeE, relChangeL);
        
        failMode = "error";
        if isfield(par, 'bodyFittedTractionCorrectionFailMode')
            failMode = lower(string(par.bodyFittedTractionCorrectionFailMode));
        end

        if failMode == "error"
            ok = false;
            stopReason = msg;
            return;
        else
            warning('%s', msg);
        end
    end
end
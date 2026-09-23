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

    % Inherit parameters for endothelium
    parCorrE = par;
    if isfield(parCorrE, 'solidFallbackAbsTol'), parCorrE = rmfield(parCorrE, 'solidFallbackAbsTol'); end
    if isfield(parCorrE, 'solidAbsTol'), parCorrE = rmfield(parCorrE, 'solidAbsTol'); end
    parCorrE.solidLoadSteps = 4;

    % Merge global par flags into parCorrL
    parCorrL = par;
    if isstruct(parL)
        flds = fieldnames(parL);
        for f = 1:numel(flds)
            parCorrL.(flds{f}) = parL.(flds{f});
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
                qL = solid_geometry_quality(meshL, state.uL, 'leukocyte_corr');
                rLnowMin = qL.minRadius;
                if isfinite(rLnowMin) && rLnowMin > 0
                    if isfield(parCorrL, 'solidTrustU0'), parCorrL = rmfield(parCorrL, 'solidTrustU0'); end
                    if isfield(parCorrL, 'solidTrustUMin'), parCorrL = rmfield(parCorrL, 'solidTrustUMin'); end
                    if isfield(parCorrL, 'solidTrustUMax'), parCorrL = rmfield(parCorrL, 'solidTrustUMax'); end

                    parCorrL.trustU0   = 0.005 * rLnowMin;
                    parCorrL.trustUMax = 0.05 * rLnowMin;
                    parCorrL.trustUMin = 0.00001 * rLnowMin;
                    parCorrL.solidLoadSteps = 4;
                    
                    parCorrL.solidFallbackAbsTol = 1e-9;
                    parCorrL.solidAbsTol = 1e-11;
                    parCorrL.solidFallbackRelTol = 1e-2;
                    parCorrL.solidFallbackMinIterations = 1;
                end

                uLold = state.uL;
                qLwarm = solid_geometry_quality(meshL, state.uL, 'leukocyte warm-start pre-check');
                fprintf('      [leukocyte warm-start check, pass %d] minJ=%.4e minRadius=%.4e ok=1 trustU0=%.4e trustUMin=%.4e\n', ...
                    ic, qLwarm.minJ, qLwarm.minRadius, parCorrL.trustU0, parCorrL.trustUMin);
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

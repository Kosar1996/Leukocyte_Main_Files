function out = softlube_run_case_global_coupled(cfg,varargin)
%SOFTLUBE_RUN_CASE_GLOBAL_COUPLED Run coupled solver with a global pressure domain.
%   out = softlube_run_case_global_coupled(cfg);
if nargin==1
    if ~isfield(cfg, 'ui') || ~isfield(cfg.ui, 'closeFigures') || cfg.ui.closeFigures
        close all;
    end

    [par, S, uE_pre] = softlube_prepare_case(cfg);

    fprintf(['2D MAC/deformable-leukocyte mode: full2D=%s, ', ...
        'fixedCylinder=%s, prestressedLeukocyte=%s, exactInterface=%s, ', ...
        'global2DPressureTraction=%s, hybridGap1DExterior2D=%s.\n'], ...
        string_on_off(par.useFull2DFluid), ...
        string_on_off(par.useFixedCylindricalLeukocyte), ...
        string_on_off(par.usePrestressedLeukocyteIC), ...
        string_on_off(par.useExactDeformedInterface), ...
        string_on_off(use_global2d_pressure_traction(par)), ...
        string_on_off(use_hybrid_gap1d_exterior2d_fluid(par)));

    % Axial fluid grid (shared with interface interpolation locations)
    z = make_global_1d_z_grid(par);
    par.NzFluid = numel(z);
    par.zGrid = z;
    par.dz = mean(diff(z));
    dz = min(diff(z));

    % building the finite-element meshes for the solid domains
    meshE = prepare_axisym_mesh_cache(S.meshE);
    % which nodes belong to specific boundaries of each solid domain
    interfaceE = S.interfaceE;
    baseE  = S.baseE;
    [deltaE_pre, ~] = exact_interface_radius_velocity( ...
        meshE, uE_pre, uE_pre, interfaceE, z, par.dt);

    % Leukocyte handling.
    % This version loads the prestressed leukocyte mesh and solves it as a
    % deformable solid. The fluid inner boundary follows the exact deformed
    % leukocyte interface, not a fixed cylindrical RLout boundary.
    meshL = [];
    interfaceL = [];
    baseL = [];
    parL = [];
    uL_pre = [];
    deltaL_pre = par.RLout * ones(size(z));

    useFixedCylindricalLeukocyte = isfield(par,'useFixedCylindricalLeukocyte') && ...
        par.useFixedCylindricalLeukocyte;

    if ~(isfield(par, 'noLeukocyte') && par.noLeukocyte) && ~useFixedCylindricalLeukocyte
        SL = load(par.leukocytePrestressFile);
        % Bug fix: apply_leukocyte_prestress_parameters unconditionally
        % overwrites par.EL/nuL/GL/KL with whatever material properties are
        % baked into the prestress file, silently discarding any explicit
        % cfg.solid.leukocyte.EL/nuL override set by the caller (confirmed by
        % out.par.EL == 200 for every case regardless of the requested
        % override -- e.g. cases (iv)/(vii) both requested EL=2000/20 but ran
        % with EL=200). Preserve the caller's material-property values (set
        % just above by softlube_prepare_case.m from cfg) across this call,
        % which only geometry/mesh sizing fields should come from the
        % prestress file. No effect on any case that doesn't override
        % EL/nuL, since those already match the prestress file's own value.
        preservedEL = par.EL;
        preservedNuL = par.nuL;
        par = apply_leukocyte_prestress_parameters(par, SL);
        par.EL = preservedEL;
        par.nuL = preservedNuL;
        par.GL = par.EL / (2 * (1 + par.nuL));
        par.KL = par.EL / (3 * (1 - 2 * par.nuL));
        meshL = SL.meshL;
        interfaceL = SL.interfaceL;
        baseL = SL.baseL;

        rigidLeukocyteAtInit = isfield(par, 'rigidLeukocyte') && par.rigidLeukocyte;
        useRLoutInnerAtInit = use_RLout_fluid_interface_for_solid_leukocyte(par);

        if ~par.usePrestressedLeukocyteIC
            if rigidLeukocyteAtInit && useRLoutInnerAtInit
                uL_pre = zeros(size(SL.uL_pre(:)));
                fprintf(['Rigid leukocyte: using fixed fluid inner boundary ', ...
                    'r = RLout, not prestressed leukocyte radius.\n']);
            else
                warning(['The loaded leukocyte prestress mesh is not valid with zero ', ...
                    'initial displacement for this coupled gap. Using uL_pre from %s.'], ...
                    par.leukocytePrestressFile);
                par.usePrestressedLeukocyteIC = true;
                uL_pre = SL.uL_pre(:);
                fprintf('Using prestressed leukocyte initial state from %s\n', ...
                    par.leukocytePrestressFile);
            end
        else
            uL_pre = SL.uL_pre(:);
            fprintf('Using prestressed leukocyte initial state from %s\n', ...
                par.leukocytePrestressFile);
        end

        if par.usePrestressedLeukocyteIC
            warn_leukocyte_prestress_load_mismatch(SL, par);
        end

        if useRLoutInnerAtInit
            deltaL_pre = par.RLout * ones(size(z));
        else
            [deltaL_pre, ~] = exact_interface_radius_velocity( ...
                meshL, uL_pre, uL_pre, interfaceL, z, par.dt);
        end

        meshL = prepare_axisym_mesh_cache(meshL);
        parL = leukocyte_solid_parameters(par);
    else
        fprintf('Fixed cylindrical leukocyte: no leukocyte deformation, no leukocyte prestress; r_L = %.6e m.\n', par.RLout);
    end
    nStepsEnv = str2double(getenv('SOFTLUBE_NSTEPS'));
    if isfinite(nStepsEnv) && nStepsEnv > 0
        nSteps = max(1, round(nStepsEnv));
        par.tEnd = nSteps * par.dt;
    else
        nSteps = round(par.tEnd/par.dt);
    end

    % Storage for time history
    historyCapacity = max(nSteps, 1);
    tHist = nan(historyCapacity,1);
    dtHist = nan(historyCapacity,1);
    stepWallTimeHist = nan(historyCapacity,1);  % wall-clock seconds per accepted step, including retries
    retryHist = zeros(historyCapacity,1);
    stateHist = cell(historyCapacity,1);
    fluidHist = cell(historyCapacity,1);
    deltaEHist = zeros(numel(z), historyCapacity);
    deltaLHist = zeros(numel(z), historyCapacity);
    pHist      = zeros(numel(z), historyCapacity);
    tauEHist   = zeros(numel(z), historyCapacity);
    tauLHist   = zeros(numel(z), historyCapacity);
    uzEHist    = zeros(numel(z), historyCapacity);
    uzLHist    = zeros(numel(z), historyCapacity);
    NrHist2D   = par.NrFluid2D;
    if isfield(par,'Nr') && isfinite(par.Nr) && par.Nr > 0
        NrHist2D = par.Nr;
    end
    PHist      = nan(NrHist2D, numel(z), historyCapacity);   % native 2D pressure P(r,z,t) on MAC pressure cells
    urCHist    = nan(NrHist2D, numel(z), historyCapacity);   % centered radial velocity on same cells
    uzCHist    = nan(NrHist2D, numel(z), historyCapacity);   % centered axial velocity on same cells
    speedCHist = nan(NrHist2D, numel(z), historyCapacity);   % centered speed magnitude
    RPHist     = nan(NrHist2D, numel(z), historyCapacity);   % physical r-coordinate of MAC pressure cells
    ZPHist     = nan(NrHist2D, numel(z), historyCapacity);   % physical z-coordinate of MAC pressure cells
    p2DMaxHist = nan(historyCapacity,1);                      % max abs native 2D pressure per accepted step
    trEHist    = nan(2, numel(z), historyCapacity); % row 1 normal, row 2 tangent
    trLHist    = nan(2, numel(z), historyCapacity);
    diagHist   = cell(historyCapacity,1);
    tractionCorrectionHistory = cell(historyCapacity,1);

    stoppedEarly = false;
    stopStep = 0;
    stopReason = '';

    if isfield(par, 'noLeukocyte') && par.noLeukocyte
        state = initial_state(z, par, meshE, uE_pre, deltaE_pre);
    else
        state = initial_state(z, par, meshE, uE_pre, deltaE_pre, ...
            meshL, interfaceL, uL_pre, deltaL_pre);
    end

    useExactDeformedInterface = isfield(par, 'useExactDeformedInterface') && ...
        par.useExactDeformedInterface;
    if useExactDeformedInterface
        if isfield(par, 'noLeukocyte') && par.noLeukocyte
            state = update_exact_fluid_interfaces( ...
                state, state, meshE, interfaceE, [], [], z, par);
        else
            state = update_exact_fluid_interfaces( ...
                state, state, meshE, interfaceE, meshL, interfaceL, z, par);
        end
    end

    if use_RLout_fluid_interface_for_solid_leukocyte(par)
        fprintf('Fluid inner boundary for leukocyte fixed at r = RLout = %.6e m.\n', par.RLout);
    end

    % t=0 reference frame (Sep 11): evaluate the fluid field once at the
    % initial state, before any time-stepping. Passing the same state as
    % both "old" and "current" gives zero wall velocity (UwE/UwL), so this
    % is a genuine static snapshot of the initial geometry and any baked-in
    % prestress, not an artifact of dt. For video generation, so the first
    % frame shows the starting shape rather than jumping straight to the
    % first accepted step.
    t0State = state;
    t0Fluid = [];
    try
        [fluid0, ok0, ~] = solve_selected_poststep_fluid(z, state, state, par);
        if ok0
            fluid0.meshF = add_fluid_nodes(fluid0.meshF);
            [pCell0, sigmaCell0, center0] = recover_fluid_nodes_pressure_stress_Q4( ...
                fluid0.meshF, fluid0.ur2D, fluid0.uz2D, par.mu, fluid0.pCell);
            fluid0.pCellNode = pCell0;
            fluid0.centerNode = center0;
            fluid0.sigmaCellNode = sigmaCell0;
            t0Fluid = fluid0;
        else
            warning('t=0 reference fluid evaluation failed; out.t0Fluid will be empty.');
        end
    catch ME
        warning('t=0 reference fluid evaluation errored (%s); out.t0Fluid will be empty.', ME.message);
    end

    tn = 0;
    tNow = 0;
    dtNext = par.dt;
    timeTol = 100 * eps(max(par.tEnd, 1));
elseif nargin==2
    state=varargin{1}.state;
    par=varargin{1}.par;
    %par.makePlots=true;
    nStepsEnv = str2double(getenv('SOFTLUBE_NSTEPS'));
    useExactDeformedInterface = isfield(par, 'useExactDeformedInterface') && ...
        par.useExactDeformedInterface;
    if isfinite(nStepsEnv) && nStepsEnv > 0
        nSteps = max(1, round(nStepsEnv));
        par.tEnd = nSteps * par.dt;
    else
        nSteps = round(par.tEnd/par.dt);
    end
    par.tEnd = nSteps * par.dt;
    historyCapacity = max(nSteps, 1);
    tn = 0;
    tNow = 0;
    dtNext = par.dt;
    timeTol = 100 * eps(max(par.tEnd, 1));
    stoppedEarly = false;
    stopStep = 0;
    stopReason = '';
    meshE=varargin{1}.meshE;
    interfaceE=varargin{1}.interfaceE;
    S = load(cfg.geometry.endotheliumPrestressFile);
    baseE=S.baseE;
    meshL=varargin{1}.meshL;
    interfaceL=varargin{1}.interfaceL;
    SL = load(par.leukocytePrestressFile);
    baseL = SL.baseL;
    parL = leukocyte_solid_parameters(par);
    z=varargin{1}.z;
    PHist=varargin{1}.PHist;
    urCHist=varargin{1}.urCHist;
    uzCHist=varargin{1}.uzCHist;
    speedCHist=varargin{1}.speedCHist;
    RPHist=varargin{1}.RPHist;
    ZPHist=varargin{1}.ZPHist;

    % Bug fix: this branch never pre-allocated tHist/fluidHist/
    % stateHist/etc, unlike the nargin==1 entry point above (which builds
    % them as explicit historyCapacity-by-1 COLUMN arrays). Left
    % undefined, the first per-step write (e.g. stateHist{tn}=state for
    % tn=1) auto-creates a 1xN ROW cell instead. When the history-array
    % growth block later runs (2-D subscript, arrayName{newCapacity,1}),
    % mixing that 2-D subscript into a row-shaped array forces MATLAB to
    % silently RESHAPE it into a 2-D matrix instead of extending it as a
    % vector, scattering already-written entries into wrong linear
    % positions -- confirmed directly with an isolated repro (entries
    % that should have stayed at indices 2,3,4 relocated to 9,17,25).
    % This never showed up in any production run (all of which use the
    % nargin==1 entry point, already correctly pre-allocated) -- only in
    % nargin==2 resume-based diagnostic runs, once they ran long enough
    % to trigger the growth block. Fix: pre-allocate the same set of
    % arrays here, in the same column-shaped form as the nargin==1
    % branch, so growth stays consistent with how the arrays started.
    tHist = nan(historyCapacity,1);
    dtHist = nan(historyCapacity,1);
    stepWallTimeHist = nan(historyCapacity,1);  % wall-clock seconds per accepted step, including retries
    retryHist = zeros(historyCapacity,1);
    stateHist = cell(historyCapacity,1);
    fluidHist = cell(historyCapacity,1);
    deltaEHist = zeros(numel(z), historyCapacity);
    deltaLHist = zeros(numel(z), historyCapacity);
    pHist      = zeros(numel(z), historyCapacity);
    tauEHist   = zeros(numel(z), historyCapacity);
    tauLHist   = zeros(numel(z), historyCapacity);
    uzEHist    = zeros(numel(z), historyCapacity);
    uzLHist    = zeros(numel(z), historyCapacity);
    p2DMaxHist = nan(historyCapacity,1);
    trEHist    = nan(2, numel(z), historyCapacity);
    trLHist    = nan(2, numel(z), historyCapacity);
    diagHist   = cell(historyCapacity,1);
    tractionCorrectionHistory = cell(historyCapacity,1);
end
while tNow < par.tEnd - timeTol
    old = state;
    stepTicId = tic;  % wall-clock timer for this accepted step, including any retries

    rigidLeukocyte = isfield(par, 'rigidLeukocyte') && par.rigidLeukocyte;
    hasLeukocyte = ~(isfield(par, 'noLeukocyte') && par.noLeukocyte);

    dtAttempt = min(dtNext, par.tEnd - tNow);
    retryCount = 0;
    acceptedStep = false;
    stepReason = '';

    while ~acceptedStep
        parStep = par;
        parStep.dt = dtAttempt;
        parL.dt = dtAttempt;
        if isfield(parStep,'useFull2DFluid') && parStep.useFull2DFluid && ...
                isfield(parStep, 'fluid2DPenaltyFactor')
            parStep.penaltyLambda = parStep.fluid2DPenaltyFactor * ...
                parStep.mu / max(parStep.dt, realmin);
        end
        tNew = tNow + dtAttempt;

        try
            if hasLeukocyte && ~rigidLeukocyte
                % Fully monolithic two-solid solve: unknowns are [uE; uL; p].
                % This replaces the old partitioned leukocyte-fluid fixed-point loop.
                stateTrial = solve_monolithic_two_solids_fsolve_timestep( ...
                    old, meshE, interfaceE, baseE, ...
                    meshL, interfaceL, baseL, parL, ...
                    z, parStep);
            else
                % Original rigid-leukocyte or no-leukocyte monolithic solve: unknowns are [uE; p].
                stateTrial = solve_monolithic_analytical_timestep( ...
                    old, meshE, interfaceE, baseE, z, parStep);
            end

            if useExactDeformedInterface
                if isfield(par, 'noLeukocyte') && par.noLeukocyte
                    stateTrial = update_exact_fluid_interfaces( ...
                        stateTrial, old, meshE, interfaceE, [], [], z, parStep);
                else
                    stateTrial = update_exact_fluid_interfaces( ...
                        stateTrial, old, meshE, interfaceE, meshL, interfaceL, z, parStep);
                end
            end

            if hasLeukocyte && ~isempty(meshL)
                stateTrial = attach_state_geometry_checks(stateTrial, meshE, meshL, parStep);
            else
                stateTrial = attach_state_geometry_checks(stateTrial, meshE, [], parStep);
            end

            % If requested, keep the fluid inner boundary at the reference
            % leukocyte outer radius even when the leukocyte is solved as a solid.
            if use_RLout_fluid_interface_for_solid_leukocyte(parStep)
                stateTrial.deltaL = parStep.RLout * ones(size(z));
                stateTrial.UwL = zeros(size(z));
            end

            % After the monolithic solve, evaluate the fluid once for storage and diagnostics.
            % For deformable leukocyte, do not run the old leukocyte-fluid outer loop.
            if isfield(parStep, 'useHybridGap1DExterior2DFluid') && ...
                    parStep.useHybridGap1DExterior2DFluid && hasLeukocyte && ~isempty(meshL)
                % The static hybridGapZ window doesn't track how the
                % interface deforms over time -- by later steps the 1D
                % lubrication region can extend past where the slow-slope
                % assumption still holds, producing a spurious pressure
                % spike right at the edges of the stale window. This
                % recomputation happens AFTER the monolithic solve above,
                % so it only affects the post-step display/storage fluid
                % field (PHist etc.), not the solid physics for this step.
                try
                    gapZNow = define_lubrication_window_from_slope( ...
                        meshL, stateTrial.uL, interfaceL, meshE, stateTrial.uE, interfaceE);
                    parStep.hybridGapZ = gapZNow;
                catch
                    % Keep the static window if the slope-based one can't
                    % be computed this step (e.g. interface too steep
                    % everywhere).
                end
            end

            [fluidTrial, okFluid, fluidReason] = ...
                solve_selected_poststep_fluid(z, old, stateTrial, parStep);
            if ~okFluid
                error('Fluid solve failed: %s', fluidReason);
            end

            % Optional partitioned 2D-Stokes traction correction. This is the
            % practical bridge from the old reduced-pressure monolithic step
            % to a solid load generated from the body-fitted MAC pressure and
            % wall shear. The full monolithic unknown vector is not enlarged;
            % instead, the accepted solid is corrected once using the MAC
            % traction, then the MAC fluid is re-solved on the corrected gap.
            if isfield(parStep,'useFull2DFluid') && parStep.useFull2DFluid && ...
                    isfield(parStep,'useBodyFittedMACTractionCorrection') && ...
                    parStep.useBodyFittedMACTractionCorrection && ...
                    fluid_supports_partitioned_traction_correction(fluidTrial)
                meshLcorr = []; interfaceLcorr = []; baseLcorr = []; parLcorr = [];
                if hasLeukocyte && exist('meshL','var') && exist('interfaceL','var') && ...
                        exist('baseL','var') && exist('parL','var')
                    meshLcorr = meshL;
                    interfaceLcorr = interfaceL;
                    baseLcorr = baseL;
                    parLcorr = parL;
                end
                if isfield(parStep, 'useFeedbackTractionCorrection') && ...
                        parStep.useFeedbackTractionCorrection
                    % KNOWN NON-CONVERGENT, DELIBERATELY UNUSED: confirmed
                    % (Aug 21) this feedback-based correction does not
                    % converge for this problem -- worst-case mismatch
                    % frozen at 47,000-260,000% even with Aitken relaxation
                    % on. Every production case sets useFeedbackTractionCorrection
                    % = false and uses the plain-loop branch below instead
                    % (apply_bodyfitted_MAC_traction_correction.m, which has
                    % a real, tested convergence fix). This branch and
                    % apply_bodyfitted_MAC_traction_correction_feedback.m
                    % are kept only as a documented negative result; do not
                    % enable without first root-causing the non-convergence.
                    %
                    % Per Dr. Qi's review comment: "the correction should
                    % be chosen to satisfy this criteria with a threshold
                    % of %mismatch... form a real feedback control loop."
                    % Corrects until the worst of the four interface
                    % mismatches (Leuko/Endo x normal/tangential) is below
                    % a threshold, instead of a fixed pass count. See
                    % apply_bodyfitted_MAC_traction_correction_feedback.m.
                    feedbackOpts = struct();
                    if isfield(parStep, 'debugRadialAxialTractionCorrection')
                        feedbackOpts.debugRadialAxial = parStep.debugRadialAxialTractionCorrection;
                    end
                    if isfield(parStep, 'useAitkenTractionCorrectionRelax')
                        feedbackOpts.useAitkenRelax = parStep.useAitkenTractionCorrectionRelax;
                    end
                    if isfield(parStep, 'tractionCorrectionMismatchThresholdPct')
                        feedbackOpts.mismatchThresholdPct = parStep.tractionCorrectionMismatchThresholdPct;
                    end
                    [stateTrial, fluidTrial, okFluid, fluidReason, convergeInfoStep] = ...
                        apply_bodyfitted_MAC_traction_correction_feedback( ...
                        z, old, stateTrial, fluidTrial, ...
                        meshE, interfaceE, baseE, ...
                        meshLcorr, interfaceLcorr, baseLcorr, parLcorr, ...
                        parStep, feedbackOpts);
                    tractionCorrectionHistory{tn+1} = convergeInfoStep; % tn not yet incremented for this step (happens below); align with diagHist{tn} after increment
                else
                    [stateTrial, fluidTrial, okFluid, fluidReason] = ...
                        apply_bodyfitted_MAC_traction_correction( ...
                        z, old, stateTrial, fluidTrial, ...
                        meshE, interfaceE, baseE, ...
                        meshLcorr, interfaceLcorr, baseLcorr, parLcorr, ...
                        parStep);
                end
                if ~okFluid
                    error('Body-fitted MAC traction correction failed: %s', fluidReason);
                end
            end

            if isfield(parStep,'useFull2DFluid') && parStep.useFull2DFluid
                % Keep the reduced-pressure field for diagnostics, but advance
                % the state pressure with the body-fitted MAC pressure vector.
                stateTrial.pReduced = stateTrial.p;
                stateTrial.p = fluidTrial.p;
                stateTrial.p2D = fluidTrial.p;
                stateTrial.pL2D = fluidTrial.pL;
                stateTrial.pE2D = fluidTrial.pE;
                % Keep the native 2D field in the accepted state so the next
                % adaptive step can compare the full 2D pressure field.
                if isfield(fluidTrial,'P') && ~isempty(fluidTrial.P)
                    stateTrial.P2DField = fluidTrial.P;
                end
                if isfield(fluidTrial,'urC') && ~isempty(fluidTrial.urC)
                    stateTrial.urC2DField = fluidTrial.urC;
                end
                if isfield(fluidTrial,'uzC') && ~isempty(fluidTrial.uzC)
                    stateTrial.uzC2DField = fluidTrial.uzC;
                end
                % Raw MAC face velocities (not just cell-centered), needed
                % by the optional unsteady-Stokes term (par.useUnsteadyStokes)
                % so the previous-step velocity source term lands on the
                % exact same DOF layout as the current step's momentum rows,
                % rather than going through an extra cell-center<->face
                % averaging step. Mesh topology (Nr, Nz) is fixed for the
                % whole run, so index (i,j) is the same LOGICAL face on both
                % steps (same fractional distance across the gap, same z) --
                % but the mesh is rebuilt every step from the current gap
                % geometry, so that logical face sits at a DIFFERENT PHYSICAL
                % r when the gap shape has changed. The physical face
                % positions (Rur/Ruz) are persisted here too so the consumer
                % (solve_fluid_2D_bodyfitted_MAC.m) can re-grid the old
                % velocity field onto the new mesh's physical r-positions
                % before using it as the unsteady term's previous-step value,
                % instead of assuming the two coincide.
                if isfield(fluidTrial,'ur') && ~isempty(fluidTrial.ur)
                    stateTrial.ur2DFaceField = fluidTrial.ur;
                end
                if isfield(fluidTrial,'uz') && ~isempty(fluidTrial.uz)
                    stateTrial.uz2DFaceField = fluidTrial.uz;
                end
                if isfield(par, 'useUnsteadyStokes') && par.useUnsteadyStokes
                    if isfield(fluidTrial,'meshF') && isfield(fluidTrial.meshF,'Rur')
                        stateTrial.Rur2DField = fluidTrial.meshF.Rur;
                    end
                    if isfield(fluidTrial,'meshF') && isfield(fluidTrial.meshF,'Ruz')
                        stateTrial.Ruz2DField = fluidTrial.meshF.Ruz;
                    end
                end
                if isfield(fluidTrial,'meshF') && isfield(fluidTrial.meshF,'Rp')
                    stateTrial.Rp2DField = fluidTrial.meshF.Rp;
                    stateTrial.Zp2DField = fluidTrial.meshF.Zp;
                end
                if isfield(fluidTrial,'tractionE'), stateTrial.tractionE2D = fluidTrial.tractionE; end
                if isfield(fluidTrial,'tractionL'), stateTrial.tractionL2D = fluidTrial.tractionL; end
            else
                stateTrial.p = fluidTrial.p;
                stateTrial.pReduced = fluidTrial.p;
            end
            stateTrial.tauE = fluidTrial.tauE;
            stateTrial.tauL = fluidTrial.tauL;

            [stateTrial, fluidTrial, pressureLimited, pressureLimitReason] = ...
                apply_pressure_temporal_limiter(stateTrial, fluidTrial, old, z, parStep);
            if pressureLimited && isfield(par, 'debugVerbose') && par.debugVerbose
                fprintf('   pressure limiter: %s\n', pressureLimitReason);
            end

            % Reject this time step if the newly solved pressure jumps too
            % much compared with the previously accepted pressure. The outer
            % try/catch block will reduce dtAttempt and recompute the step.
            check_pressure_jump_retry(old, fluidTrial, z, parStep, tn);

            % Optional event-detection refinement (default off, zero effect
            % unless explicitly enabled): par.gapStopFactor*par.minGap can
            % represent a physically-motivated "real contact" distance
            % (e.g. glycocalyx thickness) rather than just a numerical
            % floor, but the ordinary post-step stop check (below, after
            % this step is fully accepted) only fires AFTER a step
            % completes -- near the gap floor a single accepted step can
            % overshoot straight past the target (observed: 33nm -> 1.25nm
            % in one step), making the exact target value moot in practice.
            % This reuses the existing dt-halving retry mechanism to home
            % in on the target crossing instead of accepting the overshoot,
            % the same way it already retries on a hard solve failure.
            if isfield(par, 'refineToGapTarget') && par.refineToGapTarget && ...
                    isfield(par, 'gapStopFactor') && isfield(par, 'minGap')
                gapMinTrial = min(stateTrial.deltaE(:) - stateTrial.deltaL(:));
                gapMinOld = min(old.deltaE(:) - old.deltaL(:));
                gapTargetVal = par.gapStopFactor * par.minGap;
                tol = 0.1;
                if isfield(par, 'gapTargetTol'), tol = par.gapTargetTol; end
                if gapMinOld > gapTargetVal && gapMinTrial < (1-tol)*gapTargetVal && ...
                        dtAttempt > par.dtMin
                    dtNew = max(par.dtMin, par.dtRetryFactor * dtAttempt);
                    fprintf(['   event-detection: gap overshot target %.3e m ' ...
                        '(landed at %.3e m) -- retrying dt %.3e -> %.3e\n'], ...
                        gapTargetVal, gapMinTrial, dtAttempt, dtNew);
                    dtAttempt = dtNew;
                    retryCount = retryCount + 1;
                    acceptedStep = false;
                    continue;
                end
            end

            acceptedStep = true;
        catch ME
            stepReason = ['Monolithic solve failed: ', ME.message];

            if should_retry_time_step(stepReason, dtAttempt, retryCount, par)
                dtNew = max(par.dtMin, par.dtRetryFactor * dtAttempt);
                fprintf(['   retrying time step from t=%.6e s: ', ...
                    'dt %.3e -> %.3e after %s\n'], ...
                    tNow, dtAttempt, dtNew, compact_failure_reason(stepReason));
                dtAttempt = dtNew;
                retryCount = retryCount + 1;
                continue;
            end

            stoppedEarly = true;
            stopStep = tn;
            stopReason = stepReason;
            warning('Stopped early after step %d, attempted t = %.6e s. %s', ...
                tn, tNow + dtAttempt, stopReason);
            break;
        end
    end

    if stoppedEarly
        break;
    end

    tn = tn + 1;
    if tn > historyCapacity
        growBy = max(historyCapacity, max(nSteps, 1));
        newCapacity = historyCapacity + growBy;
        tHist(newCapacity,1) = nan;
        dtHist(newCapacity,1) = nan;
        stepWallTimeHist(newCapacity,1) = nan;
        retryHist(newCapacity,1) = 0;
        stateHist{newCapacity,1} = [];
        fluidHist{newCapacity,1} = [];
        deltaEHist(:,newCapacity) = 0;
        deltaLHist(:,newCapacity) = 0;
        pHist(:,newCapacity) = 0;
        tauEHist(:,newCapacity) = 0;
        tauLHist(:,newCapacity) = 0;
        uzEHist(:,newCapacity) = 0;
        uzLHist(:,newCapacity) = 0;
        PHist(:,:,newCapacity) = nan;
        urCHist(:,:,newCapacity) = nan;
        uzCHist(:,:,newCapacity) = nan;
        speedCHist(:,:,newCapacity) = nan;
        RPHist(:,:,newCapacity) = nan;
        ZPHist(:,:,newCapacity) = nan;
        p2DMaxHist(newCapacity,1) = nan;
        trEHist(:,:,newCapacity) = nan;
        trLHist(:,:,newCapacity) = nan;
        diagHist{newCapacity,1} = [];
        tractionCorrectionHistory{newCapacity,1} = [];
        historyCapacity = newCapacity;
    end

    state = stateTrial;
    fluid = fluidTrial;
    tNow = tNew;
    state.t = tNow;

    maxSolidChange = max(abs(state.uE - old.uE));
    if isfield(state, 'uL') && ~isempty(state.uL)
        maxSolidChange = max(maxSolidChange, max(abs(state.uL - old.uL)));
    end
    oldFluidPressure = old.p;
    if isfield(old, 'p2D') && isfield(par, 'useFull2DFluid') && par.useFull2DFluid
        oldFluidPressure = old.p2D;
    end
    maxPressureChange = max(abs(fluid.p - oldFluidPressure));
    printEvery = 1;
    if isfield(par, 'printEvery') && isfinite(par.printEvery) && par.printEvery > 0
        printEvery = par.printEvery;
    end
    doPrintStep = tn == 1 || tNow >= par.tEnd - timeTol || ...
        mod(tn, printEvery) == 0 || retryCount > 0;
    if doPrintStep
        fprintf(['Time step %d (t=%.4e, dt=%.3e): min(deltaE)=%.6e, ', ...
            'max|p|=%.6e, max|du|=%.3e, max|dp|=%.3e, retries=%d\n'], ...
            tn, tNow, dtAttempt, min(state.deltaE), max(abs(fluid.p)), ...
            maxSolidChange, maxPressureChange, retryCount);
    end

    % store full structs. Optionally strip heavy MAC arrays from the
    % cell history while preserving compact 2D arrays in PHist/urCHist/uzCHist.
    fluidStore = fluid;
    if isfield(par,'storeFull2DFluidHist') && ~par.storeFull2DFluidHist && ...
            isfield(fluidStore,'meshType') && strcmpi(fluidStore.meshType,'bodyfitted_MAC')
        % Keep compact cell-centered fields in fluidHist, but remove the
        % heavier native face-centered velocity arrays. The physical grid is
        % stored separately in RPHist/ZPHist.
        fluidStore.ur = [];
        fluidStore.uz = [];
        if isfield(fluidStore,'meshF') && isfield(fluidStore.meshF,'Rp')
            meshLight = struct();
            meshLight.Rp = fluidStore.meshF.Rp;
            meshLight.Zp = fluidStore.meshF.Zp;
            meshLight.zc = fluidStore.meshF.zc;
            meshLight.zF = fluidStore.meshF.zF;
            % add_fluid_nodes (called right below) needs Nr/Nz to index
            % Rp/Zp; keep them if present, otherwise derive them from the
            % Rp grid size so the light mesh is still self-contained.
            if isfield(fluidStore.meshF,'Nr') && isfield(fluidStore.meshF,'Nz')
                meshLight.Nr = fluidStore.meshF.Nr;
                meshLight.Nz = fluidStore.meshF.Nz;
            else
                [meshLight.Nr, meshLight.Nz] = size(fluidStore.meshF.Rp);
            end
            % These are small per-z-cell interface-position vectors (not
            % heavy face-centered arrays), but several post-run diagnostics
            % (check_interface_traction_continuity.m,
            % check_interface_stress_continuity.m,
            % check_wall_shear_approximation.m,
            % check_interface_traction_mismatch_report.m) require
            % fluid.meshF.deltaL_c/deltaE_c to locate the interface, so keep
            % them too.
            if isfield(fluidStore.meshF,'deltaL_c')
                meshLight.deltaL_c = fluidStore.meshF.deltaL_c;
            end
            if isfield(fluidStore.meshF,'deltaE_c')
                meshLight.deltaE_c = fluidStore.meshF.deltaE_c;
            end
            fluidStore.meshF = meshLight;
        end
    end

    fluidStore.meshF=add_fluid_nodes(fluidStore.meshF);
    %[pCell, sigmaCell, center] = recover_fluid_nodes_pressure_stress_Q4(fluidStore.meshF, fluidStore.ur2D, fluidStore.uz2D, par.mu);
    [pCell, sigmaCell, center] = recover_fluid_nodes_pressure_stress_Q4(fluidStore.meshF, fluidStore.ur2D, fluidStore.uz2D, par.mu, fluidStore.pCell);

    fluidStore.pCellNode=pCell;
    fluidStore.centerNode=center;
    fluidStore.sigmaCellNode=sigmaCell;
    stateHist{tn} = state;
    fluidHist{tn} = fluidStore;
    tHist(tn) = tNow;
    dtHist(tn) = dtAttempt;
    stepWallTimeHist(tn) = toc(stepTicId);
    retryHist(tn) = retryCount;
    deltaEHist(:,tn) = state.deltaE;
    deltaLHist(:,tn) = state.deltaL;
    pHist(:,tn)      = fluid.p;
    tauEHist(:,tn)   = fluid.tauE;
    tauLHist(:,tn)   = fluid.tauL;
    uzEHist(:,tn)    = fluid.uzE;
    uzLHist(:,tn)    = fluid.uzL;
    if isfield(fluid,'P') && ~isempty(fluid.P)
        [nrP,nzP] = size(fluid.P);
        nrS = min(size(PHist,1), nrP);
        nzS = min(size(PHist,2), nzP);
        PHist(1:nrS,1:nzS,tn) = fluid.P(1:nrS,1:nzS);
    end
    if isfield(fluid,'urC') && ~isempty(fluid.urC)
        [nrU,nzU] = size(fluid.urC);
        nrS = min(size(urCHist,1), nrU);
        nzS = min(size(urCHist,2), nzU);
        urCHist(1:nrS,1:nzS,tn) = fluid.urC(1:nrS,1:nzS);
    end
    if isfield(fluid,'uzC') && ~isempty(fluid.uzC)
        [nrU,nzU] = size(fluid.uzC);
        nrS = min(size(uzCHist,1), nrU);
        nzS = min(size(uzCHist,2), nzU);
        uzCHist(1:nrS,1:nzS,tn) = fluid.uzC(1:nrS,1:nzS);
    end
    if isfield(fluid,'urC') && isfield(fluid,'uzC') && ...
            ~isempty(fluid.urC) && ~isempty(fluid.uzC)
        [nrU,nzU] = size(fluid.urC);
        nrS = min(size(speedCHist,1), nrU);
        nzS = min(size(speedCHist,2), nzU);
        speedCHist(1:nrS,1:nzS,tn) = sqrt(fluid.urC(1:nrS,1:nzS).^2 + fluid.uzC(1:nrS,1:nzS).^2);
    end
    if isfield(fluid,'P') && ~isempty(fluid.P)
        Pvals = abs(fluid.P(:));
        Pvals = Pvals(isfinite(Pvals));
        if ~isempty(Pvals), p2DMaxHist(tn) = max(Pvals); end
    end
    if isfield(fluid,'meshF') && isfield(fluid.meshF,'Rp') && ~isempty(fluid.meshF.Rp)
        [nrR,nzR] = size(fluid.meshF.Rp);
        nrS = min(size(RPHist,1), nrR);
        nzS = min(size(RPHist,2), nzR);
        RPHist(1:nrS,1:nzS,tn) = fluid.meshF.Rp(1:nrS,1:nzS);
    end
    if isfield(fluid,'meshF') && isfield(fluid.meshF,'Zp') && ~isempty(fluid.meshF.Zp)
        [nrZ,nzZ] = size(fluid.meshF.Zp);
        nrS = min(size(ZPHist,1), nrZ);
        nzS = min(size(ZPHist,2), nzZ);
        ZPHist(1:nrS,1:nzS,tn) = fluid.meshF.Zp(1:nrS,1:nzS);
    end
    if isfield(fluid,'tractionE')
        trEHist(1,:,tn) = safe_interp1_same_or_resample(fluid.tractionE.z(:), fluid.tractionE.normal(:), z(:), 'fluid.tractionE.normal');
        trEHist(2,:,tn) = safe_interp1_same_or_resample(fluid.tractionE.z(:), fluid.tractionE.tangent(:), z(:), 'fluid.tractionE.tangent');
    end
    if isfield(fluid,'tractionL')
        trLHist(1,:,tn) = safe_interp1_same_or_resample(fluid.tractionL.z(:), fluid.tractionL.normal(:), z(:), 'fluid.tractionL.normal');
        trLHist(2,:,tn) = safe_interp1_same_or_resample(fluid.tractionL.z(:), fluid.tractionL.tangent(:), z(:), 'fluid.tractionL.tangent');
    end

    diag = compute_step_diagnostics(z, old, state, fluid, parStep, ...
        pressureLimited, pressureLimitReason);
    diag.dt = dtAttempt;
    diag.retries = retryCount;
    diagHist{tn} = diag;

    % Periodic checkpoint (Sep 11): saves a partial 'out' during the run
    % itself, independent of the normal end-of-run save. Without this, a
    % run that hits the wall-time limit or crashes in post-solve plotting
    % loses everything, even if hundreds of good steps were completed --
    % this happened to two full cluster runs (12h and 9h) before this fix.
    % Wrapped in try/catch: a checkpoint failure (e.g. disk full) must
    % never abort the simulation itself.
    if isfield(par, 'checkpointFile') && ~isempty(par.checkpointFile) && ...
            isfield(par, 'checkpointEvery') && par.checkpointEvery > 0 && ...
            mod(tn, par.checkpointEvery) == 0
        try
            out = struct();
            out.z = z;
            out.t = tHist(1:tn);
            out.dtHist = dtHist(1:tn);
            out.stepWallTimeHist = stepWallTimeHist(1:tn);
            out.retryHist = retryHist(1:tn);
            if exist('t0State', 'var')
                out.t0State = t0State;
                out.t0Fluid = t0Fluid;
            end
            out.state = state;
            out.stateHist = stateHist(1:tn);
            out.fluidHist = fluidHist(1:tn);
            out.deltaEHist = deltaEHist(:,1:tn);
            out.deltaLHist = deltaLHist(:,1:tn);
            out.pHist = pHist(:,1:tn);
            out.tauEHist = tauEHist(:,1:tn);
            out.tauLHist = tauLHist(:,1:tn);
            out.uzEHist = uzEHist(:,1:tn);
            out.uzLHist = uzLHist(:,1:tn);
            out.PHist = PHist(:,:,1:tn);
            out.urCHist = urCHist(:,:,1:tn);
            out.uzCHist = uzCHist(:,:,1:tn);
            out.speedCHist = speedCHist(:,:,1:tn);
            out.RPHist = RPHist(:,:,1:tn);
            out.ZPHist = ZPHist(:,:,1:tn);
            out.p2DMaxHist = p2DMaxHist(1:tn);
            out.trEHist = trEHist(:,:,1:tn);
            out.trLHist = trLHist(:,:,1:tn);
            out.diagHist = diagHist(1:tn);
            out.tractionCorrectionHistory = tractionCorrectionHistory(1:tn);
            out.par = par;
            out.meshE = meshE;
            out.interfaceE = interfaceE;
            if exist('meshL','var') && ~isempty(meshL)
                out.meshL = meshL;
            end
            if exist('interfaceL','var') && ~isempty(interfaceL)
                out.interfaceL = interfaceL;
            end
            out.stoppedEarly = false;
            out.stopStep = tn;
            out.stopReason = '';
            out.isCheckpoint = true;  % marks this as a partial, in-progress save
            save(par.checkpointFile, 'out', '-v7.3');
            clear out;
        catch MEcp
            warning('Checkpoint save failed at step %d (continuing run): %s', tn, MEcp.message);
        end
    end

    gapMinNow = diag.gapMin;
    if doPrintStep
        fprintf('   min gap after accepted step = %.6e m\n', gapMinNow);
        if isfield(par, 'diagnosticsEnabled') && par.diagnosticsEnabled
            print_step_diagnostics(diag);
        end
    end

    warn_step_diagnostics(diag, par);

    if isfield(par, 'stopAtMinGap') && par.stopAtMinGap && ...
            gapMinNow <= par.gapStopFactor * parStep.minGap
        stoppedEarly = true;
        stopStep = tn;
        stopReason = sprintf('Reached minimum gap: hmin = %.6e m', gapMinNow);
        fprintf('Reached minimum gap at step %d, t = %.6e s. hmin = %.6e m\n', ...
            tn, tNew, gapMinNow);
        break;
    end

    if isfield(par, 'enableAdaptiveTimeStep') && par.enableAdaptiveTimeStep && ...
            (retryCount > 0 || dtAttempt < par.dt)
        dtNext = min(par.dt, par.dtGrowFactor * dtAttempt);
    else
        dtNext = par.dt;
    end
end

if ~stoppedEarly
    stopStep = tn;
end

stateHist = stateHist(1:stopStep);
fluidHist = fluidHist(1:stopStep);
deltaEHist = deltaEHist(:,1:stopStep);
deltaLHist = deltaLHist(:,1:stopStep);
pHist      = pHist(:,1:stopStep);
tauEHist   = tauEHist(:,1:stopStep);
tauLHist   = tauLHist(:,1:stopStep);
uzEHist    = uzEHist(:,1:stopStep);
uzLHist    = uzLHist(:,1:stopStep);
PHist      = PHist(:,:,1:stopStep);
urCHist    = urCHist(:,:,1:stopStep);
uzCHist    = uzCHist(:,:,1:stopStep);
speedCHist = speedCHist(:,:,1:stopStep);
RPHist     = RPHist(:,:,1:stopStep);
ZPHist     = ZPHist(:,:,1:stopStep);
p2DMaxHist = p2DMaxHist(1:stopStep);
trEHist    = trEHist(:,:,1:stopStep);
trLHist    = trLHist(:,:,1:stopStep);
diagHist   = diagHist(1:stopStep);
tractionCorrectionHistory = tractionCorrectionHistory(1:stopStep);
tHist      = tHist(1:stopStep);
dtHist     = dtHist(1:stopStep);
stepWallTimeHist = stepWallTimeHist(1:stopStep);
retryHist  = retryHist(1:stopStep);
adaptiveSummary = summarize_time_step_adaptation( ...
    dtHist, retryHist, par, stoppedEarly, stopReason);
if adaptiveSummary.totalRetries > 0 || stoppedEarly
    print_adaptive_summary(adaptiveSummary);
end

out = struct();
out.z = z;
out.t = tHist;
out.dtHist = dtHist;
out.stepWallTimeHist = stepWallTimeHist;  % wall-clock seconds per accepted step (incl. retries)
out.retryHist = retryHist;
if exist('t0State', 'var')
    out.t0State = t0State;  % initial state before any time-stepping (nargin==1 entry point only)
    out.t0Fluid = t0Fluid;  % static fluid field at t=0 (zero wall velocity); empty if evaluation failed
end
out.adaptiveSummary = adaptiveSummary;
out.state = state;              % final state
out.stateHist = stateHist;      % all time steps
out.fluidHist = fluidHist;      % all time steps
out.deltaEHist = deltaEHist;
out.deltaLHist = deltaLHist;
out.pHist = pHist;
out.tauEHist = tauEHist;
out.tauLHist = tauLHist;
out.uzEHist = uzEHist;
out.uzLHist = uzLHist;
out.PHist = PHist;          % native 2D pressure field, size Nr x Nz x Nt
out.urCHist = urCHist;      % centered radial velocity field, size Nr x Nz x Nt
out.uzCHist = uzCHist;      % centered axial velocity field, size Nr x Nz x Nt
out.speedCHist = speedCHist;
out.RPHist = RPHist;        % physical r grid for PHist/urCHist/uzCHist
out.ZPHist = ZPHist;        % physical z grid for PHist/urCHist/uzCHist
out.p2DMaxHist = p2DMaxHist;
out.trEHist = trEHist;
out.trLHist = trLHist;
out.diagHist = diagHist;
out.tractionCorrectionHistory = tractionCorrectionHistory; % per-step convergeInfo from apply_bodyfitted_MAC_traction_correction_feedback, empty cell for steps that didn't use it
out.par = par;
out.meshE = meshE;
out.interfaceE = interfaceE;
if exist('meshL','var') && ~isempty(meshL)
    out.meshL = meshL;
end
if exist('interfaceL','var') && ~isempty(interfaceL)
    out.interfaceL = interfaceL;
end

out.stoppedEarly = stoppedEarly;
out.stopStep = stopStep;
out.stopReason = stopReason;
if isfield(par, 'useGlobal1DPressure') && par.useGlobal1DPressure
    out.global1D = build_global_1d_pressure_view(out.z, out.pHist, par);
    if out.stopStep >= 1
        out.global1D = add_global_1d_blank_solid_pressure_view(out, par);
    end
end
if out.stopStep >= 1 && use_global2d_pressure_traction(par)
    out.global2DPressureTractionComparison = ...
        build_global2d_pressure_traction_comparison(out);
end

% choose time step to plot
if out.stopStep < 1
    warning('Simulation failed before completing the first time step. No plots generated.');
    return;
end

nPlot = out.stopStep;
statePlot = out.stateHist{nPlot};
fluidPlot = out.fluidHist{nPlot};
if par.saveOutput
    outputFile = 'simulation_output_two_solid_full_analytical_nopre.mat';
    if isfield(par, 'outputFile') && ~isempty(par.outputFile)
        outputFile = par.outputFile;
    end
    save(outputFile,'out','-v7.3');
end

plotFinalGlobalPressureContour = isfield(par, 'useGlobal1DPressure') && ...
    par.useGlobal1DPressure && ...
    (par.makePlots || (isfield(par, 'plotFinalGlobalPressureContour') && ...
    par.plotFinalGlobalPressureContour));
if plotFinalGlobalPressureContour
    plot_global_1d_pressure_contour_with_blank_solids(out, par);
end

plotGlobal2DPressureTractionComparison = use_global2d_pressure_traction(par) && ...
    isfield(out, 'global2DPressureTractionComparison') && ...
    out.global2DPressureTractionComparison.available && ...
    (par.makePlots || (isfield(par, 'plotGlobal2DPressureTractionComparison') && ...
    par.plotGlobal2DPressureTractionComparison));
if plotGlobal2DPressureTractionComparison
    plot_global2d_pressure_traction_comparison(out, par);
end

if ~par.makePlots
    return;
end

[R, Z, Uz] = velocity_field_for_plot(z, fluidPlot, statePlot, par);
dz = z(2)-z(1);

% -------------------------------
% Pressure field
% -------------------------------
figure;
set(gca, 'FontSize', 24);
plot(z*1e6, fluidPlot.p, 'LineWidth', 1.8);
grid off;
xlabel('z [\mum]');
ylabel('Pressure [Pa]');
title(sprintf('Pressure field at t = %.4f s', out.t(nPlot)));

figure;
plot(out.t, out.pHist(round(par.NzFluid/2),:), 'LineWidth', 2);
set(gca, 'FontSize', 24);
xlabel('t [s]');
ylabel('Pressure at mid-point [Pa]');
title(sprintf('Pressure evolution at z = %.3f \\mum', out.z(round(par.NzFluid/2))*1e6));
grid off;



% -------------------------------
% Native 2D pressure and velocity fields on the body-fitted MAC grid
% -------------------------------
if isfield(out,'PHist') && ~isempty(out.PHist)
    Pplot = out.PHist(:,:,nPlot);
    if any(isfinite(Pplot(:)))
        figure;
        set(gca, 'FontSize', 24);
        contourf(out.RPHist(:,:,nPlot)*1e6, out.ZPHist(:,:,nPlot)*1e6, ...
            Pplot, 40, 'LineColor', 'none');
        colorbar;
        hold on;
        plot(statePlot.deltaL*1e6, z*1e6, 'k-', 'LineWidth', 1.2);
        plot(statePlot.deltaE*1e6, z*1e6, 'k-', 'LineWidth', 1.2);
        xlabel('r [\mum]');
        ylabel('z [\mum]');
        title(sprintf('Native 2D pressure P(r,z) at t = %.4f s', out.t(nPlot)));
    end
end

if isfield(out,'uzCHist') && ~isempty(out.uzCHist)
    UzCplot = out.uzCHist(:,:,nPlot);
    if any(isfinite(UzCplot(:)))
        figure;
        set(gca, 'FontSize', 24);
        contourf(out.RPHist(:,:,nPlot)*1e6, out.ZPHist(:,:,nPlot)*1e6, ...
            UzCplot*1e6, 40, 'LineColor', 'none');
        colorbar;
        hold on;
        plot(statePlot.deltaL*1e6, z*1e6, 'k-', 'LineWidth', 1.2);
        plot(statePlot.deltaE*1e6, z*1e6, 'k-', 'LineWidth', 1.2);
        xlabel('r [\mum]');
        ylabel('z [\mum]');
        title(sprintf('Native 2D axial velocity u_z(r,z) at t = %.4f s', out.t(nPlot)));
    end
end

if isfield(out,'speedCHist') && ~isempty(out.speedCHist)
    speedPlot = out.speedCHist(:,:,nPlot);
    if any(isfinite(speedPlot(:)))
        figure;
        set(gca, 'FontSize', 24);
        contourf(out.RPHist(:,:,nPlot)*1e6, out.ZPHist(:,:,nPlot)*1e6, ...
            speedPlot*1e6, 40, 'LineColor', 'none');
        colorbar;
        hold on;
        plot(statePlot.deltaL*1e6, z*1e6, 'k-', 'LineWidth', 1.2);
        plot(statePlot.deltaE*1e6, z*1e6, 'k-', 'LineWidth', 1.2);
        xlabel('r [\mum]');
        ylabel('z [\mum]');
        title(sprintf('Native 2D speed |u|(r,z) at t = %.4f s', out.t(nPlot)));
    end
end

% -------------------------------
% Shear stress on both interfaces
% -------------------------------
figure;
set(gca, 'FontSize', 24);
plot(z*1e6, fluidPlot.tauL, 'LineWidth', 1.8); hold on;
plot(z*1e6, -fluidPlot.tauE, 'LineWidth', 1.8);
grid off;
xlabel('z [\mum]');
ylabel('Shear stress \tau_{rz} [Pa]');
legend('Inner wall', 'Outer wall', 'Location', 'best');
title(sprintf('Wall shear stress at t = %.4f s', out.t(nPlot)));


% % -------------------------------
% % Interfacial velocities
% % -------------------------------
% figure;
% set(gca, 'FontSize', 24);
% plot(z*1e6, fluidPlot.uzL, 'LineWidth', 1.8); hold on;
% plot(z*1e6, fluidPlot.uzE, 'LineWidth', 1.8);
% grid off;
% xlabel('z [\mum]');
% ylabel('u_z at wall [m/s]');
% legend('u_z at inner wall', 'u_z at outer wall', 'Location', 'best');
% title(sprintf('Interfacial axial velocities at t = %.4f s', out.t(nPlot)));

% -------------------------------
% Velocity field contour
% -------------------------------
figure;
set(gca, 'FontSize', 24);
%contourf(R*1e6, Z*1e6, Uz, 30, 'LineColor', 'none');
contourf(R*1e6, Z*1e6, Uz*1e6, 40, 'LineColor', 'none');
colorbar;
hold on;
plot( statePlot.deltaL*1e6, z*1e6, 'k-', 'LineWidth', 1.2);
plot( statePlot.deltaE*1e6, z*1e6,'k-', 'LineWidth', 1.2);
xlabel('r [\mum]');
ylabel('z [\mum]');
title(sprintf('Axial velocity field at t = %.4f s', out.t(nPlot)));
%xlim([par.RLout*1e6 3])

% -------------------------------
% Final deformed solid meshes
% -------------------------------
figure;
hold on;
set(gca, 'FontSize', 24);
if ~(isfield(par, 'noLeukocyte') && par.noLeukocyte) && ...
        ~isempty(meshL) && isfield(statePlot,'uL') && ~isempty(statePlot.uL)
    hLmesh = plot_deformed_mesh(meshL, statePlot.uL, [0.10 0.55 0.25], 0.45);
    [rLDef, zLDef] = deformed_interface_curve(meshL, statePlot.uL, interfaceL);
    plot(rLDef*1e6, zLDef*1e6, 'Color', [0.00 0.35 0.10], 'LineWidth', 2.0);
else
    hLmesh = gobjects(0);
    plot(par.RLout*1e6*ones(size(z)), z*1e6, 'Color', [0.00 0.35 0.10], 'LineWidth', 2.0);
end
hEmesh = plot_deformed_mesh(meshE, statePlot.uE, [0.10 0.35 0.85], 0.45);
[rEDef, zEDef] = deformed_interface_curve(meshE, statePlot.uE, interfaceE);
plot(rEDef*1e6, zEDef*1e6, 'Color', [0.00 0.15 0.65], 'LineWidth', 2.0);
axis equal tight;
grid off;
xlabel('r [\mum]');
ylabel('deformed z [\mum]');
title(sprintf('Final deformed solid meshes at t = %.4f s', out.t(nPlot)));
if isempty(hLmesh)
    legend(hEmesh, 'Endothelium mesh', 'Location', 'best');
else
    legend([hLmesh, hEmesh], {'Leukocyte mesh', 'Endothelium mesh'}, 'Location', 'best');
end

% -------------------------------
% Mid-plane velocity profile
% -------------------------------
jmid = round(numel(z)/2);
figure;
set(gca, 'FontSize', 24);
plot(R(:,jmid)*1e6, Uz(:,jmid)*1e6, 'LineWidth', 2);
grid off;
ylabel('u_z [\mum/s]');
xlabel('r [\mum]');
title(sprintf('Velocity profile at z = %.3f \\mum, t = %.4f s', z(jmid)*1e6, out.t(nPlot)));
xlim([min(R(:,jmid)), max(R(:,jmid))] * 1e6)


figure;
hold on
set(gca,'FontSize',24)
stepPlot = 1;   % plot every n time steps
for n = 1:stepPlot:out.stopStep
    state_n = out.stateHist{n};
    if isfield(out.par, 'useExactDeformedInterface') && out.par.useExactDeformedInterface
        [rDef_n, zDef_n] = deformed_interface_curve(meshE, state_n.uE, interfaceE);
        plot(rDef_n*1e6, zDef_n*1e6, 'LineWidth', 1.2);
    else
        zDef_n = out.z + extract_interface_axial_displacement( ...
            meshE, state_n.uE, interfaceE, out.z);

        plot(state_n.deltaE*1e6, zDef_n*1e6, 'LineWidth', 1.2);
    end
end
xlabel('r [\mum]');
ylabel('deformed z [\mum]');
title('Endothelium shape evolution with axial drag');
grid off

figure;
hold on
set(gca,'FontSize',24)
for n = 1:stepPlot:out.stopStep
    state_n = out.stateHist{n};
    uz_n = extract_interface_axial_displacement( ...
        meshE, state_n.uE, interfaceE, out.z);

    plot(out.z*1e6, uz_n*1e6, 'LineWidth', 1.2);
end
xlabel('reference z [\mum]');
ylabel('interface u_z [\mum]');
title('Endothelium axial displacement evolution');
grid off

% -------------------------------
% Additional physical diagnostics
% -------------------------------
validSteps = 1:out.stopStep;
gapMinHist = nan(out.stopStep,1);
gapMaxHist = nan(out.stopStep,1);
volRelHist = nan(out.stopStep,1);
maxPTime   = nan(out.stopStep,1);
maxDuTime  = nan(out.stopStep,1);

for n = validSteps
    st = out.stateHist{n};
    gap_n = st.deltaE(:) - st.deltaL(:);
    gapMinHist(n) = min(gap_n);
    gapMaxHist(n) = max(gap_n);
    maxPTime(n) = max(abs(out.pHist(:,n)));

    duNow = max(abs(st.uE(:) - st.uEPrev(:)));
    if isfield(st, 'uL') && ~isempty(st.uL)
        duNow = max(duNow, max(abs(st.uL(:) - st.uLPrev(:))));
    end
    maxDuTime(n) = duNow;

    if ~isempty(out.diagHist{n}) && isfield(out.diagHist{n}, 'volumeChangeRel')
        volRelHist(n) = out.diagHist{n}.volumeChangeRel;
    end
end

figure;
set(gca,'FontSize',24)
plot(out.t*1e3, gapMinHist*1e6, 'LineWidth', 2); hold on;
plot(out.t*1e3, gapMaxHist*1e6, '--', 'LineWidth', 1.5);
grid off
xlabel('t [ms]');
ylabel('gap [\mum]');
legend('minimum gap', 'maximum gap', 'Location', 'best');
title('Gap evolution');

figure;
set(gca,'FontSize',24)
plot(out.t*1e3, maxPTime, 'LineWidth', 2);
grid off
xlabel('t [ms]');
ylabel('max |p| [Pa]');
title('Maximum axial pressure over time');

if isfield(out,'p2DMaxHist') && any(isfinite(out.p2DMaxHist))
    figure;
    set(gca,'FontSize',24)
    plot(out.t*1e3, out.p2DMaxHist, 'LineWidth', 2);
    grid off
    xlabel('t [ms]');
    ylabel('max |P(r,z)| [Pa]');
    title('Maximum native 2D pressure over time');
end

figure;
set(gca,'FontSize',24)
plot(out.t*1e3, maxDuTime*1e9, 'LineWidth', 2);
grid off
xlabel('t [ms]');
ylabel('max |\Delta u| per step [nm]');
title('Maximum solid displacement increment per step');

figure;
set(gca,'FontSize',24)
plot(out.t*1e3, volRelHist, 'LineWidth', 2);
grid off
xlabel('t [ms]');
ylabel('\Delta V/V per step');
title('Relative fluid-volume change per step');

figure;
set(gca,'FontSize',24)
plot(z*1e6, (statePlot.deltaE - statePlot.deltaL)*1e6, 'LineWidth', 2);
grid off
xlabel('z [\mum]');
ylabel('gap h [\mum]');
title(sprintf('Final gap profile at t = %.4f s', out.t(nPlot)));

figure;
set(gca,'FontSize',24)
plot(statePlot.deltaL*1e6, z*1e6, 'LineWidth', 2); hold on;
plot(statePlot.deltaE*1e6, z*1e6, 'LineWidth', 2);
grid off
xlabel('r [\mum]');
ylabel('fluid-grid z [\mum]');
legend('leukocyte interface', 'endothelium interface', 'Location', 'best');
title(sprintf('Fluid interfaces at t = %.4f s', out.t(nPlot)));

if ~isempty(meshL) && isfield(statePlot,'uL') && ~isempty(statePlot.uL)
    figure;
    hold on
    set(gca,'FontSize',24)
    stepPlot = max(1, ceil(out.stopStep/20));
    for n = 1:stepPlot:out.stopStep
        state_n = out.stateHist{n};
        [rLDef_n, zLDef_n] = deformed_interface_curve(meshL, state_n.uL, interfaceL);
        plot(rLDef_n*1e6, zLDef_n*1e6, 'LineWidth', 1.2);
    end
    xlabel('r [\mum]');
    ylabel('deformed z [\mum]');
    title('Leukocyte interface shape evolution');
    grid off

    figure;
    hold on
    set(gca,'FontSize',24)
    stepPlot = max(1, ceil(out.stopStep/20));
    for n = 1:stepPlot:out.stopStep
        state_n = out.stateHist{n};
        uzE_n = extract_interface_axial_displacement(meshE, state_n.uE, interfaceE, out.z);
        uzL_n = extract_interface_axial_displacement(meshL, state_n.uL, interfaceL, out.z);
        plot(out.z*1e6, uzE_n*1e6, 'LineWidth', 1.2);
        plot(out.z*1e6, uzL_n*1e6, '--', 'LineWidth', 1.2);
    end
    xlabel('reference/fluid z [\mum]');
    ylabel('interface u_z [\mum]');
    title('Endothelium and leukocyte axial displacement evolution');
    grid off
end



end

function [R, Z, Uz] = velocity_field_for_plot(z, fluid, state, par)
if isfield(par, 'useFull2DFluid') && par.useFull2DFluid && ...
        isfield(fluid, 'meshF') && ~isempty(fluid.meshF)

    % New body-fitted MAC Stokes solver format.
    if isfield(fluid, 'meshType') && strcmpi(fluid.meshType, 'bodyfitted_MAC') && ...
            isfield(fluid.meshF, 'Rp') && isfield(fluid.meshF, 'Zp') && ...
            isfield(fluid, 'uzC') && ~isempty(fluid.uzC)
        R = fluid.meshF.Rp;
        Z = fluid.meshF.Zp;
        Uz = fluid.uzC;
        return;
    end

    % Legacy Q4 penalty-Stokes format.
    if isfield(fluid.meshF, 'nodes') && isfield(fluid, 'uz2D') && ~isempty(fluid.uz2D)
        Nr = fluid.meshF.Nr;
        Nz = fluid.meshF.Nz;
        R = reshape(fluid.meshF.nodes(:,1), Nr, Nz);
        Z = reshape(fluid.meshF.nodes(:,2), Nr, Nz);
        Uz = reshape(fluid.uz2D(:), Nr, Nz);
        return;
    end

    R = fluid.meshF.Rp;
    Z = fluid.meshF.Zp;
    Uz = fluid.uzC;
    return;
end

[R, Z, Uz] = build_velocity_field(z, fluid.p, state.deltaL, ...
    state.deltaE, state.UwL, state.UwE, par);
end

function print_step_diagnostics(diag)
fprintf(['   diagnostics: fluidRes=%.3e, globalMass=%.3e, ', ...
    'dV/V=%.3e, fluxJump=%.3e, sourceInt=%.3e\n'], ...
    diag.fluidResidualInf, diag.globalMassResidual, ...
    diag.volumeChangeRel, diag.fluxJump, diag.sourceIntegral);
if isfield(diag, 'minJE')
    fprintf('   geometry: JEmin=%.3e, rEmin=%.3e m', ...
        diag.minJE, diag.minRadiusE);
    if isfield(diag, 'minJL') && isfinite(diag.minJL)
        fprintf(', JLmin=%.3e, rLmin=%.3e m', ...
            diag.minJL, diag.minRadiusL);
    end
    fprintf('\n');
end
end

function tf = should_retry_time_step(reason, dtAttempt, retryCount, par)
tf = false;
if ~isfield(par, 'enableAdaptiveTimeStep') || ~par.enableAdaptiveTimeStep
    return;
end

maxRetries = 6;
if isfield(par, 'maxTimeStepRetries') && isfinite(par.maxTimeStepRetries)
    maxRetries = par.maxTimeStepRetries;
end
if retryCount >= maxRetries
    return;
end

dtMin = 0;
if isfield(par, 'dtMin') && isfinite(par.dtMin)
    dtMin = par.dtMin;
end
if dtAttempt <= dtMin * (1 + 10*eps)
    return;
end

retryTokens = { ...
    'Negative or zero J', ...
    'Element inverted', ...
    'Non-positive radius', ...
    'Solid geometry guard failed', ...
    'Gap violates minGap', ...
    'violates minGap', ...
    'Fluid solve failed', ...
    'Pressure jump too large', ...
    'fsolve failed', ...
    'did not reach equilibrium', ...
    'residual too large', ...
    'line search failed', ...
    'did not converge'};

for k = 1:numel(retryTokens)
    if contains(reason, retryTokens{k})
        tf = true;
        return;
    end
end
end

function summary = summarize_time_step_adaptation(dtHist, retryHist, par, stoppedEarly, stopReason)
valid = isfinite(dtHist) & dtHist > 0;
summary = struct();
summary.enabled = isfield(par, 'enableAdaptiveTimeStep') && par.enableAdaptiveTimeStep;
summary.acceptedSteps = nnz(valid);
summary.totalRetries = sum(retryHist(valid));
summary.retrySteps = nnz(retryHist(valid) > 0);
summary.stoppedEarly = stoppedEarly;
summary.stopReason = stopReason;

if any(valid)
    summary.minDt = min(dtHist(valid));
    summary.maxDt = max(dtHist(valid));
    summary.finalDt = dtHist(find(valid, 1, 'last'));
    summary.maxRetriesInStep = max(retryHist(valid));
else
    summary.minDt = NaN;
    summary.maxDt = NaN;
    summary.finalDt = NaN;
    summary.maxRetriesInStep = 0;
end
end

function print_adaptive_summary(summary)
fprintf(['Adaptive stepping summary: accepted=%d, retrySteps=%d, ', ...
    'totalRetries=%d, minDt=%.3e, finalDt=%.3e\n'], ...
    summary.acceptedSteps, summary.retrySteps, summary.totalRetries, ...
    summary.minDt, summary.finalDt);
if summary.stoppedEarly
    fprintf('   stopped early: %s\n', compact_failure_reason(summary.stopReason));
end
end

function s = compact_failure_reason(reason)
s = regexprep(char(reason), '\s+', ' ');
maxChars = 180;
if numel(s) > maxChars
    s = [s(1:maxChars), '...'];
end
end

function state = attach_state_geometry_checks(state, meshE, meshL, par)
qE = solid_geometry_quality(meshE, state.uE, 'endothelium');
assert_solid_geometry_ok(qE, par);
state.geometryE = qE;

if ~isempty(meshL) && isfield(state, 'uL') && ~isempty(state.uL)
    qL = solid_geometry_quality(meshL, state.uL, 'leukocyte');
    assert_solid_geometry_ok(qL, par);
    state.geometryL = qL;
end
end

function warn_step_diagnostics(diag, par)
if ~isfield(par, 'diagnosticsEnabled') || ~par.diagnosticsEnabled
    return;
end

fluidTol = inf;
if isfield(par, 'diagnosticsWarnFluidResidual')
    fluidTol = par.diagnosticsWarnFluidResidual;
end
if isfinite(fluidTol) && diag.fluidResidualInf > fluidTol
    warning('Fluid diagnostic residual %.3e exceeds %.3e.', ...
        diag.fluidResidualInf, fluidTol);
end

volumeTol = inf;
if isfield(par, 'diagnosticsWarnVolumeJumpRel')
    volumeTol = par.diagnosticsWarnVolumeJumpRel;
end
if isfinite(volumeTol) && diag.volumeChangeRel > volumeTol
    warning('Relative fluid-volume change %.3e exceeds %.3e.', ...
        diag.volumeChangeRel, volumeTol);
end
end

function parL = leukocyte_solid_parameters(par)
parL = par;

% Bug fix: par.solidTrustU0/solidTrustUMax (2e-8 m / 2e-7 m)
% are set once, globally, on par -- tuned via extensive testing
% to work well for the ENDOTHELIUM, whose smallest mesh feature
% is order ~3 microns. Reused verbatim for the leukocyte, a 20 nm
% trust-region step becomes comparable to the leukocyte's OWN geometry
% once it has compressed significantly during the simulation (its
% deformed minimum radius reaches ~30 nm at the point observed, ~100x
% smaller than the endothelium) -- confirmed directly via a live
% production trace showing every single Newton trial step rejected by
% the line search, trustU collapsing geometrically to its floor within
% ~20 iterations, relNorm frozen at ~72% (a genuine, not
% false-convergence, basin-of-attraction failure specific to the
% leukocyte side of the body-fitted MAC traction correction).
%
% An earlier version of this fix rescaled against the REFERENCE
% (undeformed) mesh's minimum radius here, at one-time case-setup -- but
% the reference geometry's minimum radius (~4 um) is nowhere near as
% small as the DEFORMED radius the leukocyte actually reaches after
% compression, so that static rescale was a no-op in practice. The
% correction loop calls solve_finite_def_solid with the CURRENT deformed
% state, so the trust region must be rescaled dynamically there, against
% the current deformed geometry, not once here against the reference
% mesh. See apply_bodyfitted_MAC_traction_correction.m /
% apply_bodyfitted_MAC_traction_correction_feedback.m for the actual
% fix; this function only carries the material-property remapping.

if isfield(par, 'EL')
    parL.Ee = par.EL;
end
if isfield(par, 'nuL')
    parL.nuE = par.nuL;
end
if isfield(par, 'GL')
    parL.Ge = par.GL;
elseif isfield(parL, 'Ee') && isfield(parL, 'nuE')
    parL.Ge = parL.Ee/(2*(1+parL.nuE));
end
if isfield(par, 'KL')
    parL.Ke = par.KL;
elseif isfield(parL, 'Ee') && isfield(parL, 'nuE')
    parL.Ke = parL.Ee/(3*(1-2*parL.nuE));
end
if isfield(par, 'etaL')
    parL.etaE = par.etaL;
end
if isfield(par, 'etaBulkL')
    parL.etaBulkE = par.etaBulkL;
end
if isfield(par, 'useViscoelasticLeukocyte')
    parL.useViscoelasticEndothelium = par.useViscoelasticLeukocyte;
end
end

function warn_leukocyte_prestress_load_mismatch(SL, par)
if isfield(par, 'warnPrestressLoadMismatch') && ~par.warnPrestressLoadMismatch
    return;
end

[isMismatch, prestressLoad, runtimeLoad] = leukocyte_prestress_load_mismatch(SL, par);
if ~isMismatch
    return;
end

warning(['Leukocyte prestress load scale is %.3g Pa, while runtime ', ...
    'pIn/pOut are %.3g/%.3g Pa. Make sure the coupled initial ', ...
    'fluid/solid load is intentional, otherwise the first step may ', ...
    'mostly relax a prestress mismatch.'], ...
    prestressLoad, par.pIn, par.pOut);
end

function [isMismatch, prestressLoad, runtimeLoad] = leukocyte_prestress_load_mismatch(SL, par)
prestressLoad = 0;
if isfield(SL, 'P0') && isnumeric(SL.P0)
    p0Vals = SL.P0(:);
    p0Vals = p0Vals(isfinite(p0Vals));
    if ~isempty(p0Vals)
        prestressLoad = max(prestressLoad, max(abs(p0Vals)));
    end
end
if isfield(SL, 'trL') && isfield(SL.trL, 'normal') && isnumeric(SL.trL.normal)
    normalVals = SL.trL.normal(:);
    normalVals = normalVals(isfinite(normalVals));
    if ~isempty(normalVals)
        prestressLoad = max(prestressLoad, max(abs(normalVals)));
    end
end

runtimeLoad = max(abs([par.pIn, par.pOut]));
isMismatch = prestressLoad > max(10 * runtimeLoad, 1e-9);
end

function state = initial_state(z, par, meshE, uE_pre, deltaE_pre, ...
    meshL, interfaceL, uL_pre, deltaL_pre)
state = struct();
state.UwL = zeros(size(z));
state.UwE = zeros(size(z));

if isfield(par, 'noLeukocyte') && par.noLeukocyte
    state.deltaL = zeros(size(z));
    state.uL = [];
    state.uLPrev = [];
else
    state.uL = uL_pre(:);
    state.uLPrev = state.uL;
    if use_RLout_fluid_interface_for_solid_leukocyte(par)
        state.deltaL = par.RLout * ones(size(z));
        state.UwL = zeros(size(z));
    else
        state.deltaL = deltaL_pre(:);
    end
end
state.deltaE = deltaE_pre;
state.p = linspace(par.pIn, par.pOut, numel(z)).';
state.pReduced = state.p;
state.p2D = state.p;

state.uE = uE_pre;
state.uEPrev = uE_pre;
state.pPrev = state.p;
state.dtPrev = par.dt;
end

function h = plot_deformed_mesh(mesh, u, color, lineWidth)
rDef = mesh.nodes(:,1) + u(1:2:end);
zDef = mesh.nodes(:,2) + u(2:2:end);

closedConn = mesh.conn(:, [1 2 3 4 1]);
x = rDef(closedConn).';
y = zDef(closedConn).';
x = [x; nan(1, size(x,2))];
y = [y; nan(1, size(y,2))];

h = plot(x(:)*1e6, y(:)*1e6, 'Color', color, 'LineWidth', lineWidth);
end

function rOuter = rounded_outer_profile(z, Rout, zMin, L, Rc)
zloc = z - zMin;
rOuter = Rout * ones(size(zloc));

idxL = zloc < Rc;
xiL = Rc - zloc(idxL);
cutL = Rc - sqrt(max(Rc^2 - xiL.^2, 0));

idxR = zloc > (L - Rc);
xiR = zloc(idxR) - (L - Rc);
cutR = Rc - sqrt(max(Rc^2 - xiR.^2, 0));

rOuter(idxL) = Rout - cutL;
rOuter(idxR) = Rout - cutR;
end

function fe = kelvin_voigt_objective_element_residual_only(Xe, ue, ueOld, mesh, par)
fe = zeros(8,1);

Rnod = Xe(:,1);

rnodOld = Rnod + ueOld(1:2:end);
znodOld = Xe(:,2) + ueOld(2:2:end);

for g = 1:mesh.ngp
    xi  = mesh.gp(g,1);
    eta = mesh.gp(g,2);
    w   = mesh.gw(g);

    [N, dNdxi, ~] = q4_shape(xi, eta, 1.0);
    [~, dNdX, detJ0] = jacobian_2d(Xe, dNdxi);

    Rg = N * Rnod;
    if Rg <= 0
        error('Non-positive radius encountered in viscoelastic element.');
    end

    Fold = deformation_gradient_from_nodal(Xe, rnodOld, znodOld, N, dNdX, Rg);
    F = current_deformation_gradient_from_ue(Xe, ue, N, dNdX, Rg);

    Pvisc = objective_kelvin_voigt_piola(F, Fold, par);
    Wgp = (2*pi*Rg) * detJ0 * w;

    for a = 1:4
        dNa_dR = dNdX(a,1);
        dNa_dZ = dNdX(a,2);
        Na     = N(a);

        fe(2*a-1) = fe(2*a-1) + ...
            (Pvisc(1,1)*dNa_dR + Pvisc(1,3)*dNa_dZ + Pvisc(2,2)*(Na/Rg)) * Wgp;

        fe(2*a) = fe(2*a) + ...
            (Pvisc(3,1)*dNa_dR + Pvisc(3,3)*dNa_dZ) * Wgp;
    end
end
end

function F = current_deformation_gradient_from_ue(Xe, ue, N, dNdX, Rg)
Rnod = Xe(:,1);
Znod = Xe(:,2);
rnod = Rnod + ue(1:2:end);
znod = Znod + ue(2:2:end);
F = deformation_gradient_from_nodal(Xe, rnod, znod, N, dNdX, Rg);
end

% ========================================================================
% FULL-2D AXISYMMETRIC STOKES FLUID MODULE
% ========================================================================

function tf = fluid_supports_partitioned_traction_correction(fluid)
tf = false;
if ~isstruct(fluid) || ~isfield(fluid, 'meshType')
    return;
end
meshType = char(fluid.meshType);
hasTraction = isfield(fluid, 'tractionE') && isfield(fluid, 'tractionL');
tf = strcmpi(meshType, 'bodyfitted_MAC') || ...
    (strcmpi(meshType, 'hybrid_gap1d_exterior2d') && hasTraction);
end

function [N, dNdxi] = shape_Q4(xi, eta)
N = 0.25 * [
    (1-xi)*(1-eta);
    (1+xi)*(1-eta);
    (1+xi)*(1+eta);
    (1-xi)*(1+eta)];

dN_dxi = 0.25 * [
    -(1-eta);
    (1-eta);
    (1+eta);
    -(1+eta)];

dN_deta = 0.25 * [
    -(1-xi);
    -(1+xi);
    (1+xi);
    (1-xi)];

dNdxi = [dN_dxi, dN_deta];
end

function cmp = build_global2d_pressure_traction_comparison(out)
cmp = struct('available', false, ...
    'message', 'global 2D pressure traction data are not available');

if ~isfield(out, 'state') || ~isfield(out.state, 'pEGlobal2D') || ...
        ~isfield(out.state, 'pLGlobal2D')
    return;
end

z = out.z(:);
p1D = out.state.p(:);
pE = out.state.pEGlobal2D(:);
pL = out.state.pLGlobal2D(:);
if ~(numel(p1D) == numel(z) && numel(pE) == numel(z) && numel(pL) == numel(z))
    cmp.message = 'pressure vectors do not match the axial grid length';
    return;
end

dE = pE - p1D;
dL = pL - p1D;
scale = max([max(abs(p1D)), max(abs(pE)), max(abs(pL)), 1]);

cmp = struct();
cmp.available = true;
cmp.z = z;
cmp.p1D = p1D;
cmp.pEGlobal2D = pE;
cmp.pLGlobal2D = pL;
cmp.diffE = dE;
cmp.diffL = dL;
cmp.maxAbsDiffE = max(abs(dE));
cmp.maxAbsDiffL = max(abs(dL));
cmp.rmsDiffE = sqrt(mean(dE.^2));
cmp.rmsDiffL = sqrt(mean(dL.^2));
cmp.maxRelDiffE = cmp.maxAbsDiffE / scale;
cmp.maxRelDiffL = cmp.maxAbsDiffL / scale;
cmp.note = ['p1D is the reduced coupled pressure; pEGlobal2D and ', ...
    'pLGlobal2D are the final projected 2D pressures sampled from the ', ...
    'fluid side of the endothelium and leukocyte interfaces.'];
end

function plot_global2d_pressure_traction_comparison(out, par) %#ok<INUSD>
if ~isfield(out, 'global2DPressureTractionComparison') || ...
        ~out.global2DPressureTractionComparison.available
    return;
end

cmp = out.global2DPressureTractionComparison;
zUm = cmp.z * 1e6;

figure;
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
set(ax1, 'FontSize', 20);
plot(ax1, cmp.p1D, zUm, 'k-', 'LineWidth', 2.0); hold(ax1, 'on');
plot(ax1, cmp.pEGlobal2D, zUm, '-', 'Color', [0.00 0.15 0.65], 'LineWidth', 1.8);
plot(ax1, cmp.pLGlobal2D, zUm, '-', 'Color', [0.00 0.45 0.15], 'LineWidth', 1.8);
grid(ax1, 'on');
xlabel(ax1, 'pressure [Pa]');
ylabel(ax1, 'z [\mum]');
title(ax1, 'Final normal traction pressure');
legend(ax1, {'1D p(z)', 'endothelium 2D sample', 'leukocyte 2D sample'}, ...
    'Location', 'best');

ax2 = nexttile;
set(ax2, 'FontSize', 20);
plot(ax2, cmp.diffE, zUm, '-', 'Color', [0.00 0.15 0.65], 'LineWidth', 1.8); hold(ax2, 'on');
plot(ax2, cmp.diffL, zUm, '-', 'Color', [0.00 0.45 0.15], 'LineWidth', 1.8);
xline(ax2, 0, 'k--', 'LineWidth', 1.0);
grid(ax2, 'on');
xlabel(ax2, '2D sample - 1D p [Pa]');
ylabel(ax2, 'z [\mum]');
title(ax2, sprintf('Difference: max E %.2g Pa, L %.2g Pa', ...
    cmp.maxAbsDiffE, cmp.maxAbsDiffL));
legend(ax2, {'endothelium - 1D', 'leukocyte - 1D'}, 'Location', 'best');

linkaxes([ax1 ax2], 'y');
ylim(ax1, [min(zUm), max(zUm)]);
end

function global1D = build_global_1d_pressure_view(z, pHist, par)
% A light radial lift for visualization/BC bookkeeping. The coupled pressure
% solve remains 1D in z; this field is not a replacement for a 2D solve.
rMax = global_1d_outer_radius(par);
Nr = 161;
if isfield(par, 'global1DPlotNr') && isfinite(par.global1DPlotNr) && par.global1DPlotNr >= 3
    Nr = round(par.global1DPlotNr);
end
r = linspace(0, rMax, Nr).';
radialShape = 1 - (r / rMax).^2;

global1D = struct();
global1D.r = r;
global1D.z = z(:);
global1D.radialShape = radialShape;
global1D.pressureBC = struct( ...
    'zMinPressure', par.pIn, ...
    'zMaxPressure', par.pOut, ...
    'rMaxPressure', 0, ...
    'axisCondition', 'dP/dr = 0');
global1D.note = ['PfullFinal is a radial lift of the 1D pressure: ', ...
    'P(r,z)=p(z)*(1-(r/rMax)^2). It enforces P(rMax)=0 and ', ...
    'dP/dr at r=0, but it is not a true 2D pressure solve.'];

if isempty(pHist) || size(pHist,2) < 1
    global1D.PfullFinal = zeros(numel(r), numel(z));
else
    global1D.PfullFinal = radialShape * pHist(:,end).';
end
end

function global1D = add_global_1d_blank_solid_pressure_view(out, par)
% Store the final pressure field with deformed solid cells removed from view.
global1D = out.global1D;
if ~isfield(global1D, 'PfullFinal') || isempty(global1D.PfullFinal) || ...
        ~isfield(global1D, 'r') || ~isfield(global1D, 'z')
    return;
end

statePlot = out.state;
if isfield(out, 'stateHist') && isfield(out, 'stopStep') && ...
        out.stopStep >= 1 && numel(out.stateHist) >= out.stopStep && ...
        ~isempty(out.stateHist{out.stopStep})
    statePlot = out.stateHist{out.stopStep};
end

solidMask = false(size(global1D.PfullFinal));
if isfield(out, 'meshL') && isfield(statePlot, 'uL') && ...
        ~isempty(out.meshL) && ~isempty(statePlot.uL) && ...
        ~(isfield(par, 'noLeukocyte') && par.noLeukocyte)
    solidMask = solidMask | deformed_solid_mask_on_grid( ...
        out.meshL, statePlot.uL, global1D.r, global1D.z);
end
if isfield(out, 'meshE') && isfield(statePlot, 'uE') && ...
        ~isempty(out.meshE) && ~isempty(statePlot.uE)
    solidMask = solidMask | deformed_solid_mask_on_grid( ...
        out.meshE, statePlot.uE, global1D.r, global1D.z);
end

Pblank = global1D.PfullFinal;
Pblank(solidMask) = NaN;
global1D.solidMaskFinal = solidMask;
global1D.PblankFinal = Pblank;
global1D.blankNote = ['PblankFinal is PfullFinal with points inside the ', ...
    'deformed leukocyte/endothelium solid meshes set to NaN for plotting.'];
end

% ========================================================================
% BODY-FITTED MAC CYLINDRICAL STOKES FLUID MODULE
% Imported from mac_bodyfitted_2D_stokes_rigid_leukocyte_pressureBC.m
% ========================================================================



%% ============================================================
% Diagnostics and post-processing
% ============================================================


%% ============================================================
% Geometry helper
% ============================================================
function mesh = build_rounded_leukocyte_mesh(par)
zMinL = par.zMin;
zMaxL = par.zMax;
if isfield(par, 'zMinL')
    zMinL = par.zMinL;
end
if isfield(par, 'zMaxL')
    zMaxL = par.zMaxL;
end

LzL = zMaxL - zMinL;
if isfield(par, 'LzL')
    LzL = par.LzL;
end

NzSolidL = par.NzSolid;
if isfield(par, 'NzSolidL')
    NzSolidL = par.NzSolidL;
end

zvec = linspace(zMinL, zMaxL, NzSolidL).';
snapTargets = [zMinL; zMaxL; par.zMin; par.zMax];
tolZ = max(100 * eps(max(abs([zvec; snapTargets; LzL]))), ...
    1e-10 * max(abs(LzL), realmin));
for q = snapTargets.'
    zvec(abs(zvec - q) <= tolZ) = q;
end
zvec = unique([zvec; snapTargets]);
nzL = numel(zvec);
rInner = par.RLin * ones(size(zvec));
rOuter = rounded_outer_profile(zvec, par.RLout, zMinL, LzL, par.Rc);

nodes = zeros(nzL * par.NrL, 2);
s = linspace(0, 1, par.NrL);
beta = 2;
sBias = 1 - (1 - s).^beta;

for j = 1:nzL
    rline = rInner(j) + (rOuter(j) - rInner(j)) * sBias;
    for i = 1:par.NrL
        id = sub2ind([nzL, par.NrL], j, i);
        nodes(id,:) = [rline(i), zvec(j)];
    end
end

conn = zeros((par.NrL - 1) * (nzL - 1), 4);
e = 0;
for j = 1:nzL-1
    for i = 1:par.NrL-1
        n1 = sub2ind([nzL, par.NrL], j,   i);
        n2 = sub2ind([nzL, par.NrL], j,   i+1);
        n3 = sub2ind([nzL, par.NrL], j+1, i+1);
        n4 = sub2ind([nzL, par.NrL], j+1, i);
        e = e + 1;
        conn(e,:) = [n1 n2 n3 n4];
    end
end

gp1 = [-1, 1] / sqrt(3);
gw1 = [1, 1];
[g1, g2] = meshgrid(gp1, gp1);
[w1, w2] = meshgrid(gw1, gw1);

mesh = struct();
mesh.nodes = nodes;
mesh.conn = conn;
mesh.nelem = size(conn,1);
mesh.ngp = numel(g1);
mesh.gp = [g1(:), g2(:)];
mesh.gw = w1(:) .* w2(:);
end

function fe = kelvin_voigt_element_residual_only(Xe, ue, ueOld, mesh, par)
if isfield(par, 'useObjectiveKelvinVoigt') && par.useObjectiveKelvinVoigt
    fe = kelvin_voigt_objective_element_residual_only(Xe, ue, ueOld, mesh, par);
    return;
end

fe = zeros(8,1);

Rnod = Xe(:,1);

rnodOld = Rnod + ueOld(1:2:end);
znodOld = Xe(:,2) + ueOld(2:2:end);

for g = 1:mesh.ngp
    xi  = mesh.gp(g,1);
    eta = mesh.gp(g,2);
    w   = mesh.gw(g);

    [N, dNdxi, ~] = q4_shape(xi, eta, 1.0);
    [~, dNdX, detJ0] = jacobian_2d(Xe, dNdxi);

    Rg = N * Rnod;
    if Rg <= 0
        error('Non-positive radius encountered in viscoelastic element.');
    end

    Fold = deformation_gradient_from_nodal(Xe, rnodOld, znodOld, N, dNdX, Rg);

    F = current_deformation_gradient_from_ue(Xe, ue, N, dNdX, Rg);
    Pvisc = par.etaE * (F - Fold) / par.dt;
    Wgp = (2*pi*Rg) * detJ0 * w;

    for a = 1:4
        dNa_dR = dNdX(a,1);
        dNa_dZ = dNdX(a,2);
        Na     = N(a);

        fe(2*a-1) = fe(2*a-1) + ...
            (Pvisc(1,1)*dNa_dR + Pvisc(1,3)*dNa_dZ + Pvisc(2,2)*(Na/Rg)) * Wgp;

        fe(2*a) = fe(2*a) + ...
            (Pvisc(3,1)*dNa_dR + Pvisc(3,3)*dNa_dZ) * Wgp;
    end
end
end

function [stateCur, fluid, ok, stopReason] = ...
    fluid_state_for_endothelium_u(u, pGuess, old, meshE, interfaceE, z, par)

N = numel(z);

[deltaE, UwE] = monolithic_interface_kinematics( ...
    meshE, u, old.uE, interfaceE, z, par);

stateCur = old;
stateCur.uE = u;
stateCur.deltaE = deltaE;
stateCur.deltaL = old.deltaL;
stateCur.p = pGuess;
stateCur.UwE = UwE;
stateCur.UwL = par.UwL * ones(N,1);

[fluid, ok, stopReason] = solve_fluid_reynolds_slip(z, old, stateCur, par);

if ok
    stateCur.p = fluid.p;
    stateCur.tauE = fluid.tauE;
    stateCur.tauL = fluid.tauL;
end
end

function [resNorm, refNorm, Rsolid] = solid_residual_norm( ...
    mesh, u, traction, interfaceNodes, baseNodes, par)

ndof = size(mesh.nodes,1)*2;
[fixDofs, fixVals] = solid_support_conditions(baseNodes, par.supportE);
free = setdiff((1:ndof).', unique(fixDofs(:)));

u(fixDofs) = fixVals;

Fext = zeros(ndof,1);
[Fext, ~] = apply_interface_traction(mesh, u, Fext, interfaceNodes, traction);
Fint = assemble_finite_def_internal_force_only(mesh, u, par);

R = Fint - Fext;
Rsolid = R(free);

resNorm = norm(Rsolid, inf);
refNorm = max([norm(Fext(free), inf), norm(Fint(free), inf), 1e-14]);
end

function [pWall, tauWall, uzWall] = wall_traction_from_stress(meshF, u, pCell, sigmaCell, side)
Nr = meshF.Nr;
Nz = meshF.Nz;

pWall = zeros(Nz,1);
tauWall = zeros(Nz,1);
uzWall = zeros(Nz,1);

for j = 1:Nz
    if strcmpi(side, 'L')
        node = meshF.nodeId(1,j);
        iCell = 1;
    else
        node = meshF.nodeId(Nr,j);
        iCell = Nr-1;
    end

    if j == 1
        jCell = 1;
    elseif j == Nz
        jCell = Nz-1;
    else
        jCell = j-1;
    end

    e = (jCell-1)*(Nr-1) + iCell;

    pWall(j) = pCell(e);
    tauWall(j) = sigmaCell(e,4); % sigma_rz
    uzWall(j) = u(2*node);
end
end

function Q = compute_axisym_flux_from_velocity(meshF, u)
Nr = meshF.Nr;
Nz = meshF.Nz;
Q = zeros(Nz-1,1);

for j = 1:Nz-1
    qsum = 0;
    for i = 1:Nr-1
        nA = meshF.nodeId(i,j);
        nB = meshF.nodeId(i+1,j);
        nC = meshF.nodeId(i+1,j+1);
        nD = meshF.nodeId(i,j+1);

        rmean = mean(meshF.nodes([nA nB nC nD],1));
        uzmean = mean([u(2*nA), u(2*nB), u(2*nC), u(2*nD)]);
        dr = abs(meshF.nodes(nB,1) - meshF.nodes(nA,1));

        % Axisymmetric physical volume flux through a z-cross-section.
        qsum = qsum + 2*pi*rmean*uzmean*dr;
    end
    Q(j) = qsum;
end
end

function h = rounded_gap_profile2(z, H0, L, Rc)
h = H0 * ones(size(z));

% left rounded corner: z in [0, Rc]
idxL = z < Rc;
xiL  = Rc - z(idxL);   % xiL in [0, Rc]
riseL = Rc - sqrt(Rc^2 - xiL.^2);

% right rounded corner: z in [L-Rc, L]
idxR = z > (L - Rc);
xiR  = z(idxR) - (L - Rc);   % xiR in [0, Rc]
riseR = Rc - sqrt(Rc^2 - xiR.^2);

h(idxL) = H0 - riseL;
h(idxR) = H0 - riseR;
end

function ids = find_interface_nodes(mesh, whichSide)
zvals = unique(mesh.nodes(:,2));
nz = numel(zvals);
nnode = size(mesh.nodes,1);
nr = nnode / nz;

if abs(nr - round(nr)) > 1e-12
    error('Cannot infer structured mesh dimensions.');
end

nr = round(nr);

if strcmpi(whichSide,'outer')
    i = nr;
elseif strcmpi(whichSide,'inner')
    i = 1;
else
    error('unknown side');
end

j = (1:nz).';
ids = sub2ind([nz,nr], j, i*ones(nz,1));
end

function delta = extract_deformed_interface_radius(mesh, u, interfaceNodes, zq)
[rDef, zDef] = deformed_interface_curve(mesh, u, interfaceNodes);
delta = interp_curve_values(zDef, rDef, zq);
end

function rL = rounded_leukocyte_profile(z, R0, L, Rc)

rL = R0*ones(size(z));

if Rc <= 0
    return;
end

idxL = z < Rc;
xiL = Rc - z(idxL);
riseL = Rc - sqrt(max(Rc^2 - xiL.^2,0));

idxR = z > (L - Rc);
xiR = z(idxR) - (L - Rc);
riseR = Rc - sqrt(max(Rc^2 - xiR.^2,0));

rL(idxL) = R0 - riseL;
rL(idxR) = R0 - riseR;
end

function [fluid, ok, stopReason, meshF] = solve_fluid_2D_stokes_penalty(z, old, state, par)
% Backward-compatible alias. New code should call solve_fluid_2D_bodyfitted_MAC.
[fluid, ok, stopReason, meshF] = solve_fluid_2D_bodyfitted_MAC(z, old, state, par);
end

function meshF = build_gap_q4_mesh(z, rl, re, Nr)
Nz = numel(z);
nodeId = zeros(Nr,Nz);
nodes = zeros(Nr*Nz,2);

id = 0;
eta = linspace(0,1,Nr).';

for j = 1:Nz
    rcol = rl(j) + eta * (re(j)-rl(j));
    for i = 1:Nr
        id = id + 1;
        nodeId(i,j) = id;
        nodes(id,:) = [rcol(i), z(j)];
    end
end

elems = zeros((Nr-1)*(Nz-1),4);
e = 0;
for j = 1:Nz-1
    for i = 1:Nr-1
        e = e + 1;
        n1 = nodeId(i,j);
        n2 = nodeId(i+1,j);
        n3 = nodeId(i+1,j+1);
        n4 = nodeId(i,j+1);
        elems(e,:) = [n1 n2 n3 n4];
    end
end

meshF = struct();
meshF.nodes = nodes;
meshF.elems = elems;
meshF.nodeId = nodeId;
meshF.Nr = Nr;
meshF.Nz = Nz;
meshF.z = z(:);
end

function mesh = build_rect_mesh(r0, r1, z0, z1, nr, nz)
[R,Z] = meshgrid(linspace(r0,r1,nr), linspace(z0,z1,nz));
nodes = [R(:), Z(:)];
conn = zeros((nr-1)*(nz-1),4);
e = 0;
for j=1:nz-1
    for i=1:nr-1
        n1 = sub2ind([nz,nr], j,   i  );   % lower-left
        n2 = sub2ind([nz,nr], j,   i+1);   % lower-right
        n3 = sub2ind([nz,nr], j+1, i+1);   % upper-right
        n4 = sub2ind([nz,nr], j+1, i  );   % upper-left
        e = e + 1;
        conn(e,:) = [n1 n2 n3 n4];
    end
end
gp1=[-1, 1]/sqrt(3);
gw1 = [1, 1];
[g1,g2] = meshgrid(gp1,gp1);
[w1,w2] = meshgrid(gw1,gw1);
mesh = struct();
mesh.nodes = nodes;
mesh.conn = conn;          % element connectivity
mesh.nelem = size(conn,1); % number of elements
mesh.ngp = numel(g1);
mesh.gp = [g1(:), g2(:)];  % Gauss points
mesh.gw = w1(:).*w2(:);    % weights
end

function uProp = solid_newton_proposal(mesh, u, traction, interfaceNodes, baseNodes, par)
% Residual-reducing solid Newton proposal.
%
% The previous contact march used the raw Newton vector directly.  Near the
% small-gap pressure spike that vector can contain non-equilibrium interface
% modes, which then show up as a wavy surface after the outer gap limiter.  This
% proposal now behaves like the solid Newton line search: it only returns a
% candidate displacement that reduces the solid residual for the prescribed
% traction and preserves element orientation.

ndof = size(mesh.nodes,1)*2;
[fixDofs, fixVals] = solid_support_conditions(baseNodes, par.supportE);
free = setdiff((1:ndof).', unique(fixDofs(:)));

u(fixDofs) = fixVals;

Fext = zeros(ndof,1);
[Fext, Kext] = apply_interface_traction(mesh, u, Fext, interfaceNodes, traction);

[Fint, Ktan] = assemble_finite_def_axisym(mesh, u, par);

R = Fint - Fext;
Ktot = Ktan - Kext;
Rf = R(free);
Kff = Ktot(free,free);
res0 = norm(Rf, inf);

du = zeros(ndof,1);
du(free) = -Kff \ Rf;

if ~all(isfinite(du))
    error('Solid Newton proposal produced non-finite values.');
end

alpha = 1.0;
accepted = false;
uBest = u;
resBest = res0;

for ls = 1:par.lineSearchMax
    uTrial = u + alpha * du;
    uTrial(fixDofs) = fixVals;

    try
        FextTrial = zeros(ndof,1);
        [FextTrial, ~] = apply_interface_traction(mesh, uTrial, FextTrial, interfaceNodes, traction);
        FintTrial = assemble_finite_def_internal_force_only(mesh, uTrial, par);
        Rtrial = FintTrial - FextTrial;
        resTrial = norm(Rtrial(free), inf);

        if resTrial < resBest
            uBest = uTrial;
            resBest = resTrial;
        end

        if resTrial < res0
            accepted = true;
            break;
        end

    catch ME
        if contains(ME.message, 'Negative or zero J') || ...
                contains(ME.message, 'Non-positive radius') || ...
                contains(ME.message, 'Element inverted')
            % Reject and reduce alpha.
        else
            rethrow(ME);
        end
    end

    alpha = 0.5 * alpha;
end

if ~accepted
    if resBest < res0
        uProp = uBest;
    else
        error('Solid Newton proposal line search failed to reduce residual.');
    end
else
    uProp = uTrial;
end
end

function [K, f] = assemble_axisym_stokes_penalty_Q4(meshF, mu, lambda)
nodes = meshF.nodes;
elems = meshF.elems;
nn = size(nodes,1);
ndof = 2*nn;

g = 1/sqrt(3);
gps = [-g -g; g -g; g g; -g g];
wts = [1;1;1;1];

ne = size(elems,1);
rows = zeros(ne*64,1);
cols = zeros(ne*64,1);
vals = zeros(ne*64,1);
ptr = 1;

f = zeros(ndof,1);

for e = 1:ne
    conn = elems(e,:);
    xe = nodes(conn,1);
    ze = nodes(conn,2);
    Ke = zeros(8,8);

    for q = 1:4
        xi = gps(q,1);
        eta = gps(q,2);
        wt = wts(q);

        [N, dNdxi] = shape_Q4(xi, eta);

        J = dNdxi.' * [xe ze];
        detJ = det(J);
        if detJ <= 0
            error('Fluid mesh element has non-positive Jacobian.');
        end

        dNdx = dNdxi / J;
        dNdr = dNdx(:,1);
        dNdz = dNdx(:,2);

        r = N.' * xe;
        if r <= 0
            error('Fluid mesh has non-positive radius.');
        end

        % Axisymmetric velocity-gradient/strain operator.
        B = zeros(4,8);
        Div = zeros(1,8);

        for a = 1:4
            ia = 2*a-1;
            iz = 2*a;

            B(1,ia) = dNdr(a);    % d ur / dr
            B(2,ia) = N(a)/r;     % ur / r
            B(3,iz) = dNdz(a);    % d uz / dz
            B(4,ia) = dNdz(a);    % d ur / dz
            B(4,iz) = dNdr(a);    % d uz / dr

            Div(ia) = dNdr(a) + N(a)/r;
            Div(iz) = dNdz(a);
        end

        Cvis = diag([2*mu, 2*mu, 2*mu, mu]);

        weight = 2*pi*r*detJ*wt;
        Ke = Ke + (B.'*Cvis*B + lambda*(Div.'*Div)) * weight;
    end

    dofs = zeros(8,1);
    for a = 1:4
        dofs(2*a-1) = 2*conn(a)-1;
        dofs(2*a)   = 2*conn(a);
    end

    [rr, cc] = ndgrid(dofs,dofs);
    nadd = numel(rr);
    rows(ptr:ptr+nadd-1) = rr(:);
    cols(ptr:ptr+nadd-1) = cc(:);
    vals(ptr:ptr+nadd-1) = Ke(:);
    ptr = ptr + nadd;
end

rows = rows(1:ptr-1);
cols = cols(1:ptr-1);
vals = vals(1:ptr-1);
K = sparse(rows, cols, vals, ndof, ndof);
end
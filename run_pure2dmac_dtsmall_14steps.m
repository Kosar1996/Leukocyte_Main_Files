%% RUN_INPUT_FULL2D_PRESSURE
% Input deck for the hybrid pressure calculation:
%   1D lubrication mesh in the leukocyte/endothelium gap,
%   2D body-fitted MAC/Stokes mesh in the upstream/downstream exterior.
%
% This file keeps run_input.m unchanged and starts a separate case that uses
% the mixed-dimensional fluid solve inside:
%   softlube_run_case_global_coupled(cfg)
%
% Important distinction:
%   - The gap interval uses the Reynolds/lubrication p(z) solve.
%   - The exterior intervals are solved as true 2D body-fitted MAC/Stokes
%     blocks, matched to the 1D pressure at z = -0.2 and z = 4.2 um.
%   - The final composite pressure field is stored as P(r,z) in out.PHist,
%     on the physical grid out.RPHist/out.ZPHist.
%
% To run:
%   cd('/Users/shu/Softlubrication')
%   run_input_full2D_pressure
%
% Main hybrid outputs:
%   out.PHist(:,:,end)       composite hybrid pressure P(r,z) at final step
%   out.RPHist(:,:,end)      physical r locations for P
%   out.ZPHist(:,:,end)      physical z locations for P
%   out.urCHist(:,:,end)     radial velocity at pressure-cell centers
%   out.uzCHist(:,:,end)     axial velocity at pressure-cell centers
%   out.native2D             compact final-step hybrid 2D fields
%   out.hybridMesh           final 1D-gap / 2D-exterior mesh description

clc;
clear all;

%% 1. Paths
softlubeDir = '';
addpath(softlubeDir);

%% 2. Time Controls
dt = 3e-4;           % time step size [s]
nSteps = 14;          % requested time steps
tEnd = nSteps * dt;

useAdaptiveTimeStep = true;
dtMin = dt / 1024;

%% 3. Global Domain and Pressure Boundary Conditions
% The hybrid pressure uses a 1D gap solve and 2D exterior body-fitted
% r-z blocks. The axial ends use P = 0. In the exterior MAC/Stokes blocks,
% r = 0 is the axis outside the leukocyte and r = 15 um is the exterior wall
% outside the endothelium.
zMin = -6e-6;        % -6 um
zMax = 10e-6;        % 10 um
rOuter = 4e-6;      % 15 um

pAtZMin = 0;         % P = 0 at z = -6 um [Pa]
pAtZMax = 0;         % P = 0 at z = 10 um [Pa]
axisPressureBC = 'dPdr=0';

%% 4. Fluid Mesh Resolution
% Hybrid split. The gap interval is handled by the 1D lubrication mesh; the
% exterior intervals are handled by the 2D body-fitted MAC/Stokes mesh.
gap1DWindow = [-0.2e-6, 4.2e-6];

% Resolve the 1D lubrication region finely and keep the exterior coarse.
finePressureWindow = gap1DWindow;
coarseDz = 1.0e-6;
fineDz = 0.05e-6;

% 2D MAC cells across the radial exterior fluid blocks. The radial grid is
% nonuniform below, with extra resolution in the thin-gap band r = 2..4 um.
NrExterior2D = 20;
NrFluid2D = NrExterior2D;
nrEnv = str2double(getenv('SOFTLUBE_FULL2D_NR'));
if isfinite(nrEnv) && nrEnv >= 5
    NrFluid2D = round(nrEnv);
    NrExterior2D = NrFluid2D;
end
radialFineWindow = [2e-6, 4e-6];
radialFineWeight = 4.0;

%% 5. Plotting and Output
plotNative2DPressure = true;
plotNative2DStress = true;
plotNative2DVelocity = false;
plotFirstStepHybridMesh = false;
plotGlobalDomainSchematic = false;

% Keep the old solver plot block off. The local plot at the bottom shows
% only the composite hybrid result from out.PHist.
makeLegacyPlots = false;

saveOutput = false;
outputFile = fullfile(softlubeDir, 'simulation_output_full2D_pressure.mat');

closeFiguresAtStart = true;
printEvery = 1;

%% 6. Build the Base Case
cfg = struct();
% load('out_2_combined_t40.mat');
% cfg=out.cfg;
% out=rmfield(out,'cfg');


if numel(fieldnames(struct()))==0
    % Prestressed initial geometries. Leukocyte geometry/mesh settings are read
    % inside the solver from solid_leu_P600.mat and are not repeated here.
    cfg.geometry.endotheliumPrestressFile = fullfile(softlubeDir, 'solid_endo_P300_wide.mat');
    cfg.geometry.leukocytePrestressFile = fullfile(softlubeDir, 'solid_leu_P600.mat');
    cfg.geometry.usePrestressedLeukocyteIC = true;
    cfg.geometry.useLeukocyteCenterline = true;

    cfg.fluid.zMin = zMin;
    cfg.fluid.zMax = zMax;
    cfg.fluid.mu = 1.2e-3;
    cfg.fluid.slipE = 0e-9;
    cfg.fluid.slipL = 0e-9;
    cfg.fluid.NrFluid2D = NrFluid2D;
    cfg.fluid.Nr = NrFluid2D;
    cfg.fluid.NrExterior2D = NrExterior2D;
    cfg.fluid.bodyFittedRadialFineWindow = radialFineWindow;
    cfg.fluid.bodyFittedRadialFineWeight = radialFineWeight;
    cfg.fluid.innerBoundary = 'deformed_leukocyte';
    cfg.fluid.useSlopeAwareWallKinematics = true;

    cfg.bc.pIn = pAtZMin;
    cfg.bc.pOut = pAtZMax;
    cfg.bc.UwL = 0;
    cfg.bc.UwE = 0;

    cfg.solid.noLeukocyte = false;
    cfg.solid.leukocyte.enabled = true;
    cfg.solid.leukocyte.deformable = true;
    cfg.solid.leukocyte.EL = 200;
    cfg.solid.leukocyte.nuL = 0.46;
    cfg.solid.leukocyte.useViscoelasticLeukocyte = true;
    cfg.solid.leukocyte.etaL = 1;
    cfg.solid.leukocyte.supportL = 'axis';

    cfg.solid.endothelium.Ee = 500;
    cfg.solid.endothelium.nuE = 0.46;
    cfg.solid.endothelium.useViscoelasticEndothelium = true;
    cfg.solid.endothelium.etaE = 1;

    cfg.interface.useExactDeformedInterface = true;
    cfg.interface.useExactDeformedInterfaceInMonolithic = true;
    cfg.interface.useReferenceZForTractionMapping = false;
    cfg.interface.zeroLeukocyteTractionOutsideOverlap = false;
    cfg.interface.leukocytePressureSupportZ = [zMin, zMax];
    cfg.interface.warnPrestressLoadMismatch = false;
    cfg.numerics.dt = dt;
    cfg.numerics.tEnd = tEnd;
    cfg.numerics.minGap = 0.001e-6;
    cfg.numerics.maxNewtonMono = 100;
    cfg.numerics.tolNewtonMono = 1e-6;
    cfg.numerics.monoSolidAbsTol = 1e-7;
    cfg.numerics.monoFluidAbsTol = 1e-8;
    cfg.numerics.maxFunctionEvaluationsTwoSolid = 400;
    cfg.numerics.enablePressureJumpRetry = false;
    cfg.numerics.maxPressureJumpAbs = inf;
    cfg.numerics.maxPressureJump2DAbs = inf;
    cfg.numerics.useObjectiveKelvinVoigt = true;
    cfg.numerics.useDirectFluidSolve = true;
    cfg.numerics.useFsolveMono = true;
    cfg.numerics.enableAdaptiveTimeStep = useAdaptiveTimeStep;
    cfg.numerics.dtMin = dtMin;
    cfg.numerics.dtRetryFactor = 0.5;
    cfg.numerics.dtGrowFactor = 1.25;
    cfg.numerics.maxTimeStepRetries = 18;
    cfg.numerics.fsolveDisplay = 'off';
    cfg.numerics.diagnosticsEnabled = true;

    cfg.output.saveOutput = saveOutput;
    cfg.output.outputFile = outputFile;
    cfg.output.makePlots = makeLegacyPlots;
    cfg.output.plotFinalGlobalPressureContour = false;
    cfg.output.plotGlobal2DPressureTractionComparison = false;
    cfg.output.printEvery = printEvery;
    cfg.output.storeFull2DFluidHist = false;
    cfg.output.store2DFluidArrays = true;
    cfg.output.store2DPhysicalGridHist = true;
    cfg.output.store2DSpeedHist = true;

    cfg.runtime.useEnvironment = false;
    cfg.ui.closeFigures = closeFiguresAtStart;
    cfg.parOverrides = struct();

    %% 7. Global Axial Mesh and Exterior Boundary Embedding
    % The flag name is historical. Here it is kept on to reuse the global axial
    % mesh and the outside-solid radius rules:
    %   leukocyte absent: fluid inner boundary becomes r = 0
    %   endothelium absent: fluid outer boundary becomes r = 15 um
    % The composite hybrid pressure output is out.PHist, not out.global1D.
    cfg.global1D.zDomain = [zMin, zMax];
    cfg.global1D.rDomain = [0, rOuter];
    cfg.global1D.axisPressureBC = axisPressureBC;
    cfg.global1D.outerPressureBC = 0;

    cfg.parOverrides.useGlobal1DPressure = true;
    cfg.parOverrides.global1DOuterRadius = rOuter;
    cfg.parOverrides.global1DAxisRadius = 1e-9;
    cfg.parOverrides.useGlobal1DCoarseEdgeMesh = true;
    cfg.parOverrides.global1DFineWindow = finePressureWindow;
    cfg.parOverrides.global1DCoarseDz = coarseDz;
    cfg.parOverrides.global1DFineDz = fineDz;
    cfg.parOverrides.global1DPlotNr = 161;

    %% 8. Hybrid Gap-1D / Exterior-2D Fluid Module
    cfg.fluid.useFull2DFluid = true;
    %added_temp
    cfg.fluid.useHybridGap1DExterior2DFluid = false;
    %cfg.fluid.useHybridGap1DExterior2DFluid = false;
    cfg.fluid.hybridGapZ = gap1DWindow;
    cfg.fluid.hybridExteriorMinCells = 4;
    cfg.fluid.hybridExteriorFailMode = 'warn';
    cfg.fluid.hybridUseExterior2DPressureVector = true;
    cfg.fluid.useBodyFittedMACFluid = true;
    cfg.fluid.useBodyFittedMACTractionCorrection = true;
    %added_temp: bumped from 1/0.05 to test whether the interface
    %traction mismatch shrinks with a stronger partitioned correction
    %toward the 2D MAC traction (revert to 1 / 0.05 to go back).
    % apply_bodyfitted_MAC_traction_correction_feedback.m exists in this
    % codebase and was tried here, but confirmed (Aug 21) to NOT converge
    % for this problem -- worst-case mismatch frozen at 47,000-260,000%
    % even with Aitken relaxation on. Left OFF here so this entry point
    % uses the plain loop (apply_bodyfitted_MAC_traction_correction.m),
    % which has a real, tested convergence fix (see Time_Stepping_
    % Followup_Aug_21.docx). Revisit useFeedbackTractionCorrection only
    % after the feedback mechanism's non-convergence is root-caused.
    cfg.fluid.useFeedbackTractionCorrection = false;
    cfg.fluid.useAitkenTractionCorrectionRelax = false;
    cfg.fluid.debugRadialAxialTractionCorrection = false;
    cfg.fluid.maxBodyFittedTractionCorrections = 100;
    cfg.fluid.bodyFittedTractionCorrectionRelax = 0.5;
    cfg.fluid.bodyFittedTractionCorrectionFailMode = 'warn';
    cfg.fluid.resolveFluidAfterTractionCorrection = true;
    cfg.fluid.useBodyFittedMACTractionInSolid = true;

    cfg.parOverrides.useFull2DFluid = true;
    %added_temp
    cfg.parOverrides.useHybridGap1DExterior2DFluid = false;
    %cfg.parOverrides.useHybridGap1DExterior2DFluid = false;
    cfg.parOverrides.hybridGapZ = gap1DWindow;
    cfg.parOverrides.NrExterior2D = NrExterior2D;
    cfg.parOverrides.bodyFittedRadialFineWindow = radialFineWindow;
    cfg.parOverrides.bodyFittedRadialFineWeight = radialFineWeight;
    cfg.parOverrides.hybridExteriorMinCells = 4;
    cfg.parOverrides.hybridExteriorFailMode = 'warn';
    cfg.parOverrides.hybridUseExterior2DPressureVector = true;
    cfg.parOverrides.useBodyFittedMACFluid = true;
    cfg.parOverrides.useBodyFittedMACTractionCorrection = true;
    % parOverrides wins over cfg.fluid, so this copy must match the one
    % above -- see the comment there for why these are off.
    cfg.parOverrides.useFeedbackTractionCorrection = false;
    cfg.parOverrides.useAitkenTractionCorrectionRelax = false;
    cfg.parOverrides.maxBodyFittedTractionCorrections = 100;
    cfg.parOverrides.bodyFittedTractionCorrectionRelax = 0.5;
    cfg.parOverrides.bodyFittedTractionCorrectionFailMode = 'warn';
    cfg.parOverrides.resolveFluidAfterTractionCorrection = true;
    cfg.parOverrides.useBodyFittedMACTractionInSolid = true;
    cfg.parOverrides.NrFluid2D = NrFluid2D;
    cfg.parOverrides.Nr = NrFluid2D;

    % Turn off the earlier projected-pressure bridge so this file reports the
    % hybrid pressure field cleanly.
    cfg.parOverrides.useGlobal2DPressureTraction = false;
else
    cfg.numerics.dt = dt;
    cfg.numerics.tEnd = tEnd;
    if ~isempty(out)
        out.par.dt = dt;
        out.par.tEnd = tEnd;
    end
end

%% 9. Runtime
if plotNative2DPressure || plotNative2DVelocity || ...
        plotFirstStepHybridMesh || plotGlobalDomainSchematic
    set(0, 'DefaultFigureVisible', 'on');
end

fprintf('\nRunning Softlubrication hybrid gap-1D / exterior-2D pressure case\n');
fprintf('   dt                  = %.6e s\n', dt);
fprintf('   tEnd                = %.6e s\n', tEnd);
fprintf('   requested steps     = %.0f\n', nSteps);
fprintf('   z-domain            = [%.3f, %.3f] um\n', zMin*1e6, zMax*1e6);
fprintf('   exterior radius     = %.3f um\n', rOuter*1e6);
fprintf('   1D gap window       = [%.3f, %.3f] um\n', gap1DWindow(1)*1e6, gap1DWindow(2)*1e6);
fprintf('   exterior 2D radial Nr = %d\n\n', NrExterior2D);

if exist('out')
out = softlube_run_case_global_coupled(cfg, out);
else
    out = softlube_run_case_global_coupled(cfg);
end

%% 10. Compact Final Native 2D Fields
out.native2D = collect_final_native2d(out);
out.firstStepHybridMesh = collect_hybrid_mesh_at_step( ...
    out, 1, gap1DWindow, NrExterior2D);
out.hybridMesh = collect_final_hybrid_mesh(out, gap1DWindow, NrExterior2D);

if plotFirstStepHybridMesh
    plot_hybrid_mesh_at_step(out, out.firstStepHybridMesh, 1);
end

if plotNative2DPressure
    %plot_final_native2d_pressure(out);
    % Bug fix: plot the actual last accepted step (out.stopStep), not the
    % originally-requested step count (nSteps) -- these differ whenever a
    % run stops early (non-convergence, reaching the min-gap floor, etc),
    % and indexing by nSteps crashes with an out-of-bounds error in that
    % case even though the run itself completed and saved successfully.
    plot_select_native2d_pressure(out,out.stopStep);
end

if plotNative2DStress
    % Bug fix: same out.stopStep vs nSteps issue as the pressure plot above.
    plot_select_native2d_stress(out,out.stopStep);
end

if plotNative2DVelocity
    plot_final_native2d_velocity(out);
end

if plotGlobalDomainSchematic
    out.postPlots.globalDomain = softlube_plot_global1D_domain(cfg, out);
end

%% 11. Summary
fprintf('\nHybrid gap-1D / exterior-2D pressure run finished\n');
fprintf('   stopStep = %d\n', out.stopStep);
fprintf('   final t  = %.6e s\n', out.t(end));

if isfield(out.native2D, 'P') && any(isfinite(out.native2D.P(:)))
    fprintf('   native P range = [%.6e, %.6e] Pa\n', ...
        min(out.native2D.P(:), [], 'omitnan'), ...
        max(out.native2D.P(:), [], 'omitnan'));
end

if isfield(out, 'p2DMaxHist') && any(isfinite(out.p2DMaxHist))
    fprintf('   final max |P2D| = %.6e Pa\n', out.p2DMaxHist(end));
end

out.cfg=cfg;
save('out_pure2dmac_dtsmall_14steps.mat','out');

%% Local Plot Helpers
function native2D = collect_final_native2d(out)
native2D = struct();
if ~isfield(out, 'stopStep') || out.stopStep < 1 || ...
        ~isfield(out, 'PHist') || isempty(out.PHist)
    return;
end

k = out.stopStep;
native2D.P = out.PHist(:,:,k);
native2D.R = out.RPHist(:,:,k);
native2D.Z = out.ZPHist(:,:,k);
native2D.ur = out.urCHist(:,:,k);
native2D.uz = out.uzCHist(:,:,k);
native2D.speed = out.speedCHist(:,:,k);
native2D.note = ['Hybrid pressure field: 1D lubrication in the gap ', ...
    'and body-fitted 2D MAC/Stokes in the exterior blocks. ', ...
    'Coordinates are physical r,z locations of pressure-cell centers.'];
end

function hybridMesh = collect_final_hybrid_mesh(out, gap1DWindow, NrExterior2D)
hybridMesh = struct();
if ~isfield(out, 'stopStep') || out.stopStep < 1
    return;
end
hybridMesh = collect_hybrid_mesh_at_step( ...
    out, out.stopStep, gap1DWindow, NrExterior2D);
end

function hybridMesh = collect_hybrid_mesh_at_step(out, stepIdx, gap1DWindow, NrExterior2D)
hybridMesh = struct();
if nargin < 2 || isempty(stepIdx)
    stepIdx = 1;
end
stepIdx = round(stepIdx);

if ~isfield(out, 'stopStep') || out.stopStep < stepIdx || ...
        stepIdx < 1 || ~isfield(out, 'fluidHist') || ...
        numel(out.fluidHist) < stepIdx || isempty(out.fluidHist{stepIdx})
    return;
end

statePlot = state_for_plot_at_step(out, stepIdx);

opts = struct();
opts.gapZ = gap1DWindow;
opts.NrExterior = NrExterior2D + 1;
hybridMesh = softlube_build_hybrid_gap1d_exterior2d_mesh( ...
    out, statePlot, out.fluidHist{stepIdx}, opts);
hybridMesh.stepIndex = stepIdx;
if isfield(out, 't') && numel(out.t) >= stepIdx
    hybridMesh.t = out.t(stepIdx);
end
end


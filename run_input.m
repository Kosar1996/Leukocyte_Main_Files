%% RUN_INPUT
% Input deck for the global-domain leukocyte/endothelium calculation.
%
% This script uses the coupled solver:
%   softlube_run_case_global_coupled(cfg)
%
% Model used here:
%   1. The pressure unknown is solved on the full axial domain z = [-6, 10] um.
%   2. The exterior radial box is r = 15 um.
%   3. Boundary conditions are:
%        P = 0 at z = -6 um
%        P = 0 at z = 10 um
%        P = 0 at r = 15 um
%        dP/dr = 0 at r = 0
%   4. The solid normal traction uses the projected global 2D pressure field.
%   5. The shear traction still comes from the stable 1D lubrication solve.
%
% To run:
%   cd('/Users/shu/Softlubrication')
%   run_input
%
% Main outputs:
%   out                                  complete solver output
%   out.global1D                         final lifted pressure field
%   out.state.global2DPressureTraction   final projected 2D pressure field
%   out.global2DPressureTractionComparison

clc;
clearvars;

%% 1. Paths and Reloading
% Keep all project files on the MATLAB path. Clearing these functions makes
% MATLAB use the latest edited versions instead of a cached copy.
softlubeDir = '';
addpath(softlubeDir);


%% 2. Time Controls
% Change these three lines for most runs.
dt = 1e-4;           % time step size [s]
nSteps = 40;          % number of requested time steps
tEnd = nSteps * dt;  % final time [s]

% Adaptive retry is kept on. If a step fails, the solver may temporarily
% reduce dt, so accepted dt values are stored in out.dtHist.
useAdaptiveTimeStep = true;
dtMin = dt / 1024;

%% 3. Global Pressure Domain and Boundary Conditions
% All lengths are in meters. The comments show the equivalent microns.
zMin = -6e-6;        % -6 um, inlet/far-left pressure boundary
zMax = 10e-6;        % 10 um, outlet/far-right pressure boundary
rOuter = 15e-6;      % 15 um, exterior radial pressure boundary

pAtZMin = 0;         % P = 0 at z = -6 um [Pa]
pAtZMax = 0;         % P = 0 at z = 10 um [Pa]
pAtROuter = 0;       % P = 0 at r = 15 um [Pa]
axisPressureBC = 'dPdr=0';  % symmetry at r = 0

%% 4. Fluid Mesh Resolution
% Fine axial spacing is used near the solids; crude spacing is used far away.
% The fine window should cover the leukocyte/endothelium interaction region.
finePressureWindow = [-2e-6, 6e-6];
coarseDz = 0.5e-6;
fineDz = 0.1e-6;

% Radial points used for the final contour and for the 2D pressure-traction
% projection. Larger numbers give smoother plots but cost more inside fsolve.
global1DPlotNr = 161;
global2DPressureTractionNr = 81;
global2DPressureTractionScreenLength = 0.5e-6;

%% 5. Plotting and Output
% Available plot switches for this case:
%   plotFinalGlobalPressureContour
%       Final r-z pressure contour. Solid regions are blank/white.
%       Function location: softlube_run_case_global_coupled.m
%
%   plotPressureTractionComparison
%       Final comparison between p(z), endothelium sampled pressure, and
%       leukocyte sampled pressure.
%       Function location: softlube_run_case_global_coupled.m
%
%   makeLegacyPlots
%       Older solver diagnostic figures: pressure history, velocity field,
%       mesh/displacement/gap plots, and stored full-2D fields if enabled.
%       This can open many figures, so the default is false.
%       Function location: softlube_run_case_global_coupled.m
%
%   plotGlobalDomainSchematic
%       Standalone picture of the global pressure box, axial mesh, full
%       solid shapes, and boundary conditions.
%       Function location: softlube_plot_global1D_domain.m
%
%   plotStandalone2DProjection
%       Standalone screened 2D pressure projection diagnostic. This is a
%       post-processing plot, not a replacement for the coupled pressure solve.
%       Function location: softlube_global2D_pressure_projection.m
makeLegacyPlots = true;
plotFinalGlobalPressureContour = true;
plotPressureTractionComparison = true;
plotGlobalDomainSchematic = true;
plotStandalone2DProjection = true;

domainPlotOptions = struct();
domainPlotOptions.makeFigure = true;
domainPlotOptions.showTractionSupport = true;

pressureProjectionPlotOptions = struct();
pressureProjectionPlotOptions.makePlot = true;
pressureProjectionPlotOptions.screenLength = global2DPressureTractionScreenLength;
pressureProjectionPlotOptions.nContour = 48;

saveOutput = false;
outputFile = fullfile(softlubeDir, 'simulation_output_global1D_coupled.mat');

closeFiguresAtStart = true;
printEvery = 1;

%% 6. Build the Base Case
% The input deck is written out explicitly so this file can be read on its
% own. These are the baseline geometry, material, interface, and solver
% settings before the global-domain overrides below.
cfg = struct();

% Prestressed initial geometries. The leukocyte size/mesh settings
% (RLin, RLout, zMinL, zMaxL, NrL, etc.) are read inside the solver from
% solid_leu_P600.mat, so they are intentionally not repeated here.
cfg.geometry.endotheliumPrestressFile = fullfile(softlubeDir, 'prestress_IC.mat');
cfg.geometry.leukocytePrestressFile = fullfile(softlubeDir, 'solid_leu_P600.mat');
cfg.geometry.usePrestressedLeukocyteIC = true;
cfg.geometry.useLeukocyteCenterline = true;

% Fluid properties. The local 0..4 um values are overwritten by the global
% z-domain in the next section, but the viscosity/slip/interface choices are
% still used.
cfg.fluid.zMin = 0e-6;
cfg.fluid.zMax = 4e-6;
cfg.fluid.mu = 1.2e-3;
cfg.fluid.slipE = 0e-9;
cfg.fluid.slipL = 0e-9;
cfg.fluid.NrFluid2D = 30;
cfg.fluid.innerBoundary = 'deformed_leukocyte';
cfg.fluid.useSlopeAwareWallKinematics = true;

% Boundary conditions. The pressure values are overwritten below by the
% named global BC variables pAtZMin and pAtZMax.
cfg.bc.pIn = 0;
cfg.bc.pOut = 0;
cfg.bc.UwL = 0;
cfg.bc.UwE = 0;

% Leukocyte material and support.
cfg.solid.noLeukocyte = false;
cfg.solid.leukocyte.enabled = true;
cfg.solid.leukocyte.deformable = true;
cfg.solid.leukocyte.EL = 200;
cfg.solid.leukocyte.nuL = 0.46;
cfg.solid.leukocyte.useViscoelasticLeukocyte = true;
cfg.solid.leukocyte.etaL = 1;
cfg.solid.leukocyte.supportL = 'axis';

% Endothelium material.
cfg.solid.endothelium.Ee = 500;
cfg.solid.endothelium.nuE = 0.46;
cfg.solid.endothelium.useViscoelasticEndothelium = true;
cfg.solid.endothelium.etaE = 1;

% Exact interface mapping is important here because pressure traction is
% applied on the deformed solid boundaries.
cfg.interface.useExactDeformedInterface = true;
cfg.interface.useExactDeformedInterfaceInMonolithic = true;
cfg.interface.useReferenceZForTractionMapping = false;
cfg.interface.zeroLeukocyteTractionOutsideOverlap = true;
cfg.interface.leukocytePressureSupportZ = [cfg.fluid.zMin, cfg.fluid.zMax];
cfg.interface.warnPrestressLoadMismatch = false;

% Baseline nonlinear solver controls.
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

cfg.output.storeFull2DFluidHist = false;
cfg.runtime.useEnvironment = false;
cfg.ui.closeFigures = closeFiguresAtStart;
cfg.parOverrides = struct();

%% 7. Domain and Boundary-Condition Overrides
% Important: z = 0 and z = 4 um are interior locations in this global-domain
% run. Do not impose pressure boundary conditions there.
cfg.fluid.zMin = zMin;
cfg.fluid.zMax = zMax;

cfg.global1D.zDomain = [zMin, zMax];
cfg.global1D.rDomain = [0, rOuter];
cfg.global1D.axisPressureBC = axisPressureBC;
cfg.global1D.outerPressureBC = pAtROuter;

cfg.bc.pIn = pAtZMin;
cfg.bc.pOut = pAtZMax;

%% 8. Interface and Solid-Traction Support
% The leukocyte extends beyond the original endothelium overlap. For the
% global-domain calculation, allow pressure traction over the full global
% axial domain instead of zeroing it outside z = 0..4 um.
cfg.interface.zeroLeukocyteTractionOutsideOverlap = false;
cfg.interface.leukocytePressureSupportZ = [zMin, zMax];

% Keep exact deformed-interface interpolation in the monolithic solve.
cfg.interface.useExactDeformedInterface = true;
cfg.interface.useExactDeformedInterfaceInMonolithic = true;
cfg.interface.useReferenceZForTractionMapping = false;

%% 9. Time-Stepping and Nonlinear Solver Controls
cfg.numerics.dt = dt;
cfg.numerics.tEnd = tEnd;
cfg.numerics.enableAdaptiveTimeStep = useAdaptiveTimeStep;
cfg.numerics.dtMin = dtMin;
cfg.numerics.dtRetryFactor = 0.5;
cfg.numerics.dtGrowFactor = 1.25;
cfg.numerics.maxTimeStepRetries = 18;
cfg.numerics.fsolveDisplay = 'off';
cfg.numerics.diagnosticsEnabled = true;

%% 10. Output Controls
cfg.output.saveOutput = saveOutput;
cfg.output.outputFile = outputFile;
cfg.output.makePlots = makeLegacyPlots;
cfg.output.plotFinalGlobalPressureContour = plotFinalGlobalPressureContour;
cfg.output.plotGlobal2DPressureTractionComparison = plotPressureTractionComparison;
cfg.output.printEvery = printEvery;

% The current global-projection case does not need to store full 2D MAC
% Stokes arrays, because that module is intentionally disabled below.
cfg.output.store2DFluidArrays = false;
cfg.output.store2DPhysicalGridHist = false;
cfg.output.store2DSpeedHist = false;

%% 11. Global 1D Pressure Solve
% The 1D pressure equation is still the coupled pressure unknown. It is solved
% over the full z-domain, then lifted radially for visualization/BC bookkeeping.
cfg.parOverrides.useGlobal1DPressure = true;
cfg.parOverrides.global1DOuterRadius = rOuter;
cfg.parOverrides.global1DAxisRadius = 1e-9;
cfg.parOverrides.useGlobal1DCoarseEdgeMesh = true;
cfg.parOverrides.global1DFineWindow = finePressureWindow;
cfg.parOverrides.global1DCoarseDz = coarseDz;
cfg.parOverrides.global1DFineDz = fineDz;
cfg.parOverrides.global1DPlotNr = global1DPlotNr;

%% 12. Global 2D Pressure Traction
% This is the bridge toward a full 2D pressure model. The solid normal loads
% sample a projected global r-z pressure field:
%   endothelium normal load: +pEGlobal2D
%   leukocyte normal load:   -pLGlobal2D
%
% Blend = 1 means use the projected 2D pressure completely for normal loads.
% Blend = 0 would recover the original 1D normal pressure load.
cfg.parOverrides.useGlobal2DPressureTraction = true;
cfg.parOverrides.global2DPressureTractionBlend = 1.0;
cfg.parOverrides.global2DPressureTractionNr = global2DPressureTractionNr;
cfg.parOverrides.global2DPressureTractionScreenLength = global2DPressureTractionScreenLength;

%% 13. Full 2D Stokes/MAC Module
% Keep this off for the current model. The current model uses:
%   pressure: global 1D solve plus projected global 2D traction
%   shear:    1D lubrication shear
cfg.fluid.useFull2DFluid = false;
cfg.fluid.useBodyFittedMACFluid = false;
cfg.fluid.useBodyFittedMACTractionCorrection = false;
cfg.fluid.useBodyFittedMACTractionInSolid = false;
cfg.fluid.resolveFluidAfterTractionCorrection = false;

cfg.parOverrides.useFull2DFluid = false;
cfg.parOverrides.useBodyFittedMACFluid = false;
cfg.parOverrides.useBodyFittedMACTractionCorrection = false;
cfg.parOverrides.useBodyFittedMACTractionInSolid = false;
cfg.parOverrides.resolveFluidAfterTractionCorrection = false;
cfg.parOverrides.store2DFluidArrays = false;
cfg.parOverrides.store2DPhysicalGridHist = false;
cfg.parOverrides.store2DSpeedHist = false;

%% 14. Runtime Behavior
% Set useEnvironment to true only if you want SOFTLUBE_* environment variables
% to override values from this input deck.
cfg.runtime.useEnvironment = false;
cfg.ui.closeFigures = closeFiguresAtStart;

if plotFinalGlobalPressureContour || plotPressureTractionComparison || ...
        makeLegacyPlots || plotGlobalDomainSchematic || plotStandalone2DProjection
    set(0, 'DefaultFigureVisible', 'on');
end

%% 15. Run
fprintf('\nRunning Softlubrication global pressure case\n');
fprintf('   dt                 = %.6e s\n', dt);
fprintf('   tEnd               = %.6e s\n', tEnd);
fprintf('   requested steps    = %.0f\n', nSteps);
fprintf('   z-domain           = [%.3f, %.3f] um\n', zMin*1e6, zMax*1e6);
fprintf('   r-domain           = [0, %.3f] um\n', rOuter*1e6);
fprintf('   global 2D traction = on\n\n');
fprintf('   plots: contour=%d, traction=%d, legacy=%d, domain=%d, projection=%d\n\n', ...
    plotFinalGlobalPressureContour, plotPressureTractionComparison, ...
    makeLegacyPlots, plotGlobalDomainSchematic, plotStandalone2DProjection);

out = softlube_run_case_global_coupled(cfg);

%% 16. Optional Post-Run Plots
% These two plotting functions live outside the solver. Their returned
% diagnostics are stored in out.postPlots when enabled.
postPlots = struct();

if plotGlobalDomainSchematic
    postPlots.globalDomain = softlube_plot_global1D_domain( ...
        cfg, out, domainPlotOptions);
end

if plotStandalone2DProjection
    postPlots.global2DProjection = softlube_global2D_pressure_projection( ...
        out, pressureProjectionPlotOptions);
end

if ~isempty(fieldnames(postPlots))
    out.postPlots = postPlots;
end

%% 17. Summary
fprintf('\nRun finished\n');
fprintf('   stopStep = %d\n', out.stopStep);
fprintf('   final t  = %.6e s\n', out.t(end));

if isfield(out, 'dtHist') && ~isempty(out.dtHist)
    fprintf('   accepted dt range = [%.6e, %.6e] s\n', ...
        min(out.dtHist), max(out.dtHist));
end

if isfield(out, 'global2DPressureTractionComparison') && ...
        out.global2DPressureTractionComparison.available
    cmp = out.global2DPressureTractionComparison;
    fprintf('   max |pEGlobal2D - p1D| = %.6e Pa\n', cmp.maxAbsDiffE);
    fprintf('   max |pLGlobal2D - p1D| = %.6e Pa\n', cmp.maxAbsDiffL);
end

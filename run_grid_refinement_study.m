%% RUN_GRID_REFINEMENT_STUDY
% Way-1 diagnostic sweep for the PI's pressure/stress-trace consistency
% question: reruns the run_input_full2D_pressure2.m case at increasing
% exterior-2D grid resolution, and at each resolution checks how well
% -(1/3)*trace(sigma) (the recovered fluid stress) matches the solver's
% actual pressure field p, using check_pressure_stress_consistency.m.
%
% If the mismatch shrinks under refinement, it is ordinary convergent
% post-processing discretization error (see the explanation for the
% 0.37 Pa observation) and can be reported/bounded rather than treated
% as a bug. If it does not shrink, something else needs investigating
% (e.g. the 1D lubrication-gap region, which uses an analytic velocity
% profile rather than the MAC solve and is not affected by this
% refinement in the same way).
%
% Requires build_cfg_full2D_pressure2.m and
% check_pressure_stress_consistency.m on the MATLAB path alongside the
% rest of the solver codebase.
%
% NOTE: each resolution below reruns the full nonlinear FSI timestep
% solve (fsolve). Start with the default list; NrExterior2D=80 in
% particular can be considerably slower than NrExterior2D=20. Trim the
% list if this is too slow on your machine.

clc;
clearvars;
close all;

resolutions = [10, 20, 40, 80];   % NrExterior2D values to test
plotstep = 1;                     % which accepted timestep to check

nRes = numel(resolutions);
diffMax  = nan(nRes,1);
diffRelP = nan(nRes,1);
diffRelV = nan(nRes,1);
divInfPa = nan(nRes,1);
gridDz   = nan(nRes,1);

results = struct('NrExterior2D', {}, 'diag', {});

for k = 1:nRes
    Nr = resolutions(k);
    % Halve fineDz/coarseDz each time Nr doubles, so radial and axial
    % resolution refine together rather than one lagging the other.
    fineDz = 0.05e-6 * (20 / Nr);
    coarseDz = 1.0e-6 * (20 / Nr);

    fprintf('\n=== Grid refinement step %d/%d: NrExterior2D=%d, fineDz=%.4g um ===\n', ...
        k, nRes, Nr, fineDz*1e6);

    cfg = build_cfg_full2D_pressure2( ...
        'NrExterior2D', Nr, 'fineDz', fineDz, 'coarseDz', coarseDz);

    out = softlube_run_case_global_coupled(cfg);

    if out.stopStep < plotstep
        warning('Run at NrExterior2D=%d stopped before step %d; skipping.', Nr, plotstep);
        continue;
    end

    diag = check_pressure_stress_consistency(out, plotstep);

    results(end+1) = struct('NrExterior2D', Nr, 'diag', diag); %#ok<SAGROW>
    diffMax(k)  = diag.diffMax;
    diffRelP(k) = diag.diffRelToP;
    diffRelV(k) = diag.diffRelToVisc;
    divInfPa(k) = diag.divInfAsPa;
    gridDz(k)   = fineDz;
end

%% Summary table
fprintf('\n=== Grid refinement summary ===\n');
fprintf('%12s %12s %14s %14s %14s %14s\n', ...
    'NrExt2D', 'fineDz[um]', 'diffMax[Pa]', '%% of maxP', '%% of visc', 'divInf[Pa]');
for k = 1:nRes
    if isnan(diffMax(k)), continue; end
    fprintf('%12d %12.4f %14.6e %14.4f %14.4f %14.6e\n', ...
        resolutions(k), gridDz(k)*1e6, diffMax(k), ...
        100*diffRelP(k), 100*diffRelV(k), divInfPa(k));
end

%% Convergence plot
valid = ~isnan(diffMax);
if nnz(valid) >= 2
    figure;
    loglog(gridDz(valid)*1e6, diffMax(valid), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    grid on; hold on;
    xlabel('fine axial grid spacing [\mum]');
    ylabel('max|p - (-1/3 trace(\sigma))|  [Pa]');
    title('Pressure/stress-trace consistency vs. grid resolution');

    % Reference 1st- and 2nd-order slopes for visual comparison.
    dzRef = gridDz(valid)*1e6;
    firstIdx = find(valid, 1, 'first');
    y0 = diffMax(firstIdx);
    dz0 = gridDz(firstIdx)*1e6;
    loglog(dzRef, y0*(dzRef/dz0), '--');
    loglog(dzRef, y0*(dzRef/dz0).^2, ':');
    legend({'measured', '1st order ref.', '2nd order ref.'}, 'Location', 'best');
end

save('grid_refinement_study_results.mat', 'results', 'resolutions', ...
    'diffMax', 'diffRelP', 'diffRelV', 'divInfPa', 'gridDz');

fprintf('\nSaved results to grid_refinement_study_results.mat\n');

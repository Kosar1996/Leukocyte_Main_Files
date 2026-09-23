function plot_select_native2d_stress(out, plotstep)
% PLOT_SELECT_NATIVE2D_STRESS
% Per PI request: instead of the single averaged hydrostatic-pressure
% plot (see plot_select_native2d_stress_backup.m), this now plots the
% four Cauchy stress components sigma_rr, sigma_zz, sigma_tt
% (theta-theta), sigma_rz for the fluid and both solids together, one
% component per tile, in a single 2x2 figure.
%
% REVISION HISTORY & MERGED BUG FIXES:
% -------------------------------------------------------------------------
% 1. Lines 75-98 (Explicit Viscoelastic Parameter Remapping):
%    [OLD]: parL = out.par; if isfield(out.par, 'GL') ...
%    FIXED: Built explicit parE and parL structures to ensure both endothelial 
%    (etaE) and leukocyte (etaL remapped to etaE) viscous terms are passed 
%    to recover_nodal_stress_axisym_viscoelastic.m. Without explicit etaE in parE, 
%    the endothelial recovery would silently fall back to elastic-only stress.
%
% 2. Lines 100-112 (Signature Alignment):
%    FIXED: Verified function invocation strictly matches the 5-argument 
%    signature (mesh, u, uOld, dt, par) for both solid domains.
% -------------------------------------------------------------------------
% Ensure z grid is present for axis calculations
    if ~isfield(out, 'z') && isfield(out, 'meshE') && isfield(out.meshE, 'nodes')
        out.z = out.meshE.nodes(:,2);
    end

fluid = out.fluidHist{plotstep};

% ---------------------------------------------------------
% Locate the recovered fluid total-stress array and its matching
% cell-center coordinates.
% ---------------------------------------------------------
if isfield(fluid, 'sigmaCellNode') && ~isempty(fluid.sigmaCellNode)
    fluidSigma  = fluid.sigmaCellNode;
    fluidCenter = fluid.centerNode;
elseif isfield(fluid, 'sigmaCell') && ~isempty(fluid.sigmaCell)
    fluidSigma  = fluid.sigmaCell;
    fluidCenter = fluid.cellCenter;
else
    warning(['No recovered fluid stress field was found. ', ...
        'Expected sigmaCellNode or sigmaCell.']);
    return;
end

if size(fluidSigma, 2) < 4
    warning('The fluid stress array must have four columns: rr, tt, zz, rz.');
    return;
end

if isempty(fluidCenter) || size(fluidCenter, 1) ~= size(fluidSigma, 1)
    warning(['Fluid stress values and fluid center coordinates do not ', ...
        'match; skipping fluid stress plot.']);
    fluidCenter = [];
end

if ~any(isfinite(fluidSigma(:)))
    warning('The recovered fluid stress field contains no finite values.');
    return;
end

% ---------------------------------------------------------
% Recover solid stresses once. Uses the viscoelastic recovery 
% (elastic + Kelvin-Voigt) by default so the solid stress is
% directly comparable to the fluid's total stress at the interface;
% falls back to elastic-only if the previous-step state or the accepted
% dt for this step aren't available.
% ---------------------------------------------------------
st = out.stateHist{plotstep};

% Rebuild leukocyte parameter structure (parL)
parL = out.par;
if isfield(out.par, 'GL') && isfinite(out.par.GL)
    parL.Ge = out.par.GL;
end
if isfield(out.par, 'KL') && isfinite(out.par.KL)
    parL.Ke = out.par.KL;
end

% [OLD]: if isfield(out.par, 'etaL') && isfinite(out.par.etaL)
% [OLD]:     parL.etaE = out.par.etaL;
% [OLD]: end
if isfield(out.par, 'etaL') && isfinite(out.par.etaL)
    parL.etaE = out.par.etaL;
end

% Explicit parameter structure for endothelium (parE)
parE = out.par;
if ~isfield(parE, 'etaE') && isfield(out.par, 'etaE')
    parE.etaE = out.par.etaE;
end

useViscoRecovery = isfield(st, 'uEPrev') && isfield(st, 'uLPrev') && ...
    ~isempty(st.uEPrev) && ~isempty(st.uLPrev) && ...
    isfield(out, 'dtHist') && numel(out.dtHist) >= plotstep && ...
    isfinite(out.dtHist(plotstep));

if useViscoRecovery
    dtStep = out.dtHist(plotstep);

    % [OLD]: stressE = recover_nodal_stress_axisym_viscoelastic(out.meshE, st.uE, st.uEPrev, dtStep, out.par);
    stressE = recover_nodal_stress_axisym_viscoelastic( ...
        out.meshE, st.uE, st.uEPrev, dtStep, parE);

    stressL = recover_nodal_stress_axisym_viscoelastic( ...
        out.meshL, st.uL, st.uLPrev, dtStep, parL);
else
    warning(['plot_select_native2d_stress: uEPrev/uLPrev or out.dtHist ', ...
        'not available for step %d; falling back to the elastic-only ', ...
        'solid stress. The viscous Kelvin-Voigt term is omitted, so the ', ...
        'solid-vs-fluid interface comparison in this plot will be unfair.'], ...
        plotstep);

    % [OLD]: stressE = recover_nodal_stress_axisym(out.meshE, st.uE, out.par);
    stressE = recover_nodal_stress_axisym( ...
        out.meshE, st.uE, parE);

    stressL = recover_nodal_stress_axisym( ...
        out.meshL, st.uL, parL);
end

% Fluid sigmaCellNode/sigmaCell ordering is [sigma_rr, sigma_tt, sigma_zz, sigma_rz].
% Map that onto the display order [rr, zz, tt, rz] used below.
fluidColumns = [1, 3, 2, 4];

componentTitles = { ...
    '\sigma_{rr}', ...
    '\sigma_{zz}', ...
    '\sigma_{\theta\theta}', ...
    '\sigma_{rz}'};

% Corresponding solid stress components, same [rr, zz, tt, rz] order.
solidStressE = { ...
    stressE.sigma_rr, ...
    stressE.sigma_zz, ...
    stressE.sigma_tt, ...
    stressE.sigma_rz};

solidStressL = { ...
    stressL.sigma_rr, ...
    stressL.sigma_zz, ...
    stressL.sigma_tt, ...
    stressL.sigma_rz};

% ---------------------------------------------------------
% Create four plots in one figure
% ---------------------------------------------------------
%figure;

tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for k = 1:4

    ax = nexttile;
    set(ax, 'FontSize', 16);
    hold(ax, 'on');

    % Plot the selected fluid total-stress component
    if is_hybrid_fluid_plot(out) && ~isempty(fluidCenter)
        selectedFluidStress = fluidSigma(:, fluidColumns(k));

        contour_masked_regular_field_select( ...
            ax, ...
            out, ...
            plotstep, ...
            fluidCenter, ...
            selectedFluidStress, ...
            1.0, ...
            48);
    end

    % Plot the same stress component in the endothelial solid
    plot_select_nodal_stress_contour( ...
        out.meshE, ...
        st.uE, ...
        solidStressE{k});

    % Plot the same stress component in the leukocyte solid
    plot_select_nodal_stress_contour( ...
        out.meshL, ...
        st.uL, ...
        solidStressL{k});

    % Draw the interfaces/boundaries
    draw_select_fluid_boundaries(ax, out, plotstep);

    colorbar(ax);

    xlabel(ax, 'r [\mum]');
    ylabel(ax, 'z [\mum]');

    title(ax, sprintf('%s, t = %.4g s', ...
        componentTitles{k}, out.t(plotstep)));

    axis(ax, 'equal');
    axis(ax, 'tight');

    rMaxE = max(out.meshE.nodes(:,1)) * 1e6;   % meters -> um
    xlim(ax, [0 rMaxE]);
    box(ax, 'on');
end

if useViscoRecovery
    sgtitle('Solid (elastic + Kelvin-Voigt viscous) and fluid (total) stress components');
else
    sgtitle('Solid (elastic ONLY -- viscous term unavailable) and fluid (total) stress components');
end

end

function contour_masked_regular_field_select( ...
    ax, out, plotstep, fluidCenter, fluidField, valueScale, nLevels)

native2D.S = fluidField(:);
native2D.R = fluidCenter(:,1);
native2D.Z = fluidCenter(:,2);

valid = isfinite(native2D.R) & ...
        isfinite(native2D.Z) & ...
        isfinite(native2D.S);
if nnz(valid) < 3
    warning('Not enough finite hybrid field points to plot.');
    return;
end

rVals = native2D.R(valid);
zVals = native2D.Z(valid);
fVals = valueScale * native2D.S(valid);

nrPlot = 260;
nzPlot = 260;
rGrid = linspace(min(rVals), max(rVals), nrPlot);
zGrid = linspace(min(zVals), max(zVals), nzPlot);
[RGrid, ZGrid] = meshgrid(rGrid, zGrid);

F = scatteredInterpolant(rVals(:), zVals(:), fVals(:), 'linear', 'none');
FGrid = F(RGrid, ZGrid);

statePlot = state_for_plot_at_step(out, plotstep);
solidMask = deformed_solids_mask_on_grid(out, statePlot, RGrid, ZGrid);
FGrid(solidMask) = NaN;

finiteVals = FGrid(isfinite(FGrid));
if isempty(finiteVals)
    warning('No finite hybrid field values remain after solid masking.');
    return;
end
if max(finiteVals) > min(finiteVals)
    contourf(ax, RGrid*1e6, ZGrid*1e6, FGrid, nLevels, 'LineColor', 'none');
else
    contourf(ax, RGrid*1e6, ZGrid*1e6, FGrid, 1, 'LineColor', 'none');
end
end

function solidMask = deformed_solids_mask_on_grid(out, statePlot, RGrid, ZGrid)
solidMask = false(size(RGrid));
if isfield(out, 'meshL') && isfield(statePlot, 'uL') && ...
        ~isempty(out.meshL) && ~isempty(statePlot.uL)
    solidMask = solidMask | local_deformed_solid_mask( ...
        out.meshL, statePlot.uL, RGrid, ZGrid);
end
if isfield(out, 'meshE') && isfield(statePlot, 'uE') && ...
        ~isempty(out.meshE) && ~isempty(statePlot.uE)
    solidMask = solidMask | local_deformed_solid_mask( ...
        out.meshE, statePlot.uE, RGrid, ZGrid);
end
end

function solidMask = local_deformed_solid_mask(mesh, u, RGrid, ZGrid)
solidMask = false(size(RGrid));
if isempty(mesh) || isempty(u) || ~isfield(mesh, 'nodes') || ...
        ~isfield(mesh, 'conn') || numel(u) < 2*size(mesh.nodes,1)
    return;
end

rDef = mesh.nodes(:,1) + u(1:2:end);
zDef = mesh.nodes(:,2) + u(2:2:end);
for e = 1:size(mesh.conn,1)
    ids = mesh.conn(e,:);
    rv = rDef(ids);
    zv = zDef(ids);
    inBox = RGrid >= min(rv) & RGrid <= max(rv) & ...
        ZGrid >= min(zv) & ZGrid <= max(zv);
    if any(inBox(:))
        localMask = false(size(solidMask));
        localMask(inBox) = inpolygon(RGrid(inBox), ZGrid(inBox), rv, zv);
        solidMask = solidMask | localMask;
    end
end
end

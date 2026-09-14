function plot_select_native2d_stress(out, plotstep)
% PLOT_SELECT_NATIVE2D_STRESS
% Per PI request: instead of the single averaged hydrostatic-pressure
% plot (see plot_select_native2d_stress_backup.m), this now plots the
% four Cauchy stress components sigma_rr, sigma_zz, sigma_tt
% (theta-theta), sigma_rz for the fluid and both solids together, one
% component per tile, in a single 2x2 figure.
%
% NOTE ON WHAT "STRESS" MEANS HERE:
%   - Fluid: sigmaCellNode/sigmaCell is the full Newtonian Cauchy stress
%     (pressure + viscous), i.e. a true total stress.
%   - Solids: recover_nodal_stress_axisym_viscoelastic.m evaluates the
%     neo-Hookean elastic Cauchy stress from the current deformation
%     gradient F, PLUS the Kelvin-Voigt viscous Cauchy stress
%     (par.etaE/par.etaL), built from Fdot using the previous accepted
%     step (uEPrev/uLPrev) and the actual accepted dt (out.dtHist). So
%     the solid panels are now a fair total-stress comparison against
%     the fluid's total stress at the interface.
%
%     UPDATED: this function used to call the elastic-only
%     recover_nodal_stress_axisym.m, which omitted that viscous term and
%     made any solid-vs-fluid interface comparison unfair -- especially
%     at step 1, right after starting from a prestressed initial
%     condition, when strain rates (and the viscous stress they
%     produce) can be large. Same fix as already applied in
%     check_interface_traction_continuity.m. Falls back to the
%     elastic-only recovery (with a warning) if uEPrev/uLPrev or
%     out.dtHist are not available for this step.
%
% FIX vs. the original draft of this function: the leukocyte's stress
% must be recovered with its OWN shear/bulk modulus (GL/KL, from
% EL/nuL = 200 Pa here), not the endothelium's (out.par.Ge/Ke, from
% Ee/nuE = 500 Pa). out.par only stores the endothelium-tagged
% parameters used by the monolithic solve; the leukocyte-specific parL
% used internally by the solver is not saved to out, so it is rebuilt
% below from the GL/KL fields that ARE kept on out.par. Without this,
% the leukocyte panels (which cover most of the plotted domain) are
% overstated by roughly Ee/EL = 2.5x.

fluid = out.fluidHist{plotstep};

% ---------------------------------------------------------
% Locate the recovered fluid total-stress array and its matching
% cell-center coordinates. sigmaCellNode/centerNode and sigmaCell/
% cellCenter are always produced as pairs by the solver, so keep them
% paired here rather than mixing sources.
%
% Column order (see recover_fluid_nodes_pressure_stress_Q4.m):
%   column 1 = sigma_rr
%   column 2 = sigma_tt
%   column 3 = sigma_zz
%   column 4 = sigma_rz
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
% Recover solid stresses once (see FIX note above for why the leukocyte
% needs its own parL instead of out.par). Uses the viscoelastic
% recovery (elastic + Kelvin-Voigt) by default so the solid stress is
% directly comparable to the fluid's total stress at the interface;
% falls back to elastic-only if the previous-step state or the accepted
% dt for this step aren't available.
% ---------------------------------------------------------
st = out.stateHist{plotstep};

parL = out.par;
if isfield(out.par, 'GL') && isfinite(out.par.GL)
    parL.Ge = out.par.GL;
end
if isfield(out.par, 'KL') && isfinite(out.par.KL)
    parL.Ke = out.par.KL;
end

useViscoRecovery = isfield(st, 'uEPrev') && isfield(st, 'uLPrev') && ...
    ~isempty(st.uEPrev) && ~isempty(st.uLPrev) && ...
    isfield(out, 'dtHist') && numel(out.dtHist) >= plotstep && ...
    isfinite(out.dtHist(plotstep));

if useViscoRecovery
    dtStep = out.dtHist(plotstep);
    if isfield(out.par, 'etaL') && isfinite(out.par.etaL)
        parL.etaE = out.par.etaL;
    end

    stressE = recover_nodal_stress_axisym_viscoelastic( ...
        out.meshE, st.uE, st.uEPrev, dtStep, out.par);

    stressL = recover_nodal_stress_axisym_viscoelastic( ...
        out.meshL, st.uL, st.uLPrev, dtStep, parL);
else
    warning(['plot_select_native2d_stress: uEPrev/uLPrev or out.dtHist ', ...
        'not available for step %d; falling back to the elastic-only ', ...
        'solid stress. The viscous Kelvin-Voigt term is omitted, so the ', ...
        'solid-vs-fluid interface comparison in this plot will be unfair.'], ...
        plotstep);

    stressE = recover_nodal_stress_axisym( ...
        out.meshE, st.uE, out.par);

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
figure;

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

    % than the 4 um gap region -- show the whole thing (0 to the mesh's
    % true outer radius) so the far-field zero-stress boundary condition
    % is actually visible, instead of truncating right where the
    % interesting near-field behavior is and hiding whether stress ever
    % decays to zero.
    rMaxE = max(out.meshE.nodes(:,1)) * 1e6;   % meters -> um, matching RGrid*1e6/ZGrid*1e6 used for the actual contour data above
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
% Plot the hybrid field on a regular r-z image grid and mask the deformed
% solids. This avoids drawing artificial curvilinear cells that connect the
% 1D gap strip directly to the 2D exterior reservoir across solid caps.

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

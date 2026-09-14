function plot_select_native2d_velocity(out,plotstep)
    native2D = out.native2D;
    if ~isfield(native2D, 'uz') || isempty(native2D.uz) || ...
            ~any(isfinite(native2D.uz(:)))
        warning('No finite native 2D velocity field is available to plot.');
        return;
    end

    figure;
    ax = gca;
    set(ax, 'FontSize', 22);
    if is_hybrid_fluid_plot(out)
        contour_masked_regular_field_select(ax, out, plotstep, 1e6, 48);
    else
        contourf(ax, native2D.R*1e6, native2D.Z*1e6, native2D.uz*1e6, ...
            48, 'LineColor', 'none');
    end
    colorbar;
    hold(ax, 'on');
    draw_select_solid_blanks(ax, out,plotstep);
    draw_select_fluid_boundaries(ax, out,plotstep);
    xlabel(ax, 'r [\mum]');
    ylabel(ax, 'z [\mum]');
    title(ax, sprintf('Hybrid axial velocity u_z(r,z), t = %.4g s', ...
        out.t(plotstep)));
    axis(ax, 'equal');
    axis(ax, 'tight');
    box(ax, 'on');
end

function contour_masked_regular_field_select(ax, out, plotstep, valueScale, nLevels)
% Plot the hybrid field on a regular r-z image grid and mask the deformed
% solids. This avoids drawing artificial curvilinear cells that connect the
% 1D gap strip directly to the 2D exterior reservoir across solid caps.

%native2D.P = out.PHist(:,:,plotstep);
native2D.R = out.RPHist(:,:,plotstep);
native2D.Z = out.ZPHist(:,:,plotstep);
% native2D.ur = out.urCHist(:,:,k);
native2D.uz = out.uzCHist(:,:,plotstep);
% native2D.speed = out.speedCHist(:,:,k);
    valid = isfinite(native2D.R) & isfinite(native2D.Z) & isfinite(native2D.uz);
    if nnz(valid) < 3
        warning('Not enough finite hybrid field points to plot.');
        return;
    end

    rVals = native2D.R(valid);
    zVals = native2D.Z(valid);
    fVals = valueScale * native2D.uz(valid);

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



function q = solid_geometry_quality(mesh, u, label)
    if ~isfield(mesh, 'axisymCache')
        mesh = prepare_axisym_mesh_cache(mesh);
    end
    cache = mesh.axisymCache;

    q = struct();
    q.label = label;
    q.minJ = inf;
    q.minRadius = inf;
    q.minJElement = NaN;
    q.minJGaussPoint = NaN;
    q.minRadiusElement = NaN;
    q.minRadiusGaussPoint = NaN;
    q.ok = true;

    for e = 1:cache.nelem
        dofs = cache.dofs(e,:).';
        ue = u(dofs);
        Rnod = cache.Rnod(:,e);
        rnod = Rnod + ue(1:2:end);

        for g = 1:cache.ngp
            N = cache.N(:,g,e).';
            Rg = cache.Rg0(g,e);
            rg = N * rnod;

            if rg < q.minRadius
                q.minRadius = rg;
                q.minRadiusElement = e;
                q.minRadiusGaussPoint = g;
            end

            if Rg <= 0 || rg <= 0
                q.ok = false;
                q.minJ = -inf;
                q.minJElement = e;
                q.minJGaussPoint = g;
                return;
            end

            F = current_deformation_gradient_from_cached(cache, e, ue, g);
            J = det(F);
            if J < q.minJ
                q.minJ = J;
                q.minJElement = e;
                q.minJGaussPoint = g;
            end
            if ~isfinite(J) || J <= 0
                q.ok = false;
                return;
            end
        end
    end
end
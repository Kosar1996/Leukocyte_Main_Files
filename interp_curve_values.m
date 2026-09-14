function vq = interp_curve_values(zNodes, vNodes, zq)
    zNodes = zNodes(:);
    vNodes = vNodes(:);
    zq = zq(:);

    [zs, idx] = sort(zNodes);
    vs = vNodes(idx);
    [zu, ~, ic] = unique(zs);
    if numel(zu) < numel(zs)
        vs = accumarray(ic, vs, [], @mean);
    end

    if numel(zu) == 1
        vq = vs(1) * ones(size(zq));
    else
        vq = safe_interp1_same_or_resample(zu, vs, zq, 'interp_curve_values');
    end
end

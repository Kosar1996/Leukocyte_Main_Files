function [Fint, K] = assemble_finite_def_axisym(mesh, u, par)

    ndof = size(mesh.nodes,1)*2;
    Fint = zeros(ndof,1);
    useCache = isfield(mesh, 'axisymCache');
    if useCache
        cache = mesh.axisymCache;
        iK = cache.iK;
        jK = cache.jK;
        vK = zeros(size(iK));
    else
        nnzLocal = mesh.nelem * 64;
        iK = zeros(nnzLocal,1);
        jK = zeros(nnzLocal,1);
        vK = zeros(nnzLocal,1);
        ptr = 1;
    end

    for e = 1:mesh.nelem
        if useCache
            dofs = cache.dofs(e,:).';
            [fe, Ke] = finite_def_element_residual_tangent_cached( ...
                cache, e, u(dofs), par);
        else
            conn = mesh.conn(e,:);
            Xe   = mesh.nodes(conn,:);
            dofs = reshape([2*conn-1; 2*conn], [], 1);
            [fe, Ke] = finite_def_element_residual_tangent(Xe, u(dofs), mesh, par);
        end

        Fint(dofs) = Fint(dofs) + fe;
        if useCache
            loc = (64*(e-1)+1):(64*e);
        else
            [ii, jj] = ndgrid(dofs, dofs);
            loc = ptr:(ptr + 63);
            iK(loc) = ii(:);
            jK(loc) = jj(:);
            ptr = ptr + 64;
        end
        vK(loc) = Ke(:);
    end

    K = sparse(iK, jK, vK, ndof, ndof);
end

function [fe, Ke] = finite_def_element_residual_tangent_cached(cache, e, ue, par)
    fe = zeros(8,1);
    Ke = zeros(8,8);

    Rnod = cache.Rnod(:,e);
    Znod = cache.Znod(:,e);

    rnod = Rnod + ue(1:2:end);
    znod = Znod + ue(2:2:end);

    I3 = eye(3);

    % LINE-SEARCH STABILITY FIX: Use unscaled physical dt for rate evaluation
    dt_phys = par.dt;

    for g = 1:cache.ngp
        N = cache.N(:,g,e).';
        dNdX = cache.dNdX(:,:,g,e);
        detJ0 = cache.detJ0(g,e);
        Rg = cache.Rg0(g,e);
        w = cache.gw(g);

        rg = N * rnod;

        if Rg <= 0 || rg <= 0
            error('Non-positive radius encountered in finite-deformation element.');
        end

        drdR = dNdX(:,1).' * rnod;
        drdZ = dNdX(:,2).' * rnod;
        dzdR = dNdX(:,1).' * znod;
        dzdZ = dNdX(:,2).' * znod;

        F = [drdR,   0,    drdZ;
               0,   rg/Rg, 0;
             dzdR,   0,    dzdZ];

        J = det(F);
        if J <= 0
            error('Negative or zero J encountered. Element inverted.');
        end

        Finv  = F \ I3;
        FinvT = Finv.';
        B = F * F.';
        trB = B(1,1) + B(2,2) + B(3,3);
        devB = B - (trB/3)*I3;

        % KINEMATIC SCALING FIX: Neo-Hookean J^(-2/3)
        aIso = J^(-2/3);
        T = par.Ge * aIso * devB + par.Ke * (J - 1) * I3;
        P = J * T * FinvT;
        Wgp = (2*pi*Rg) * detJ0 * w;

        for a = 1:4
            dNa_dR = dNdX(a,1);
            dNa_dZ = dNdX(a,2);
            Na     = N(a);

            fe(2*a-1) = fe(2*a-1) + ...
                ( P(1,1)*dNa_dR + P(1,3)*dNa_dZ + P(2,2)*(Na/Rg) ) * Wgp;

            fe(2*a) = fe(2*a) + ...
                ( P(3,1)*dNa_dR + P(3,3)*dNa_dZ ) * Wgp;
        end

        for alpha = 1:8
            dF = local_dF_from_dof(alpha, N, dNdX, Rg);

            % TANGENT TRACE FIX: In-lined scalar product (eliminates sum(sum(...)))
            trFinv_dF = dF(1,1)*Finv(1,1) + dF(1,3)*Finv(3,1) + ...
                        dF(2,2)*Finv(2,2) + dF(3,1)*Finv(1,3) + dF(3,3)*Finv(3,3);

            dJ = J * trFinv_dF;
            dB = dF * F.' + F * dF.';
            trdB = dB(1,1) + dB(2,2) + dB(3,3);
            dDevB = dB - (trdB/3)*I3;

            % KINEMATIC SCALING FIX: Derivative of J^(-2/3)
            daIso = -(2/3) * aIso * trFinv_dF;

            dT = par.Ge * ( daIso * devB + aIso * dDevB ) ...
               + par.Ke * dJ * I3;
            dFinvT = -FinvT * dF.' * FinvT;
            dP = dJ * T * FinvT + J * dT * FinvT + J * T * dFinvT;

            for a = 1:4
                dNa_dR = dNdX(a,1);
                dNa_dZ = dNdX(a,2);
                Na     = N(a);

                Ke(2*a-1, alpha) = Ke(2*a-1, alpha) + ...
                    ( dP(1,1)*dNa_dR + dP(1,3)*dNa_dZ + dP(2,2)*(Na/Rg) ) * Wgp;

                Ke(2*a, alpha) = Ke(2*a, alpha) + ...
                    ( dP(3,1)*dNa_dR + dP(3,3)*dNa_dZ ) * Wgp;
            end
        end
    end
end

function [fe, Ke] = finite_def_element_residual_tangent(Xe, ue, mesh, par)

    fe = zeros(8,1);
    Ke = zeros(8,8);

    Rnod = Xe(:,1);
    Znod = Xe(:,2);

    rnod = Rnod + ue(1:2:end);
    znod = Znod + ue(2:2:end);

    I3 = eye(3);

    % LINE-SEARCH STABILITY FIX: Use unscaled physical dt for rate evaluation
    dt_phys = par.dt;

    for g = 1:mesh.ngp
        xi  = mesh.gp(g,1);
        eta = mesh.gp(g,2);
        w   = mesh.gw(g);

        [N, dNdxi, ~] = q4_shape(xi, eta, 1.0);
        [~, dNdX, detJ0] = jacobian_2d(Xe, dNdxi);

        Rg = N * Rnod;
        rg = N * rnod;

        if Rg <= 0 || rg <= 0
            error('Non-positive radius encountered in finite-deformation element.');
        end

        drdR = dNdX(:,1).' * rnod;
        drdZ = dNdX(:,2).' * rnod;
        dzdR = dNdX(:,1).' * znod;
        dzdZ = dNdX(:,2).' * znod;

        F = [drdR,   0,    drdZ;
               0,   rg/Rg, 0;
             dzdR,   0,    dzdZ];

        J = det(F);
        if J <= 0
            error('Negative or zero J encountered. Element inverted.');
        end

        Finv  = F \ I3;
        FinvT = Finv.';
        B = F * F.';
        trB = B(1,1) + B(2,2) + B(3,3);
        devB = B - (trB/3)*I3;

        % KINEMATIC SCALING FIX: Neo-Hookean J^(-2/3)
        aIso = J^(-2/3);

        T = par.Ge * aIso * devB + par.Ke * (J - 1) * I3;
        P = J * T * FinvT;
        Wgp = (2*pi*Rg) * detJ0 * w;

        for a = 1:4
            dNa_dR = dNdX(a,1);
            dNa_dZ = dNdX(a,2);
            Na     = N(a);

            fe(2*a-1) = fe(2*a-1) + ...
                ( P(1,1)*dNa_dR + P(1,3)*dNa_dZ + P(2,2)*(Na/Rg) ) * Wgp;

            fe(2*a) = fe(2*a) + ...
                ( P(3,1)*dNa_dR + P(3,3)*dNa_dZ ) * Wgp;
        end

        for alpha = 1:8
            dF = local_dF_from_dof(alpha, N, dNdX, Rg);

            % TANGENT TRACE FIX: In-lined scalar product (eliminates sum(sum(...)))
            trFinv_dF = dF(1,1)*Finv(1,1) + dF(1,3)*Finv(3,1) + ...
                        dF(2,2)*Finv(2,2) + dF(3,1)*Finv(1,3) + dF(3,3)*Finv(3,3);

            dJ = J * trFinv_dF;

            dB = dF * F.' + F * dF.';
            trdB = dB(1,1) + dB(2,2) + dB(3,3);
            dDevB = dB - (trdB/3)*I3;

            % KINEMATIC SCALING FIX: Derivative of J^(-2/3)
            daIso = -(2/3) * aIso * trFinv_dF;

            dT = par.Ge * ( daIso * devB + aIso * dDevB ) ...
               + par.Ke * dJ * I3;

            dFinvT = -FinvT * dF.' * FinvT;
            dP = dJ * T * FinvT + J * dT * FinvT + J * T * dFinvT;

            for a = 1:4
                dNa_dR = dNdX(a,1);
                dNa_dZ = dNdX(a,2);
                Na     = N(a);

                Ke(2*a-1, alpha) = Ke(2*a-1, alpha) + ...
                    ( dP(1,1)*dNa_dR + dP(1,3)*dNa_dZ + dP(2,2)*(Na/Rg) ) * Wgp;

                Ke(2*a, alpha) = Ke(2*a, alpha) + ...
                    ( dP(3,1)*dNa_dR + dP(3,3)*dNa_dZ ) * Wgp;
            end
        end
    end
end

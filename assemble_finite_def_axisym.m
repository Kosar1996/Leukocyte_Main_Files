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

        aIso = J^(-5/3);
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

            trFinv_dF = sum(sum(Finv.' .* dF));
            dJ = J * trFinv_dF;
            dB = dF * F.' + F * dF.';
            trdB = dB(1,1) + dB(2,2) + dB(3,3);
            dDevB = dB - (trdB/3)*I3;
            daIso = -(5/3) * aIso * trFinv_dF;

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
% Consistent analytical tangent for the axisymmetric finite-deformation Q4 element.
% This replaces the old finite-difference tangent.
%
% Unknown ordering in ue:
%   ue = [u_r1; u_z1; u_r2; u_z2; u_r3; u_z3; u_r4; u_z4]

    fe = zeros(8,1);
    Ke = zeros(8,8);

    Rnod = Xe(:,1);
    Znod = Xe(:,2);

    rnod = Rnod + ue(1:2:end);
    znod = Znod + ue(2:2:end);

    I3 = eye(3);

    for g = 1:mesh.ngp
        xi  = mesh.gp(g,1);
        eta = mesh.gp(g,2);
        w   = mesh.gw(g);

        [N, dNdxi, ~] = q4_shape(xi, eta, 1.0);
        [~, dNdX, detJ0] = jacobian_2d(Xe, dNdxi);

        % Reference/current radii at GP
        Rg = N * Rnod;
        rg = N * rnod;

        if Rg <= 0 || rg <= 0
            error('Non-positive radius encountered in finite-deformation element.');
        end

        % Current deformation gradient ingredients
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

        aIso = J^(-5/3);

        % Cauchy stress
        T = par.Ge * aIso * devB + par.Ke * (J - 1) * I3;

        % First Piola
        P = J * T * FinvT;

        % Constant GP weight in reference configuration
        Wgp = (2*pi*Rg) * detJ0 * w;

        % ---- residual contribution ----
        for a = 1:4
            dNa_dR = dNdX(a,1);
            dNa_dZ = dNdX(a,2);
            Na     = N(a);

            fe(2*a-1) = fe(2*a-1) + ...
                ( P(1,1)*dNa_dR + P(1,3)*dNa_dZ + P(2,2)*(Na/Rg) ) * Wgp;

            fe(2*a) = fe(2*a) + ...
                ( P(3,1)*dNa_dR + P(3,3)*dNa_dZ ) * Wgp;
        end

        % ---- consistent analytical tangent ----
        for alpha = 1:8
            dF = local_dF_from_dof(alpha, N, dNdX, Rg);

            % variations of kinematics
                trFinv_dF = sum(sum(Finv.' .* dF));
                dJ = J * trFinv_dF;

                dB = dF * F.' + F * dF.';
                trdB = dB(1,1) + dB(2,2) + dB(3,3);
                dDevB = dB - (trdB/3)*I3;

                daIso = -(5/3) * aIso * trFinv_dF;

            % variation of Cauchy stress
            dT = par.Ge * ( daIso * devB + aIso * dDevB ) ...
               + par.Ke * dJ * I3;

            % variation of F^{-T}
            dFinvT = -FinvT * dF.' * FinvT;

            % variation of First Piola
            dP = dJ * T * FinvT + J * dT * FinvT + J * T * dFinvT;

            % assemble tangent column alpha
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

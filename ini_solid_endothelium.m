clc;clear;close all;
% =========================
par = struct();

par.REin     = 1.5e-6;       % endothelium inner radius [m]
par.REout    = 15e-6;       % endothelium outer radius [m]

par.Rc      = 0.5e-6;

par.zMin    = 0.0;          % [m]
par.zMax    = 4e-6;        % [m]
par.Lz      = par.zMax - par.zMin;

par.NrE    = 30;            % radial elements, endothelium
par.NzSolid = 41;           % axial resolution for solid meshes

par.Ee  = 500;      % Young's modulus [Pa]
par.nuE = 0.46;      % Poisson ratio

par.Ge = par.Ee/(2*(1+par.nuE));
par.Ke = par.Ee/(3*(1-2*par.nuE));

par.pIn     = 0;          % inlet pressure [Pa]
par.pOut    = 0;            % outlet pressure [Pa]

par.supportE = 'roller';

par.newtonMaxItSolid = 100;
par.lineSearchMax    = 50;
par.newtonTolSolid = 1e-9;
par.nLoadStepsSolid = 40;

% building the finite-element meshes for the two solid domains
meshE = build_rounded_endothelium_mesh(par);
interfaceE = find_interface_nodes(meshE, 'inner');  % vector of node indices of surface facing the fluid gap receives fluid traction
baseE  = find_interface_nodes(meshE, 'outer');      % substrate/anchored side gets boundary conditions (fixed or supported)

par.NzFluid = par.NzSolid;
% Axial fluid grid (shared with interface interpolation locations)
z = linspace(par.zMin, par.zMax, par.NzFluid).';

state = initial_state(z, par, meshE);

P0 = 400;                      % Pa
tau0 = 0;                      % Pa

zLoad = z;   % use the existing axial grid

P_edge = 400;      % highest pressure at two ends [Pa]
P_mid  = 100;      % lowest pressure at z = 2 um [Pa]

zCenter = 2e-6;    % location of minimum pressure [m]

d = zLoad - zCenter;
dMax = max(abs(d));

loadProfile = (d ./ dMax).^2;

% trE.normal  = P_mid + (P_edge - P_mid) * loadProfile;

trE.normal  = P0 * ones(size(zLoad));
trE.tangent = tau0 * zeros(size(zLoad));

uEnew = solve_finite_def_endothelium(meshE, state.uE, trE, interfaceE, baseE, par);

% postprocess interface shape
deltaE = extract_interface_radius(meshE, uEnew, interfaceE, zLoad);
uzE    = extract_interface_axial_displacement(meshE, uEnew, interfaceE, zLoad);

figure;
plot_deformed_axisym_mesh(meshE, uEnew, 1);   % scale = 1 means exact deformed position
xlabel('r [\mum]');
ylabel('z [\mum]');
axis equal;
set(gca,'FontSize',14);
box on;
xlim([0 15])

% plot radial position of inner surface
figure;
plot(zLoad*1e6, deltaE*1e6, 'LineWidth', 2);
xlabel('z [\mum]');
ylabel('deformed inner radius [\mum]');
title('Endothelium inner surface');

surfE = compute_surface_traction(meshE, uEnew, interfaceE, 'inner', par);

% Recover nodal stresses on endothelium
stressE = recover_nodal_stress_axisym(meshE, uEnew, par);

% Plot contour of nodal sigma_rr
plot_nodal_stress_contour(meshE, uEnew, stressE.sigma_rr, ...
    '\sigma_{rr} [Pa]', 'Endothelium nodal stress contour: \sigma_{rr}');

% Plot contour of nodal sigma_zz
plot_nodal_stress_contour(meshE, uEnew, stressE.sigma_zz, ...
    '\sigma_{zz} [Pa]', 'Endothelium nodal stress contour: \sigma_{zz}');

% Plot contour of nodal sigma_rz
plot_nodal_stress_contour(meshE, uEnew, stressE.sigma_rz, ...
    '\sigma_{rz} [Pa]', 'Endothelium nodal stress contour: \sigma_{rz}');

% Plot contour of nodal von Mises stress
plot_nodal_stress_contour(meshE, uEnew, stressE.vonMises, ...
    '\sigma_{vM} [Pa]', 'Endothelium nodal stress contour: von Mises');

% Plot contour of nodal sigma_thetatheta
plot_nodal_stress_contour(meshE, uEnew, stressE.sigma_tt, ...
    '\sigma_{\theta\theta} [Pa]', ...
    'Endothelium nodal stress contour: \sigma_{\theta\theta}');


uE_pre = uEnew;
deltaE_pre = extract_interface_radius(meshE, uE_pre, interfaceE, z);

% save('prestress_IC.mat', 'par', 'uE_pre', 'deltaE_pre', 'z', 'meshE', 'interfaceE', 'baseE' );

function state = initial_state(z, par, meshE)
    state = struct();
    state.UwL = zeros(size(z));
    state.UwE = zeros(size(z));

    state.deltaL = zeros(size(z));
    state.deltaE = rounded_gap_profile(z, par.REin, par.Lz, par.Rc);

    state.p = linspace(par.pIn, par.pOut, numel(z)).';

    state.uE = zeros(size(meshE.nodes,1)*2,1);
end


function h = rounded_gap_profile(z, H0, L, Rc)
    h = H0 * ones(size(z));

    % left rounded corner: z in [0, Rc]
    idxL = z < Rc;
    xiL  = Rc - z(idxL);   % xiL in [0, Rc]
    riseL = Rc - sqrt(Rc^2 - xiL.^2);

    % right rounded corner: z in [L-Rc, L]
    idxR = z > (L - Rc);
    xiR  = z(idxR) - (L - Rc);   % xiR in [0, Rc]
    riseR = Rc - sqrt(Rc^2 - xiR.^2);

    h(idxL) = H0 + riseL;
    h(idxR) = H0 + riseR;
end


function ids = find_interface_nodes(mesh, whichSide)

    % infer structured-grid dimensions
    zvals = unique(mesh.nodes(:,2));
    nz = numel(zvals);
    nnode = size(mesh.nodes,1);
    nr = nnode / nz;

    if abs(nr - round(nr)) > 1e-12
        error('Cannot infer structured mesh dimensions.');
    end

    nr = round(nr);

    switch lower(whichSide)
        case 'inner'
            i = 1;      % first radial column
        case 'outer'
            i = nr;     % last radial column
        otherwise
            error('unknown side');
    end

    j = (1:nz).';
    ids = sub2ind([nz,nr], j, i*ones(nz,1));
end


function uzInt = extract_interface_axial_displacement(mesh, u, interfaceNodes, zq)
    zn = mesh.nodes(interfaceNodes,2); % z-locations of interface nodes
    uz = u(2*interfaceNodes);
    [zs, idx] = sort(zn);
    uzs = uz(idx);
    uzInt = interp1(zs, uzs, zq, 'linear', 'extrap');
end


function [fixDofs, fixVals] = solid_support_conditions(baseNodes, supportType)
    switch lower(supportType)
        case 'clamped'
            fixDofs = sort([2*baseNodes(:)-1; 2*baseNodes(:)]);
            fixVals = zeros(size(fixDofs));

        case 'roller'
            fixDofs = sort(2*baseNodes(:)-1);
            n0 = baseNodes(round(end/2));
            fixDofs = sort([fixDofs; 2*n0]);
            fixVals = zeros(size(fixDofs));

        otherwise
            error('Unknown support type: %s', supportType);
    end
end

function delta = extract_interface_radius(mesh, u, interfaceNodes, zq)

    zn = mesh.nodes(interfaceNodes,2);
    rn = mesh.nodes(interfaceNodes,1);
    un = u(2*interfaceNodes-1);

    [zs, idx] = sort(zn);
    rs = rn(idx) + un(idx);

    delta = interp1(zs, rs, zq, 'linear', 'extrap');
end

function zq = traction_to_z(traction, zNodes)
    N = numel(traction.normal);
    z0 = linspace(min(zNodes), max(zNodes), N).';
    zq.normal  = interp1(z0, traction.normal,  zNodes, 'linear', 'extrap');
    zq.tangent = interp1(z0, traction.tangent, zNodes, 'linear', 'extrap');
end


function [N, dNdxi, w] = q4_shape(xi, eta, gw)
    N = 0.25*[(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
    dNdxi = 0.25*[ -(1-eta), -(1-xi);
                    +(1-eta), -(1+xi);
                    +(1+eta), +(1+xi);
                    -(1+eta), +(1-xi) ];
    w = gw;
end


function [Jm, dNdx, detJ] = jacobian_2d(Xe, dNdxi)
    Jm = Xe.' * dNdxi;
    detJ = det(Jm);
    if detJ <= 0
        error('Non-positive element Jacobian.');
    end
    dNdx = dNdxi / Jm;
end


function uNew = solve_finite_def_endothelium(mesh, uOld, traction, interfaceNodes, baseNodes, par)

    ndof = size(mesh.nodes,1)*2;
    u = uOld;

    [fixDofs, fixVals] = solid_support_conditions(baseNodes, par.supportE);
    free = setdiff((1:ndof).', unique(fixDofs(:)));
    u(fixDofs) = fixVals;

    nLoadSteps = par.nLoadStepsSolid;
    tractionFull = traction;

    for loadStep = 1:nLoadSteps
        scale = loadStep/nLoadSteps;
        traction.normal = scale * tractionFull.normal;
        traction.tangent = scale * tractionFull.tangent;
        fprintf('Load step %d/%d\n', loadStep, nLoadSteps);

        for it = 1:par.newtonMaxItSolid

            Fext = zeros(ndof,1);
            [Fext, Kext] = apply_interface_traction(mesh, u, Fext, interfaceNodes, traction);
        
            [Fint, Ktan] = assemble_finite_def_axisym(mesh, u, par);
        
            R = Fint - Fext;
            Ktot = Ktan - Kext;
        
            Rf  = R(free);
            Kff = Ktot(free, free);
            fprintf('   it=%02d rcond(Kff)=%.3e\n', it, rcond(full(Kff)));

            resNorm = norm(Rf, inf);
            refNorm = max(norm(Fext(free), inf), 1e-14);
            relRes = resNorm/refNorm;
            fprintf('   relRes=%.3e\n', relRes);
            if relRes < par.newtonTolSolid
                break;
            end

            du_free = -Kff \ Rf;

            fprintf('   max|du_free|=%.3e, min|du_free|=%.3e\n', max(abs(du_free)), min(abs(du_free)));

            alpha = 1.0;
            accepted = false;

            for ls = 1:par.lineSearchMax
                uTrial = u;
                uTrial(free) = uTrial(free) + alpha*du_free;
                uTrial(fixDofs) = fixVals;
        
                try
                    FextTrial = zeros(ndof,1);
                    [FextTrial, ~] = apply_interface_traction(mesh, uTrial, FextTrial, interfaceNodes, traction);
        
                    FintTrial = assemble_finite_def_internal_force_only(mesh, uTrial, par);
                    Rtrial = FintTrial - FextTrial;
        
                    if norm(Rtrial(free), inf) < resNorm
                        u = uTrial;
                        accepted = true;
                        break;
                    end
        
                catch ME
                    if contains(ME.message, 'Negative or zero J') || ...
                       contains(ME.message, 'Non-positive radius') || ...
                       contains(ME.message, 'Element inverted')
                        % Bad trial step: reject and reduce alpha
                    else
                        rethrow(ME);
                    end
                end
        
                alpha = 0.5 * alpha;
            end
        
            if ~accepted
                warning('Line search failed: no non-inverted residual-reducing step found.');
                uNew = u;
                return;
            end
        end

        if relRes >= par.newtonTolSolid
            warning('Finite-deformation endothelium solve hit max iterations at load step %d.', loadStep);
            uNew = u;
            return;
        end
    end

    uNew = u;
end

function [Fint, K] = assemble_finite_def_axisym(mesh, u, par)

    ndof = size(mesh.nodes,1)*2;
    Fint = zeros(ndof,1);
    K    = sparse(ndof, ndof);

    for e = 1:mesh.nelem
        conn = mesh.conn(e,:);
        Xe   = mesh.nodes(conn,:);
        dofs = reshape([2*conn-1; 2*conn], [], 1);
        ue   = u(dofs);

        [fe, Ke] = finite_def_element_residual_tangent(Xe, ue, mesh, par);

        Fint(dofs) = Fint(dofs) + fe;
        K(dofs,dofs) = K(dofs,dofs) + Ke;
    end
end

function Fint = assemble_finite_def_internal_force_only(mesh, u, par)

    ndof = size(mesh.nodes,1)*2;
    Fint = zeros(ndof,1);

    for e = 1:mesh.nelem
        conn = mesh.conn(e,:);
        Xe   = mesh.nodes(conn,:);
        dofs = reshape([2*conn-1; 2*conn], [], 1);
        ue   = u(dofs);

        fe = finite_def_element_residual_only(Xe, ue, mesh, par);
        Fint(dofs) = Fint(dofs) + fe;
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

        Finv  = inv(F);
        FinvT = Finv.';
        B = F * F.';
        trB = trace(B);
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
            dJ = J * trace(Finv * dF);

            dB = dF * F.' + F * dF.';
            trdB = trace(dB);
            dDevB = dB - (trdB/3)*I3;

            daIso = -(5/3) * aIso * trace(Finv * dF);

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


function dF = local_dF_from_dof(alpha, N, dNdX, Rg)
% Returns the exact variation dF associated with one local DOF alpha.
%
% Local DOF numbering:
%   1 = ur1, 2 = uz1, 3 = ur2, 4 = uz2, 5 = ur3, 6 = uz3, 7 = ur4, 8 = uz4

    aNode = ceil(alpha/2);
    isRadial = mod(alpha,2)==1;

    dF = zeros(3,3);

    if isRadial
        % Only r-components vary
        d_r_g   = N(aNode);
        d_drdR  = dNdX(aNode,1);
        d_drdZ  = dNdX(aNode,2);

        dF(1,1) = d_drdR;
        dF(1,3) = d_drdZ;
        dF(2,2) = d_r_g / Rg;

    else
        % Only z-components vary
        d_dzdR = dNdX(aNode,1);
        d_dzdZ = dNdX(aNode,2);

        dF(3,1) = d_dzdR;
        dF(3,3) = d_dzdZ;
    end
end

function fe = finite_def_element_residual_only(Xe, ue, mesh, par)

    fe = zeros(8,1);

    Rnod = Xe(:,1);
    Znod = Xe(:,2);

    rnod = Rnod + ue(1:2:end);
    znod = Znod + ue(2:2:end);

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

        B = F * F.';
        I = eye(3);

        % Exact constitutive law requested by user
        T = par.Ge * J^(-5/3) * ( B - (trace(B)/3)*I ) ...
          + par.Ke * (J - 1) * I;

        P = J * T / (F.');

        for a = 1:4
            dNa_dR = dNdX(a,1);
            dNa_dZ = dNdX(a,2);
            Na     = N(a);

            fe(2*a-1) = fe(2*a-1) + ...
                ( P(1,1)*dNa_dR + P(1,3)*dNa_dZ + P(2,2)*(Na/Rg) ) ...
                * (2*pi*Rg) * detJ0 * w;

            fe(2*a) = fe(2*a) + ...
                ( P(3,1)*dNa_dR + P(3,3)*dNa_dZ ) ...
                * (2*pi*Rg) * detJ0 * w;
        end
    end
end


function [F, Kext] = apply_interface_traction(mesh, u, F, interfaceNodes, traction)
% APPLY_INTERFACE_TRACTION
% Adds distributed follower traction on the deformed interface to the
% global external force vector F, and returns the consistent tangent Kext.
%
% traction.normal  : positive in local outward normal direction
% traction.tangent : positive in local tangent direction (+increasing z
%                    when interface is undeformed and ordered by z)
%
% Inputs:
%   mesh, u, F, interfaceNodes, traction
%
% Outputs:
%   F    : updated global external force vector
%   Kext : global consistent tangent of the follower traction load
%
% Notes:
% - This is a deformed-configuration line load:
%       f_e = \int N^T t ds * 2*pi*r
% - Both the line Jacobian ds and the local basis (n_hat,t_hat) depend on u.
% - The traction magnitudes tn, tt are assumed prescribed functions of z only;
%   their geometric directions follow the deformed interface.

    ndof = size(mesh.nodes,1) * 2;
    Kext = sparse(ndof, ndof);

    zn = mesh.nodes(interfaceNodes,2);
    [zs, idx] = sort(zn);
    interface = interfaceNodes(idx);

    zq = traction_to_z(traction, zs);
    tr_n = zq.normal;
    tr_t = zq.tangent;

    % 2-point Gauss rule on each line segment
    xi_gp = [-1, 1] / sqrt(3);
    w_gp  = [1, 1];

    for k = 1:numel(interface)-1
        n1 = interface(k);
        n2 = interface(k+1);

        dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2];

        % reference coordinates
        r1 = mesh.nodes(n1,1);  z1 = mesh.nodes(n1,2);
        r2 = mesh.nodes(n2,1);  z2 = mesh.nodes(n2,2);

        % current coordinates
        ur1 = u(2*n1-1); uz1 = u(2*n1);
        ur2 = u(2*n2-1); uz2 = u(2*n2);

        x1 = [r1 + ur1; z1 + uz1];
        x2 = [r2 + ur2; z2 + uz2];

        dx = x2 - x1;
        L  = norm(dx);

        if L <= 0
            continue;
        end

        t_hat = dx / L;
        n_hat = [ t_hat(2); -t_hat(1) ];

        % nodal traction magnitudes on this segment
        tn_nodes = [tr_n(k);   tr_n(k+1)];
        tt_nodes = [tr_t(k);   tr_t(k+1)];

        fe = zeros(4,1);
        ke = zeros(4,4);

        for g = 1:2
            xi = xi_gp(g);
            wg = w_gp(g);

            N1 = 0.5 * (1 - xi);
            N2 = 0.5 * (1 + xi);

            Nline = [N1, N2];

            % scalar radius at GP in current configuration
            r_gp = N1 * x1(1) + N2 * x2(1);

            % traction magnitudes interpolated from nodal values
            tn_gp = Nline * tn_nodes;
            tt_gp = Nline * tt_nodes;

            % local traction vector in current basis
            tvec = tn_gp * n_hat + tt_gp * t_hat;   % 2x1

            % shape matrix for line element
            Nmat = [N1 0  N2 0;
                    0  N1 0  N2];                  % 2x4

            Jline = L / 2;
            fac   = (2*pi*r_gp) * Jline * wg;

            % force contribution
            fe = fe + (Nmat.' * tvec) * fac;

            % ---------------------------------------------------------
            % consistent tangent of follower load
            % ---------------------------------------------------------
            %
            % x_xi = d x / d xi = 0.5*(x2-x1) = 0.5*dx
            % L = |dx|
            % t_hat = dx/L
            % n_hat = R * t_hat, with R = [0 1; -1 0]
            %
            % Variation wrt nodal displacement vector q=[ur1 uz1 ur2 uz2]^T:
            %   delta dx = Bdx * delta q
            %   Bdx = [ -1  0  1  0
            %            0 -1  0  1 ]
            %
            %   delta L = t_hat^T delta dx
            %   delta t_hat = (I - t_hat t_hat^T) delta dx / L
            %   delta n_hat = R delta t_hat
            %   delta r_gp = Br delta q, Br = [N1 0 N2 0]
            %
            % Then:
            %   delta tvec = tn_gp delta n_hat + tt_gp delta t_hat
            %   delta fac  = 2*pi*( Jline delta r_gp + r_gp delta Jline )*wg
            %   delta Jline = 0.5 delta L
            %
            %   delta f = N^T (delta tvec) fac + N^T tvec delta fac
            %

            Bdx = [-1  0  1  0;
                    0 -1  0  1];                  % 2x4

            Br  = [N1 0 N2 0];                    % 1x4

            I2 = eye(2);
            Ptan = I2 - (t_hat * t_hat.');        % projector normal to tangent
            R90 = [0 1; -1 0];

            % delta t_hat = At * delta q
            At = (Ptan / L) * Bdx;                % 2x4

            % delta n_hat = An * delta q
            An = R90 * At;                        % 2x4

            % delta Jline = AJ * delta q
            % Jline = L/2, delta L = t_hat^T delta dx
            AJ = 0.5 * (t_hat.' * Bdx);           % 1x4

            % delta fac = Afac * delta q
            Afac = 2*pi * wg * ( Jline * Br + r_gp * AJ );  % 1x4

            % delta tvec = Avec * delta q
            Avec = tn_gp * An + tt_gp * At;       % 2x4

            % consistent tangent contribution:
            % ke(:,j) = d(fe)/dq_j
            ke = ke + (Nmat.' * Avec) * fac + (Nmat.' * tvec) * Afac;
        end

        F(dofs) = F(dofs) + fe;
        Kext(dofs,dofs) = Kext(dofs,dofs) + ke;
    end
end


function surf = compute_surface_traction(mesh, u, interfaceNodes, whichSide, par)
% Computes solid-surface traction from the Cauchy stress on a deformed
% axisymmetric boundary.
%
% whichSide = 'inner' or 'outer'
%
% Outputs:
%   surf.z   : axial locations of quadrature points
%   surf.r   : radial locations of quadrature points
%   surf.tn  : normal traction scalar [Pa]
%   surf.tt  : tangential traction scalar [Pa]
%   surf.tr  : radial component of traction vector [Pa]
%   surf.tz  : axial component of traction vector [Pa]
%   surf.sigma_rr, sigma_rz, sigma_zz : Cauchy stress components [Pa]
%   surf.Fr_total, surf.Fz_total      : net force components [N]
%   surf.Fn_total, surf.Ft_total      : integrated normal/tangential forces [N]

    zn = mesh.nodes(interfaceNodes,2);
    [~, idx] = sort(zn);
    interface = interfaceNodes(idx);

    xi_gp = [-1, 1] / sqrt(3);
    w_gp  = [1, 1];

    z_all  = [];
    r_all  = [];
    tn_all = [];
    tt_all = [];
    tr_all = [];
    tz_all = [];
    srr_all = [];
    srz_all = [];
    szz_all = [];

    Fr_total = 0;
    Fz_total = 0;
    Fn_total = 0;
    Ft_total = 0;

    for k = 1:numel(interface)-1
        n1 = interface(k);
        n2 = interface(k+1);

        e = find_boundary_element(mesh, n1, n2, whichSide);
        if isempty(e)
            continue;
        end

        conn = mesh.conn(e,:);
        Xe   = mesh.nodes(conn,:);
        dofs = reshape([2*conn-1; 2*conn], [], 1);
        ue   = u(dofs);

        % current nodal positions
        xcurr = Xe;
        xcurr(:,1) = xcurr(:,1) + ue(1:2:end);
        xcurr(:,2) = xcurr(:,2) + ue(2:2:end);

        % choose boundary side in parent coordinates
        switch lower(whichSide)
            case 'inner'
                xi_fixed = -1;
            case 'outer'
                xi_fixed = +1;
            otherwise
                error('whichSide must be ''inner'' or ''outer''.');
        end

        for g = 1:2
            eta = xi_gp(g);
            wg  = w_gp(g);

            [N, dNdxi, ~] = q4_shape(xi_fixed, eta, 1.0);
            [~, dNdX, ~]  = jacobian_2d(Xe, dNdxi);

            % current position at GP
            rg = N * xcurr(:,1);
            zg = N * xcurr(:,2);

            % tangent along boundary direction eta
            dN_deta = 0.25 * [ -(1-xi_fixed);
                               -(1+xi_fixed);
                               +(1+xi_fixed);
                               +(1-xi_fixed) ];
            dx_deta = [dN_deta.' * xcurr(:,1);
                       dN_deta.' * xcurr(:,2)];

            ds_deta = norm(dx_deta);
            if ds_deta <= 0
                continue;
            end

            s_hat = dx_deta / ds_deta;

            % base normal
            n_hat_base = [ s_hat(2);
                          -s_hat(1) ];

            % outward normal from the solid
            if strcmpi(whichSide, 'inner')
                n_hat = -n_hat_base;
            else
                n_hat = n_hat_base;
            end

            % stress at this GP
            sigma = cauchy_stress_at_qp(Xe, ue, N, dNdX, par);

            srr = sigma(1,1);
            srz = sigma(1,3);
            szz = sigma(3,3);

            % traction vector in (r,z)
            tvec = [srr, srz;
                    srz, szz] * n_hat;

            % tangential unit vector
            t_hat = s_hat;

            tn = dot(tvec, n_hat);
            tt = dot(tvec, t_hat);

            % axisymmetric surface measure:
            % dA = 2*pi*r*ds = 2*pi*r*(ds/deta)*deta
            % with 2-pt Gauss on eta in [-1,1]
            dA = 2*pi*rg*ds_deta*wg;

            Fr_total = Fr_total + tvec(1) * dA;
            Fz_total = Fz_total + tvec(2) * dA;
            Fn_total = Fn_total + tn      * dA;
            Ft_total = Ft_total + tt      * dA;

            z_all   = [z_all; zg];
            r_all   = [r_all; rg];
            tn_all  = [tn_all; tn];
            tt_all  = [tt_all; tt];
            tr_all  = [tr_all; tvec(1)];
            tz_all  = [tz_all; tvec(2)];
            srr_all = [srr_all; srr];
            srz_all = [srz_all; srz];
            szz_all = [szz_all; szz];
        end
    end

    [z_all, ord] = sort(z_all);

    surf = struct();
    surf.z = z_all;
    surf.r = r_all(ord);

    surf.tn = tn_all(ord);
    surf.tt = tt_all(ord);
    surf.tr = tr_all(ord);
    surf.tz = tz_all(ord);

    surf.sigma_rr = srr_all(ord);
    surf.sigma_rz = srz_all(ord);
    surf.sigma_zz = szz_all(ord);

    surf.Fr_total = Fr_total;
    surf.Fz_total = Fz_total;
    surf.Fn_total = Fn_total;
    surf.Ft_total = Ft_total;
end


function e = find_boundary_element(mesh, n1, n2, whichSide)
% Finds the element sharing boundary segment (n1,n2) on the chosen side.

    e = [];

    for ii = 1:size(mesh.conn,1)
        conn = mesh.conn(ii,:);

        switch lower(whichSide)
            case 'inner'
                edge = conn([1 4]);
            case 'outer'
                edge = conn([2 3]);
            otherwise
                error('whichSide must be ''inner'' or ''outer''.');
        end

        if isequal(edge(:), [n1; n2]) || isequal(edge(:), [n2; n1])
            e = ii;
            return;
        end
    end
end


function sigma = cauchy_stress_at_qp(Xe, ue, N, dNdX, par)
% Returns 3x3 Cauchy stress tensor at one quadrature point.

    Rnod = Xe(:,1);
    Znod = Xe(:,2);

    rnod = Rnod + ue(1:2:end);
    znod = Znod + ue(2:2:end);

    Rg = N * Rnod;
    rg = N * rnod;

    if Rg <= 0 || rg <= 0
        error('Non-positive radius encountered while post-processing stress.');
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
        error('Negative or zero J encountered while post-processing stress.');
    end

    B = F * F.';
    I = eye(3);

    sigma = par.Ge * J^(-5/3) * ( B - (trace(B)/3)*I ) ...
          + par.Ke * (J - 1) * I;
end


function stress = recover_nodal_stress_axisym(mesh, u, par)
% Recover nodal Cauchy stresses by averaging element-center stresses
% from all connected elements.

    nnode = size(mesh.nodes,1);

    sigma_rr_sum = zeros(nnode,1);
    sigma_zz_sum = zeros(nnode,1);
    sigma_rz_sum = zeros(nnode,1);
    sigma_tt_sum = zeros(nnode,1);
    vm_sum       = zeros(nnode,1);
    count        = zeros(nnode,1);

    % Use element-center stress as representative element stress
    xi = 0.0;
    eta = 0.0;
    [N, dNdxi, ~] = q4_shape(xi, eta, 1.0);

    for e = 1:mesh.nelem
        conn = mesh.conn(e,:);
        Xe   = mesh.nodes(conn,:);
        dofs = reshape([2*conn-1; 2*conn], [], 1);
        ue   = u(dofs);

        [~, dNdX, ~] = jacobian_2d(Xe, dNdxi);
        sigma = cauchy_stress_at_qp(Xe, ue, N, dNdX, par);

        srr = sigma(1,1);
        stt = sigma(2,2);
        szz = sigma(3,3);
        srz = sigma(1,3);

        % axisymmetric 3D von Mises stress
        svm = sqrt(0.5*((srr-stt)^2 + (stt-szz)^2 + (szz-srr)^2 + 6*srz^2));

        for a = 1:4
            node = conn(a);
            sigma_rr_sum(node) = sigma_rr_sum(node) + srr;
            sigma_tt_sum(node) = sigma_tt_sum(node) + stt;
            sigma_zz_sum(node) = sigma_zz_sum(node) + szz;
            sigma_rz_sum(node) = sigma_rz_sum(node) + srz;
            vm_sum(node)       = vm_sum(node) + svm;
            count(node)        = count(node) + 1;
        end
    end

    stress = struct();
    stress.sigma_rr = sigma_rr_sum ./ max(count,1);
    stress.sigma_tt = sigma_tt_sum ./ max(count,1);
    stress.sigma_zz = sigma_zz_sum ./ max(count,1);
    stress.sigma_rz = sigma_rz_sum ./ max(count,1);
    stress.vonMises = vm_sum ./ max(count,1);
end





function plot_nodal_stress_contour(mesh, u, nodalStress, cbarLabel, plotTitle)
% Plot nodal stress contour on deformed mesh.
% Convention: x-axis = r, y-axis = z.

    nodes = mesh.nodes;
    conn  = mesh.conn;

    ur = u(1:2:end);
    uz = u(2:2:end);

    rDef = nodes(:,1) + ur;
    zDef = nodes(:,2) + uz;
    
    figure;

    patch('Faces', conn, ...
          'Vertices', [rDef*1e6, zDef*1e6], ...
          'FaceVertexCData', nodalStress, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');

    cb = colorbar;
    ylabel(cb, cbarLabel, 'Interpreter', 'tex');

    xlabel('r [\mum]');
    ylabel('z [\mum]');
    title(plotTitle, 'Interpreter', 'tex');
    xlim([0 15])
    ylim([-1 5])

    box on;
    set(gca,'FontSize',14);
end

function mesh = build_rounded_endothelium_mesh(par)

    zvec = linspace(par.zMin, par.zMax, par.NzSolid).';
    rInner = rounded_gap_profile(zvec, par.REin, par.Lz, par.Rc);
    rOuter = par.REout * ones(size(zvec));

    nodes = zeros(par.NzSolid * par.NrE, 2);

    % for uniform mesh
    % for j = 1:par.NzSolid
    %     rline = linspace(rInner(j), rOuter(j), par.NrE);
    %     for i = 1:par.NrE
    %         id = sub2ind([par.NzSolid, par.NrE], j, i);
    %         nodes(id,:) = [rline(i), zvec(j)];
    %     end
    % end


    % for finner mesh near inner wall
    s = linspace(0,1,par.NrE);
    beta = 2;          % beta > 1 gives finer mesh near inner wall
    sBias = s.^beta;

    for j = 1:par.NzSolid
        rline = rInner(j) + (rOuter(j)-rInner(j))*sBias;

        for i = 1:par.NrE
            id = sub2ind([par.NzSolid, par.NrE], j, i);
            nodes(id,:) = [rline(i), zvec(j)];
        end
    end


    conn = zeros((par.NrE-1)*(par.NzSolid-1),4);
    e = 0;

    for j = 1:par.NzSolid-1
        for i = 1:par.NrE-1
            n1 = sub2ind([par.NzSolid, par.NrE], j,   i);
            n2 = sub2ind([par.NzSolid, par.NrE], j,   i+1);
            n3 = sub2ind([par.NzSolid, par.NrE], j+1, i+1);
            n4 = sub2ind([par.NzSolid, par.NrE], j+1, i);
            e = e + 1;
            conn(e,:) = [n1 n2 n3 n4];
        end
    end

    gp1 = [-1, 1]/sqrt(3);
    gw1 = [1, 1];
    [g1,g2] = meshgrid(gp1,gp1);
    [w1,w2] = meshgrid(gw1,gw1);

    mesh = struct();
    mesh.nodes = nodes;
    mesh.conn  = conn;
    mesh.nelem = size(conn,1);
    mesh.ngp   = numel(g1);
    mesh.gp    = [g1(:), g2(:)];
    mesh.gw    = w1(:).*w2(:);
end


function plot_deformed_axisym_mesh(mesh, u, scale)

if nargin < 3
    scale = 1;
end

r = mesh.nodes(:,1);
z = mesh.nodes(:,2);

ur = u(1:2:end);
uz = u(2:2:end);

rDef = r + scale*ur;
zDef = z + scale*uz;


patch('Faces', mesh.conn, ...
      'Vertices', [rDef*1e6, zDef*1e6], ...
      'FaceColor','none','EdgeColor','k');

axis equal;
xlabel('r [\mum]');
ylabel('z [\mum]');
box on;

end

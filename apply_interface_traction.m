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

    if isfield(traction, 'z')
        zn = mesh.nodes(interfaceNodes,2) + u(2*interfaceNodes(:));
    else
        zn = mesh.nodes(interfaceNodes,2);
    end
    [zs, idx] = sort(zn);
    interface = interfaceNodes(idx);

    if isfield(traction, 'z')
        tr_n = interp_curve_values(traction.z, traction.normal, zs);
        tr_t = interp_curve_values(traction.z, traction.tangent, zs);
    else
        zq = traction_to_z(traction, zs);
        tr_n = zq.normal;
        tr_t = zq.tangent;
    end
    [tr_n, tr_t] = apply_traction_support_window(traction, zs, tr_n, tr_t);

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

function zq = traction_to_z(traction, zNodes)
    N = numel(traction.normal);
    z0 = linspace(min(zNodes), max(zNodes), N).';
    zq.normal  = safe_interp1_same_or_resample(z0, traction.normal,  zNodes, 'traction.normal');
    zq.tangent = safe_interp1_same_or_resample(z0, traction.tangent, zNodes, 'traction.tangent');
end

function [tr_n, tr_t] = apply_traction_support_window(traction, zs, tr_n, tr_t)
    useWindow = isfield(traction, 'outsideZero') && traction.outsideZero;
    if ~useWindow
        return;
    end

    if isfield(traction, 'supportInterval') && ...
            numel(traction.supportInterval) == 2
        support = sort(traction.supportInterval(:));
    elseif isfield(traction, 'z') && ~isempty(traction.z)
        support = [min(traction.z(:)); max(traction.z(:))];
    elseif isempty(zs)
        support = [0; 0];
    else
        support = [min(zs(:)); max(zs(:))];
    end

    outside = zs(:) < support(1) | zs(:) > support(2);
    tr_n(outside) = 0;
    tr_t(outside) = 0;
end
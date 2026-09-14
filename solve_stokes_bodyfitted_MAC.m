%% ============================================================
% Body-fitted MAC Stokes solver
% ============================================================
function fluid = solve_stokes_bodyfitted_MAC(mesh, state, bc, par, tNow)

    Nr = mesh.Nr;
    Nz = mesh.Nz;
    dz = mesh.dz;
    mu = par.mu;

    % Optional unsteady-Stokes term: rho*du/dt added to the momentum
    % equations (backward Euler in the already-available par.dt), instead
    % of the pure steady Stokes balance -mu*lap(u) + grad(p) = 0 used by
    % default. Off unless explicitly enabled -- when off, massCoef=0 and
    % every line below that uses it is a no-op, so default behavior is
    % byte-for-byte unchanged. rho is not otherwise a codebase parameter
    % (this solver has no density anywhere else); default 1000 kg/m^3
    % (water/plasma-like, consistent with the existing mu=1.2e-3 Pa*s) if
    % not explicitly set on par.
    useUnsteady = isfield(par, 'useUnsteadyStokes') && par.useUnsteadyStokes && ...
        isfield(par, 'dt') && isfinite(par.dt) && par.dt > 0;
    massCoef = 0;
    if useUnsteady
        rho = 1000;
        if isfield(par, 'rho') && isfinite(par.rho) && par.rho > 0
            rho = par.rho;
        end
        massCoef = rho / par.dt;
    end
    haveUrPrev = useUnsteady && isfield(bc, 'urPrev') && isequal(size(bc.urPrev), [Nr+1, Nz]);
    haveUzPrev = useUnsteady && isfield(bc, 'uzPrev') && isequal(size(bc.uzPrev), [Nr, Nz+1]);
    % CAUTION for future callers: bc.urPrev/uzPrev are currently only
    % populated by solve_fluid_2D_bodyfitted_MAC.m (Pure/Full 2D MAC mode).
    % solve_fluid_hybrid_gap1d_exterior2d.m and
    % solve_fluid_hybrid_gap1d_exterior2d_withnodes.m call this function
    % directly and do not set them. Enabling par.useUnsteadyStokes for
    % Hybrid mode without also wiring that plumbing through will NOT error
    % -- it silently falls back to treating every step as starting from
    % rest (haveUrPrev/haveUzPrev false), which is safe but not a correct
    % time-history-aware unsteady solve. Wire the same bc.urPrev/uzPrev
    % pattern into those two callers before relying on this flag there.

    % Unknown numbering
    pId = reshape(1:Nr*Nz, Nr, Nz);
    Np = Nr*Nz;

    % ur unknowns: interior radial faces only, i=2:Nr
    urId = zeros(Nr+1,Nz);
    urList = [];
    count = 0;
    for j = 1:Nz
        for i = 2:Nr
            count = count + 1;
            urId(i,j) = count;
            urList(end+1,1) = sub2ind([Nr+1,Nz],i,j); %#ok<AGROW>
        end
    end
    Nur = count;

    % uz unknowns: all axial faces
    uzId = reshape(1:(Nr*(Nz+1)), Nr, Nz+1);
    Nuz = Nr*(Nz+1);

    Nu = Nur + Nuz;
    Ntot = Nu + Np;

    ii = [];
    jj = [];
    vv = [];
    rhs = zeros(Ntot,1);

    function add(row,col,val)
        if col ~= 0 && val ~= 0 && isfinite(val)
            ii(end+1,1) = row; %#ok<AGROW>
            jj(end+1,1) = col; %#ok<AGROW>
            vv(end+1,1) = val; %#ok<AGROW>
        end
    end

    function val = ur_known(i,j)
        if i == 1
            val = bc.urL(j);
        elseif i == Nr+1
            val = bc.urE(j);
        else
            val = 0;
        end
    end

    function val = uz_wall(side,jFace)
        zq = mesh.zF(jFace);
        switch side
            case 'L'
                val = interp1(state.zc, bc.uzL, zq, 'linear','extrap');
            case 'E'
                val = interp1(state.zc, bc.uzE, zq, 'linear','extrap');
        end
    end

    %% -------------------------------
    % Radial momentum, ur on interior radial faces
    % -------------------------------
    for a = 1:Nur
        lin = urList(a);
        [i,j] = ind2sub([Nr+1,Nz], lin);

        row = a;
        rFaces = mesh.Rur(:,j);
        r = max(rFaces(i), 1e-30);
        drM = max(rFaces(i) - rFaces(i-1), 1e-30);
        drP = max(rFaces(i+1) - rFaces(i), 1e-30);
        drCV = max(0.5 * (drM + drP), 1e-30);

        % Cylindrical radial diffusion coefficients at radial face
        rM = max(0.5 * (rFaces(i-1) + rFaces(i)), 1e-30);
        rP = max(0.5 * (rFaces(i) + rFaces(i+1)), 1e-30);
        cM = rM/(r*drM*drCV);
        cP = rP/(r*drP*drCV);
        cZ = 1/dz^2;

        center = cM + cP + 2*cZ + 1/r^2;

        add(row,row,mu*center);

        if massCoef > 0
            add(row,row,massCoef);
            if haveUrPrev
                rhs(row) = rhs(row) + massCoef * bc.urPrev(i,j);
            end
        end

        % radial minus neighbor
        if i-1 >= 2
            add(row, urId(i-1,j), -mu*cM);
        else
            rhs(row) = rhs(row) + mu*cM*ur_known(1,j);
        end

        % radial plus neighbor
        if i+1 <= Nr
            add(row, urId(i+1,j), -mu*cP);
        else
            rhs(row) = rhs(row) + mu*cP*ur_known(Nr+1,j);
        end

        % axial neighbors: zero-gradient at open z ends
        if j > 1
            add(row, urId(i,j-1), -mu*cZ);
        else
            center = center - cZ; %#ok<NASGU>
            add(row,row,-mu*cZ);
        end

        if j < Nz
            add(row, urId(i,j+1), -mu*cZ);
        else
            add(row,row,-mu*cZ);
        end

        % Pressure gradient dp/dr across the radial face
        idL = pId(i-1,j);
        idR = pId(i,j);
        drPCells = max(mesh.Rp(i,j) - mesh.Rp(i-1,j), 1e-30);
        add(row, Nu + idR,  1/drPCells);
        add(row, Nu + idL, -1/drPCells);
    end

    %% -------------------------------
    % Axial momentum, uz at axial faces
    % -------------------------------
    for j = 1:Nz+1
        for i = 1:Nr
            localUz = uzId(i,j);
            row = Nur + localUz;

            rCells = mesh.Ruz(:,j);
            rFaces = mesh.Rzf(:,j);
            r = max(rCells(i), 1e-30);

            % Use local nonuniform radial spacing at the axial face.
            drCell = max(rFaces(i+1) - rFaces(i), 1e-30);
            if i > 1
                drM = max(rCells(i) - rCells(i-1), 1e-30);
            else
                drM = max(rCells(i) - rFaces(i), 1e-30);
            end
            if i < Nr
                drP = max(rCells(i+1) - rCells(i), 1e-30);
            else
                drP = max(rFaces(i+1) - rCells(i), 1e-30);
            end

            rM = max(rFaces(i), 1e-30);
            rP = max(rFaces(i+1), 1e-30);
            cM = rM/(r*drM*drCell);
            cP = rP/(r*drP*drCell);
            cZ = 1/dz^2;

            center = cM + cP;

            % radial minus neighbor or inner wall no-slip
            if i > 1
                add(row, Nur + uzId(i-1,j), -mu*cM);
            else
                rhs(row) = rhs(row) + mu*cM*uz_wall('L',j);
            end

            % radial plus neighbor or outer wall no-slip
            if i < Nr
                add(row, Nur + uzId(i+1,j), -mu*cP);
            else
                rhs(row) = rhs(row) + mu*cP*uz_wall('E',j);
            end

            % axial diffusion neighbors; zero-gradient at open ends
            if j > 1
                add(row, Nur + uzId(i,j-1), -mu*cZ);
                center = center + cZ;
            end
            if j < Nz+1
                add(row, Nur + uzId(i,j+1), -mu*cZ);
                center = center + cZ;
            end

            add(row, Nur + localUz, mu*center);

            if massCoef > 0
                add(row, Nur + localUz, massCoef);
                if haveUzPrev
                    rhs(row) = rhs(row) + massCoef * bc.uzPrev(i,j);
                end
            end

            % Pressure gradient dp/dz
            if j == 1
                pIn = prescribed_pressure(par,'in',mesh.Ruz(i,j),tNow);
                idT = pId(i,1);
                add(row, Nu + idT,  2/dz);
                rhs(row) = rhs(row) + (2/dz)*pIn;
            elseif j == Nz+1
                pOut = prescribed_pressure(par,'out',mesh.Ruz(i,j),tNow);
                idB = pId(i,Nz);
                add(row, Nu + idB, -2/dz);
                rhs(row) = rhs(row) - (2/dz)*pOut;
            else
                idB = pId(i,j-1);
                idT = pId(i,j);
                add(row, Nu + idT,  1/dz);
                add(row, Nu + idB, -1/dz);
            end
        end
    end

    %% -------------------------------
    % Continuity / pressure equations
    % -------------------------------
    for j = 1:Nz
        for i = 1:Nr
            pLocal = pId(i,j);
            row = Nu + pLocal;

            % Conservative finite-volume continuity per radian:
            %
            %   Fr_R - Fr_L + Fz_T - Fz_B = 0
            %
            % where
            %   Fr = dz * r * u_r
            %   Fz = A_z * u_z
            %   A_z = 0.5*(r_outer^2 - r_inner^2)

            rL = mesh.Rur(i,j);
            rR = mesh.Rur(i+1,j);

            % Cell volume per radian
            V = 0.5*(rR^2 - rL^2)*dz;
            V = max(V, 1e-30);

            % Radial flux coefficients
            cL = -dz*rL/V;
            cR =  dz*rR/V;

            if i == 1
                rhs(row) = rhs(row) - cL*bc.urL(j);
            else
                add(row, urId(i,j), cL);
            end

            if i == Nr
                rhs(row) = rhs(row) - cR*bc.urE(j);
            else
                add(row, urId(i+1,j), cR);
            end

            % Axial-face areas for this eta-cell
            rB_L = mesh.Rzf(i,  j);
            rB_R = mesh.Rzf(i+1,j);
            rT_L = mesh.Rzf(i,  j+1);
            rT_R = mesh.Rzf(i+1,j+1);

            AB = 0.5*(rB_R^2 - rB_L^2);
            AT = 0.5*(rT_R^2 - rT_L^2);

            % Axial flux coefficients
            add(row, Nur + uzId(i,j),   -AB/V);
            add(row, Nur + uzId(i,j+1),  AT/V);

            if par.pressurePenalty > 0
                add(row, Nu + pLocal, -par.pressurePenalty);
            end
        end
    end

    %% -------------------------------
    % Assemble and solve
    % -------------------------------
    A = sparse(ii,jj,vv,Ntot,Ntot);

    rowScale = 1 ./ max(sum(abs(A),2), 1);
    S = spdiags(rowScale,0,Ntot,Ntot);

    x = (S*A) \ (S*rhs);

    linRes = norm(A*x - rhs, inf) / max(norm(rhs, inf), 1);

    % Diagnostic only, default off: reports the row-scaled system's
    % estimated condition number (1-norm, via condest -- exact cond() is
    % too expensive for a matrix this size). Added Sep 3 to directly test
    % Issue 2 (fluid matrix conditioning near the gap floor) rather than
    % inferring it only through Issue 3's side effects. No effect on the
    % solve itself; purely an optional extra field on the returned struct.
    condNumber = NaN;
    if isfield(par, 'reportMatrixConditionNumber') && par.reportMatrixConditionNumber
        condNumber = condest(S*A);
    end

    %% -------------------------------
    % Recover arrays
    % -------------------------------
    P = reshape(x(Nu+1:Nu+Np), Nr, Nz);

    ur = nan(Nr+1,Nz);
    ur(1,:) = bc.urL(:).';
    ur(Nr+1,:) = bc.urE(:).';
    for a = 1:Nur
        lin = urList(a);
        [i,j] = ind2sub([Nr+1,Nz], lin);
        ur(i,j) = x(a);
    end

    uz = reshape(x(Nur+1:Nur+Nuz), Nr, Nz+1);

    urC = 0.5*(ur(1:Nr,:) + ur(2:Nr+1,:));
    uzC = 0.5*(uz(:,1:Nz) + uz(:,2:Nz+1));

    div = compute_bodyfitted_divergence(mesh, ur, uz);

    contMask = true(Nr,Nz);
    divInf = max(abs(div(contMask)),[],'omitnan');

    fluid = struct();
    fluid.P = P;
    fluid.ur = ur;
    fluid.uz = uz;
    fluid.urC = urC;
    fluid.uzC = uzC;

    fluid.div = div;
    fluid.divInf = divInf;
    fluid.linRes = linRes;
    fluid.condNumber = condNumber;
    fluid.Nunknown = Ntot;

    fluid.contMask = contMask;
    fluid.pBC = false(Nr,Nz);
end

function pval = prescribed_pressure(par, side, r, t) 

    switch lower(side)
        case 'in'
            if isfield(par,'pInFun') && ~isempty(par.pInFun)
                pval = par.pInFun(r,t);
            else
                pval = par.pIn + 0*r;
            end

        case 'out'
            if isfield(par,'pOutFun') && ~isempty(par.pOutFun)
                pval = par.pOutFun(r,t);
            else
                pval = par.pOut + 0*r;
            end

        otherwise
            error('Unknown prescribed-pressure side.');
    end
end

function div = compute_bodyfitted_divergence(mesh, ur, uz)

    Nr = mesh.Nr;
    Nz = mesh.Nz;
    dz = mesh.dz;

    div = nan(Nr,Nz);

    for j = 1:Nz
        for i = 1:Nr

            rL = mesh.Rur(i,j);
            rR = mesh.Rur(i+1,j);

            V = 0.5*(rR^2 - rL^2)*dz;
            V = max(V, 1e-30);

            Fr = dz*(rR*ur(i+1,j) - rL*ur(i,j));

            rB_L = mesh.Rzf(i,  j);
            rB_R = mesh.Rzf(i+1,j);
            rT_L = mesh.Rzf(i,  j+1);
            rT_R = mesh.Rzf(i+1,j+1);

            AB = 0.5*(rB_R^2 - rB_L^2);
            AT = 0.5*(rT_R^2 - rT_L^2);

            Fz = AT*uz(i,j+1) - AB*uz(i,j);

            div(i,j) = (Fr + Fz)/V;
        end
    end
end
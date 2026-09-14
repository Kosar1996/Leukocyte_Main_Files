function [fluid, ok, stopReason, meshF] = solve_fluid_2D_bodyfitted_MAC(z, old, state, par)
%SOLVE_FLUID_2D_BODYFITTED_MAC
% Body-fitted MAC wrapper that replaces the old Q4 penalty-Stokes post-solve
% with the body-fitted MAC cylindrical Stokes solver.
%
% The calling code is unchanged: it still requests a staged 2D fluid solve
% after each accepted solid step. Internally this now uses
%
%     r(eta,z,t) = deltaL(z,t) + eta*(deltaE(z,t)-deltaL(z,t))
%
% and pressure ghost values at the physical inlet/outlet, as in
% mac_bodyfitted_2D_stokes_rigid_leukocyte_pressureBC.m.
%
% Output fields are kept compatible with the old solver:
%   fluid.p, pL, pE, Q, tauL, tauE, uzL, uzE, gap
% plus body-fitted fields:
%   fluid.P, ur, uz, urC, uzC, meshF, divInf, linRes.

    ok = true;
    stopReason = '';
    meshF = [];

    zOut = z(:);
    NzOut = numel(zOut);

    if NzOut < 2
        ok = false;
        stopReason = 'Body-fitted MAC Stokes fluid requires at least two axial points.';
        fluid = empty_fluid_2D_return(NzOut, nan(NzOut,1));
        return;
    end

    rlOut = state.deltaL(:);
    reOut = state.deltaE(:);

    useRLoutInner = use_RLout_fluid_interface_for_solid_leukocyte(par);
    if useRLoutInner
        rlOut = par.RLout * ones(NzOut,1);
    end

    if numel(rlOut) ~= NzOut || numel(reOut) ~= NzOut
        ok = false;
        stopReason = 'Body-fitted MAC Stokes fluid: state.deltaL/deltaE size does not match z.';
        fluid = empty_fluid_2D_return(NzOut, nan(NzOut,1));
        return;
    end

    hOut = reOut - rlOut;
    if any(hOut <= par.minGap)
        [hmin, imin] = min(hOut);
        ok = false;
        stopReason = sprintf('Body-fitted MAC Stokes fluid: gap below minGap at z=%.6e, h=%.6e.', ...
            zOut(imin), hmin);
        fluid = empty_fluid_2D_return(NzOut, hOut);
        return;
    end

    % Use a cell-centered body-fitted Stokes grid over the physical axial
    % domain. The solid/interface state is interpolated from the coupling
    % grid zOut onto these fluid pressure-cell centers.
    zMin = par.zMin;
    zMax = par.zMax;
    if ~(isfinite(zMin) && isfinite(zMax) && zMax > zMin)
        zMin = zOut(1);
        zMax = zOut(end);
    end

    NzMAC = NzOut;
    zF = linspace(zMin, zMax, NzMAC+1).';
    zc = 0.5*(zF(1:end-1) + zF(2:end));

    parMAC = par;
    if isfield(par, 'NrFluid2D') && ~isempty(par.NrFluid2D)
        parMAC.Nr = par.NrFluid2D;
    elseif ~isfield(parMAC, 'Nr') || isempty(parMAC.Nr)
        parMAC.Nr = 16;
    end
    parMAC.Nz = NzMAC;
    parMAC.NrFluid2D = parMAC.Nr;

    % continuity matrix becomes severely ill-conditioned as the gap
    % narrows toward the minGap floor. solve_stokes_bodyfitted_MAC.m
    % already supports an artificial-compressibility pressure-stabilization
    % term (par.pressurePenalty > 0 adds -pressurePenalty*p to each
    % continuity row), but this wrapper unconditionally zeroed it, so no
    % run could ever exercise it. Overridable via par so default behavior
    % (0, off) is unchanged for configs that never asked for it.
    parMAC.pressurePenalty = 0;
    if isfield(par, 'fluidPressurePenalty') && isfinite(par.fluidPressurePenalty)
        parMAC.pressurePenalty = par.fluidPressurePenalty;
    end
    if ~isfield(parMAC, 'pIn') || isempty(parMAC.pIn),  parMAC.pIn = 0;  end
    if ~isfield(parMAC, 'pOut') || isempty(parMAC.pOut), parMAC.pOut = 0; end
    parMAC.pInFun  = @(r,t) 0*r + parMAC.pIn;
    parMAC.pOutFun = @(r,t) 0*r + parMAC.pOut;

    oldDeltaL = old.deltaL(:);
    oldDeltaE = old.deltaE(:);
    if useRLoutInner
        oldDeltaL = par.RLout * ones(NzOut,1);
    end

    stateMAC = struct();
    stateMAC.zc = zc;
    stateMAC.zF = zF;
    stateMAC.deltaL = safe_interp1_same_or_resample(zOut, rlOut, zc, 'state.deltaL');
    stateMAC.deltaE = safe_interp1_same_or_resample(zOut, reOut, zc, 'state.deltaE');

    oldMAC = struct();
    oldMAC.deltaL = safe_interp1_same_or_resample(zOut, oldDeltaL, zc, 'old.deltaL');
    oldMAC.deltaE = safe_interp1_same_or_resample(zOut, oldDeltaE, zc, 'old.deltaE');

    if any(stateMAC.deltaE - stateMAC.deltaL <= par.minGap)
        hMAC = stateMAC.deltaE - stateMAC.deltaL;
        [hmin, imin] = min(hMAC);
        ok = false;
        stopReason = sprintf('Body-fitted MAC Stokes fluid: interpolated gap below minGap at z=%.6e, h=%.6e.', ...
            zc(imin), hmin);
        fluid = empty_fluid_2D_return(NzOut, hOut);
        return;
    end

    if useRLoutInner
        stateMAC.UwL = zeros(NzMAC,1);
    elseif isfield(state, 'UwL') && numel(state.UwL) == NzOut
        stateMAC.UwL = safe_interp1_same_or_resample(zOut, state.UwL(:), zc, 'state.UwL');
    else
        stateMAC.UwL = zeros(NzMAC,1);
    end
    if isfield(state, 'UwE') && numel(state.UwE) == NzOut
        stateMAC.UwE = safe_interp1_same_or_resample(zOut, state.UwE(:), zc, 'state.UwE');
    else
        stateMAC.UwE = zeros(NzMAC,1);
    end

    bc = struct();
    % Axial wall velocities are carried from monolithic_interface_kinematics,
    % i.e. from solid axial displacement increments. Fall back to prescribed
    % wall speeds only when the state has no solid-derived Uw fields.
    bc.uzL = stateMAC.UwL;
    bc.uzE = stateMAC.UwE;

    graphUrL = (stateMAC.deltaL - oldMAC.deltaL) / par.dt;
    graphUrE = (stateMAC.deltaE - oldMAC.deltaE) / par.dt;
    useSlopeKinematics = ~isfield(par, 'useSlopeAwareWallKinematics') || ...
        par.useSlopeAwareWallKinematics;
    if useSlopeKinematics
        slopeL = curve_slope_1d(stateMAC.zc, stateMAC.deltaL);
        slopeE = curve_slope_1d(stateMAC.zc, stateMAC.deltaE);
        bc.urL = graphUrL + slopeL .* bc.uzL;
        bc.urE = graphUrE + slopeE .* bc.uzE;
    else
        bc.urL = graphUrL;
        bc.urE = graphUrE;
    end

    try
        meshMAC = build_body_fitted_gap_mesh(stateMAC, parMAC);
        tNow = 0;
        if isfield(state, 't') && isfinite(state.t)
            tNow = state.t;
        end

        % Optional unsteady-Stokes term (par.useUnsteadyStokes): pass the
        % previous step's MAC face velocities through as bc.urPrev/bc.uzPrev,
        % on the same (Nr+1,Nz)/(Nr,Nz+1) DOF layout the current mesh uses.
        %
        % The mesh is rebuilt every step from the current gap geometry, so
        % index (i,j) is the same LOGICAL face on both steps (same eta, same
        % z -- Nr/Nz and the axial grid zc/zF are fixed for the whole run)
        % but NOT generally the same PHYSICAL r: when the gap shape changes
        % between steps, the radial position of a given logical face moves.
        % Reusing old.ur2DFaceField(i,j) directly as "u at this physical
        % point, last step" is therefore only exact for a non-deforming
        % mesh -- not the case here, especially in the near-contact regime
        % where the gap changes fastest. Since z itself never moves, this is
        % corrected with a 1-D radial re-grid at each fixed z-column: the old
        % velocity (known at the old mesh's physical r-positions) is
        % interpolated onto the current mesh's physical r-positions before
        % being used in the unsteady term. See interp_facefield_radial_to_mesh.m.
        %
        % No data on the very first step (or on any size mismatch, e.g. an
        % older saved state without the Rur/Ruz fields) is treated as
        % starting from rest, which is the standard, defensible initial
        % condition for a first-order-in-time discretization's first step.
        if isfield(par, 'useUnsteadyStokes') && par.useUnsteadyStokes
            if isfield(old, 'ur2DFaceField') && ~isempty(old.ur2DFaceField) && ...
                    isequal(size(old.ur2DFaceField), [parMAC.Nr+1, parMAC.Nz]) && ...
                    isfield(old, 'Rur2DField') && ~isempty(old.Rur2DField) && ...
                    isequal(size(old.Rur2DField), [parMAC.Nr+1, parMAC.Nz])
                bc.urPrev = interp_facefield_radial_to_mesh(old.ur2DFaceField, old.Rur2DField, meshMAC.Rur);
            else
                bc.urPrev = zeros(parMAC.Nr+1, parMAC.Nz);
            end
            if isfield(old, 'uz2DFaceField') && ~isempty(old.uz2DFaceField) && ...
                    isequal(size(old.uz2DFaceField), [parMAC.Nr, parMAC.Nz+1]) && ...
                    isfield(old, 'Ruz2DField') && ~isempty(old.Ruz2DField) && ...
                    isequal(size(old.Ruz2DField), [parMAC.Nr, parMAC.Nz+1])
                bc.uzPrev = interp_facefield_radial_to_mesh(old.uz2DFaceField, old.Ruz2DField, meshMAC.Ruz);
            else
                bc.uzPrev = zeros(parMAC.Nr, parMAC.Nz+1);
            end
        end

        fluidMAC = solve_stokes_bodyfitted_MAC(meshMAC, stateMAC, bc, parMAC, tNow);
        fluidMAC.bc = bc;
        [tauLMAC, tauEMAC] = estimate_wall_shear_bodyfitted(meshMAC, fluidMAC, stateMAC, parMAC);
    catch ME
        ok = false;
        stopReason = ['Body-fitted MAC Stokes fluid failed: ', ME.message];
        fluid = empty_fluid_2D_return(NzOut, hOut);
        return;
    end

    pMidMAC = extract_midgap_pressure_line_bodyfitted(fluidMAC);
    pLMAC = fluidMAC.P(1,:).';
    pEMAC = fluidMAC.P(end,:).';
    uzLMAC = fluidMAC.uz(1,:).';
    uzEMAC = fluidMAC.uz(end,:).';

    % Convert the MAC result back to the original coupling-grid format.
    pOutVec    = safe_interp1_same_or_resample(zc, pMidMAC(:), zOut, 'pMidMAC');
    pLOutVec   = safe_interp1_same_or_resample(zc, pLMAC(:),    zOut, 'pLMAC');
    pEOutVec   = safe_interp1_same_or_resample(zc, pEMAC(:),    zOut, 'pEMAC');

    % Compatibility with the original coupled code: fluid.p is stored on
    % the old axial grid, whose first/last entries are treated as prescribed
    % inlet/outlet pressures in diagnostics and pressure limiting. The native
    % MAC cell-center pressure remains available as fluid.P.
    pOutVec(1) = parMAC.pIn;       pOutVec(end) = parMAC.pOut;
    pLOutVec(1) = parMAC.pIn;      pLOutVec(end) = parMAC.pOut;
    pEOutVec(1) = parMAC.pIn;      pEOutVec(end) = parMAC.pOut;

    tauLOutVec = safe_interp1_same_or_resample(zc, tauLMAC(:),  zOut, 'tauLMAC');
    tauEOutVec = safe_interp1_same_or_resample(zc, tauEMAC(:),  zOut, 'tauEMAC');
    uzLOutVec  = safe_interp1_same_or_resample(zc, uzLMAC(:),   zOut, 'uzLMAC');
    uzEOutVec  = safe_interp1_same_or_resample(zc, uzEMAC(:),   zOut, 'uzEMAC');

    % The legacy diagnostics expect Q on intervals between zOut nodes.
    % Compute per-radian flux Q = int r*u_z dr at MAC axial faces and then
    % sample it at the original interval midpoints.
    Qfaces = compute_bodyfitted_axisym_flux_faces(meshMAC, fluidMAC);
    zQfaces = meshMAC.zF(:);
    zQout = 0.5*(zOut(1:end-1) + zOut(2:end));
    Qout = safe_interp1_same_or_resample(zQfaces, Qfaces(:), zQout, 'Qfaces');

    fluid = struct();
    fluid.meshType = 'bodyfitted_MAC';
    fluid.leukocyteInnerBoundary = ternary_local(useRLoutInner, 'RLout_reference', 'state_deltaL');
    fluid.usesRLoutFluidInterfaceForSolidLeukocyte = useRLoutInner;
    fluid.meshF = meshMAC;
    meshF = meshMAC;

    % Compatibility fields used elsewhere in the coupled code.
    fluid.p = pOutVec(:);
    fluid.pL = pLOutVec(:);
    fluid.pE = pEOutVec(:);
    fluid.Q = Qout(:);
    fluid.tauL = tauLOutVec(:);
    fluid.tauE = tauEOutVec(:);
    fluid.uzL = uzLOutVec(:);
    fluid.uzE = uzEOutVec(:);
    fluid.gap = hOut(:);

    % Native MAC fields for plotting/debugging.
    fluid.P = fluidMAC.P;
    fluid.ur = fluidMAC.ur;
    fluid.uz = fluidMAC.uz;
    fluid.urC = fluidMAC.urC;
    fluid.uzC = fluidMAC.uzC;
    fluid.divInf = fluidMAC.divInf;
    fluid.linRes = fluidMAC.linRes;
    fluid.condNumber = fluidMAC.condNumber;
    fluid.Nunknown = fluidMAC.Nunknown;
    fluid.bc = bc;
    [fluid.tractionL, fluid.tractionE] = compute_bodyfitted_wall_traction(meshMAC, fluidMAC, parMAC);
    fluid.pAvg = radial_average_pressure_bodyfitted(fluidMAC.P, meshMAC);
    fluid.pMid = pMidMAC(:);

    % Legacy field names, retained so old plotting branches do not error.
    fluid.ur2D = fluidMAC.urC(:);
    fluid.uz2D = fluidMAC.uzC(:);
    fluid.pCell = fluidMAC.P(:);
    fluid.sigmaCell = [];
    fluid.cellCenter = [meshMAC.Rp(:), meshMAC.Zp(:)];
end



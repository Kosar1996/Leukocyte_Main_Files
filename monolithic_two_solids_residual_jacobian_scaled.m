function [R, J] = monolithic_two_solids_residual_jacobian_scaled( ...
    y, old, meshE, interfaceE, baseE, ...
    meshL, interfaceL, baseL, parL, ...
    z, par, freeE, fixE, valsE, freeL, fixL, valsL, ...
    JuE, JuL, Jp, solidTargetE, solidTargetL, fluidTarget)
% Semi-analytical two-solid residual/Jacobian.
%
% Unknown vector:
%   y = [uE_free/JuE; uL_free/JuL; p_internal/Jp]
%
% The solid residuals and solid tangents are analytical. The fluid derivatives
% with respect to local wall radii/velocities are computed by cheap local
% finite differences of local_flux_and_shear, not by finite-differencing every
% global solid DOF. This avoids the thousands of expensive solid residual
% evaluations that killed the full-DOF residual-only fsolve test.

    N = numel(z);
    nE = numel(freeE);
    nL = numel(freeL);

    try
        [uE, uL, p] = unpack_two_solid_y( ...
            y, old, freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp, par);
        assert_solid_geometry_ok( ...
            solid_geometry_quality(meshE, uE, 'endothelium fsolve iterate'), par);
        assert_solid_geometry_ok( ...
            solid_geometry_quality(meshL, uL, 'leukocyte fsolve iterate'), par);

        [RE, RL, RF, JEE, JEL, JEp, JLE, JLL, JLp, JFE, JFL, JFp] = ...
            monolithic_two_solids_residual_jacobian_unscaled( ...
            uE, uL, p, old, meshE, interfaceE, baseE, ...
            meshL, interfaceL, baseL, parL, z, par, freeE, freeL);

        R = [
            RE / solidTargetE
            RL / solidTargetL
            RF / fluidTarget
        ];

        if nargout > 1
            J = [
                (JEE*JuE)/solidTargetE, (JEL*JuL)/solidTargetE, (JEp*Jp)/solidTargetE
                (JLE*JuE)/solidTargetL, (JLL*JuL)/solidTargetL, (JLp*Jp)/solidTargetL
                (JFE*JuE)/fluidTarget,  (JFL*JuL)/fluidTarget,  (JFp*Jp)/fluidTarget
            ];
        end

        if any(~isfinite(R))
            R = 1e12 * ones(nE+nL+N-2,1);
            if nargout > 1
                J = speye(numel(R), numel(y));
            end
        end
    catch ME
        if contains(ME.message, 'Gap violates minGap') || ...
           contains(ME.message, 'Negative or zero J') || ...
           contains(ME.message, 'Non-positive radius') || ...
           contains(ME.message, 'Solid geometry guard failed') || ...
           contains(ME.message, 'Element inverted')
            R = 1e12 * ones(nE+nL+N-2,1);
            if nargout > 1
                J = speye(numel(R), numel(y));
            end
        else
            rethrow(ME);
        end
    end
end

function [RE, RL, RF, JEE, JEL, JEp, JLE, JLL, JLp, JFE, JFL, JFp] = ...
    monolithic_two_solids_residual_jacobian_unscaled( ...
    uE, uL, p, old, meshE, interfaceE, baseE, ... %#ok<INUSD>
    meshL, interfaceL, baseL, parL, z, par, freeE, freeL) %#ok<INUSD>

    ndofE = size(meshE.nodes,1) * 2;
    ndofL = size(meshL.nodes,1) * 2;
    N = numel(z);

    p(1) = par.pIn;
    p(end) = par.pOut;

    [deltaE, UwE, HrE, HUE] = monolithic_interface_kinematics( ...
        meshE, uE, old.uE, interfaceE, z, par);
    [deltaL, UwL, HrL, HUL] = monolithic_interface_kinematics( ...
        meshL, uL, old.uL, interfaceL, z, par);

    gap = deltaE - deltaL;
    % Tested Sep 3 (Issue 1/2/3 shared mechanism) and DISPROVEN: this
    % guard turns any trial Newton iterate with gap <= par.minGap into a
    % discontinuous 1e12 penalty residual, which looked like a plausible
    % contributor to the repeated exitflag=-2 events at the pinned floor.
    % Tried widening the guard (an overridable par.minGapNewtonFloor set
    % below par.minGap, giving Newton's trial iterates room to explore
    % past the official floor) to see if a softer boundary let the solve
    % converge more cleanly. Tested directly against a real near-floor
    % state: instead of helping, every trial iterate near the new floor
    % produced literally infinite residuals (RE=Inf, RL=Inf, RF=Inf), and
    % the outer dt-retry mechanism could not escape it even after four
    % successive dt halvings (down to 1/16 of the original). Consistent
    % with the local_flux_and_shear.m finding the same day: part of what
    % "blows up" near the floor is genuine physical divergence (thin-film
    % shear ~1/gap, pressure gradient ~1/gap^3), not just numerical noise
    % a wider exploration zone can safely ride through. This guard needs
    % to stay a hard boundary at par.minGap; a real fix has to address the
    % underlying near-floor physics (see local_flux_and_shear.m), not
    % relax where Newton is allowed to look.
    if any(gap <= par.minGap)
        error('Gap violates minGap.');
    end

    [Q, ~, ~, tauL, tauE, ~, ~] = ...
        local_flux_and_shear(z, p, deltaL, deltaE, UwL, UwE, par);

    if isfield(par, 'useTwoSolidFullAnalyticalJacobian') && par.useTwoSolidFullAnalyticalJacobian
        fs = local_flux_shear_sensitivities_analytical_full(z, p, deltaL, deltaE, UwL, UwE, par);
    else
        fs = local_flux_shear_sensitivities_fd_full(z, p, deltaL, deltaE, UwL, UwE, par, Q, tauL, tauE);
    end

    [pLoadE, pLoadL, global2DJacE, global2DJacL] = ...
        global2d_pressure_traction_loads( ...
        z, p, deltaE, deltaL, meshE, uE, meshL, uL, par);

    % ------------------------------------------------------------
    % Endothelium residual and same-solid tangent.
    % Fluid-on-endothelium traction: normal = +p, tangent = -tauE.
    %
    % Sign convention independently re-derived and confirmed: traction
    % on a solid from the fluid is
    % t = T*n, n = the SOLID's own outward normal (pointing INTO the
    % fluid). Endothelium is the OUTER body enclosing the fluid gap, so
    % its outward normal points in the -r direction (n_r=-1). With fluid
    % stress T_rr=-p, T_rz=tau (tau as computed in local_flux_and_shear.m,
    % i.e. mu*du_z/dr at that radius): t_r = T_rr*n_r = -p*(-1) = +p,
    % t_z = T_rz*n_r = tau*(-1) = -tau. Matches trE.normal=+pLoadE,
    % trE.tangent=-tauE exactly -- confirmed correct, not a bug.
    % ------------------------------------------------------------
    trE.normal = pLoadE;
    trE.tangent = -tauE;
    if use_exact_interface_in_monolithic(par)
        trE.z = z;
    end

    FextE = zeros(ndofE,1);
    [FextE, KextE, BnE, BtE] = ...
        apply_interface_traction_sensitivity(meshE, uE, FextE, interfaceE, trE);

    [FintE, KintE] = assemble_finite_def_axisym(meshE, uE, par);
    [FviscE, KviscE] = assemble_axisym_kelvin_voigt_viscous(meshE, uE, old.uE, par);
    FintE = FintE + FviscE;
    KintE = KintE + KviscE;

    REfull = FintE - FextE;
    RE = REfull(freeE);

    % ------------------------------------------------------------
    % Leukocyte residual and same-solid tangent.
    % Fluid-on-leukocyte traction: normal = -p, tangent = +tauL.
    %
    % Sign convention independently re-derived and confirmed: leukocyte
    % is the INNER body, so its
    % outward normal (into the fluid) points in the +r direction
    % (n_r=+1). t_r = T_rr*n_r = -p*(+1) = -p, t_z = T_rz*n_r = tau*(+1)
    % = +tau. Matches trL.normal=-pLoadL, trL.tangent=+tauL exactly --
    % opposite sign from the endothelium's, as required by the two
    % bodies having opposite outward normals across the same gap.
    % Confirmed correct, not a bug.
    % ------------------------------------------------------------
    trL.normal = -pLoadL;
    trL.tangent = tauL;
    if use_exact_interface_in_monolithic(par)
        trL.z = z;
    end
    trL = apply_leukocyte_traction_support(trL, par);

    FextL = zeros(ndofL,1);
    [FextL, KextL, BnL, BtL] = ...
        apply_interface_traction_sensitivity(meshL, uL, FextL, interfaceL, trL);

    [FintL, KintL] = assemble_finite_def_axisym(meshL, uL, parL);
    [FviscL, KviscL] = assemble_axisym_kelvin_voigt_viscous(meshL, uL, old.uL, parL);
    FintL = FintL + FviscL;
    KintL = KintL + KviscL;

    RLfull = FintL - FextL;
    RL = RLfull(freeL);

    % ------------------------------------------------------------
    % Fluid residual.
    % ------------------------------------------------------------
    re    = deltaE(:);
    rl    = deltaL(:);
    reOld = old.deltaE(:);
    rlOld = old.deltaL(:);

    A    = 0.5 * (re.^2    - rl.^2);
    Aold = 0.5 * (reOld.^2 - rlOld.^2);

    Ssrc = -par.SsrcFactor * (A - Aold) / par.dt;

    dzControl = global_1d_control_lengths(z);
    RF = zeros(N-2,1);
    for i = 2:N-1
        RF(i-1) = (Q(i) - Q(i-1))/dzControl(i) - Ssrc(i);
    end

    % Divergence matrix mapping face fluxes Q(1:N-1) to interior residuals.
    Dq = sparse(N-2, N-1);
    for i = 2:N-1
        row = i-1;
        Dq(row,i)   =  1/dzControl(i);
        Dq(row,i-1) = -1/dzControl(i);
    end

    KFp_all = Dq * fs.dQdp;
    KFre = Dq * fs.dQdre;
    KFrl = Dq * fs.dQdrl;
    KFUwE = Dq * fs.dQdUwE;
    KFUwL = Dq * fs.dQdUwL;

    for i = 2:N-1
        row = i-1;
        % R = div(Q) - Ssrc, Ssrc = -(A-Aold)/dt.
        % dR/dre = +re/dt; dR/drl = -rl/dt.
        KFre(row,i) = KFre(row,i) + par.SsrcFactor * re(i) / par.dt;
        KFrl(row,i) = KFrl(row,i) - par.SsrcFactor * rl(i) / par.dt;
    end

    JFEfull = KFre*HrE + KFUwE*HUE;
    JFLfull = KFrl*HrL + KFUwL*HUL;
    JFE = JFEfull(:,freeE);
    JFL = JFLfull(:,freeL);
    JFp = KFp_all(:,2:end-1);

    % ------------------------------------------------------------
    % Solid Jacobian blocks.
    % RE = FintE - FextE(p, -tauE)
    % RL = FintL - FextL(-p, +tauL)
    % ------------------------------------------------------------
    % Tested Sep 3 (Issue 2) and DISPROVEN: fs.dtauEdre/dtauEdrl/dtauLdre/
    % dtauLdrl (shear sensitivity to radius) grow as ~1/gap^2 approaching
    % par.minGap (measured: ~230 at gap=1 micron to ~2.4e8 at the 1nm
    % floor; confirmed genuine -- matches an independent hand-derived
    % analytical formula to ~5e-5 relative, not an FD artifact or a bug).
    % Tried capping these four sensitivity matrices before they enter the
    % Jacobian (overridable par.dTauDRadiusCap, a damped-Newton tradeoff,
    % same category as this codebase's existing solidTrustU0/
    % solidTrustUMax caps -- not a correctness fix, since there was
    % nothing wrong to fix). Tested directly against a real pinned-floor
    % state at two cap magnitudes 100x apart (1e5 and 1e3, the latter
    % only ~4x above the normal away-from-floor scale): both produced
    % results IDENTICAL to the uncapped baseline in every respect
    % (iterations=25, scaledRes=1.084e-02, RE=2.354e-10, RL=1.084e-09,
    % every single step, no deviation). This term is not the binding
    % constraint on Newton's convergence at this fixed point, regardless
    % of magnitude -- reverted. The actual bottleneck remains
    % unidentified; do not assume this specific coupling term without
    % re-deriving further.
    dtauE_duE = fs.dtauEdre*HrE + fs.dtauEdUwE*HUE;
    dtauE_duL = fs.dtauEdrl*HrL + fs.dtauEdUwL*HUL;
    dtauL_duE = fs.dtauLdre*HrE + fs.dtauLdUwE*HUE;
    dtauL_duL = fs.dtauLdrl*HrL + fs.dtauLdUwL*HUL;

    % Endothelium tangent: t_t = -tauE, so -Bt*d(t_t) = +Bt*d(tauE).
    JEEfull = KintE - KextE + BtE*dtauE_duE;
    JELfull =              + BtE*dtauE_duL;
    JEpfull = -BnE*global2DJacE(:,2:end-1) + BtE*fs.dtauEdp(:,2:end-1);

    % Leukocyte tangent: t_n = -p, t_t = +tauL.
    JLEfull =              - BtL*dtauL_duE;
    JLLfull = KintL - KextL - BtL*dtauL_duL;
    JLpfull =  BnL*global2DJacL(:,2:end-1) - BtL*fs.dtauLdp(:,2:end-1);

    JEE = JEEfull(freeE, freeE);
    JEL = JELfull(freeE, freeL);
    JEp = JEpfull(freeE, :);

    JLE = JLEfull(freeL, freeE);
    JLL = JLLfull(freeL, freeL);
    JLp = JLpfull(freeL, :);
end
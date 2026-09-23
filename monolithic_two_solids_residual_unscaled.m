function [RE, RL, RF] = monolithic_two_solids_residual_unscaled( ...
    uE, uL, p, old, meshE, interfaceE, baseE, ... %#ok<INUSD>
    meshL, interfaceL, baseL, parL, z, par, freeE, freeL) %#ok<INUSD>

    ndofE = size(meshE.nodes,1) * 2;
    ndofL = size(meshL.nodes,1) * 2;
    N = numel(z);

    p(1) = par.pIn;
    p(end) = par.pOut;

    [deltaE, UwE] = monolithic_interface_kinematics_value_only( ...
        meshE, uE, old.uE, interfaceE, z, par);
    [deltaL, UwL] = monolithic_interface_kinematics_value_only( ...
        meshL, uL, old.uL, interfaceL, z, par);

    gap = deltaE - deltaL;
    if any(gap <= par.minGap)
        error('Gap violates minGap.');
    end

    [Q, ~, ~, tauL, tauE, ~, ~] = ...
        local_flux_and_shear(z, p, deltaL, deltaE, UwL, UwE, par);

    [pLoadE, pLoadL] = global2d_pressure_traction_loads( ...
        z, p, deltaE, deltaL, meshE, uE, meshL, uL, par);

    % Endothelium residual. Fluid-on-endothelium traction: +p normal, -tauE tangent.
    trE.normal = pLoadE;
    trE.tangent = -tauE;
    if use_exact_interface_in_monolithic(par)
        trE.z = z;
    end

    FextE = zeros(ndofE,1);
    [FextE, ~] = apply_interface_traction(meshE, uE, FextE, interfaceE, trE);

    FintE = assemble_finite_def_internal_force_only(meshE, uE, par) + ...
        assemble_axisym_kelvin_voigt_viscous_force_only(meshE, uE, old.uE, par);

    REfull = FintE - FextE;
    RE = REfull(freeE);

    % Leukocyte residual. Fluid-on-leukocyte traction: -p normal, +tauL tangent.
    trL.normal = -pLoadL;
    trL.tangent = tauL;
    if use_exact_interface_in_monolithic(par)
        trL.z = z;
    end
    trL = apply_leukocyte_traction_support(trL, par);

    FextL = zeros(ndofL,1);
    [FextL, ~] = apply_interface_traction(meshL, uL, FextL, interfaceL, trL);

    FintL = assemble_finite_def_internal_force_only(meshL, uL, parL) + ...
        assemble_axisym_kelvin_voigt_viscous_force_only(meshL, uL, old.uL, parL);

    RLfull = FintL - FextL;
    RL = RLfull(freeL);

    % Fluid mass residual using both current interfaces.
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
end

function [rAtZ, wzAtZ] = exact_interface_radius_velocity( ...
    mesh, uNew, uOld, interfaceNodes, zq, dt)

    [rNew, zNew, uzNew] = deformed_interface_curve(mesh, uNew, interfaceNodes);
    [~,    zOld, uzOld] = deformed_interface_curve(mesh, uOld, interfaceNodes);

    rAtZ = interp_curve_values(zNew, rNew, zq);
    uzNewAtZ = interp_curve_values(zNew, uzNew, zq);
    uzOldAtZ = interp_curve_values(zOld, uzOld, zq);
    wzAtZ = (uzNewAtZ - uzOldAtZ) / dt;
end

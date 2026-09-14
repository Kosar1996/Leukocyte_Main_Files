function [delta, Uw] = monolithic_interface_kinematics_value_only( ...
    mesh, uNew, uOld, interfaceNodes, z, par)
% Value-only wrapper around monolithic_interface_kinematics.
% Works for either endothelium or leukocyte.

    [delta, Uw] = monolithic_interface_kinematics( ...
        mesh, uNew, uOld, interfaceNodes, z, par);
end

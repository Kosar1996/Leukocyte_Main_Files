function [rAtZ, wzAtZ] = exact_interface_radius_velocity( ...
    mesh, uNew, uOld, interfaceNodes, zq, dt)
% EXACT_INTERFACE_RADIUS_VELOCITY
% Position-aware interpolation of deformed interface radius and mesh velocity.
%
% REVISION HISTORY & MERGED BUG FIXES:
% -------------------------------------------------------------------------
% 1. Lines 24-28 (Unified Eulerian Reference Configuration):
%    [OLD]: uzOldAtZ = interp_curve_values(zOld, uzOld, zq);
%    Evaluating uOld at zOld created a spatial reference mismatch relative to 
%    uNew evaluated at zNew. This generated artificial advective velocity spikes 
%    (10^6 um/s) along deforming boundaries.
%    FIXED: Evaluate both current and old displacements at the unified current 
%    deformed spatial coordinates (zNew) before computing wz = (uNew - uOld)/dt.
%
% 2. Lines 33-38 (Endpoint Boundary Profile Clamping Safeguard):
%    [OLD]: rAtZ = interp_curve_values(zNew, rNew, zq);
%    Evaluating query points (zq) outside the active node range [min(zNew), max(zNew)] 
%    with un-clamped linear extrapolation generated artificial radial bulges 
%    at domain ends (z = zMin and z = zMax).
%    FIXED: Clamped interpolated radius values (rAtZ) at boundary limits 
%    min(zNew) and max(zNew) to enforce flat physical profile continuity.
% -------------------------------------------------------------------------

    % Extract deformed coordinates and displacements
    [rNew, zNew, uzNew] = deformed_interface_curve(mesh, uNew, interfaceNodes);
    [~,    ~,    uzOld] = deformed_interface_curve(mesh, uOld, interfaceNodes);

    % Deformed radius along the current interface curve
    rAtZ = interp_curve_values(zNew, rNew, zq);
    
    % Endpoint clamping safeguard to prevent artificial boundary bulges near zMin and zMax
    rAtZ(zq < min(zNew)) = rNew(1);
    rAtZ(zq > max(zNew)) = rNew(end);

    % Interpolate current displacement field
    uzNewAtZ = interp_curve_values(zNew, uzNew, zq);

    % [OLD]: uzOldAtZ = interp_curve_values(zOld, uzOld, zq);
    % FIXED: Map previous displacement onto current spatial configuration (zNew)
    uzOldAtZ = interp_curve_values(zNew, uzOld, zq);

    % Material mesh interface velocity
    wzAtZ = (uzNewAtZ - uzOldAtZ) / dt;
end
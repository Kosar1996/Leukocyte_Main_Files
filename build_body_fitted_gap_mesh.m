% ========================================================================
% BODY-FITTED MAC CYLINDRICAL STOKES FLUID MODULE
% Imported from mac_bodyfitted_2D_stokes_rigid_leukocyte_pressureBC.m
% ========================================================================
function mesh = build_body_fitted_gap_mesh(state, par)

    Nr = par.Nr;
    Nz = par.Nz;

    zc = state.zc(:).';
    zF = state.zF(:).';

    dL_c = state.deltaL(:).';
    dE_c = state.deltaE(:).';
    h_c = dE_c - dL_c;

    dL_f = interp1(state.zc, state.deltaL, zF, 'linear','extrap').';
    dE_f = interp1(state.zc, state.deltaE, zF, 'linear','extrap').';
    h_f = dE_f - dL_f;

    if any(h_c <= 0) || any(h_f <= 0)
        error('Invalid body-fitted mesh: non-positive gap detected.');
    end

    [Rur, etaF_c] = build_radial_face_grid_bodyfitted(dL_c, dE_c, Nr, par);
    [Rzf, etaF_f] = build_radial_face_grid_bodyfitted(dL_f, dE_f, Nr, par);
    Rp = radial_cell_centers_from_faces(Rur);
    Ruz = radial_cell_centers_from_faces(Rzf);
    etaC = radial_cell_centers_from_faces(etaF_c);

    % Pressure-cell centers and radial-velocity faces at z centers
    Zp = repmat(zc,Nr,1);
    Zur = repmat(zc,Nr+1,1);

    % Axial-velocity faces at z faces
    Zuz = repmat(zF,Nr,1);

    % Radial cell faces at axial faces. These are used to compute the
    % conservative axial flux areas in the continuity equation.

    mesh = struct();
    mesh.Nr = Nr;
    mesh.Nz = Nz;
    mesh.zc = zc(:);
    mesh.zF = zF(:);
    mesh.dz = state.zF(2)-state.zF(1);

    mesh.etaF = etaF_c;
    mesh.etaF_zFace = etaF_f;
    mesh.etaC = etaC;

    mesh.deltaL_c = dL_c(:);
    mesh.deltaE_c = dE_c(:);
    mesh.h_c = h_c(:);

    mesh.deltaL_f = dL_f(:);
    mesh.deltaE_f = dE_f(:);
    mesh.h_f = h_f(:);

    mesh.Rp = Rp;
    mesh.Rur = Rur;
    mesh.Ruz = Ruz;
    mesh.Rzf = Rzf;
    mesh.Zp = Zp;
    mesh.Zur = Zur;
    mesh.Zuz = Zuz;
end

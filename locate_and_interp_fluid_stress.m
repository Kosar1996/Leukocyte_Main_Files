function [sigma, p, rz_actual, xieta, ij] = locate_and_interp_fluid_stress(mraw, meshF, ur2D, uz2D, mu, plist, rq, zq, i_guess, j_guess)
%LOCATE_AND_INTERP_FLUID_STRESS
% Robust point location + position-aware Q4 interpolation. Starts from a
% guessed (i,j) element and, if the Newton solve for local (xi,eta) wants
% to leave the [-1,1]x[-1,1] element, walks to the neighboring element in
% that direction and retries -- standard point-location-by-walking, needed
% because the body-fitted mesh is curved (sloped walls), so a simple
% single-column bracketing search can pick an adjacent-but-wrong element
% near sloped boundaries.

Nr = mraw.Nr; Nz = mraw.Nz;
i = i_guess; j = j_guess;

for walk = 1:10
    e = (j-1)*(Nr-1) + i;
    conn = meshF.elems(e,:);
    xe = meshF.nodes(conn,1);
    ze = meshF.nodes(conn,2);

    xi = 0; eta = 0;
    for it = 1:50
        [N, dNdxi] = shape_Q4(xi, eta);
        r_cur = N.' * xe;
        z_cur = N.' * ze;
        resid = [rq - r_cur; zq - z_cur];
        if norm(resid) < 1e-13
            break;
        end
        J = [xe ze].' * dNdxi;
        dxieta = J \ resid;
        xi = xi + dxieta(1);
        eta = eta + dxieta(2);
    end

    tol = 1e-6;
    if xi >= -1-tol && xi <= 1+tol && eta >= -1-tol && eta <= 1+tol
        xi = min(max(xi,-1),1);
        eta = min(max(eta,-1),1);
        break;
    end

    % walk to the neighboring element in whichever direction we overshot
    moved = false;
    if xi < -1 && i > 1, i = i-1; moved = true;
    elseif xi > 1 && i < Nr-1, i = i+1; moved = true;
    end
    if eta < -1 && j > 1, j = j-1; moved = true;
    elseif eta > 1 && j < Nz-1, j = j+1; moved = true;
    end
    if ~moved
        xi = min(max(xi,-1),1);
        eta = min(max(eta,-1),1);
        break;
    end
end

ij = [i, j];
xieta = [xi, eta];

[N, dNdxi] = shape_Q4(xi, eta);
J = [xe ze].' * dNdxi;
dNdx = dNdxi / J;
dNdr = dNdx(:,1);
dNdz = dNdx(:,2);

r_actual = N.' * xe;
z_actual = N.' * ze;
rz_actual = [r_actual, z_actual];

pc = N.' * plist(conn);

ur = zeros(4,1); uz = zeros(4,1);
for a = 1:4
    ur(a) = ur2D(conn(a));
    uz(a) = uz2D(conn(a));
end

durdr = dNdr.' * ur;
duzdz = dNdz.' * uz;
durdz = dNdz.' * ur;
duzdr = dNdr.' * uz;
ur_over_r = (N.' * ur) / r_actual;

srr = -pc + 2*mu*durdr;
stt = -pc + 2*mu*ur_over_r;
szz = -pc + 2*mu*duzdz;
srz = mu*(durdz + duzdr);

sigma = [srr, stt, szz, srz];
p = pc;
end

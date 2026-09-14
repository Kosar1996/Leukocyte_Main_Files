function [leuko, endo] = compute_interface_traction_mismatch( ...
    meshE, uE, uEPrev, parE, meshL, uL, uLPrev, parL, dtStep, fluid, opts)
%COMPUTE_INTERFACE_TRACTION_MISMATCH
% Core computation extracted from check_interface_traction_mismatch_report.m
% so it can be called on LIVE in-progress state (e.g. from inside a
% partitioned traction-correction loop, once per pass) as well as on a
% saved out.stateHist{k}/out.fluidHist{k} snapshot. No behavior change from
% the original -- same stress recovery, same interpolation, same
% normal/tangential rotation, same arc-length parametrization.
%
% Usage:
%   [leuko, endo] = compute_interface_traction_mismatch( ...
%       meshE, uE, uEPrev, par, meshL, uL, uLPrev, parL, dtStep, fluid, opts)
%
% Inputs:
%   meshE, uE, uEPrev, parE   endothelium mesh, current/previous displacement,
%                             and its material par struct (needs .mu, .Ge,
%                             .Ke, .etaE)
%   meshL, uL, uLPrev, parL   same for the leukocyte (parL should already
%                             have Ge=GL, Ke=KL, etaE=etaL remapped in)
%   dtStep                    accepted dt for the Kelvin-Voigt viscous term
%   fluid                     fluid state struct with .meshF, .ur2D, .uz2D,
%                             .pCell (same as out.fluidHist{k})
%   opts (optional): opts.nQuery (61), opts.epsFrac (0.1), opts.trimFrac
%       (0.05), opts.slopeThreshold (0.15) -- same defaults/meaning as
%       check_interface_traction_mismatch_report.m
%
% Outputs: leuko, endo -- same fields as report.leuko/report.endo in
% check_interface_traction_mismatch_report.m (.s, .en, .et, .z, .rWall,
% .slope, .stats.maxAbsEn/medianAbsEn/p90AbsEn/maxAbsEt/... /locMaxEn/locMaxEt).

if nargin < 11 || isempty(opts), opts = struct(); end
if ~isfield(opts,'nQuery')   || isempty(opts.nQuery),   opts.nQuery   = 61;   end
if ~isfield(opts,'epsFrac')  || isempty(opts.epsFrac),  opts.epsFrac  = 0.1;  end
if ~isfield(opts,'trimFrac') || isempty(opts.trimFrac), opts.trimFrac = 0.05; end

mu = parE.mu;
meshF = add_fluid_nodes(fluid.meshF);
mraw = fluid.meshF;
% NOTE: fluid stress at arbitrary interface query points is now evaluated
% via locate_and_interp_fluid_stress.m (position-aware bilinear
% quadrilateral interpolation, resolving the true local coordinates of
% the query point within its containing cell) instead of the old
% recover_fluid_nodes_pressure_stress_Q4.m + scatteredInterpolant-between-
% centroids approach. The old approach always evaluated at a fixed
% element centroid regardless of query position -- confirmed bug, fixed
% here. See locate_and_interp_fluid_stress.m for the verified fix.

stressE = recover_nodal_stress_axisym_viscoelastic(meshE, uE, uEPrev, dtStep, parE);
stressL = recover_nodal_stress_axisym_viscoelastic(meshL, uL, uLPrev, dtStep, parL);

rE = meshE.nodes(:,1) + uE(1:2:end);
zE = meshE.nodes(:,2) + uE(2:2:end);
rL = meshL.nodes(:,1) + uL(1:2:end);
zL = meshL.nodes(:,2) + uL(2:2:end);

interpF = @(x,y,v) scatteredInterpolant(x(:), y(:), v(:), 'linear', 'none');
FLrr = interpF(rL, zL, stressL.sigma_rr);
FLzz = interpF(rL, zL, stressL.sigma_zz);
FLrz = interpF(rL, zL, stressL.sigma_rz);
FErr = interpF(rE, zE, stressE.sigma_rr);
FEzz = interpF(rE, zE, stressE.sigma_zz);
FErz = interpF(rE, zE, stressE.sigma_rz);

mesh = fluid.meshF;
if ~all(isfield(mesh, {'zc','deltaL_c','deltaE_c'}))
    error('compute_interface_traction_mismatch: needs fluid.meshF.zc/deltaL_c/deltaE_c.');
end

% Query range must respect BOTH solids' actual deformed z-extent, not just
% the leukocyte's -- otherwise points outside the endothelium's coverage
% return NaN from its scatteredInterpolant (convex-hull extrapolation is
% 'none'), silently dropping a large fraction of endothelium sample points
% (found via check_weak_force_balance.m: 40/61 points were NaN because this
% only checked zL, never zE).
zLo = max([min(mesh.zc), min(zL), min(zE)]);
zHi = min([max(mesh.zc), max(zL), max(zE)]);
span = zHi - zLo;
zLoQ = zLo + opts.trimFrac * span;
zHiQ = zHi - opts.trimFrac * span;
zQuery = linspace(zLoQ, zHiQ, opts.nQuery);

deltaLofz = @(z) interp1(mesh.zc, mesh.deltaL_c, z, 'linear', 'extrap');
deltaEofz = @(z) interp1(mesh.zc, mesh.deltaE_c, z, 'linear', 'extrap');

dz_fd = 1e-8;
slopeLofz = @(z) (deltaLofz(z+dz_fd) - deltaLofz(z-dz_fd)) / (2*dz_fd);
slopeEofz = @(z) (deltaEofz(z+dz_fd) - deltaEofz(z-dz_fd)) / (2*dz_fd);

nQ = numel(zQuery);
enL = nan(nQ,1); etL = nan(nQ,1);
enE = nan(nQ,1); etE = nan(nQ,1);
slopeLv = nan(nQ,1); slopeEv = nan(nQ,1);
tnFluidL = nan(nQ,1); ttFluidL = nan(nQ,1);
tnFluidE = nan(nQ,1); ttFluidE = nan(nQ,1);

for k = 1:nQ
    z = zQuery(k);
    rLwall = deltaLofz(z);
    rEwall = deltaEofz(z);
    gap = rEwall - rLwall;
    eps_ = opts.epsFrac * gap;

    slopeL = slopeLofz(z);
    slopeLv(k) = slopeL;
    nL = [1, -slopeL] / norm([1, -slopeL]);

    slopeE = slopeEofz(z);
    slopeEv(k) = slopeE;
    nE = [1, -slopeE] / norm([1, -slopeE]);

    rq  = rLwall + eps_;
    [igL,jgL] = initial_guess_ij_local(mraw, rq, z);
    sigFvec = locate_and_interp_fluid_stress(mraw, meshF, fluid.ur2D, fluid.uz2D, mu, fluid.pCell, rq, z, igL, jgL);
    sigF = [sigFvec(1), sigFvec(4); sigFvec(4), sigFvec(3)];
    rq2 = rLwall - eps_;
    sigS = [FLrr(rq2,z), FLrz(rq2,z); FLrz(rq2,z), FLzz(rq2,z)];
    [tnF, ttF] = traction_components(sigF, nL);
    [tnS, ttS] = traction_components(sigS, nL);
    enL(k) = tnS - tnF;
    etL(k) = ttS - ttF;
    tnFluidL(k) = tnF; ttFluidL(k) = ttF;

    rq  = rEwall - eps_;
    [igE,jgE] = initial_guess_ij_local(mraw, rq, z);
    sigFvec = locate_and_interp_fluid_stress(mraw, meshF, fluid.ur2D, fluid.uz2D, mu, fluid.pCell, rq, z, igE, jgE);
    sigF = [sigFvec(1), sigFvec(4); sigFvec(4), sigFvec(3)];
    rq2 = rEwall + eps_;
    sigS = [FErr(rq2,z), FErz(rq2,z); FErz(rq2,z), FEzz(rq2,z)];
    [tnF, ttF] = traction_components(sigF, nE);
    [tnS, ttS] = traction_components(sigS, nE);
    enE(k) = tnS - tnF;
    etE(k) = ttS - ttF;
    tnFluidE(k) = tnF; ttFluidE(k) = ttF;
end

sL = arc_length_from_slope(zQuery, slopeLv);
sE = arc_length_from_slope(zQuery, slopeEv);

leuko = pack_side(sL, enL, etL, zQuery(:), deltaLofz(zQuery(:)), slopeLv, tnFluidL, ttFluidL);
endo  = pack_side(sE, enE, etE, zQuery(:), deltaEofz(zQuery(:)), slopeEv, tnFluidE, ttFluidE);

end

function side = pack_side(s, en, et, z, rWall, slope, tnFluid, ttFluid)
side = struct();
side.s = s(:);
side.en = en(:);
side.et = et(:);
side.z = z(:);
side.rWall = rWall(:);
side.slope = slope(:);
side.tnFluid = tnFluid(:);
side.ttFluid = ttFluid(:);

absEn = abs(en(isfinite(en)));
absEt = abs(et(isfinite(et)));

side.stats = struct();
if isempty(absEn)
    side.stats.maxAbsEn = NaN; side.stats.medianAbsEn = NaN; side.stats.p90AbsEn = NaN;
else
    side.stats.maxAbsEn = max(absEn);
    side.stats.medianAbsEn = median(absEn);
    side.stats.p90AbsEn = prctile(absEn, 90);
end
if isempty(absEt)
    side.stats.maxAbsEt = NaN; side.stats.medianAbsEt = NaN; side.stats.p90AbsEt = NaN;
else
    side.stats.maxAbsEt = max(absEt);
    side.stats.medianAbsEt = median(absEt);
    side.stats.p90AbsEt = prctile(absEt, 90);
end

[~, iEn] = max(abs(en(:)));
[~, iEt] = max(abs(et(:)));
side.stats.locMaxEn = [side.rWall(iEn), side.z(iEn)];
side.stats.locMaxEt = [side.rWall(iEt), side.z(iEt)];


% relative to that point's own fluid traction component -- not a global
% reference scale. %mismatch = |solid - fluid| / |fluid| * 100, normal and
% tangential each normalized against their own local fluid component.
pctEn = 100 * abs(en(:)) ./ abs(tnFluid(:));
pctEt = 100 * abs(et(:)) ./ abs(ttFluid(:));
side.pctEn = pctEn;
side.pctEt = pctEt;

validEn = pctEn(isfinite(pctEn));
validEt = pctEt(isfinite(pctEt));
if isempty(validEn)
    side.stats.pctMaxEn = NaN; side.stats.pctP90En = NaN; side.stats.pctMedianEn = NaN;
else
    side.stats.pctMaxEn = max(validEn);
    side.stats.pctP90En = prctile(validEn, 90);
    side.stats.pctMedianEn = median(validEn);
end
if isempty(validEt)
    side.stats.pctMaxEt = NaN; side.stats.pctP90Et = NaN; side.stats.pctMedianEt = NaN;
else
    side.stats.pctMaxEt = max(validEt);
    side.stats.pctP90Et = prctile(validEt, 90);
    side.stats.pctMedianEt = median(validEt);
end
end

function s = arc_length_from_slope(zQuery, slope)
ds = sqrt(1 + slope(:).^2);
dz = diff(zQuery(:));
segMean = 0.5 * (ds(1:end-1) + ds(2:end));
s = [0; cumsum(segMean .* dz)];
end

function [tn, tt] = traction_components(sig, n)
n = n(:);
that = [-n(2); n(1)];
t = sig * n;
tn = t.' * n;
tt = t.' * that;
end

function [i,j] = initial_guess_ij_local(mraw, rq, z)
% Structured-grid bracketing to seed locate_and_interp_fluid_stress's
% point-location walk. Only needs to be approximately right -- the walk
% corrects it if the true containing element is a neighbor.
Nr = mraw.Nr; Nz = mraw.Nz;
Zcol = mraw.Zp(1,:);
[~,j0] = min(abs(Zcol - z)); j = min(max(j0,1),Nz-1);
if Zcol(j) > z && j>1, j=j-1; end
Rcol = mraw.Rp(:,j);
i0 = find(Rcol <= rq, 1, 'last'); if isempty(i0), i0=1; end
i = min(max(i0,1),Nr-1);
end

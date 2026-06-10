% VERIFY_PROJECTION_PORT Check the minvtran->projinv port on R2026a.
% Exercises every ported path with real Shasta data:
%   1. mstruct2projcrs on the real topo projection (eqaconicstd)
%   2. include_vars_melt path: domain center -> lon -> timezone
%   3. TopoSunAngle 'projected' case (per-pixel lat/lon -> sun angles)
%   4. load_coarse_topo -> reprojectRaster geographic->projected
% Cross-check the lat/lon numbers against pyproj (independent implementation).
topofile = 'D:\code\ucsb\data\ParBal\Shasta\inputs\ShastaTopography.h5';
ltfile   = 'D:\code\ucsb\data\ParBal\Shasta\GLDAS\GLDAS_topo.h5';

fprintf('minvtran exists: %d (0 expected on R2026a)\n', exist('minvtran','file'));

% 1. build the projcrs
[slope,hdr] = GetTopography(topofile,'slope');
p = mstruct2projcrs(hdr.ProjectionStructure);
fprintf('projcrs built: %s\n', p.Name);

% 2. include_vars_melt path (timezone longitude)
x = mean(hdr.RasterReference.XWorldLimits);
y = mean(hdr.RasterReference.YWorldLimits);
[lat,lon] = projinv(p,x,y);
fprintf('center: x=%.2f y=%.2f -> lat=%.6f lon=%.6f\n', x, y, lat, lon);
fprintf('timezone: %g h (expect 8; pipeline negates -> UTC-8, US Pacific)\n', timezone(lon));
[clat,clon] = projinv(p, hdr.RasterReference.XWorldLimits([1 2 2 1]), ...
    hdr.RasterReference.YWorldLimits([1 1 2 2]));
fprintf('corners: lat %.4f..%.4f, lon %.4f..%.4f\n', ...
    min(clat), max(clat), min(clon), max(clon));

% 3. TopoSunAngle, projected case (hits the ported per-pixel projinv)
topo.slope  = slope;
topo.hdr    = hdr;
topo.aspect = GetTopography(topofile,'aspect');
topo.topofile = topofile;
T = TopoSunAngle(datenum([2019 4 15 20 0 0]), topo); % 20Z = 13:00 PDT
fprintf('TopoSunAngle: mu0 %.3f..%.3f, mu %.3f..%.3f (sun up, expect ~0.8 max)\n', ...
    min(T.mu0(:)), max(T.mu0(:)), min(T.mu(:)), max(T.mu(:)));

% 4. load_coarse_topo -> reprojectRaster (geographic in, projected out)
topo.dem = GetTopography(topofile,'elevation');
ct = load_coarse_topo(ltfile, {'Z','aspect','slope'}, topo);
fprintf('coarse Z reprojected to fine grid: %s, NaNs %d/%d\n', ...
    mat2str(size(ct.Z)), nnz(isnan(ct.Z)), numel(ct.Z));
fprintf('Zdiff: %.0f..%.0f m (summit expect ~ +2300)\n', ...
    min(ct.Zdiff(:)), max(ct.Zdiff(:)));
disp('PORT VERIFY DONE');

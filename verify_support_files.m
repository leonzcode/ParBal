% VERIFY_SUPPORT_FILES Structural check of the two built Stage-1 support files.
% Reads them exactly as ParBal does (include_vars_melt.m / load_coarse_topo.m).
lcfile = 'E:\ucsb\data\ParBal\Shasta\inputs\Shasta_landcover.h5';
ltfile = 'E:\ucsb\data\ParBal\Shasta\GLDAS_topo.h5';

% --- landcover, as include_vars_melt.m:85-86 reads it ---
LandCover.Z  = h5read(lcfile,'/Grid/Z');
LandCover.cc = h5read(lcfile,'/Grid/cc');
fprintf('landcover: Z %s %s, cc %s %s\n', mat2str(size(LandCover.Z)), ...
    class(LandCover.Z), mat2str(size(LandCover.cc)), class(LandCover.cc));
FOREST = makeFOREST(LandCover);
fprintf('FOREST: cc %s, tau values %s, dec=%d con=%d bare=%d\n', ...
    mat2str(size(FOREST.cc)), mat2str(unique(FOREST.tau)'), ...
    nnz(FOREST.type.num_val==1), nnz(FOREST.type.num_val==2), ...
    nnz(FOREST.type.num_val==0));

% --- coarse topo, as load_coarse_topo.m:34-40 reads it ---
Z      = h5read(ltfile,'/Grid/Z');
slope  = h5read(ltfile,'/Grid/slope');
aspect = h5read(ltfile,'/Grid/aspect');
R      = h5readatt(ltfile,'/Grid','ReferencingMatrix');
fprintf('coarse topo: Z %s %s, slope %s, aspect %s\n', mat2str(size(Z)), ...
    class(Z), mat2str(size(slope)), mat2str(size(aspect)));
fprintf('RefMatrix:\n'); disp(R);
% world coords from refmat math: [x y] = [row col 1]*R
fprintf('cell (1,1) center lon/lat = %.3f %.3f (expect -122.875 41.875)\n', [1 1 1]*R);
fprintf('cell (4,6) center lon/lat = %.3f %.3f (expect -121.625 41.125)\n', [4 6 1]*R);
fprintf('Z range %.0f..%.0f m\n', min(Z(:)), max(Z(:)));
disp('VERIFY OK');

function run_shasta_stage1(days)
% RUN_SHASTA_STAGE1  Driver for ParBal Stage 1 (downscale_energy) on Shasta.
% Writes one YYYYMMDD.mat (containing M + SWE) per processed day into outdir.
% All inputs present and R2026a port done as of 2026-06-05.
%
%   run_shasta_stage1        runs the full water year (days 1:365, ~40 min)
%   run_shasta_stage1(197)   runs just day 197 = Apr 15 2019 (~40 s; for demos)

addpath('D:\code\ucsb\ParBal');

% ---- paths (all present) ----
sFile        = 'D:\code\ucsb\data\ParBal\Shasta\inputs\fsca_Shasta_2019.h5';
topofile     = 'D:\code\ucsb\data\ParBal\Shasta\inputs\ShastaTopography.h5';
landcoverfile= 'D:\code\ucsb\data\ParBal\Shasta\inputs\Shasta_landcover.h5';   % built by build_shasta_support_files.py
ldas_dir     = 'D:\code\ucsb\data\ParBal\Shasta\GLDAS';                        % .nc4 in \YYYY\DDD\
ldas_topo    = 'D:\code\ucsb\data\ParBal\Shasta\GLDAS\GLDAS_topo.h5';          % built by build_shasta_support_files.py; MUST live in a dir whose path contains 'GLDAS' (include_vars_melt.m:107 keys the variable list off it)
outdir       = 'D:\code\ucsb\data\ParBal\Shasta\energy';                       % output .mat folder

fast_flag    = true;     % only solve for melt energy M (what reconstruction needs)
LDASOnlyFlag = true;     % GLDAS-only: CERES/MERRA ignored
metvars_flag = true;     % also save hourly Ta/direct/diffuse/Lin/albedo/... (diagnostics vs UCSB forcings)
if nargin < 1
    days = 1:365;        % day index into the fSCA cube (1 = Oct 1 2018; 197 = Apr 15 2019)
end

% ---- pre-flight: report any missing inputs instead of a cryptic crash ----
need = {sFile,'fSCA'; topofile,'topo'; landcoverfile,'land cover'; ldas_topo,'GLDAS topo'};
missing = false;
for i = 1:size(need,1)
    if exist(need{i,1},'file')~=2
        fprintf('MISSING (%s): %s\n', need{i,2}, need{i,1}); missing = true;
    end
end
if exist(ldas_dir,'dir')~=7
    fprintf('MISSING (GLDAS dir): %s\n', ldas_dir); missing = true;
end
if missing
    error('run_shasta_stage1:inputs','Acquire the inputs above before running Stage 1.');
end
if ~exist(outdir,'dir'), mkdir(outdir); end

% ---- run ----
for d = days
    downscale_energy(d, sFile, topofile, landcoverfile, ...
        ldas_dir, ldas_topo, '', '', '', '', fast_flag, outdir, LDASOnlyFlag, metvars_flag);
end
fprintf('Stage 1 done -> %s\n', outdir);
end

function run_shasta_stage2(rFile)
% RUN_SHASTA_STAGE2  Driver for ParBal Stage 2 (reconstructSWE) on Shasta.
% Reads the daily YYYYMMDD.mat files written by Stage 1 and reconstructs the
% SWE/melt cube. Run AFTER run_shasta_stage1 has populated the energy folder.
%
%   run_shasta_stage2                runs with the default output file below
%   run_shasta_stage2('...demo.h5')  writes to a different file (e.g. for demos,
%                                    so existing results are not overwritten)

addpath('D:\code\ucsb\ParBal');

poolsize   = 4;                                                       % parallel workers
energy_dir = 'D:\code\ucsb\data\ParBal\Shasta\energy';                         % .mat files from Stage 1
sFile      = 'D:\code\ucsb\data\ParBal\Shasta\inputs\fsca_Shasta_2019.h5';     % the fSCA cube
if nargin < 1
    rFile = 'D:\code\ucsb\data\ParBal\Shasta\outputs\my_reconstruction_Shasta_2019.h5';
end

% ---- pre-flight: make sure Stage 1 actually ran ----
if exist(sFile,'file')~=2,  error('fSCA missing: %s', sFile); end
if exist(energy_dir,'dir')~=7, error('energy dir missing (run Stage 1 first): %s', energy_dir); end
nmat = numel(dir(fullfile(energy_dir,'*.mat')));
if nmat==0, error('no .mat files in %s (run Stage 1 first)', energy_dir); end
fprintf('found %d daily energy files in %s\n', nmat, energy_dir);
[od,~,~] = fileparts(rFile); if ~exist(od,'dir'), mkdir(od); end

% ---- run ----
% A/B tested 2026-06-05 vs the official reconstruction (Apr 1 2019 SWE):
%   without canopycoverfile: corr 0.978, bias  -3.0%, RMSE 136 mm   <- matches official
%   with    canopycoverfile: corr 0.975, bias +26.2%, RMSE 237 mm
% i.e. UCSB evidently did NOT apply the viewable-gap correction fsca/(1-cc)
% for the Shasta product, so the default here is without. To enable:
%   ccfile = 'D:\code\ucsb\data\ParBal\Shasta\inputs\Shasta_landcover.h5';
%   reconstructSWE(poolsize, energy_dir, sFile, rFile, 'canopycoverfile', ccfile);
reconstructSWE(poolsize, energy_dir, sFile, rFile);
fprintf('Stage 2 done -> %s\n', rFile);
end

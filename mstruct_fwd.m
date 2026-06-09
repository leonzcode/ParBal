function [x,y] = mstruct_fwd(mstruct,lat,lon)
%MSTRUCT_FWD mfwdtran replacement that works before and after R2026a.
%   [x,y] = MSTRUCT_FWD(mstruct,lat,lon) forward-projects lat/lon to map
%   x/y. Uses legacy mfwdtran where it still works; on R2026a (where
%   mfwdtran is a stub that only throws "has been removed" - note exist()
%   still returns 2 for it) falls back to projfwd via mstruct2projcrs.
try
    [x,y] = mfwdtran(mstruct,lat,lon);
catch
    [x,y] = projfwd(mstruct2projcrs(mstruct),lat,lon);
end
end

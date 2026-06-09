function [lat,lon] = mstruct_inv(mstruct,x,y)
%MSTRUCT_INV minvtran replacement that works before and after R2026a.
%   [lat,lon] = MSTRUCT_INV(mstruct,x,y) inverse-projects map x/y to
%   lat/lon. Uses legacy minvtran where it still works; on R2026a (where
%   minvtran is a stub that only throws "has been removed" - note exist()
%   still returns 2 for it) falls back to projinv via mstruct2projcrs.
try
    [lat,lon] = minvtran(mstruct,x,y);
catch
    [lat,lon] = projinv(mstruct2projcrs(mstruct),x,y);
end
end

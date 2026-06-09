function R = makerefmat(x11,y11,dx,dy)
%MAKEREFMAT Construct an affine referencing matrix (R2026a shim).
%   R = MAKEREFMAT(x11,y11,dx,dy) reproduces the 4-scalar form of the
%   Mapping Toolbox function removed in R2026a. ParBal passes referencing
%   matrices (not raster reference objects) throughout, so the original
%   semantics are kept rather than converting every caller.
%
%   (x11,y11) is the world position of the center of pixel (1,1); dx and
%   dy are the column and row spacing. Mapping: [x y] = [row col 1] * R.
%
%   On MATLAB releases that still ship makerefmat, this file shadows the
%   toolbox version with identical results for the 4-scalar form (the only
%   form ParBal uses).
narginchk(4,4)
assert(isscalar(x11) && isscalar(y11) && isscalar(dx) && isscalar(dy), ...
    'only the 4-scalar form of makerefmat is supported by this shim')
R = [0 dy; dx 0; x11-dx, y11-dy];
end

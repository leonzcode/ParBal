function p = mstruct2projcrs(mstruct)
%MSTRUCT2PROJCRS Convert a legacy projection structure to a projcrs object.
%   p = MSTRUCT2PROJCRS(mstruct) builds the projcrs equivalent of a legacy
%   mstruct (defaultm) projection structure, e.g. the ProjectionStructure
%   that GetCoordinateInfo assembles from ParBal HDF5 attributes.
%
%   R2026a removed minvtran/mfwdtran (with the rest of the mstruct
%   transform functions). Replacements:
%       [lat,lon] = minvtran(mstruct,x,y)  ->  [lat,lon] = projinv(p,x,y)
%       [x,y] = mfwdtran(mstruct,lat,lon)  ->  [x,y] = projfwd(p,lat,lon)
%
%   Supports the projections that occur in ParBal data; extend the switch
%   below with another WKT template if a new one shows up.

% ellipsoid from mstruct.geoid = [semimajor_axis eccentricity]
a = mstruct.geoid(1);
ecc = mstruct.geoid(2);
if ecc == 0
    rf = 0; % WKT convention for a sphere
else
    rf = 1/(1 - sqrt(1 - ecc^2)); % inverse flattening
end

if numel(mstruct.origin) >= 3 && mstruct.origin(3) ~= 0
    error('mstruct2projcrs:rotation','nonzero grid rotation not supported');
end
fe   = mstruct.falseeasting;
fn   = mstruct.falsenorthing;
lat0 = mstruct.origin(1);
lon0 = mstruct.origin(2);

geogcs = sprintf(['GEOGCS["GCS_from_mstruct",DATUM["D_from_mstruct",', ...
    'SPHEROID["from_mstruct",%.10g,%.12g]],PRIMEM["Greenwich",0.0],', ...
    'UNIT["Degree",0.0174532925199433]]'], a, rf);

switch mstruct.mapprojection
    case {'eqaconicstd','eqaconic'} % (Albers) Equal-Area Conic, standard parallels
        sp = mstruct.mapparallels;
        if isscalar(sp), sp = [sp sp]; end
        proj = sprintf(['PROJECTION["Albers"],', ...
            'PARAMETER["False_Easting",%.10g],PARAMETER["False_Northing",%.10g],', ...
            'PARAMETER["Central_Meridian",%.10g],', ...
            'PARAMETER["Standard_Parallel_1",%.10g],', ...
            'PARAMETER["Standard_Parallel_2",%.10g],', ...
            'PARAMETER["Latitude_Of_Origin",%.10g]'], ...
            fe, fn, lon0, sp(1), sp(2), lat0);
    case 'tranmerc' % Transverse Mercator (e.g. geotiff2mstruct of UTM tiffs)
        proj = sprintf(['PROJECTION["Transverse_Mercator"],', ...
            'PARAMETER["False_Easting",%.10g],PARAMETER["False_Northing",%.10g],', ...
            'PARAMETER["Central_Meridian",%.10g],PARAMETER["Scale_Factor",%.10g],', ...
            'PARAMETER["Latitude_Of_Origin",%.10g]'], ...
            fe, fn, lon0, mstruct.scalefactor, lat0);
    otherwise
        error('mstruct2projcrs:unsupported', ...
            'projection ''%s'' not implemented yet - add a WKT template for it in mstruct2projcrs.m', ...
            mstruct.mapprojection);
end

wkt = sprintf('PROJCS["%s_from_mstruct",%s,%s,UNIT["Meter",1.0]]', ...
    mstruct.mapprojection, geogcs, proj);
p = projcrs(wkt);
end

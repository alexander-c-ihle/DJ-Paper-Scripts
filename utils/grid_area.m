% compute area of pixels on a lat/lon grid

function [da,dx,dy] = grid_area(lon,lat,dlon,dlat)

% expand to 2d
[LON,LAT] = meshgrid(lon,lat);

% radius of Earth
a = 6370.661745873249136e3; 

% get lengths
dy = (a*pi/180)*dlat;
dx = (a*cos(LAT*pi/180)*pi/180)*dlon;

% area
da = dy*dx;
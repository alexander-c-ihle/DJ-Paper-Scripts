function [mu,var,n] = bin2d(x,y,d,derr,X,Y)

% bins data on a 2-d grid with [X,Y] coordinates
% output:
%   mu, the mean of observations at each grid point
%   var, the variance of observations at each grid point
%   n, the number of observations at each grid point
% input:
%   X, Y, are 2-d matrices defining the grid produced by meshgrid
%   x, y, are 1-d objects with the x and y coordinates of the data
%   d (also a 1-d object) is the original data 
%   derr is the 1 std. dev. estimate

% total number of observations
nobs = length(d);

% fix points outside of grid
x(x<X(1))=X(1); x(x>X(end))=X(end);
y(y<Y(1))=Y(1); y(y>Y(end))=Y(end);

% grid size
[ny,nx,nz] = size(X);
m = prod(size(X));

% bin indices in the horizontal
indx = zeros(ny,nx);
indx(:) = 1:m;
Q.indx = interp2(X,Y,indx,x,y,'nearest');

% make binning operator
ikeep = find(~isnan(Q.indx));
BIN = sparse(Q.indx(ikeep),ikeep,ones(length(ikeep),1),m,length(Q.indx));
Q.BIN = BIN;

% set up variables to receive binned data
mu = zeros(ny,nx);
n = zeros(ny,nx);
var = zeros(ny,nx);

% bin the data
n(:) = Q.BIN*ones(nobs,1);
mu(:) = Q.BIN*d./n(:);
var(:) = Q.BIN*(d-(Q.BIN'*mu(:))).^2./n(:) +...
    Q.BIN*derr.^2./n(:);
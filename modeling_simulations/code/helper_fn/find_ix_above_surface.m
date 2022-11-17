% find surface normal corresponding to shortest distance per point
% estimate min distance from each point to surface

function [above_surface,N_nearest,min_dist,min_ix] = ...
    find_ix_above_surface(Xq,Yq,Zq,x,y,z)

x = x(:);
y = y(:);
z = z(:);

% get normal vectors to surface
[Nx,Ny,Nz] = surfnorm(Xq,Yq,Zq);

all_dist = sqrt((x' - Xq(:)).^2 + (y' - Yq(:)).^2 + ...
    (z' - Zq(:)).^2);
[min_dist,min_ix] = min(all_dist,[],1);

N_nearest = [Nx(min_ix)',Ny(min_ix)',Nz(min_ix)'];

% dot product
above_surface = (diag(N_nearest*([x,y,z] - ...
    [Xq(min_ix)',Yq(min_ix)',Zq(min_ix)'])') > 0)';

end
function [col,d] = trace_single_ray_2d(x,y,xcur,ycur,a,b,isDisp)

x0 = x(1);
x1 = x(end);
y0 = y(1);
y1 = y(end);

nx = length(x)-1;
ny = length(y)-1;

dx = (x1-x0)/nx;
dy = (y1-y0)/ny;

% Use the parametrisation of line to get the y-coordinates of
% intersections with x = k, i.e. x constant.
tx = (x - xcur)/a;
yx = b*tx + ycur;

% Use the parametrisation of line to get the x-coordinates of
% intersections with y = k, i.e. y constant.
ty = (y - ycur)/b;
xy = a*ty + xcur;

% Illustration of the rays
if isDisp
    
    plot(x,yx,'-','color',[220 0 0]/255,'linewidth',1.5)
    plot(xy,y,'-','color',[220 0 0]/255,'linewidth',1.5)
    
    %set(gca,'Xticklabel',{})
    %set(gca,'Yticklabel',{})
    pause(isDisp)
end

% Minimum and maximum times
if isempty(tx)
    mintx = nan;
    maxtx = nan;
else
    mintx = min(tx(1),tx(end));
    maxtx = max(tx(1),tx(end));
end
if isempty(ty)
    minty = nan;
    maxty = nan;
else
    minty = min(ty(1),ty(end));
    maxty = max(ty(1),ty(end));
end
tmin = max([mintx, minty]);
tmax = min([maxtx, maxty]);

% Collect the intersection times and coordinates.
t = [tx; ty];
xxy = [x; xy];
yxy = [yx; y];
I_inside = t>=tmin & t<=tmax;

% Only keep times inside
t = t(I_inside);
xxy = xxy(I_inside);
yxy = yxy(I_inside);

% Sort the coordinates according to intersection time.
% Possible speed-up: Do merge sort since tx and ty are already sorted.
[t, I] = sort(t);
xxy = xxy(I);
yxy = yxy(I);

% Skip double points. Modified to check only difference in t, instead of
% both x and y directions. This check can lead to small differences from
% old implementation.
%I = (abs(diff(xxy)) <= 1e-10 & abs(diff(yxy)) <= 1e-10);
%I = (abs(diff(xxy)) <= 1e-10 & abs(diff(yxy)) <= 1e-10);
I = diff(t) <= 1e-10;
xxy(I) = [];
yxy(I) = [];

% Calculate the lengths within cells.
% Square-root is most expensive operation. Possible speed-up ideas:
% Approximate square-root by fewer power series terms. Or do a basis change
% by rotating the array to the direction of the ray, so that path length is
% calculable along single dimension.
d = sqrt(diff(xxy,1,1).^2 + diff(yxy,1,1).^2);

% If the ray is on the boundary of the box in the top or to the
% right the ray does not by definition lie with in a valid cell.
if ~((b == 0 && (abs(ycur - y0) < 1e-15 ...
        || abs(ycur - y1) < 1e-15)) || ...
        (a == 0 && (abs(xcur - x0) < 1e-15 || ...
        abs(xcur - x1) < 1e-15))   )
    
    % Calculates the midpoints of the line within the cells.
    xm = 0.5*(xxy(1:end-1)+xxy(2:end)) - x0;
    ym = 0.5*(yxy(1:end-1)+yxy(2:end)) - y0;
    
    % Translate the midpoint coordinates to index.
    col = floor(xm/dx)*ny + (ny - floor(ym/dy));
end


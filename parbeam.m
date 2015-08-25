function A = parbeam(conf)
%PARBEAM Build 2D parallel beam tomography system matrix
%
% x0,x1,nx,y0,y1,ny,r0,r1,nr,theta,isDisp
%
%   A = parbeam(conf)
%
% This function builds a sparse matrix corresponding to a discretized
% forward operator for a 2D parallel beam computed tomography setup. The
% discretization is done using the line intersection method, also sometimes
% referred to as Siddon's method.
%
% The single input is a struct with some or all of the fields described
% below. If any field is unset, a default value is set.
%
% Input fields of the struct CONF:
%   x           2-entry vector with left and right bounds for physical
%               object, in physical units. First entry must be smaller.
%               Default: [-1,1].
%   y           2-entry vector with bottom and top bounds for physical
%               object, in physical units. First entry must be smaller.
%               Default: [-1,1].
%   dims        2-entry vector with integers stating the number of
%               discretization intervals (pixel divisions) in the x and y
%               directions, respectively.
%               Default: [100,100].
%   r           2-entry vector with left and right bounds for detector when
%               aligned with projection direction parallel to y-axis.
%               Physical units, first entry must be smaller.
%               Default: [-1,1].
%   nr          Integer stating the number of discretization intervals of
%               the detector.
%               Default: 100.
%   theta       Vector containing the angles of the projection directions
%               in degrees measured in the positive direction from the
%               y-axis.
%               Default: 0:1:179.
%               each dimesion, such that the domain consists of N^2 cells.
%   theta       Vector containing the angles in degrees. Default: theta =
%               0:1:179.
%   isDisp      Non-negative scaler. If non-zero, an illustration of the
%               x-ray paths is shown. The numeric value gives the time in
%               seconds to pause between showing individual rays.
%               Default: 0.
%
% Output:
%   A           CT system matrix with prod(conf.dims) columns and
%               length(theta)*conf.nr rows. Each row contains the length of
%               a given ray's intersection with all pixels, including zeros
%               at non-intersected pixels.

% Jakob Sauer Joergensen, DTU Compute, 2014-01-31.
%
% Modified from paralleltomo from package AIRtools.
% Jakob Heide Joergensen, Maria Saxild-Hansen and Per Christian Hansen,
% June 21, 2011, DTU Informatics.
%
% Reference: A. C. Kak and M. Slaney, Principles of Computerized
% Tomographic Imaging, SIAM, Philadelphia, 2001.

% Unpack the fields of the conf input struct to variables.
[x0,x1,nx,y0,y1,ny,r0,r1,nr,theta,isDisp] = unpack_conf(conf);

% Detector bin width and pixel side lengths in x and y directions.
dr = abs(r1-r0)/nr;
dy = abs(y1-y0)/ny;
dx = abs(x1-x0)/nx;

% Precompute halves for convenience.
dr2 = 0.5*dr;
dx2 = 0.5*dx;
dy2 = 0.5*dy;

% Define the number of angles.
na = length(theta);

% The starting values both the x and the y coordinates.
xinit = linspace(r0+dr2,r1-dr2,nr)';
yinit = zeros(nr,1);

% The intersection lines.
x = (x0:dx:x1)';
y = (y0:dy:y1)';

% Initialize vectors that contains the row numbers, the column numbers and
% the values for creating the matrix A effiecently.
rows = zeros((nx+ny)*na*nr,1);
cols = rows;
vals = rows;
idxend = 0;

% Prepare for illustration
if isDisp
    AA = rand(nx,ny);
    figure
end

% Translation length is the diagonal
tlen = sqrt((x1-x0)^2+(y1-y0)^2);

% Loop over the chosen angles.
for i = 1:na
    
    % Illustration of the domain
    if isDisp
        clf
        pause(isDisp)
        imagesc(x0+dx2:dx:x1-dx2,y0+dy2:dy:y1-dy2,AA'), colormap gray,
        hold on
        axis xy
        axis equal
        axis([x0-dx x1+dx y0-dy y1+dy])
    end
    
    % All the starting points for the current angle.
    xtheta = cosd(theta(i))*xinit-sind(theta(i))*yinit;
    ytheta = sind(theta(i))*xinit+cosd(theta(i))*yinit;
    
    % Translate starting points outside array. This is a new relative to
    % the version with explicit matrix. This is to have a only positive
    % intersection times in order to remove negative times before the
    % expensive sort operation. Translate by adding multiple of rotated
    % vector (0,tlen).
    xtheta = xtheta - sind(theta(i))*tlen;
    ytheta = ytheta + cosd(theta(i))*tlen;
    
    % The direction vector for all the rays corresponding to the current
    % angle.
    a = sind(theta(i));
    b = -cosd(theta(i));
    
    % Loop over the rays.
    for j = 1:nr
        
        % Intersections of single ray
        [col,d] = trace_single_ray_2d(x,y,xtheta(j),ytheta(j), a, b,isDisp);
        numvals = numel(d);
        
        % Create the indices to store the values to vector for
        % later creation of A matrix.
        idxstart = idxend + 1;
        idxend = idxstart + numvals - 1;
        idx = idxstart:idxend;
        
        % Store row numbers, column numbers and values.
        rows(idx) = (i-1)*nr + j;
        cols(idx) = col;
        vals(idx) = d;
        
    end
end

% Truncate excess zeros.
rows = rows(1:idxend);
cols = cols(1:idxend);
vals = vals(1:idxend);

% Create sparse matrix A from the stored values.
A = sparse(rows,cols,vals,nr*na,nx*ny);



function [x0,x1,nx,y0,y1,ny,r0,r1,nr,theta,isDisp] = unpack_conf(conf)

if ~isfield(conf,'x')
    conf.x = [-1,1];
end
if ~isfield(conf,'y')
    conf.y = [-1,1];
end
if ~isfield(conf,'dims')
    conf.dims = [100,100];
end
if ~isfield(conf,'r')
    conf.r = [-1,1];
end
if ~isfield(conf,'nr')
    conf.nr = 100;
end
if ~isfield(conf,'theta')
    conf.theta = 0:179;
end
if ~isfield(conf,'isDisp')
    conf.isDisp = false;
end

x0 = conf.x(1);
x1 = conf.x(2);
y0 = conf.y(1);
y1 = conf.y(2);
nx = conf.dims(1);
ny = conf.dims(2);
r0 = conf.r(1);
r1 = conf.r(2);
nr = conf.nr;
theta = conf.theta;
isDisp = conf.isDisp;
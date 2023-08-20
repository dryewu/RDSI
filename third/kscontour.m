function H = kscontour(X, varargin)
% Kernel-density-smoothed 2D density plot.
%
% USAGE
%   kscontour(X);
%
%     OR
%
%   kscontour(X, 'ArgName', ArgValue);
% 
%     OR
%
%   H = kscontour(___);
%
% INPUTS
%   X - Datapoints in 2D space (Npoints-by-2 matrix).  If passed a cell
%       array of Xs, will superimpose contours for each.
%
% OPTIONAL INPUTS
%   Edges - edges of histogram bins (cell array w/ 2 vectors of edges)
%   Nbins - number of histogram bins (if Edges not defined)
%   Alpha - alpha of contours
%   Color - color for plot. RGB triplet or string, one of: 'blue',
%       'orange', 'green', 'yellow', 'darkblue', 'vermillion', 'pink',
%       'black', or 'gray'
%   Pad - padding to use around points (proportion of range of data)
%   Nlevels - number of contour levels (default=6)
%   Sigmas - to plot contours at sigma values, pass vector of sigma values
%       at which to plot contours.
%   ShowPoints - whether or not to also plot the datapoints (default=false)
%   PointSize - size of points if showing raw points
%   PrcLo - lower percentile for bound on density estimate (default=2)
%   PrcHi - upper percentile for bound on density estimate (default=98)
%
% OUTPUT
%   H - handle to the contour object(s)
%       
% EXAMPLES
%
%   Plot a density
%     P1 = mvnrnd([0 0], [1 0;0 1], 100);
%     kscontour(P1);
%
%   Plot two densities, specifying the colors
%     P1 = mvnrnd([0 0], [1 0;0 1], 100);
%     P2 = mvnrnd([2 2], [1 0;0 1], 100);
%     kscontour(P1, 'Color', 'green');
%     kscontour(P2, 'Color', 'pink');
%
%   Plot several overlaid densities (with auto-color!)
%     P1 = mvnrnd([0 0], [1 0;0 1], 100);
%     P2 = mvnrnd([2 2], [1 0;0 1], 100);
%     P3 = mvnrnd([-2 2], [1 .7;.7 1], 100);
%     Hs = kscontour({P1, P2, P3});
%     legend(Hs, {'Thing1', 'Thing2', 'Thing3'})
%
%   Plot points on top
%     P1 = mvnrnd([0 0], [1 0;0 1], 100);
%     kscontour(P1, 'ShowPoints', true);
%
% Oct 2018
% Brendan Hasz
% haszx010@umn.edu

% Colors
black = [0 0 0]; %#ok<*NASGU>
blue = [0.35 0.7 0.9];
orange = [0.9,0.6,0];
green = [0 0.6 0.5];
yellow = [0.95 0.9 0.25];
darkblue = [0 0.45 0.7];
vermillion = [0.8, 0.4, 0];
pink = [0.8 0.6 0.7];
gray = 0.7*ones(1,3);
grey = gray;
colors = [blue; orange; green; yellow; pink; darkblue; vermillion; gray];

% Parse optional arguments
p = inputParser;
iscolor = @(x) (isvector(x) && length(x)==3) || ischar(x);
def_PointSize = max(1, min(30, 200/sqrt(length(X))));
addParameter(p, 'Edges', [], @(x) isvector(x) || isempty(x));
addParameter(p, 'Nbins', 100, @isscalar);
addParameter(p, 'Alpha', 0.3, @isscalar);
addParameter(p, 'Color', blue, iscolor);
addParameter(p, 'Pad', 0.5, @isscalar);
addParameter(p, 'Nlevels', 5, @isscalar);
addParameter(p, 'Sigmas', [], @(x) isvector(x) || isempty(x));
addParameter(p, 'ShowPoints', false, @islogical);
addParameter(p, 'PointSize', def_PointSize, @isscalar);
addParameter(p, 'PrcLo', 2, @isscalar);
addParameter(p, 'PrcHi', 98, @isscalar);
parse(p, varargin{:});
Edges = p.Results.Edges;
Nbins = p.Results.Nbins;
Alpha = p.Results.Alpha;
Color = p.Results.Color;
Pad = p.Results.Pad;
Nlevels = p.Results.Nlevels;
Sigmas = p.Results.Sigmas;
ShowPoints = p.Results.ShowPoints;
PointSize = p.Results.PointSize;
PrcLo = p.Results.PrcLo;
PrcHi = p.Results.PrcHi;

% Plot multiple densities if X is cell array
if iscell(X)
    H = nan(1, length(X)); %handles for each
    for iX = 1:length(X)
        H(iX) = kscontour(X{iX}, 'Color', colors(mod(iX-1,8)+1,:), ...
                          'Edges', Edges, 'Nbins', Nbins, ...
                          'Alpha', Alpha, 'Pad', Pad, ...
                          'Nlevels', Nlevels, 'Sigmas', Sigmas);
    end
    return;
end

% Set color
if ischar(Color)
    try 
        eval(['Color=' Color ';']);
    catch err
        error([Color ' is not a valid color string'])
    end
end

% Compute edges if not provided
if isempty(Edges)
    min1 = prctile(X(:,1), PrcLo);
    max1 = prctile(X(:,1), PrcHi);
    dx1 = max1-min1;
    E1 = linspace(min1-dx1*Pad, max1+dx1*Pad, Nbins);
    min2 = prctile(X(:,2), PrcLo);
    max2 = prctile(X(:,2), PrcHi);
    dx2 = max2-min2;
    E2 = linspace(min2-dx2*Pad, max2+dx2*Pad, Nbins);
end

% Compute KDE density
[xi, yi] = meshgrid(E1, E2);
[R, ~] = ksdensity(X, [xi(:), yi(:)]);
Z = reshape(R, length(E1), length(E2));

% Compute levels
if ~isempty(Sigmas) %set contours to sigma values
    oZ = sort(Z(:), 'descend');
    cZ = cumsum(oZ);
    cZ = cZ/cZ(end);
    levels = nan(length(Sigmas)+2, 1);
    levels(1) = oZ(end);
    levels(end) = oZ(1);
    for iS = 1:length(Sigmas)
        ts = 1-2*(1-SigmaToPerc(Sigmas(iS)));
        ti = find(cZ>ts, 1); 
        levels(iS+1) = oZ(ti);
    end
else %just use evenly-split levels
    levels = linspace(min(Z(:)), max(Z(:)), Nlevels+2);
end

% Get contour data
f = figure('visible', 'off'); %just want the contour data, not a plot!
[C, ~] = contourf(E1, E2, Z, levels);
close(f);

% Plot each contour
iC = 2;
tC = 0;
while iC<size(C,2)
    tn = C(2, iC-1);
    ix = iC:(iC+tn-1); %inds for this contour's values
    H = patch(C(1,ix), C(2,ix), max(0, (Color*256-2*tC)/256), ...
              'FaceAlpha', Alpha, 'EdgeColor', 'none');
    hold on
    iC = iC + tn + 1;
    tC = tC + 1;
end

% Show raw points?
if ShowPoints
    scatter(X(:,1), X(:,2), PointSize, 'filled' , ...
            'MarkerFaceColor', Color);
end

end

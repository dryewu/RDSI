function varargout = ploterr(X, Y, varargin)
% Plot data with shaded error bounds.  Removes NaN values.
%
% USAGE
%   H = ploterr(X, Y, 'ArgName', ArgValue);
% 
%     OR
%
%   H = ploterr(X, Y, E, 'ArgName', ArgValue);
%
%     OR
%
%   H = ploterr(X, Y, Epos, Eneg, 'ArgName', ArgValue);
%
%     OR
%
%   [H, M, Epos, Eneg] = ploterr(___);
%
% INPUTS
%   X - X values of samples in Y (Ntimepoints-length vector).  If Group 
%       (see below) is set to true, then treat X as a grouping variable
%       instead.
%   Y - sample values (Nsamples-by-Ntimepoints). Will plot error across the
%       1st dimension of Y.  If passing (X, Y, E), Y is same size as X. If
%       using Group=true, Y is same size as X.
%
% OPTIONAL INPUTS
%   E - error to plot in shaded region (X, Y, and E must be same size)
%   Epos - error to plot in shaded region above mean line (same size as X)
%   Eneg - error to plot in shaded region below mean line (same size as X)
%   Group - if set to false, treat X as timepoints (default), but if set to
%       true, treat X as a grouping variable.  If Group==true then X must
%       be same size as Y.
%   Error - type of error to use. 'sem' for standard error of the mean 
%       (default) or 'std' for standard deviation.  'prc' for median and
%       1sigma (68.3% percentile) (use PercTile to adjust bounds).
%   Labels - labels for groups when Group==true
%   Color - color for the lines and shaded error area.  Can pass a vector
%       of length 3 (an RGB triplet) or a string with a color name, 
%       one of: 'black', 'blue', 'orange', 'green', 'yellow', 'darkblue', 
%       'vermillion', 'pink', 'gray', 'grey'
%   LineWidth - width of the mean line
%   Alpha - transparency of the shaded error area (scalar from 0 to 1)
%   ShowPoints - whether or not to show individual points. Default is false
%   PointColor - color for the individual points if shown.  Set in the same
%       way as the Color option above. Default=gray
%   MarkerSize - size for the individual points if shown. Default=5
%   LineSpec - line specification string for the average line (e.g. 'o--')
%   PercTile - percentile bounds when using `prc` as Error (default = 95)
%
% OUTPUT
%   H - handle to the mean line
%   M - mean values
%   E - error which was plotted
%   Epos - positive error which was plotted
%   Eneg - negative error which was plotted
%
% EXAMPLES
%
%   Plot a line with error:
%     x = linspace(0, 1, 10);
%     y = 2*x+0.2*randn(size(x));
%     e = 0.3+0.3*rand(size(x));
%     ploterr(x, y, e)
%
%   Plot a line with separate positive and negative error:
%     Epos = 0.3+0.3*rand(size(x));
%     Eneg = 0.6+0.3*rand(size(x));
%     ploterr(x, y, Epos, Eneg)
%
%   Plot grouped data mean and SEM:
%     x = ceil(3*rand(20,1)); %grouping variable
%     y = randn(size(x));
%     y(x==1) = y(x==1) + 1;
%     y(x==2) = y(x==2) - 1;
%     ploterr(x, y, 'Group', true, 'Labels', {'A', 'B', 'C'})
%
%   Plot continuous-time data mean and SEM:
%     Nt = 10; %number of timepoints
%     Ns = 5;  %number of samples/trials
%     x = linspace(0, 1, Nt); %timepoints
%     y = 2*repmat(x, Ns, 1)+0.5*randn(Ns,Nt); %Nsamples-by-Ntimepoints
%     ploterr(x, y)
%
%   Show individual points
%     ploterr(x, y, 'ShowPoints', true)
%
%   Plot standard deviation instead of SEM
%     ploterr(x, y, 'Error', 'std')
%
%   Plot median and 68% interval (1 sigma of true distribution)
%     ploterr(x, y, 'Error', 'prc')
%
%   Plot median and 95% interval
%     ploterr(x, y, 'Error', 'prc', 'PercTile', 95)
%
%   Adjust transparency of shaded area, line width, colors, and linespec
%     ploterr(x, y, ...
%       'Color', 'orange', 'LineWidth', 5, 'LineSpec', '--', ...
%       'Alpha', 0.5, 'ShowPoints', true, 'PointColor', 'pink', ...
%       'MarkerSize', 20)
%
%   Plot two lines with a legend.
%     h1 = ploterr(x, y);
%     h2 = ploterr(x, -y, 'Color', 'orange');
%     legend([h1 h2], {'Line 1', 'Line 2'})
%
%   Plot several lines, auto-generating new colors for each:
%     x = linspace(0, 1, 10);
%     for iL = 1:5
%       y = iL+randn(length(x));
%       ploterr(x, y, 'Color', iL)
%     end
%       
% Apr 2018
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
gray = 0.5*ones(1,3);
grey = gray;
colors = [blue; orange; green; yellow; pink; darkblue; vermillion; gray];

% Were we passed (x, y, error, varargin)?
if ~isempty(varargin) && ~ischar(varargin{1})
    E = varargin{1};
    varargin = varargin(2:end);
else
    E = [];
end

% Were we passed (x, y, errorPositive, errorNegative, varargin)?
if ~isempty(varargin) && ~ischar(varargin{1})
    Eneg = varargin{1};
    varargin = varargin(2:end);
else
    Eneg = E;
end

% Parse optional arguments
p = inputParser;
iscolor = @(x) (isvector(x) && length(x)==3) || ischar(x) || isscalar(x);
addParameter(p, 'Group', false, @islogical);
addParameter(p, 'Error', 'sem', @ischar);
addParameter(p, 'PercTile', 100*erf(1/sqrt(2)), @isscalar);
addParameter(p, 'Labels', {}, @iscell);
addParameter(p, 'ShowPoints', false, @islogical);
addParameter(p, 'PointColor', gray, iscolor);
addParameter(p, 'MarkerSize', 5, @isscalar);
addParameter(p, 'Color', blue, iscolor);
addParameter(p, 'LineWidth', 2, @isscalar);
addParameter(p, 'Alpha', 0.3, @isscalar);
addParameter(p, 'LineSpec', '-', @ischar);
parse(p, varargin{:});
Group = p.Results.Group;
Error = p.Results.Error;
PercTile = p.Results.PercTile;
Labels = p.Results.Labels;
ShowPoints = p.Results.ShowPoints;
PointColor = p.Results.PointColor;
MarkerSize = p.Results.MarkerSize;
Color = p.Results.Color;
LineWidth = p.Results.LineWidth;
Alpha = p.Results.Alpha;
LineSpec = p.Results.LineSpec;

% Set color
if ischar(Color)
    try 
        eval(['Color=' Color ';']);
    catch err
        error([Color ' is not a valid color string'])
    end
end
if ischar(PointColor)
    try
        eval(['PointColor=' PointColor ';']);
    catch err
        error([PointColor ' is not a valid color string'])
    end
end
if isnumeric(Color) && isscalar(Color)
    Color = colors(round(Color), :);
end
if isnumeric(PointColor) && isscalar(PointColor)
    PointColor = colors(round(PointColor), :);
end
    
% Plot just (X, Y, E) or (X, Y, Epos, Eneg)
if ~isempty(E)
    assert(all(size(X)==size(Y)), ... 
        'X and Y must be same size when using Error');
    assert(all(size(X)==size(E)), ... 
        'X and E must be same size when using Error');
    fill([X(:); flipud(X(:))], ... plot error bounds
        [Y(:)+E(:); flipud(Y(:)-Eneg(:))], ...
        Color, 'EdgeColor', 'none', 'facealpha', Alpha)
    hold on
    H = plot(X(:), Y(:), LineSpec, ... plot mean line
             'Color', Color, 'LineWidth', LineWidth);
    errP = E;
    errN = Eneg;
    
% Plot with X as grouping variable, not timepoints
elseif Group
    uX = unique(X); %#ok<*UNRCH>
    assert(all(size(X)==size(Y)), ... 
        'X and Y must be same size when using Group');
    m = nan(1,length(uX));
    errP = nan(1,length(uX));
    errN = nan(1,length(uX));
    for iX = 1:length(uX)
        tY = Y(X==uX(iX)); %iX'ths group Y values
        if strcmp(Error, 'prc') %use median + confidence intervals
            m(iX) = nanmedian(tY); %median of the samples
            errP(iX) = prctile(tY, 100-(100-PercTile)/2) - m(iX);
            errN(iX) = m(iX) - prctile(tY, (100-PercTile)/2);
        else
            m(iX) = nanmean(tY); %mean of the samples
            if strcmp(Error, 'std') %use standard deviation
                errP(iX) = nanstd(tY);
                errN(iX) = errP(iX);
            elseif strcmp(Error, 'sem') %use standard error of the mean
                errP(iX) = nanstd(tY)/sqrt(size(tY,1));
                errN(iX) = errP(iX);
            else
                error(['Unimplemented Error type ' Error])
            end
        end
        if ShowPoints
            if iX==length(uX)
                tX_s = uX(iX)-uX(iX-1);
            else
                tX_s = uX(iX+1)-uX(iX);
            end
            plot(uX(iX)+0.5*tX_s*(betarnd(3,3,length(tY),1)-0.5), tY, ...
                '.', 'Color', PointColor, 'MarkerSize', MarkerSize);
            hold on
        end
    end
    fill([uX(:); flipud(uX(:))], ... plot error bounds
        [m-errN fliplr(m+errP)], ...
        Color, 'EdgeColor', 'none', 'facealpha', Alpha)
    hold on
    H = plot(uX, m, LineSpec, ... plot mean line
             'Color', Color, 'LineWidth', LineWidth);
    xticks(uX)
    if ~isempty(Labels)
        xticklabels(Labels)
    end
    
% Plot with X as timepoints, not grouping variable
else 
    if strcmp(Error, 'prc') %use median + confidence intervals
        m = nanmedian(Y); %median of the samples
        errP = prctile(Y, 100-(100-PercTile)/2) - m;
        errN = m - prctile(Y, (100-PercTile)/2);
    else
        m = nanmean(Y); %mean of the samples
        if strcmp(Error, 'std') %use standard deviation
            errP = nanstd(Y);
            errN = errP;
        elseif strcmp(Error, 'sem') %use standard error of the mean
            errP = nanstd(Y)/sqrt(size(Y,1));
            errN = errP;
        else
            error(['Unimplemented Error type ' Error])
        end
    end
    if ShowPoints
        tX_s = diff(X(:));
        tX_s = [tX_s; tX_s(end)];
        xx = X(:)+0.5*tX_s.*(betarnd(3,3,length(X),size(Y,1))-0.5);
        plot(xx, Y', ...
            '.', 'Color', PointColor, 'MarkerSize', MarkerSize);
        hold on
    end
    errP(isnan(errP)) = 0; %set NaN error to 0
    errN(isnan(errN)) = 0;
    k = ~isnan(m); %skip plotting nan samples
    m = m(k);
    X = X(k);
    errP = errP(k); 
    errN = errN(k);
    fill([X(:); flipud(X(:))], ... plot error bounds
        [m-errN fliplr(m+errP)], ...
        Color, 'EdgeColor', 'none', 'facealpha', Alpha)
    hold on
    H = plot(X, m, LineSpec, ... plot mean line
             'Color', Color, 'LineWidth', LineWidth);
    if ~isempty(Labels)
        warning('Cannot use Labels when Group==false')
    end
end

% Return handle to line, means, and errors, if requested
if nargout>0
    varargout{1} = H;
end
if nargout>1
    varargout{2} = m;
end
if nargout>2
    varargout{3} = errP;
end
if nargout>3
    varargout{4} = errN;
end

end

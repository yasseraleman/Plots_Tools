function varargout = my_confidence_interval_plot(dataX,dataY,varargin)
%
% Syntax :
%  [FigID,envTop,envBot] = my_confidence_interval_plot(dataY,varargin)
%
% This script creates a confidence interval plot for a given data.
% 3D matrix is asumed as Nsubjects x Npoints x Ngroups.
%
% Input Parameters:
%       dataY                 : Data matrix
%
%
% Output Parameters:
%      FigID                  : Figure handle
%      allStats               : Stats from comparison between data columns.
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% July 22th 2014
% Version $1.0


% bidData = 1;
if nargin ~=0
%     dataSize = size(dataY);
%     dimData = length(dataSize);
%     if dimData > 2
%         bidData = 0;
%     end
    dispLeg = 0; % Display Legend
    
    % Boundary Method
    boundMeth ='std';
    
    % Parameters
    transpVal = .8; % Opacity
    
    % Boolean variables
    reglinePlot = 1; % Plot Regression Line
    pointsPlot = 1;  % Plot Data Points
    jitFill = 1;     % Fill jitter markers
    
    % Figure Properties
    param.figcolor = 'white'; % Color
    param.figurevisible = 'on'; % Visible
    
    set(0,'units','centimeters');
    cm_screen = get(0,'screensize');
    figPosition = [1 1 cm_screen(3)-2 cm_screen(4)-2*cm_screen(3)/cm_screen(4)]; % Position
    
    figUnits = 'centimeters'; % Units
    AxColor = [0.95  0.95 0.95]; % Axis Color
    
    % Axis properties
    xlabelcad = 'X Axis'; % Labels
    ylabelcad = 'Y Axis'; % Labels
    titlecad = ' '; % Title
end

%% ====================== Checking input parameters ===================== %
% deal with the input arguments
if nargin<1 % the indispensable input arguments are not provided
    error('One input is mandatory');
else
    if numel(varargin)>0 % optional input arguments are provided
        while ~isempty(varargin)
            if numel(varargin)<2
                error('You need to provide optional input arguments as ''ParameterName''-''ParameterValue'' pairs.');
            end
            switch varargin{1}
                case 'figType'
                    figType = varargin{2}; % Figure Type: 'subplot', 'indivplot' or 'multiplot'
                case 'transpVal' % Opacity value
                    transpVal=varargin{2};
                case 'figcolor'
                    param.figcolor=varargin{2};
                case 'figurevisible'
                    param.figurevisible=varargin{2};
                case 'figUnits'
                    figUnits=varargin{2};
                case 'figPosition'
                    figPosition=varargin{2};
                case 'FigID'
                    FigID=varargin{2};
                case 'AxID'
                    AxID=varargin{2};
                case 'AxColor'
                    AxColor=varargin{2};
                case 'boundMeth'
                    boundMeth=varargin{2};
                case 'methLabels'
                    methLabels = varargin{2};
                case 'xlabelcad'
                    xlabelcad=varargin{2};
                case 'ylabelcad'
                    ylabelcad=varargin{2};
                case 'titlecad'
                    titlecad=varargin{2};
                case 'dispLeg'
                    dispLeg=varargin{2};
                case 'pointsPlot'
                    pointsPlot=varargin{2};
                case 'pointMarkers'
                    pointMarkers = varargin{2};
                case 'lineStyles'
                    lineStyles   = varargin{2};
                case 'pointColors'
                    pointColors  = varargin{2};
                otherwise
                    error('Unexpected ''ParameterName'' input: %s\n',varargin{1});
            end
            varargin(1:2)=[]; % this pair of optional input arguments has been dealt with -- remove...
        end
    end
end
if ischar(param.figcolor)
    switch param.figcolor
        case 'black'
            param.figcolor = [ 0 0 0];
        case 'white'
            param.figcolor = [ 1 1 1];
    end
end

if exist('AxColor','var')
    if ~isempty(AxColor)&length(AxColor) == 3
        if ~isempty(find(AxColor> 1))
            AxColor = AxColor/255;
        end
    end
end

if ~exist('FigID','var')
    FigID = figure('numbertitle','off','name','New Figure','Color',param.figcolor, 'units',figUnits,'Position',figPosition,'Visible',param.figurevisible,'InvertHardcopy','off');
else
    if ~isempty(FigID)
        figure(FigID);
    else
        FigID = figure('numbertitle','off','name','New Figure','Color',param.figcolor, 'units',figUnits,'Position',figPosition,'Visible',param.figurevisible,'InvertHardcopy','off');
    end
end
if exist('AxID','var')
    if ~isempty(AxID)
        set(FigID,'CurrentAxes',AxID);
        %                 cla;
    else
        clf;
    end
else
    AxID = gca;
end

if ~exist('pointMarkers','var')
    pointMarkers = strvcat('d','o','s','^','d','p','>','v');
end
if ~exist('lineStyles','var')
    lineStyles = strvcat('-',':','--','-.');
end

C = [...
    240,163,255;0,117,220;153,63,0;76,0,92;25,25,25;0,92,49;43,206,72;...
    255,204,153;128,128,128;148,255,181;143,124,0;157,204,0;194,0,136;...
    0,51,128;255,164,5;255,168,187;66,102,0;255,0,16;94,241,242;0,153,143;...
    224,255,102;116,10,255;153,0,0;255,255,128;255,255,0;255,80,5;250 0 0]/255;
Colors = C([27 2 15 20  8 16 1 3 4 5 6 7 9 10 11 12 13 14 17 18 19 21 22 23 24 25 26],:);

Xmin =min(dataX(:));
Xmax =max(dataX(:));
cont = 0;
if ~iscell(dataY);

    valY2plot = mean(dataY);
    switch boundMeth
        case 'std'
            stdValues = std(dataY,[],1);
            valY2plot = [valY2plot;valY2plot - stdValues;valY2plot + stdValues];
        case 'percent'
            valY2plot = [valY2plot;prctile(dataY,[2.5 97.5],1)];
        otherwise
            
    end
    cont = cont + 1;
    if ~exist('methLabels','var')
        methLabels{cont} = 'Method 01';
    end
    tempHandle = plot_ci(dataX,valY2plot', 'PatchColor', [0,102,51]/255, 'PatchAlpha', 0.1, ...
        'MainLineWidth', 2, 'MainLineColor', [0,117,220]/255, ...
        'LineWidth', 1, 'LineStyle','-', 'LineColor', [0,117,220]/255);
    handPlots(cont) = tempHandle.Patch;
    tempEnv = valY2plot;
        
    Ymin =min(valY2plot(:));
    Ymax =max(valY2plot(:));
    
else
    Nplots = length(dataY);
    if ~exist('methLabels','var')
        methLabels = cellstr([repmat('Method ',[Nplots 1]) num2str([1:Nplots]','%.3d')]); % Labels for Methods
    else
        if length(methLabels) < Nplots
            clear methLabels;
            methLabels = cellstr([repmat('Method ',[Nplots 1]) num2str([1:Nplots]','%.3d')]); % Labels for Methods
        end
    end
    
    for i = 1:length(dataY)
        cont = cont + 1;
        tempVar = dataY{i};
%         tempVar = squeeze(dataY(:,:,i));
        valY2plot = mean(tempVar);
        switch boundMeth
            case 'std'
                stdValues = std(tempVar,[],1);
                valY2plot = [valY2plot;valY2plot - stdValues;valY2plot + stdValues];
            case 'percent'
                valY2plot = [valY2plot;prctile(tempVar,[2.5 97.5],1)];
            otherwise
                
        end
        hold on;
        tempHandle= plot_ci(dataX,valY2plot', 'PatchColor', Colors(i,:), 'PatchAlpha', 0.1, ...
            'MainLineWidth', 3, 'MainLineColor', Colors(i,:), ...
            'LineWidth', 1, 'LineStyle','--', 'LineColor', Colors(i,:));
        handPlots(cont) = tempHandle.Patch;
        if i == 1
            Ymin =min(valY2plot(:));
            Ymax =max(valY2plot(:));
            tempEnv = valY2plot;
        else
            Ymin =min([Ymin;valY2plot(:)]);
            Ymax =max([Ymax;valY2plot(:)]);
            tempEnv = [tempEnv;valY2plot];
        end
    end
end
envTop = max(tempEnv);
envBot = min(tempEnv);

set(0,'units','centimeters');
cm_screen = get(0,'screensize');

oldUnits = AxID.Units;
AxID.Units = 'centimeters';
axPos = AxID.Position;
AxID.Units = oldUnits;

oldUnits = FigID.Units;
FigID.Units = 'centimeters';
figPos = FigID.Position;
FigID.Units = oldUnits;


refVal = 55*prod(axPos(3:4))/prod(figPos(3:4));
%                 refVal = min(axPos(3:4)); % Reference Value for the scaling process
fsaxis = refVal*7/10;   % Axis Ticks fontsize
fslabels = refVal*8/10; % Axis Labels fontsize
fstitle = refVal*9/10;  % Axis Title fontsize


bLine = refVal*1.5/14;
markSize = refVal*2.5/6;
medSize = refVal*4/14;


Xrange =  Xmax - Xmin;
Yrange =  Ymax - Ymin;


AxID.XLim = [Xmin Xmax+0*Xrange]; % X Limits
% AxID.XLim = [Xmin-0.05*Xrange Xmax+0.05*Xrange]; % X Limits

AxID.YLim = [Ymin-0.05*Yrange Ymax+0.05*Yrange]; % Y Limits
AxID.FontSize = fsaxis;

AxID.Color = AxColor;
AxID.XColor = [1 1 1]-param.figcolor;
AxID.YColor = [1 1 1]-param.figcolor;
AxID.ZColor = [1 1 1]-param.figcolor;

AxID.Title.String = [''];
AxID.Title.Visible = 'on';
AxID.Title.Color = [1 1 1] -param.figcolor;
AxID.Title.FontSize = fslabels*0.8;
AxID.Title.FontWeight = 'normal';

AxID.XLabel.String = xlabelcad;
AxID.XLabel.Visible = 'on';
AxID.XLabel.Color = [1 1 1] -param.figcolor;
AxID.XLabel.FontSize = fslabels;

AxID.YLabel.String = ylabelcad;
AxID.YLabel.Visible = 'on';
AxID.YLabel.Color = [1 1 1] -param.figcolor;
AxID.YLabel.FontSize = fslabels;

hTitle = suptitle(titlecad);
hTitle.FontName = 'Arial';
hTitle.Visible = 'on';
hTitle.Color = [1 1 1] -param.figcolor;
hTitle.FontSize = fstitle;
grid on;

% Displaying Legend
if dispLeg
    hleg = legend( handPlots, methLabels,'Box', 'on','Color', 'none','FontSize',fslabels*0.7, 'FontName','Helvetica','Color',param.figcolor);
end
hold off;

%% ====================== End of Main Program =========================== %

% Outputs
varargout{1} = FigID;
varargout{2} = envTop;
varargout{3} = envBot;

return;
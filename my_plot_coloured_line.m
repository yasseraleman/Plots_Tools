function varargout = my_plot_coloured_line(Coords,CData,varargin);
%
% Syntax :
%  [strlines] = my_regression_groups_plot(Coords,CData,varargin)
%
% This script plots coloured lines. 
%
% Input Parameters:
%       Coords                 : Data matrix
%       CData                  : Color data for each line point
%
% Output Parameters:
%      strlines                : Lines handles.
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% July 22th 2014
% Version $1.0


if nargin ~=0
    if size(Coords,2) > 3
        error('');
        return;
    end
    colMap = 'jet';
    lineWidth =3;
    lineStyle = '-';
    % Figure Properties
    param.figcolor = 'white'; % Color
    param.figurevisible = 'on'; % Visible
    
    set(0,'units','centimeters');
    cm_screen = get(0,'screensize');
    figPosition = [1 1 cm_screen(3)-2 cm_screen(4)-2*cm_screen(3)/cm_screen(4)]; % Position
    
    figUnits = 'centimeters'; % Units
    
end

%% ====================== Checking input parameters ===================== %
% deal with the input arguments
if nargin<2 % the indispensable input arguments are not provided
    error('Two inputs are mandatory');
else
    if numel(varargin)>0 % optional input arguments are provided
        while ~isempty(varargin)
            if numel(varargin)<2
                error('You need to provide optional input arguments as ''ParameterName''-''ParameterValue'' pairs.');
            end
            switch varargin{1}
                case 'colMap' % Colormap
                    colMap=varargin{2};
                case 'FigID'
                    FigID=varargin{2};
                case 'Linewidth'
                    lineWidth=varargin{2};
                case 'Linestyle'
                    lineStyle=varargin{2};
                case 'AxID'
                    AxID=varargin{2};
                case 'figcolor'
                    param.figcolor=varargin{2};
                case 'figurevisible'
                    param.figurevisible=varargin{2};
                case 'figUnits'
                    figUnits=varargin{2};
                case 'figPosition'
                    figPosition=varargin{2};
                otherwise
                    error('Unexpected ''ParameterName'' input: %s\n',varargin{1});
            end
            varargin(1:2)=[]; % this pair of optional input arguments has been dealt with -- remove...
        end
    end
end

%% ===================== End of Input Parameters ======================== %
%% ========================= Main Program =============================== %

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
    end
else
    AxID = gca;
end
Npoints = size(Coords,1);
poss = [[1:Npoints-1]' [2:Npoints]'];

if size(CData,1) == 1
    ColorsLine = CData; 
end
if size(CData,2) == 1
    ColorsLine = Val2colors(CData,colMap);
elseif size(CData,2) == 3
     ColorsLine = CData; 
    if ~isempty(find(CData>1))
        ColorsLine = CData/255;
    end
end

switch size(Coords,2)
    case 1
        X = [1:Npoints]';
        Y = Coords(:,1);

    case 2
        X = Coords(:,1);
        Y = Coords(:,2);
    case 3
        X = Coords(:,1);
        Y = Coords(:,2);
        Z = Coords(:,3);
        
    otherwise
end
if size(ColorsLine,1) == 1
    switch size(Coords,2)
        case {1,2}
            strlines = plot(X,Y,'Linestyle',lineStyle,'Color',ColorsLine,'Linewidth',lineWidth);
        case 3
            strlines = plot3(X,Y,Z,'Linestyle',lineStyle,'Color',ColorsLine,'Linewidth',lineWidth);
    end
else
    ColorsLine = (ColorsLine(1:end-1,:) + ColorsLine(2:end,:))/2;
    hold on;
    switch size(Coords,2)
        case {1,2}
            for k = 1:size(ColorsLine,1);
                strlines(k) = plot(X(poss(k,:)),Y(poss(k,:)),'Linestyle',lineStyle,'Color',ColorsLine(k,:),'Linewidth',lineWidth);
            end
        case 3
            for k = 1:size(ColorsLine,1);
                strlines(k) = plot3(X(poss(k,:)),Y(poss(k,:)),Z(poss(k,:)),'Linestyle',lineStyle,'Color',ColorsLine(k,:),'Linewidth',lineWidth);
            end
    end
end

%% =================== End of Main Program ============================== %
% Outputs
varargou{1} = strlines;
return
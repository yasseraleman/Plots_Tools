function varargout = my_regression_groups_plot(dataY,varargin)
%
% Syntax :
%  [FigID,outliersMat] = my_regression_groups_plot(dataY,varargin)
%
% This script creates a regression plots for a given data.
% 2D matrix is asumed as Nsubjects x Nmethods
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


bidData = 0;
if nargin ~=0
    dataSize = size(dataY);
%     groups = ones(dataSize(1),1);
    dimData = length(dataSize);
    if dimData ~= 2
        error('This script only works for 2D data');
        return;
    end
    methLabels = cellstr([repmat('Method ',[size(dataY,2) 1]) num2str([1:size(dataY,2)]','%.3d')]); % Labels for Methods
    dispLeg = 0; % Display Legend
    
    figType = 'subplot'; % Plot Type: subplot:   All regressions in the same figure but different subplots.
    %            multiplot: All regressions in the same figure and same subplot.
    %            indivplot: All regressions in the individual figures.
    refColumn = 1;  % Reference column for multiple regressions in the same graph
    allStats = '';
    
%     pointMarkers = strvcat('d','o','s','^','d','p','>','v');
%     lineStyles = strvcat('-',':','--','-.');
%     C = [...
%         240,163,255;0,117,220;153,63,0;76,0,92;25,25,25;0,92,49;43,206,72;...
%         255,204,153;128,128,128;148,255,181;143,124,0;157,204,0;194,0,136;...
%         0,51,128;255,164,5;255,168,187;66,102,0;255,0,16;94,241,242;0,153,143;...
%         224,255,102;116,10,255;153,0,0;255,255,128;255,255,0;255,80,5;250 0 0]/255;
%     pointColors = C([27 2 15 20  8 16 1 3 4 5 6 7 9 10 11 12 13 14 17 18 19 21 22 23 24 25 26],:);
    
    
    % Parameters
    transpVal = .8; % Opacity
    ellip_confPerc = 0; % Ellipese Confidence Interval
    
    
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
                case 'structsLabels'
                    structsLabels=varargin{2};
                case 'FigID'
                    FigID=varargin{2};
                case 'AxID'
                    AxID=varargin{2};
                case 'AxColor'
                    AxColor=varargin{2};
                case 'reglinePlot'
                    reglinePlot=varargin{2};
                case 'pointsPlot'
                    pointsPlot=varargin{2};
                case 'groups'
                   groups = varargin{2};
                case 'xlabelcad'
                    xlabelcad=varargin{2};
                case 'ylabelcad'
                    ylabelcad=varargin{2};
                case 'titlecad'
                    titlecad=varargin{2};
                case 'dispLeg'
                    dispLeg=varargin{2};
                case 'methLabels'
                    methLabels=varargin{2};
                case 'refColumn'
                    refColumn = varargin{2};
                case 'pointMarkers'
                    pointMarkers = varargin{2};
                case 'lineStyles'
                    lineStyles   = varargin{2};
                case 'pointColors'
                    pointColors  = varargin{2};
                case 'drawEllipse'
                    ellip_confPerc  = varargin{2};
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
[Xcorrpair,Ycorrpair] = find(triu(ones(dataSize(2),dataSize(2)),1));
Nregress = length(Xcorrpair);

nc = 1;
nr = 1;

switch figType
    case {'subplot','multiplot'}
        % Detect figure and axis
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
        
        switch figType
            case 'subplot'
                % Detecting number of regressions
                nc = ceil(sqrt(Nregress));
                nr = ceil(Nregress/nc);
                
                if ~exist('pointMarkers','var')
                    pointMarkers = 'o';
                end
                if ~exist('pointColors','var')
                      pointColors = [ 0,153,153]/255;
%                     pointColors = [ 21    50    140]/255;
                else
                    if ~isempty(find(pointColors > 1));
                        pointColors = pointColors/255;
                    end
                end
                
                for i = 1:Nregress
                    axIds(i) = subplot(nr,nc,i);
                    
                    % Data
                    X = dataY(:,Xcorrpair(i));
                    Y = dataY(:,Ycorrpair(i));
                    
                    % Labels
                    xLabelString = methLabels{Xcorrpair(i)};
                    yLabelString = methLabels{Ycorrpair(i)};
                    
                    % Removing outliers
                    plotOp = 0;
                    noutliers = 2;
%                     [X,Y,rSquares,outliers_idx] = regoutliers(X,Y,noutliers,plotOp);
%                     X(outliers_idx) = [];
%                     Y(outliers_idx) = [];
                    
                    % Regression
                    [B,BINT,R,RINT,STATS] = regress(Y,[ones(length(X),1) X]);
                    
                    % Correlation
                    [r,p] = corr(X,Y,'type','Spearman');
                    
                    if i == 1
                        oldUnits = axIds(i).Units;
                        axIds(i).Units = 'centimeters';
                        axPos = axIds(i).OuterPosition;
                        axIds(i).Units = oldUnits;
                        
                        refVal = min(axPos(3:4)); % Reference Value for the scaling process
                        fsaxis = refVal*7/7;   % Axis Ticks fontsize
                        fslabels = refVal*8/7; % Axis Labels fontsize
                        fstitle = refVal*10/7;  % Axis Title fontsize
                        
                        %                     fsaxis = 12;   % Axis Ticks fontsize
                        %                     fslabels = 14; % Axis Labels fontsize
                        %                     fstitle = 16;  % Axis Title fontsize
                        
                        
                        bLine = refVal*1.5/10;
                        markSize = 2*refVal*3/15;
                        medSize = refVal*4/15;
                    end
                    if pointsPlot
                        h = plot(X,Y,'Markersize',markSize,'MarkerEdgeColor','none','LineStyle','none','Marker',pointMarkers,'MarkerFaceColor',pointColors);
                    end
                    
                    Xrange =  max(X) - min(X);
                    Yrange =  max(Y) - min(Y);
                    
                    Xline = [min(X)-0.03*Xrange  max(X)+0.03*Xrange];
                    Yline = B(1)+B(2)*Xline;
                    hold on;
                    
                    if reglinePlot
                        h1 = plot(Xline,Yline,'Color',[0 0 0]/255,'LineWidth',bLine,'LineStyle','-');
                    end
                    
                    %                     titlecad = ['r = ' sprintf('%1.2f',r) ': B = ' sprintf('%2.2f',B(2)) ': B0 = ' sprintf('%2.2f',B(1)) ': R^2 = ' sprintf('%2.2f',STATS(1)) ];
                    
                    axIds(i).XLim = [min(X)-0.05*Xrange max(X)+0.05*Xrange]; % X Limits
                    axIds(i).YLim = [min(Y)-0.05*Yrange max(Y)+0.05*Yrange]; % Y Limits
                    axIds(i).FontSize = fsaxis;
                    
                    axIds(i).Color = AxColor;
                    axIds(i).XColor = [1 1 1]-param.figcolor;
                    axIds(i).YColor = [1 1 1]-param.figcolor;
                    axIds(i).ZColor = [1 1 1]-param.figcolor;
                    
                    axIds(i).Title.String = ['r = ' sprintf('%1.2f',r)  '(pval = ' sprintf('%1.2f',p) ') : B = ' sprintf('%2.2f',B(2)) ': B0 = ' sprintf('%2.2f',B(1)) ': R^2 = ' sprintf('%2.2f',STATS(1)) ];
                    axIds(i).Title.Visible = 'on';
                    axIds(i).Title.Color = [1 1 1] -param.figcolor;
                    axIds(i).Title.FontSize = fslabels*0.8;
                    axIds(i).Title.FontWeight = 'normal';
                    
                    axIds(i).XLabel.String = xLabelString;
                    axIds(i).XLabel.Visible = 'on';
                    axIds(i).XLabel.Color = [1 1 1] -param.figcolor;
                    axIds(i).XLabel.FontSize = fslabels;
                    
                    axIds(i).YLabel.String = yLabelString;
                    axIds(i).YLabel.Visible = 'on';
                    axIds(i).YLabel.Color = [1 1 1] -param.figcolor;
                    axIds(i).YLabel.FontSize = fslabels;
                    grid on;
                    
                end
                hTitle = suptitle(titlecad);
                hTitle.FontName = 'Arial';
                hTitle.Visible = 'on';
                hTitle.Color = [1 1 1] -param.figcolor;
                hTitle.FontSize = fstitle*1.1;
            case 'multiplot'
                if ~isempty(AxID)
                    set(FigID,'CurrentAxes',AxID);
%                     cla;
                else
                    clf;
                end
                
                if exist('groups','var')
                    groupsIds = unique(int8(groups));
                    Nregress = length(groupsIds);
                else
                    Nregress = dataSize(2)-1;
                    groups = ones(dataSize(1),1);
                    groupsIds = 1;
                end
                if ~exist('pointMarkers','var')
                    pointMarkers = strvcat('d','o','s','^','d','p','>','v');
                end
                if ~exist('lineStyles','var')
                    lineStyles = strvcat('-',':','--','-.');
                end
                
                if ~exist('pointColors','var')
                    C = [...
                        240,163,255;0,117,220;153,63,0;76,0,92;25,25,25;0,92,49;43,206,72;...
                        255,204,153;128,128,128;148,255,181;143,124,0;157,204,0;194,0,136;...
                        0,51,128;255,164,5;255,168,187;66,102,0;255,0,16;94,241,242;0,153,143;...
                        224,255,102;116,10,255;153,0,0;255,255,128;255,255,0;255,80,5;250 30 30]/255;
                    pointColors = C([27 2 15 20  8 16 1 3 4 5 6 7 9 10 11 12 13 14 17 18 19 21 22 23 24 25 26],:);
                    
                else
                    if ~isempty(find(pointColors > 1));
                        pointColors = pointColors/255;
                    end
                end
                pointMarkers = repmat(pointMarkers,[Nregress 1]);
                lineStyles   = repmat(lineStyles,[Nregress 1]);
                pointColors  = repmat(pointColors,[Nregress 1]);
                
                
                Ncolor = size(pointColors,1);
               
                
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
                fstitle = refVal*10/10;  % Axis Title fontsize

                
                bLine = refVal*1.5/14;
                markSize = refVal*2.5/6;
                medSize = refVal*4/14;
                
                dataOrig = dataY;
                Xorig = dataY(:,refColumn);
                dataY(:,refColumn) = [];
                
                % Removing 0s
                ind2rem = find(Xorig == 0);
                Xorig(ind2rem) = [];
                dataY(ind2rem,:) = [];
                groups(ind2rem,:) = [];
                
                                
                % Removing NaNs
                ind2rem = find(isnan(Xorig));
                Xorig(ind2rem) = [];
                dataY(ind2rem,:) = [];
                groups(ind2rem,:) = [];
                
                
                methLabelsOrig = methLabels;
                if exist('xlabelcad','var')
                    xLabelString = xlabelcad;
                else
                    xLabelString = methLabels{refColumn};
                end
                methX = methLabels{refColumn};
                methLabels(refColumn) = [];
                
                allYs = 0;
                allXs = 0;
                cont = 0;
                legText = '';
                allNames = '';
                for i = 1:Nregress
                    % Data
                    ind2plot = find(groups == groupsIds(i));
                    if length(unique(groupsIds))>1
                        Y = dataY(ind2plot,1);
                        yLabelString = methLabels{1};
                    else
                        Y = dataY(:,i);
                        yLabelString = methLabels{i};
                        
                    end
                    Xorig(isnan(Xorig)) = [];Y(isnan(Xorig)) = [];
                    Xorig(isnan(Y)) = [];Y(isnan(Y)) = [];
                    methY = methLabels{i};
                    % Labels
                    if i == 1
                        allMethodsString = yLabelString;
                    else
                        allMethodsString = [allMethodsString ', ' yLabelString];
                    end
                    if exist('ylabelcad','var')
                        allMethodsString = ylabelcad;
                    end
                    % Removing outliers
                    plotOp = 0;
                    noutliers = 1;
                    [X,Y,rSquares,outliers_idx] = regoutliers(Xorig(ind2plot),Y,noutliers,plotOp);
                    X(outliers_idx) = [];
                    Y(outliers_idx) = [];
                    
                    % Regression
                    [B,BINT,R,RINT,RSTATS] = regress(Y,[ones(length(X),1) X]);
                    
                    % Correlation
                    [r,rpval] = corr(X,Y,'type','Spearman');
                    
                  
                    
                    
                    % Regression Line
                    Xrange =  max(X) - min(X);
                    Xline = [min(X)-0.03*Xrange  max(X)+0.03*Xrange];
                    Yline = B(1)+B(2)*Xline;
                    hold on;
                    if reglinePlot
                        h = plot(Xline,Yline,'Color',pointColors(i,:),'LineWidth',bLine*1.5,'LineStyle',lineStyles(i,:));
                        
                        %                   cont = cont + 1;
                        %                   handPlots(cont)  = h;
                        %                   legText = strvcat(legText,['Line: ' xLabelString ' vs ' yLabelString ]);
                    end
                    
                    
                    hold on;
                    cont = cont + 1;
                    handPlots(cont) = plot(X,Y,'Marker',pointMarkers(i),'MarkerFaceColor',pointColors(i,:),'MarkerEdgeColor',[0 0 0],'LineStyle','none', 'Markersize',markSize);
                    legText = strvcat(legText,['Plot: ' xLabelString ' vs ' yLabelString ]);

                    
                    %                     titlecad = ['r = ' sprintf('%1.2f',r) ': B = ' sprintf('%2.2f',B(2)) ': B0 = ' sprintf('%2.2f',B(1)) ': R^2 = ' sprintf('%2.2f',STATS(1)) ];
                    
                    
                    if ellip_confPerc
                        [meanX, meanY, stdX, stdY, a0, Xe, Ye] =  EEP(X,Y,ellip_confPerc,100);
                        hold on;
                        ph = patch(Xe,Ye,Ye*0,'FaceColor',pointColors(i,:),'EdgeColor',[[1 1 1]-param.figcolor]*.70,'FaceAlpha',0.1);
                        allYs = [allYs;Y;Ye(:)];
                        allXs = [allXs;X;Xe(:)];
                    else
                        allYs = [allYs;Y];
                        allXs = [allXs;X];
                    end
                    
                    
                    
                    % Regression
                    [H,tpval,CI,TSTATS] = ttest(X,Y,'tail','both');
                    
                    % Cohen's D
                    cohens_D = cohensD(X,Y,'related');
                    
                    
                    corrrValues(i,1)    = r;
                    corrpValues(i,1)    = rpval;
                    slopesValues(i,1)    = B(2);
                    interceptValues(i,1) = B(1);
                    RsquareValues(i,1)   = RSTATS(1);
                    tValues(i,1)         = TSTATS.tstat;
                    ttestpvalValues(i,1) = tpval;
                    cohensDValues(i,1)   = cohens_D;
                    
                   compName = [methX ' vs ' methY];
                        
                   allNames = strvcat(allNames,compName);
                    
                end
                
                
                divColumns = repmat(' , ',[size(allNames,1)+1 1]);
                parRight = repmat(' (',[size(allNames,1) 1]);
                parLeft  = repmat(' )',[size(allNames,1) 1]);
                allStats = [strvcat('Comparisons',allNames) divColumns ...
                            strvcat('Pearson Corr (pvalue)',[num2str(corrrValues) parRight num2str(corrpValues) parLeft ]) divColumns ...
                            strvcat('Regress Slopes',num2str(slopesValues)) divColumns ...
                            strvcat('Regress Intercept',num2str(interceptValues)) divColumns ...
                            strvcat('Regress R-square',num2str(RsquareValues)) divColumns ...
                            strvcat('T-test Tvalue (pvalue, cohensD)',[num2str(tValues) parRight num2str(ttestpvalValues) divColumns(1:size(cohensDValues,1),:) num2str(cohensDValues) parLeft ])];
                
                
                
                allYs(1) = [];
                allXs(1) = [];
                
                Xrange =  max(allXs) - min(allXs);
                Yrange =  max(allYs) - min(allYs);
                
                AxID.XLim = [min(allXs)-0.05*Xrange max(allXs)+0.05*Xrange]; % X Limits
                AxID.YLim = [min(allYs)-0.05*Yrange max(allYs)+0.05*Yrange]; % Y Limits
                AxID.FontSize = fsaxis;
                if reglinePlot
%                     h = plot([ min(X)-0.03*Xrange; max(X)+0.03*Xrange],[min(allYs)-0.03*Yrange max(allYs)+0.03*Yrange],'Color',[1 1 1]-param.figcolor,'LineWidth',bLine*1.5,'LineStyle','-');
                end
                
                AxID.Color = AxColor;
                AxID.XColor = [1 1 1]-param.figcolor;
                AxID.YColor = [1 1 1]-param.figcolor;
                AxID.ZColor = [1 1 1]-param.figcolor;
                
                AxID.Title.String = titlecad;
                AxID.Title.Visible = 'on';
                AxID.Title.Color = [1 1 1] -param.figcolor;
                AxID.Title.FontSize = fslabels*1.1;
                AxID.Title.FontWeight = 'normal';
                
                AxID.XLabel.String = xLabelString;
                AxID.XLabel.Visible = 'on';
                AxID.XLabel.Color = [1 1 1] -param.figcolor;
                AxID.XLabel.FontSize = fslabels;
                
                AxID.YLabel.String = allMethodsString;
                AxID.YLabel.Visible = 'on';
                AxID.YLabel.Color = [1 1 1] -param.figcolor;
                AxID.YLabel.FontSize = fslabels;
                
                grid on;

% % % % % %                 % Computing Stats
% % % % % %                 allNames = '';
% % % % % %                 Np = size(dataOrig,2);
% % % % % %                 cont = 0;
% % % % % %                 for j = 1:Np-1
% % % % % %                     methX = methLabelsOrig{j};
% % % % % %                     for z = j+1:Np
% % % % % %                         methY = methLabelsOrig{z};
% % % % % %                         cont = cont + 1;
% % % % % %                         
% % % % % %                         
% % % % % %                         X = dataOrig(:,j);
% % % % % %                         Y = dataOrig(:,z);
% % % % % %                         ind2rem = logical(isnan(Y) + isnan(X));
% % % % % %                         X(ind2rem) = [];
% % % % % %                         Y(ind2rem) = [];
% % % % % %                         
% % % % % %                         
% % % % % %                         
% % % % % %                         % Removing outliers
% % % % % %                         plotOp = 0;
% % % % % %                         noutliers = 1;
% % % % % %                         [X,Y,rSquares,outliers_idx] = regoutliers(X,Y,noutliers,plotOp);
% % % % % %                         X(outliers_idx) = [];
% % % % % %                         Y(outliers_idx) = [];
% % % % % %                         
% % % % % %                         % Regression
% % % % % %                         [B,BINT,R,RINT,RSTATS] = regress(Y,[ones(length(X),1) X]);
% % % % % %                         
% % % % % %                         % Correlation
% % % % % %                         [r,rpval] = corr(X,Y);
% % % % % %                         
% % % % % %                         % Regression
% % % % % %                         [H,tpval,CI,TSTATS] = ttest(X,Y,'tail','both');
% % % % % %                         
% % % % % %                         % Cohen's D
% % % % % %                         cohens_D = cohensD(X,Y,'related');
% % % % % %                         
% % % % % %                         
% % % % % %                         corrrValues(cont,1)    = r;
% % % % % %                         corrpValues(cont,1)    = rpval;
% % % % % %                         slopesValues(cont,1)    = B(2);
% % % % % %                         interceptValues(cont,1) = B(1);
% % % % % %                         RsquareValues(cont,1)   = RSTATS(1);
% % % % % %                         tValues(cont,1)         = TSTATS.tstat;
% % % % % %                         ttestpvalValues(cont,1) = tpval;
% % % % % %                         cohensDValues(cont,1)   = cohens_D;
% % % % % %                         
% % % % % %                         
% % % % % %                         compName = [methX ' vs ' methY];
% % % % % %                         
% % % % % %                         allNames = strvcat(allNames,compName);
% % % % % % %                         ['r = ' sprintf('%1.2f',r) ...
% % % % % % %                             ', Slope = ' sprintf('%2.2f',B(2)) ...
% % % % % % %                             ', Intercept = ' sprintf('%2.2f',B(1)) ...
% % % % % % %                             ', R^2 = ' sprintf('%2.2f',RSTATS(1)) ...
% % % % % % %                             ', t-value = ' sprintf('%2.2f',TSTATS.tstat) ...
% % % % % % %                             ', p-value = ' sprintf('%2.2f',pval)];
% % % % % %                     end
% % % % % %                 end
end
% % % % % %                         divColumns = repmat(' , ',[size(allNames,1)+1 1]);
% % % % % %                 parRight = repmat(' (',[size(allNames,1) 1]);
% % % % % %                 parLeft  = repmat(' )',[size(allNames,1) 1]);
% % % % % %                 allStats = [strvcat('Comparisons',allNames) divColumns ...
% % % % % %                             strvcat('Pearson Corr (pvalue)',[num2str(corrrValues) parRight num2str(corrpValues) parLeft ]) divColumns ...
% % % % % %                             strvcat('Regress Slopes',num2str(slopesValues)) divColumns ...
% % % % % %                             strvcat('Regress Intercept',num2str(interceptValues)) divColumns ...
% % % % % %                             strvcat('Regress R-square',num2str(RsquareValues)) divColumns ...
% % % % % %                             strvcat('T-test Tvalue (pvalue, cohensD)',[num2str(tValues) parRight num2str(ttestpvalValues) divColumns(1:size(cohensDValues,1),:) num2str(cohensDValues) parLeft ])];
                        
        % Displaying Legend
        if dispLeg
            hleg = legend( handPlots, cellstr(legText),'Box', 'on','Color', 'none','FontSize',fslabels*0.7, 'FontName','Arial');
        end
        hold off;
        
    case 'indivplot'
        if ~exist('pointMarkers','var')
            pointMarkers = 'o';
        end
        if ~exist('pointColors','var')
            pointColors = [ 230,0,0]/255;
%             pointColors = [ 21    50    140]/255;
        else
            if ~isempty(find(pointColors > 1));
                pointColors = pointColors/255;
            end
        end
        for i = 1:Nregress
            
            if ~exist('FigID','var')
                 FigID(i) = figure('numbertitle','off','name','New Figure','Color',param.figcolor, 'units',figUnits,'Position',figPosition,'Visible',param.figurevisible,'InvertHardcopy','off');
            else
                if ~isempty(FigID)
%                     figure(FigID);
                else
                 FigID(i) = figure('numbertitle','off','name','New Figure','Color',param.figcolor, 'units',figUnits,'Position',figPosition,'Visible',param.figurevisible,'InvertHardcopy','off');
                end
            end
            axIds = gca;
            
            % Data
            X = dataY(:,Xcorrpair(i));
            Y = dataY(:,Ycorrpair(i));
            indnanX = find(isnan(X));
            indnanY = find(isnan(Y));

            X(unique([indnanX;indnanY])) = [];
            Y(unique([indnanX;indnanY])) = [];

            % Labels
            xLabelString = methLabels{Xcorrpair(i)};
            yLabelString = methLabels{Ycorrpair(i)};
            
            % Removing outliers
            plotOp = 0;
            noutliers = 2;
            [X,Y,rSquares,outliers_idx] = regoutliers(X,Y,noutliers,plotOp);
            X(outliers_idx) = [];
            Y(outliers_idx) = [];
            
            % Regression
            [B,BINT,R,RINT,STATS] = regress(Y,[ones(length(X),1) X]);
            
            % Correlation
            [r,p] = corr(X,Y,'type','Spearman');
            
            set(0,'units','centimeters');
            cm_screen = get(0,'screensize');
            
            oldUnits = axIds.Units;
            axIds.Units = 'centimeters';
            axPos = axIds.Position;
            axIds.Units = oldUnits;
            
            oldUnits = FigID(i).Units;
            FigID(i).Units = 'centimeters';
            figPos = FigID(i).Position;
            FigID(i).Units = oldUnits;
                
% %             refVal = min(axPos(3:4)); % Reference Value for the scaling process
% %             fsaxis = refVal*7/5;   % Axis Ticks fontsize
% %             fslabels = refVal*8/5; % Axis Labels fontsize
% %             fstitle = refVal*10/5;  % Axis Title fontsize
% %            
% %             
% %             bLine = refVal*1.5/10;
% %             markSize = 2*refVal*3/10;
% %             medSize = refVal*4/10;
% %             

                

            
                refVal = 55*prod(axPos(3:4))/prod(figPos(3:4))*2;
%                 refVal = min(axPos(3:4)); % Reference Value for the scaling process
                fsaxis = refVal*7/10;   % Axis Ticks fontsize
                fslabels = refVal*8/10; % Axis Labels fontsize
                fstitle = refVal*10/10;  % Axis Title fontsize

                
                bLine = refVal*1.5/14;
                markSize = refVal*2.5/6;
                medSize = refVal*4/14;
            
            

           hLine = plot(X,Y,'Markersize',markSize,'MarkerEdgeColor',[0.2 0.2 0.2],'LineStyle','none','Marker', deblank(pointMarkers), 'MarkerFaceColor', pointColors,'Linewidth',1.8); 

           drawnow;
           hMarkers = hLine.MarkerHandle;
           hMarkers.FaceColorData = uint8(255*[pointColors';transpVal]); 
           
           
            
            %             h = scatter(X,Y,markSize*5,pointColors,'Marker', deblank(pointMarkers),'MarkerEdgeColor',[0.2 0.2 0.2]); 

%             cFace = uint8(255*[pointColors transpVal])';
%             cEdge = uint8(255*[0.2 0.2 0.2 .8])';
% %              set(h.MarkerHandle,'FaceColorBinding','interpolated', 'FaceColorData',cFace)
% %             plotMarkers = h.MarkerHandle;
%              h.MarkerHandle.FaceColorBinding = 'interpolated';
%              h.MarkerHandle.FaceColorData =cFace;
%              h.MarkerHandle.EdgeColorData =cEdge;
%             plotMarkers.FaceColorData = cFace;
%             plotMarkers.EdgeColorData = cEdge;
            
            % keep transparent
%             drawnow();
            
            
            
% %             plotMarkers = h.MarkerHandle;
% %             plotMarkers.FaceColorType = 'truecoloralpha';
% %             tempVar = floor(255*[pointColors transpVal]');
% %             plotMarkers.FaceColorData = uint8(tempVar);
% %             
% %             plotMarkers.EdgeColorData = uint8(floor(255*[0.2 0.2 0.2 .8]))';

            
            
            Xrange =  max(X) - min(X);
            Yrange =  max(Y) - min(Y);
            
            Xline = [min(X)-0.03*Xrange  max(X)+0.03*Xrange];
            Yline = B(1)+B(2)*Xline;
            hold on;
            
            if reglinePlot
                h1 = plot(Xline,Yline,'Color',[0 0 0]/255,'LineWidth',bLine,'LineStyle','-');
            end
            %                     titlecad = ['r = ' sprintf('%1.2f',r) ': B = ' sprintf('%2.2f',B(2)) ': B0 = ' sprintf('%2.2f',B(1)) ': R^2 = ' sprintf('%2.2f',STATS(1)) ];
            
            axIds.XLim = [min(X)-0.05*Xrange max(X)+0.05*Xrange]; % X Limits
            axIds.YLim = [min(Y)-0.05*Yrange max(Y)+0.05*Yrange]; % Y Limits
            axIds.FontSize = fsaxis;
            
            axIds.Color = AxColor;
            axIds.XColor = [1 1 1]-param.figcolor;
            axIds.YColor = [1 1 1]-param.figcolor;
            axIds.ZColor = [1 1 1]-param.figcolor;
            
%             axIds.Title.String = ['r = ' sprintf('%1.2f',r) '(pval = ' sprintf('%1.2f',p) ') : B = ' sprintf('%2.2f',B(2)) ': B0 = ' sprintf('%2.2f',B(1)) ': R^2 = ' sprintf('%2.2f',STATS(1)) ];
            axIds.Title.String = ['r = ' sprintf('%1.2f',r) '(pval = ' sprintf('%1.2f',p) ') :  R^2 = ' sprintf('%2.2f',STATS(1)) ];

            axIds.Title.Visible = 'on';
            axIds.Title.Color = [1 1 1] -param.figcolor;
            axIds.Title.FontSize = fslabels*0.8;
            axIds.Title.FontWeight = 'normal';
            
            axIds.XLabel.String = xLabelString;
            axIds.XLabel.Visible = 'on';
            axIds.XLabel.Color = [1 1 1] -param.figcolor;
            axIds.XLabel.FontSize = fslabels;
            
            axIds.YLabel.String = yLabelString;
            axIds.YLabel.Visible = 'on';
            axIds.YLabel.Color = [1 1 1] -param.figcolor;
            axIds.YLabel.FontSize = fslabels;
            
%             hTitle = suptitle(titlecad);
%             hTitle.FontName = 'Arial';
%             hTitle.Visible = 'on';
%             hTitle.Color = [1 1 1] -param.figcolor;
%             hTitle.FontSize = fstitle*1.1;
            grid on;
%             allStats = ['r = ' sprintf('%1.2f',r) ', pval = ' sprintf('%1.2f',p) ', B = ' sprintf('%2.2f',B(2)) ', B0 = ' sprintf('%2.2f',B(1)) ', R^2 = ' sprintf('%2.2f',STATS(1)) ];
            allStats = [sprintf('%1.2f',r) ', ' sprintf('%1.2f',p) ', ' sprintf('%2.2f',B(2)) ', ' sprintf('%2.2f',B(1)) ', ' sprintf('%2.2f',STATS(1)) ];

        end
end

% % % % % % %% ================= Assigning the same properties for all graphs ====================== %
% % % % % % axIds = FigID.Children;
% % % % % % xNumticks = length(axIds(1).XTick);
% % % % % % yNumticks = length(axIds(1).YTick);
% % % % % % for i = 2:length(axIds)
% % % % % %     axIds(i).Title.FontName     = axIds(1).Title.FontName;
% % % % % %     axIds(i).Title.FontSize     = axIds(1).Title.FontSize;
% % % % % %     axIds(i).Title.FontWeight   = axIds(1).Title.FontWeight;
% % % % % %     axIds(i).Title.Color        = axIds(1).Title.Color;
% % % % % %     axIds(i).Title.FontWeight   = axIds(1).Title.FontWeight;
% % % % % %     
% % % % % %     %% XLabel Properties
% % % % % %     % Font Properties XLabel
% % % % % %     axIds(i).XLabel.FontName    = axIds(1).XLabel.FontName;
% % % % % %     axIds(i).XLabel.FontSize    = axIds(1).XLabel.FontSize; 
% % % % % %     axIds(i).XLabel.FontWeight  = axIds(1).XLabel.FontWeight;
% % % % % %     axIds(i).XLabel.Color       = axIds(1).XLabel.Color;
% % % % % %     axIds(i).XLabel.FontWeight  = axIds(1).XLabel.FontWeight;
% % % % % %     
% % % % % %     %% YLabel Properties
% % % % % %     % Font Properties XLabel
% % % % % %     axIds(i).YLabel.FontName    = axIds(1).YLabel.FontName;
% % % % % %     axIds(i).YLabel.FontSize    = axIds(1).YLabel.FontSize; 
% % % % % %     axIds(i).YLabel.FontWeight  = axIds(1).YLabel.FontWeight;
% % % % % %     axIds(i).YLabel.Color       = axIds(1).YLabel.Color;
% % % % % %     axIds(i).YLabel.FontWeight  = axIds(1).YLabel.FontWeight;
% % % % % %      
% % % % % %     
% % % % % %     
% % % % % %     xLimits = axIds(i).XLim;
% % % % % %     values = linspace(xLimits(1),xLimits(2),xNumticks);
% % % % % %     axIds(i).XTickLabel = cellstr(sprintf('%.2f\n',values(:)));
% % % % % %     axIds(i).XTick = values(:);
% % % % % %     
% % % % % %     yLimits = axIds(i).YLim;
% % % % % %     values = linspace(yLimits(1),yLimits(2),yNumticks);
% % % % % %     axIds(i).YTickLabel = cellstr(sprintf('%.2f\n',values(:)));
% % % % % %     axIds(i).YTick = values(:);
% % % % % %     
% % % % % % end


%% ====================== End of Main Program =========================== %

% Outputs
varargout{1} = FigID;
varargout{2} = allStats;
return;
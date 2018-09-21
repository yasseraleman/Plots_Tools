function varargout = my_jitter_plus_box_groups_plot(dataY,varargin)
%
% Syntax :
%  [FigID,outliersMat] = my_jitter_plus_box_groups_plot(dataY,varargin)
%
% This script creates a group boxplot over a jitter distribution for a given data.
% 3D matrix is asumed as Nstruc x Nmethods x Nsubjects
% 2D matrix is asumed as Nsubjects x Nstruc or Nmethods
%
%
% Input Parameters:
%       dataY                 : Data matrix
%
%       charVal               : Map Filename
%
% Output Parameters:
%      FigID                  : Figure handle
%      outliersMat            : Boolean matrix with the same size of dataY
%                               and 1s representing data outliers.
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
    dimData = length(dataSize);
    if dimData == 2
        if dataSize(1) == 1
            dataY = dataY(:);
        end
        dataYt(:,1,:) = dataY';
        dataY = dataYt; clear dataYt;
        bidData = 1; % Bidimensional Data
        methLabels = cellstr([repmat('Method ',[size(dataY,1) 1]) num2str([1:size(dataY,1)]','%.3d')]); % Labels for Methods
        structsLabels = methLabels; % Labels for Groups
        dispLeg = 0; % Display Legend
        
    elseif dimData == 3
        % Default Labels
        structsLabels = cellstr([repmat('Group ',[size(dataY,1) 1]) num2str([1:size(dataY,1)]','%.3d')]); % Labels for Groups
        methLabels = cellstr([repmat('Method ',[size(dataY,2) 1]) num2str([1:size(dataY,2)]','%.3d')]); % Labels for Methods    
        dispLeg = 1; % Display Legend
    
    elseif dimData > 3 
       error('This script only works for 1D, 2D or 3D data');
       return;
    end
    
    % Parameters
    grpSepIndex = 4; % Separation between groups
    jittPercent = 50; % Jitter Percentage
    transpVal = .5; % Opacity
    
    % Boolean variables
    boxPlot = 1; % Plot Boxplots
    pointsPlot = 1;  % Plot Data Points
    jitFill = 1;     % Fill jitter markers
    
    % Figure Properties
    param.figcolor = 'white'; % Color
    param.figurevisible = 'on'; % Visible
    figPosition = [1 1 30 30]; % Position
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
                case 'grpSepIndex' % Group Separation
                    grpSepIndex=varargin{2};
                case 'jittPercent' % Jitter compression percent
                    jittPercent=varargin{2};
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
                case 'figUnits'
                    figUnits=varargin{2};
                case 'structsLabels'
                    structsLabels=varargin{2};
                case 'FigID'
                    FigID=varargin{2};
                case 'AxID'
                    AxID=varargin{2};
                case 'AxColor'
                    AxColor=varargin{2};
                case 'boxPlot'
                    boxPlot=varargin{2};
                case 'pointsPlot'
                    pointsPlot=varargin{2};
                case 'jitFill'
                    jitFill=varargin{2};
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
                otherwise
                    error('Unexpected ''ParameterName'' input: %s\n',varargin{1});
            end
            varargin(1:2)=[]; % this pair of optional input arguments has been dealt with -- remove...
        end
    end
end
% Detect figure and axis
if ~exist('FigID','var')
    switch param.figcolor
        case 'black'
            colordef black; FigID = figure('numbertitle','off','name','New Figure','Color',[0 0 0], 'units',figUnits,'Position',figPosition,'Visible',param.figurevisible);
        case 'white'
            colordef white; FigID = figure('numbertitle','off','name','New Figure','Color',[1 1 1], 'units',figUnits,'Position',figPosition,'Visible',param.figurevisible,'InvertHardcopy','off');
    end
else
    if ~isempty(FigID)
%         figure(FigID);
    else
        figure;
        FigID=gcf;
    end
end
if exist('AxID','var')
    if ~isempty(AxID)
        set(FigID,'CurrentAxes',AxID);
        cla;
    else
        clf;
    end
end

if exist('AxColor','var')
     if ~isempty(AxColor)&length(AxColor) == 3
         ha = gca;
         if ~isempty(find(AxColor> 1))
             AxColor = AxColor/255;
         end
         ha.Color = AxColor;
     end
end

hold on;

%% =================== End of checking input parameters ================= %
AxID = gca;
oldUnits = get(AxID,'units');
set(AxID,'units','centimeters');
axPos = get(AxID,'outerposition');
set(AxID,'units',oldUnits);
refVal = min(axPos(3:4)); % Reference Value for the scaling process
fsaxis = refVal*7/10;   % Axis Ticks fontsize
fslabels = refVal*8/10; % Axis Labels fontsize
fstitle = refVal*10/10;  % Axis Title fontsize


fsaxis = 12;   % Axis Ticks fontsize
fslabels = 14; % Axis Labels fontsize
fstitle = 16;  % Axis Title fontsize


bLine = refVal*1.5/20;
markSize = 2*refVal*3/20;
medSize = refVal*4/20;
%% ========================== Main Program ============================== %

% Defining Colors and Markers
C = [...
    240,163,255;0,117,220;153,63,0;76,0,92;25,25,25;0,92,49;43,206,72;...
    255,204,153;128,128,128;148,255,181;143,124,0;157,204,0;194,0,136;...
    0,51,128;255,164,5;255,168,187;66,102,0;255,0,16;94,241,242;0,153,143;...
    224,255,102;116,10,255;153,0,0;255,255,128;255,255,0;255,80,5;250 0 0]/255;
Colors = C([2 27 15 20  8 16 1 3 4 5 6 7 9 10 11 12 13 14 17 18 19 21 22 23 24 25 26],:);
Ncolor = size(Colors,1);

% Colors = [240 0 0; 90 90 220 ; 44 170 156]/255;


markers = strvcat('d','o','s','^','d','p','>','v');
Nmarkers = size(markers,1);


% Detecting Groups, Methods and Subjects
[N_structs, N_methods, N_subj]=size(dataY);
re = floor(max([N_methods N_structs])/Ncolor); Colors = repmat(Colors,[re+1 1]);

re = floor(max([N_methods N_structs])/Nmarkers); markers = repmat(markers,[re+1 1]);

% Building X Coordinates
Xcoords = [1:N_structs*N_methods];
Xcoords = reshape(Xcoords,[N_methods N_structs ])';

% Separations between Groups
groupSep = grpSepIndex*repmat([0:size(Xcoords,1)-1]',[1 size(Xcoords,2)]);
Xcoords = Xcoords + groupSep;

% Jitter Distribution between Groups
jittDist = jittPercent*(rand(N_structs, N_methods, N_subj)-0.5)/100;

% Adding Jitter to X coordinates
Xcoords = repmat(Xcoords,[1 1 size(jittDist,3)]) + jittDist;

outliersMat = zeros(size(dataY));
if bidData
%     numrep = ceil(size(methLabels,1)/size(markers,1));
%     markers = repmat(markers,[numrep 1]);
    structsLabels = methLabels;
    for grp=1:N_structs
        Xpos = squeeze(Xcoords(grp,1,:))';
        Ypos = squeeze(dataY(grp,1,:))';
        hold on;
        if pointsPlot
            hplot = plot(Xpos(:),Ypos(:),'Marker',markers(grp),'MarkerFaceColor',Colors(grp,:),'MarkerEdgeColor',[0 0 0],'LineStyle','none', 'Markersize',markSize*1.2);
            
            % Legend handles
            handPlots(grp) = hplot;
        end
        if boxPlot
            % Plotting BoxPlot
            hbox = boxplot(gca,Ypos(isfinite(Ypos)),'positions',mean(Xpos), 'widths',jittPercent/100);
            set(gca,'XTickLabel',{' '});
            
            % Finding Outliers
            h = flipud(hbox(7)); % flip order of handles
            e = eps(max(Ypos(:)));
            xout =  get( h, 'XData' );
            yout =  get( h, 'YData' );
            if  ~isnan(yout)&~isempty(yout)
                outliersMat(grp,1,:) = ismember(Ypos,yout);
            end
            % Manipulating Boxplots atributes
            set(hbox(1),'Linewidth',bLine,'Color',[0 0 0],'LineStyle','-.'); % Changing Dot Lines
            set(hbox(2),'Linewidth',bLine,'Color',[0 0 0],'LineStyle','-.'); % Changing Dot Lines
            set(hbox(3),'Linewidth',bLine,'Color',[0 0 0],'LineStyle','-'); % Changing Top Lines
            set(hbox(4),'Linewidth',bLine,'Color',[0 0 0],'LineStyle','-'); % Changing Top Lines
            %if boxPlot&&~pointsPlot
            set(hbox(5),'Linewidth',bLine,'Color',[0 0 0]); % Changing Rectangle Line width
            hpatch = patch(get(hbox(5),'XData'),get(hbox(5),'YData'),Colors(grp,:),'FaceAlpha',.3);
            %end
            set(hbox(6),'Linewidth',medSize,'Color',[1 0 0]); % Changing Median Line properties
            set(hbox(7),'Marker',markers(grp),'MarkerFaceColor',Colors(grp,:),'LineStyle','none', 'Markersize',2); % Changing Outliers properties
            if ~pointsPlot
                handPlots(grp) = hpatch;
            end
            
        end
    end
    outliersMat = squeeze(outliersMat)';
else
    for grp=1:N_methods
        Xpos = squeeze(Xcoords(:,grp,:))';
        Ypos = squeeze(dataY(:,grp,:))';
        hold on;
        if pointsPlot
            if jitFill
                hplot = plot(Xpos(:),Ypos(:),'Marker',markers(grp),'MarkerFaceColor',Colors(grp,:),'MarkerEdgeColor',[0 0 0],'LineStyle','none', 'Markersize',markSize*1.2);
            else
                hplot = plot(Xpos(:),Ypos(:),'Marker',markers(grp),'MarkerFaceColor','none','MarkerEdgeColor',Colors(grp,:),'LineStyle','none', 'Markersize',markSize);
                
            end
            
            % Legend handles
            handPlots(grp) = hplot;
        end
        if boxPlot
            % Plotting BoxPlot
            
% % % % %             hbox = boxplot(gca,Ypos,'positions',mean(Xpos), 'widths',jittPercent/100);
% % % % %             set(gca,'XTickLabel',{' '});
            
            % Manipulating Boxplots atributes
            for strs = 1:N_structs
                hbox = boxplot(gca,Ypos(isfinite(Ypos(:,strs)),strs),'positions',mean(Xpos(:,strs)), 'widths',jittPercent/100);
                
                
                % Finding Outliers
                h = flipud(hbox(7)); % flip order of handles
                xout =  get( h, 'XData' );
                yout =  get( h, 'YData' );
                if  ~isnan(yout)&~isempty(yout)
                    outliersMat(strs,grp,:) = ismember(Ypos(:,strs),yout);
                end
                
                
                
                set(gca,'XTickLabel',{' '});
                set(hbox(1),'Linewidth',bLine,'Color',[0 0 0],'LineStyle','-.'); % Changing Dot Lines
                set(hbox(2),'Linewidth',bLine,'Color',[0 0 0],'LineStyle','-.'); % Changing Dot Lines
                set(hbox(3),'Linewidth',bLine,'Color',[0 0 0],'LineStyle','-'); % Changing Top Lines
                set(hbox(4),'Linewidth',bLine,'Color',[0 0 0],'LineStyle','-'); % Changing Top Lines
                %if boxPlot&&~pointsPlot
                set(hbox(5),'Linewidth',bLine,'Color',[0 0 0]); % Changing Rectangle Line width
                hpatch = patch(get(hbox(5),'XData'),get(hbox(5),'YData'),Colors(grp,:),'FaceAlpha',.3);
                %end
                set(hbox(6),'Linewidth',medSize,'Color',[1 0 0]*0); % Changing Median Line properties
                set(hbox(7),'Marker',markers(grp),'MarkerFaceColor',Colors(grp,:),'LineStyle','none', 'Markersize',2); % Changing Outliers properties
            end
            % Box Plot Legend
            if ~pointsPlot
                handPlots(grp) = hpatch;
            end
        end
        clear hplot;
    end
end

% % % % % % Setting Opacity
% % % % % if pointsPlot
% % % % %     for i = 1:length(handPlots);
% % % % %         try
% % % % %             plotMarkers = handPlots(i).MarkerHandle;
% % % % % %             set(plotMarkers,'EdgeColorData',uint8(255*[0.7 0.7 0.7 .8])');
% % % % % %             set(plotMarkers,'FaceColorData',uint8(255*[Colors(i,:) transpVal])');
% % % % %             
% % % % %             plotMarkers.EdgeColorData = uint8(255*round([0.7 0.7 0.7 .8]))';
% % % % %             plotMarkers.FaceColorData = uint8(255*round([Colors(i,:) transpVal]))';
% % % % %             
% % % % %         end
% % % % %     end
% % % % % end


set(gca,'XTickLabel',{' '});

% Tunning the figure
set(gca,'XLim',[min(Xcoords(:))-jittPercent/100 max(Xcoords(:))+jittPercent/100]); % adjusts the x axis to the plot
usData = dataY(isfinite(dataY)); % Removing Infinite values
gapval = std(usData(:))/4;
if min(usData(:))>=0&~isempty(find(get(gca,'YLim') <0))
    set(gca,'YLim',[0  max(usData(:))+gapval ]); % adjusts the y axis to the plot
else
    set(gca,'YLim',[min(usData(:))-gapval  max(usData(:))+gapval ]); % adjusts the y axis to the plot
end
box on;grid on;

% set the x tick labels
set(gca,'XTick',mean(mean(Xcoords,3),2)','XTickLabel',structsLabels);

% Computing Automatic FontSize

% fsaxis
% fslabels
% fstitle =
set(gca,'FontSize',fsaxis,'FontName','Arial'); % adjusts the y axis to the plot

% Setting Axis Labels

hx = xlabel(xlabelcad,'FontSize',fslabels, 'FontName','Arial');
hx.Position(2) = -axPos(3);
hy =  ylabel(ylabelcad,'FontSize',fslabels, 'FontName','Arial');
%hy.Position(1) = -2;
title(titlecad,'FontSize',fstitle, 'FontName','Arial');

% Displaying Legend
if dispLeg
    hleg = legend( handPlots, methLabels,'Box', 'on','Color', 'none','FontSize',fslabels, 'FontName','Arial');
end
hold off;

%% ====================== End of Main Program =========================== %


hx = xlabel(xlabelcad,'FontSize',16, 'FontName','Arial');
hy =  ylabel(ylabelcad,'FontSize',16, 'FontName','Arial');
set(gca,'FontSize',16,'FontName','Arial'); % adjusts the y axis to the plot

% Outputs
varargout{1} = FigID;
varargout{2} = outliersMat;
return;
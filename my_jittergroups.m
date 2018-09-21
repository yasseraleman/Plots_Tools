function [jitter_xtick,hb,he]=my_jittergroups(Y,varargin)
%
% ERRORBAR_GROUPS: grouped bar plot with error bars, customizable.
%
% ERRORBAR_GROUPS(mean_vals,std_vals) produces a bar plot of matrix
% mean_vals, whose first dimension indicates the number of bars per group
% and the second dimension indicates the number of groups. The bars are
% grouped together, similar to the plot produced by BAR(Y,'grouped').
%
% It then overlays an error bar plot of std_vals, the size of which
% must match that of mean_vals. Error bars are assumed to be symmetric
% around the values of the bar plot, similar to the plot produced by
% ERRORBAR(Y,E). Note that, if all the values of std_vals are set to
% zero, the function skips plotting the error bars.
%
% ERRORBAR_GROUPS(mean_vals,errorbar_lower,errorbar_upper) allows
% the lower and upper bounds of the error bars to be asymmetric around the
% values of the bar plot, similar to the plot produced by ERRORBAR(Y,L,U).
%
% Note that it is impossible to input the X coordinates of the bars using
% this function.
%
% jitter_xtick = ERRORBAR_GROUPS(...) returns the X coordinates for the center
% of each group of bars.
%
% [...,hb,he] = ERRORBAR_GROUPS(...) returns the handles to the bars
% (produced by BAR) and the error bars (produced by ERRORBAR).
%
% ERRORBAR_GROUPS(...,'ParameterName',ParameterValue) allows customizing
% the display.
%
% 'jitter_width': scalar determining bar width. Must be between 0 and 1. A
% jitter_width value of 1 causes all the bars to touch each other. This
% parameter is identical to the optional width parameter of BAR. Default
% value: 0.9.
%
% 'errorjitter_width': scalar determining error bar width, expressed as a
% fraction of the width of the bars themselves. Default value: 0.75.
%
% 'jitter_colors': N-by-3 matrix determing the RGB values for bar colors. You
% must provide at least as many colors as there are groups in your plot.
% Default value: the function DISTINGUISHABLE_COLORS is used (cf.
% References).
%
% 'errorjitter_colors': N-by-3 matrix determining the RGB values for error bar
% colors. You must provide at least as many colors as there are groups in
% your plot. Default value: the error bar colors are set to black.
%
% 'optional_bar_arguments': cell array containing any 'PropertyName' -
% 'PropertyValue' input argument pair that you would like to pass on to
% BAR. Default value: no further input arguments are passed to BAR.
%
% 'optional_errorbar_arguments': cell array containing any 'PropertyName' -
% 'PropertyValue' input argument pair that you would like to pass on to
% ERRORBAR. Default value: {'LineStyle','none','Marker','none'}
% (causes ERRORBAR to plot no line between error bars, and no marker at the
% center of each error bar)
%
% 'jitter_names': cell string array containing labels to apply to each group
% of bars. You must provide at least as many labels as there are bars in
% your plot. Default value: the bars are numbered consecutively, starting
% from 1.
%
% 'FigID': handle to the figure where the bar plot will be plotted.
% Default value: a new figure is created.
%
% 'AxID': handle to the axes where the bar plot will be plotted. Default
% value: a new axes is created.
%
% Examples:
%
% Basic usage:
%   mean_vals=rand(3,8)/2+0.5;
%   std_vals=rand(3,8)/8;
%   errorbar_groups(mean_vals,std_vals);
%
% Set the lower bound of the error bars to 0, effectively plotting only the
% upper bound:
%   mean_vals=rand(4,6)/2+0.5;
%   errorbar_lower=zeros(size(mean_vals));
%   errorbar_upper=rand(4,6)/8;
%   errorbar_groups(mean_vals,errorbar_lower,errorbar_upper);
%
% When plotting fewer groups and bars, the plot might look better with
% thinner bars and error bars. This also shows how to input custom names
% for the groups of bars:
%   mean_vals=rand(2,4)/2+0.5;
%   std_vals=rand(2,4)/8;
%   errorbar_groups(mean_vals,std_vals, ...
%       'jitter_width',0.75,'errorjitter_width',0.5, ...
%       'jitter_names',{'A','B','C','D'});
%
% Here is how to pass optional input arguments to BAR and ERRORBAR:
%   mean_vals=rand(2,4)/2+0.5;
%   std_vals=rand(2,4)/8;
%   errorbar_groups(mean_vals,std_vals, ...
%       'jitter_width',0.75,'errorjitter_width',0.5, ...
%       'optional_bar_arguments',{'LineWidth',1.5}, ...
%       'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',1.5});
%
% References:
% This function calls the function DISTINGUISHABLE_COLORS, by Timothy E.
% Holy. http://www.mathworks.com/matlabcentral/fileexchange/29702
%
% See also BAR, ERRORBAR, DISTINGUISHABLE_COLORS.

% Author of this function: pierre.megevand@gmail.com
%
% change log:
% 2014/07/16 -- first uploaded to the MATLAB File Exchange. Does not
% do anything with the 'grp_names' optional input yet.
% 2014/07/18 -- The x axis is now adjusted to fit the plot.
% 2014/07/28 -- With David Groppe's help, added optional figure and axes
% handles.
% 2014/08/19 -- Simplified how default parameters are defined. Added option
% to skip plotting the error bars altogether.

% init defaults for parameters

Colors = [[12 73 153;7 96 32;126 18 33;210 113 26]/255; 1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
boolface = 1;

if length(size(Y)) == 3
    mean_vals = mean(Y,3)';
    std_vals = std(Y,[],3)';
elseif length(size(Y)) < 3
    mean_vals = mean(Y);
    std_vals = std(Y);
elseif length(size(Y)) > 3
    error('')
    return
end

[N_grps,N_jitters]=size(mean_vals);
jitter_width=0.9;
errorjitter_width=0.75;
jitter_colors = Colors(N_grps,:); % get some distinguishable colors for the bars!
errorjitter_colors=zeros(N_grps,3); % default errorbar color is black

jitter_names=strtrim(cellstr(num2str((1:N_jitters)')));
grp_names=strtrim(cellstr(num2str((1:N_grps)')));
FigID=[];
AxID=[];

% deal with the input arguments
if nargin<1 % the indispensable input arguments are not provided
    error('You need to provide at least the data');
else
    errorbar_lower=std_vals;
    errorbar_upper=std_vals;
    if any(size(mean_vals)~=size(std_vals)) % the indispensable input arguments must have the exact same size
        error('The size of ''mean_vals'' and ''std_vals'' must be the same.');
    else
        if numel(varargin)>0 % optional input arguments are provided
            while ~isempty(varargin)
                if numel(varargin)<2
                    error('You need to provide optional input arguments as ''ParameterName''-''ParameterValue'' pairs.');
                end
                switch varargin{1}
                    case 'jitter_width'
                        jitter_width=varargin{2};
                    case 'errorjitter_width'
                        errorjitter_width=varargin{2};
                    case 'jitter_colors'
                        jitter_colors=varargin{2};
                    case 'errorjitter_colors'
                        errorjitter_colors=varargin{2};
                    case 'optional_jitter_arguments'
                        optional_jitter_arguments=varargin{2};
                    case 'optional_errorjitter_arguments'
                        optional_errorjitter_arguments=varargin{2};
                    case 'jitter_names'
                        jitter_names=varargin{2};
                    case 'grp_names'
                        grp_names=varargin{2};
                    case 'FigID'
                        FigID=varargin{2};
                    case 'AxID'
                        AxID=varargin{2};
                    case 'boolface'
                        boolface=varargin{2};
                    otherwise
                        error('Unexpected ''ParameterName'' input: %s\n',varargin{1});
                end
                varargin(1:2)=[]; % this pair of optional input arguments has been dealt with -- remove...
            end
        end
    end
end

% init group width and bar shift
shift_span=(1-jitter_width)*(N_grps-1);
jitter_shift=linspace(-shift_span/2,+shift_span/2,N_grps);

% compute position of group x ticks
% jitter_xtick=N_grps/2+0.5:N_grps:N_jitters*N_grps-N_grps/2+0.5;
jitter_separation = shift_span;
jitter_xtick=(N_grps)/2+0.5:(N_grps+jitter_separation):N_jitters*(N_grps+jitter_separation)-(N_grps+jitter_separation)/2+0.5;



% init handles vectors
hb=zeros(N_grps,1);
he=zeros(N_grps,1);

% init figure
if ~isempty(FigID)
    figure(FigID);
else
    figure;
    FigID=gcf;
end
if ~isempty(AxID)
    set(FigID,'CurrentAxes',AxID);
    cla;
else
    clf;
end
hold on;


% plot the bars themselves
markers = strvcat('o','s','^','d','p','>','v');
for grp=1:N_grps
    Pos = (grp:N_grps+jitter_separation:N_jitters*(N_grps+jitter_separation)-(N_grps+jitter_separation-grp))-jitter_shift(grp);
    % % %     for jit = 1:N_jitters
    % % %         Ytemp = nonzeros(squeeze(Y(jit,grp,:)));
    % % %
    % % %         % Detecting outliers
    % % %         indoutliers = find(abs(Ytemp - (mean(Ytemp)+std(Ytemp))) > 2*std(Ytemp));
    % % %         Ytemp(indoutliers) = [];
    % % %
    % % %         jitterPos = Pos(jit) + (2*(rand(size(Ytemp))-0.5)*0.2);
    % % %         if sum(Ytemp) > 0
    % % %             if boolface
    % % %                 h = plot(jitterPos,Ytemp,'Marker',markers(grp),'MarkerFaceColor',Colors(grp,:), optional_jitter_arguments{:});
    % % %             else
    % % %                 h = plot(jitterPos,Ytemp,'Marker',markers(grp),'MarkerEdgeColor',Colors(grp,:), 'Markersize',4,'LineStyle','none');
    % % %             end
    % % %             hold on;
    % % %             %boxplot(gca,Ytemp,'positions',Pos(jit));
    % % %             %         hl = line([min(jitterPos) max(jitterPos) ],[mean(Ytemp) mean(Ytemp)],'Color',[1 1 0],optional_errorjitter_arguments{:}); % Males Mean
    % % %             %         he(grp)=errorbar(Pos(jit), mean(Ytemp), std(Ytemp), std(Ytemp),optional_errorjitter_arguments{:});
    % % %         end
    % % %
    % % %     end
    
    AA = repmat(Pos,[size(Y,3) 1]) + (2*(rand(size(squeeze(Y(:,grp,:))'))-0.5)*0.2);
    if boolface
        h = plot(AA',squeeze(Y(:,grp,:)),'Marker',markers(grp),'MarkerFaceColor',Colors(grp,:), optional_jitter_arguments{:})
    else
        h = plot(AA',squeeze(Y(:,grp,:)),'Marker',markers(grp),'MarkerEdgeColor',Colors(grp,:), 'Markersize',4,'LineStyle','none');
    end
    test = 1;
    
    
    %boxplot(gca,squeeze(Y(:,grp,:))','positions',Pos,'plotstyle','compact');
    
    
end

% set the x tick labels
set(gca,'XTick',jitter_xtick,'XTickLabel',jitter_names);
set(gca,'FontName','Arial');

% cosmetic fine-tuning of the figure
set(gca,'XLim',[0 jitter_xtick(end)+jitter_xtick(1)]); % adjusts the x axis to the plot
gapval = std(Y(:));
set(gca,'YLim',[min(Y(:))-gapval  max(Y(:))+gapval ]); % adjusts the y axis to the plot
set(gca,'FontSize',10); % adjusts the y axis to the plot

hold off;

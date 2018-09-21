function varargout = Plot_Plane(varargin);
%
% Syntax :
%  hf = Plot_Plane(planEq, centPoint ,scalVar, opts);
%
% This script plots a plane at a given point ([Xc Yc Zc]). The plane equation must be
% specified as a 1x4 vector ([A B C D]);
%
% Input Parameters:
%       planEq                  : Plane equation must be specified as a 
%                                 1x4 vector ([A B C D]);
%       centPoint               : Coordinates to center of the plane ([Xc
%                                 Yc Zc]). If the coordinate center is not 
%                                 specified the the origin will be set as 
%                                 coordinate center.
%       scalVar                 : Scale factor. If the scale factor is not
%                                 specified then the scale factor will be
%                                 equal to 1. 
%       opts                    : Options
%                                 opts.alpha: Transparency of the patch
%                                 opts.nfigure: Boolean variable to create
%                                 or not a new figure.
%
% Output Parameters:
%       SulcMetrics             : Sulcus metrics
%        Surfo                  : Reparametrized Sulcus
%        SurfL                  : Curves: Topline, Bottom line, Length,
%                                 Depth and width curves
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 13th 2014
% Version $1.0

%% =========================== Input parameters  =========================%

if nargin < 1
    error('The plane equation is mandatory');
    return
elseif nargin == 1
    planEq = varargin{1};
    centPoint = [0 0 0];
    scalVar = 1;
    opts.alpha = .4;
    opts.nfigure = 0;
elseif nargin == 2
    planEq = varargin{1};
    centPoint = varargin{2};
    scalVar = 1;
    opts.alpha = .4;
    opts.nfigure = 0;
elseif nargin == 3
    planEq = varargin{1};
    centPoint = varargin{2};
    scalVar = varargin{3};
    opts.alpha = .4;
    opts.nfigure = 0;
elseif nargin == 4
    planEq = varargin{1};
    centPoint = varargin{2};
    scalVar = varargin{3};
    opts = varargin{4};
    if ~isfield(opts,'alpha')
        opts.alpha = 0.4;
    end
    if ~isfield(opts,'nfigure')
        opts.alpha = 0;
    end
elseif nargin > 4
    error('Please enter only four inputs as maximum');
    return
end

%% ======================= End of Input parameters  ======================%

%% ============================== Main Program  ==========================%

% % % % % % % nn = planEq(1:3)./repmat(sqrt(sum((planEq(1:3)').^2))',[1 3]); %figure; Plot_Normals(Vert, nn,1);
% % % % % % % betas = acos(nn(:,3)./sqrt((nn(:,1).^2+nn(:,3).^2)));
% % % % % % % gammas = acos(sqrt((nn(:,1).^2+nn(:,3).^2))./ones(size(nn,1),1));
% % % % % % % aa = cross(nn,repmat([0 0 1],[size(nn,1) 1]),2);indf = find(aa(:,1)>0);gammas(indf) =-1*gammas(indf); betas(indf) =-1*betas(indf);
% % % % % % % 
% % % % % % % X = [-1 1 1 -1]*scalVar;
% % % % % % % Y = [-1 -1 1 1]*scalVar;
% % % % % % % X = X(:);Y = Y(:);Z = X(:)*0;
% % % % % % % 
% % % % % % % x = gammas;
% % % % % % % y = betas;
% % % % % % % z = 0;
% % % % % % % Rx = [1 0 0 0; 0 cos(x) -sin(x) 0; 0 sin(x) cos(x) 0; 0 0 0 1];
% % % % % % % Ry = [cos(y) 0 sin(y) 0; 0 1 0 0; -sin(y) 0 cos(y) 0; 0 0 0 1];
% % % % % % % Rz = [cos(z) -sin(z) 0 0; sin(z) cos(z) 0 0; 0 0 1 0; 0 0 0 1];
% % % % % % % Mat = Rx*Ry*Rz;
% % % % % % % %     t = [X Y Z ones(size(X,1),1)]*Mat';t = t(:,1:3);
% % % % % % % t = [X Y Z ones(size(X,1),1)]*Mat';t = t(:,1:3);
% % % % % % % X1 =  t(:,1)+centPoint(1);
% % % % % % % Y1 =  t(:,2)+centPoint(2);
% % % % % % % Z1 =  t(:,3)+centPoint(3);
% % % % % % % 
% % % % % % % if opts.nfigure
% % % % % % %     colordef black;h = figure('numbertitle','off','Color',[0 0 0],'Position',[0 0 1200 900],'name',['IBASPM Plane Ploting...:  ']);hold on;  grid on; box on; % Creating new figure
% % % % % % %     plot3(centPoint(1),centPoint(2),centPoint(3),'.b','Markersize',30); % patch plotting
% % % % % % %     strsurf = patch(X1,Y1,Z1,'b'); % patch plotting
% % % % % % %     varargout{1} = strsurf;
% % % % % % % else
% % % % % % %     hold on;
% % % % % % %     strsurf = patch(X1,Y1,Z1,'b'); % patch plotting
% % % % % % %     plot3(centPoint(1),centPoint(2),centPoint(3),'.b','Markersize',30); % patch plotting
% % % % % % %     varargout{1} = strsurf;
% % % % % % % end
% % % % % % % set(strsurf,'FaceAlpha',opts.alpha);

X = [-1 1]';
Y = [-1 1]';
Z = [-1 1]';

%[X,Y,Z] = sphere(10); % Creating the unit sphere
x = X(:)*scalVar + centPoint(1); % scaling coordinates
y = Y(:)*scalVar + centPoint(2); % scaling coordinates
z = Z(:)*scalVar + centPoint(3); % scaling coordinates

% Plane equation
A = planEq(1);
B = planEq(2);
C = planEq(3);
D = planEq(4);

% Detecting Limits
xLim = [min(x) max(x)];
zLim = [min(z) max(z)];
[X,Z] = meshgrid(xLim,zLim);
if B~=0
    Y = (A * X + C * Z + D)/ (-B+eps); % Computing the Z coordinate
else
    
end
reOrder = [1 2  4 3];
if opts.nfigure
    colordef black;h = figure('numbertitle','off','Color',[0 0 0],'Position',[0 0 1200 900],'name',['IBASPM Plane Ploting...:  ']);hold on;  grid on; box on; % Creating new figure
    plot3(centPoint(1),centPoint(2),centPoint(3),'.b','Markersize',30); % patch plotting
    strsurf = patch(X(reOrder),Y(reOrder),Z(reOrder),'b'); % patch plotting
    varargout{1} = h;
else
    hold on;
    strsurf = patch(X(reOrder),Y(reOrder),Z(reOrder),'b'); % patch plotting
    plot3(centPoint(1),centPoint(2),centPoint(3),'.b','Markersize',30); % patch plotting
    varargout{1} = '';
end
set(strsurf,'FaceAlpha',opts.alpha);
%% =========================== End of Main Program   =====================%
return;
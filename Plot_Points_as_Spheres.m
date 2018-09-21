function h = Plot_Points_as_Spheres(Coord, varargin);
%
% Syntax :
%     h = Plot_Points_as_Spheres(Coord, varargin);
%
% This function plots points as spheres.
%
% Input Parameters:
%        Coord                          : Point Coordinates: (Nx3) matrix
%        Radius                         : Spheres radius: (Nx1) vector. If
%                                         N = 1 then the radius will be
%                                         same for all the spheres.
%        Color                          : Spheres colors: (Nx3) Matrix. If
%                                         N = 1 then the color will be
%                                         same for all the spheres.
%        FaceAlpha                      : Spheres transparency: (Nx1) vector. If
%                                         N = 1 then the transparency will be
%                                         same for all the spheres.
%       NumSphFaces                     : Number of faces.
%
% Output Parameters:
%         h                             : Spheres surface handles
%
%
% Examples: h = Plot_Points_as_Spheres(Coord); 
%           h = Plot_Points_as_Spheres(Coord,'Radius',3,'Color',[.7 .7 .7]); 
% 
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% November 13th 2014
% Version $1.0



%% ====================== Checking input parameters ===================== %
if nargin<1 % the indispensable input arguments are not provided
    error('One input is mandatory');
else
    % Parameters
    Radius = 1;            % Sphere Radius
    NumSphFaces = 30;      % Number of faces
    Color = [.7 .7 .7];    % Faces Colors
    FaceAlpha = 1;         % Face Alpha
end

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
                case 'Radius' % Sphere Radius
                    Radius=varargin{2};
                case 'Color' % Face Color
                    Color=varargin{2};
                case 'FaceAlpha' % Face Alpha: Transparency Value
                    FaceAlpha=varargin{2};
                case 'NumSphFaces'
                    NumSphFaces=varargin{2};
                case 'FigID'
                    FigID=varargin{2};
                otherwise
                    error('Unexpected ''ParameterName'' input: %s\n',varargin{1});
            end
            varargin(1:2)=[]; % this pair of optional input arguments has been dealt with -- remove...
        end
    end
end

if length(Radius) == 1
    Radius = repmat(Radius, size(Coord,1), 1);
end

if size(Color,1) == 1
    Color = repmat(Color, size(Coord,1), 1);
end

if length(FaceAlpha) == 1
    FaceAlpha = repmat(FaceAlpha, size(Coord,1), 1);
end

%% ================== End of Checking input parameters ================= %%

%% ============================ Main Program =========================== %%

[SX,SY,SZ] = sphere(NumSphFaces);
washeld = ishold();
hold on;
x = Coord(:,1);
y = Coord(:,2);
z = Coord(:,3);
for K = 1 : length(x)
    h(K) = surf(SX*Radius(K) + x(K), SY*Radius(K) + y(K), SZ*Radius(K) + z(K));
    h(K).EdgeColor = 'none';
    h(K).FaceColor = Color(K,:);
    h(K).FaceAlpha = FaceAlpha(K);
end
axis image;
lighting phong;

if ~washeld
    hold off
end
%% ======================= End of Main Program ========================= %%

return
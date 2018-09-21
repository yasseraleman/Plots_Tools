function hf = Plot_Circular_Graph_Org(CMatrix, param,CharFile);
%
% Syntax :
% Plot_Circular_Graph(CMatrix);
%
% This script plots a "bullseye" according to a specified connectivity
% matrix. Different connection colors means different networks.
%
% Input Parameters:
%       CMatrix                 : Connectivity Matrix. It can be a textfile
%       CharFile                : Text file with structures characteristics
%                               (cortical thickness, area, girification index)
%       param                   : Plotting parameters
%                                param.linethick - Lines thickness
%                                param.connecthr - Connectivity threshold
%
%                                Rings parameters
%                                param.N - Patches resolution
%                                param.rho - Interior and exterior rings radius
%                                param.tht - Rings theta limits
%                                param.tht0  =0;
%                                param.label =0;
%                                param.lines =[1];
%                                param.ringwidth - Ring Width
%                                param.intringwidth - Inter ring distance
%                                param.blacklines - Boolean variable to specified colored or black lines
%                                param.figcolor - Figure color
%                                param.textspace - Text space for multiple colorbars
%
% Output Parameters:
%      hf                       : Figure Handle
%
% See also: Examples proposed at line 40.
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% November 13th 2013
% Version $1.0

% --------------------- Some Examples ------------------------------------
temp = which('Plot_Circular_Graph.m');
hf = '';
atype = 'aparc+aseg';
% Textfile - Connectivity Matrix and Characteristics File
% [pth, nm, ext] = fileparts(temp);                         % Uncomment this line to see this example
% CharFile = [pth filesep 'Thick+Area.txt'];                % Uncomment this line to see this example
% CMatrix = [pth filesep 'Pvalue_Matrix-aparc+aseg.txt'];   % Uncomment this line to see this example

% -------------------- Just Connectivity Matrices -------------------------
% ----- No 1 Single connectivity matrix
%  CMatrix = zeros(20,20);      % Uncomment this line to see this example
%  CMatrix(1,[2 3]) = 0.5;    % Uncomment this line to see this example
%  CMatrix(4,[5 8]) = 1;      % Uncomment this line to see this example
%  CMatrix([2 4 5 6],8) = 1;  % Uncomment this line to see this example
%  C = (CMatrix+CMatrix')/2;  % Uncomment this line to see this example

% ----- No 2 Single connectivity matrix but with a predefined structure organization
% CMatrix = zeros(68,68);  % Uncomment this line to see this example
% CMatrix(33,56) = 1;      % Uncomment this line to see this example
% CMatrix(28,54) = 1;      % Uncomment this line to see this example
% CMatrix(54,28) = 1;      % Uncomment this line to see this example
% CMatrix(26,54) = 1;      % Uncomment this line to see this example
% CMatrix(54,26) = 1;      % Uncomment this line to see this example
% CMatrix(34,68) = 1;      % Uncomment this line to see this example
% CMatrix(68,34) = 1;      % Uncomment this line to see this example

% -------------------- Just Connectivity Matrices -------------------------
% ----- No 3 Single connectivity matrix
%  CMatrix = rand(35,35);     % Uncomment this line to see this example
%  ind = find(CMatrix < 0.990);
%  CMatrix(ind)= 0;
%  ind = find(CMatrix);
% a = randi(142,70,1);
% CMatrix(ind(a)) = 0;
% CMatrix = (CMatrix+CMatrix')/2;  % Uncomment this line to see this example

% Textfile - Connectivity Matrix Only
% CMatrix = [pth filesep '0007-20091002-Connectivity_Matrix_FAwACP-aparc+aseg.txt'];  % Uncomment this line to see this example


% --------------------- End of Some Examples ------------------------------
%% ====================== Freesurfer parcellation ========================%
% Subcortical Structures
SubCL = [10:13 17:18 26];
SubCR = [ 49:54 58];

% Selecting Frontal regions
FroIds = sort([1028 1003 1027 1018 1019 1020 1012 1014 1024 1017 1032 1026 1002]);
FroIdsR = FroIds + 1000;

% Selecting Temporal regions
TempIds = sort([1009 1015 1030 1001 1007 1034 1006 1033 1016]);
TempIdsR = TempIds + 1000;

% Selecting Parietal regions
ParIds = sort([1029 1008 1031 1022 1025 1010 1023]);
ParIdsR = ParIds + 1000;

% Selecting Occipital regions
OccIds = sort([1011 1013 1005 1021]);
OccIdsR = OccIds + 1000;

% Selecting Insula regions
InsIds = [1035];
InsIdsR = [2035];
OldOrg = sort([SubCL SubCR FroIds FroIdsR ParIds ParIdsR TempIds TempIdsR OccIds OccIdsR InsIds InsIdsR]);
%% ====================== End of Freesurfer parcellation =================%

%% ===================== Checking input parameters ====================== %
% close all;

% =================== Detecting Connectivity type ====================== %

if ischar(CMatrix)
    Connect = Read_Connectivity_Matrix(CMatrix,0);
    Connect.Matrix = triu(Connect.Matrix);
    C = (Connect.Matrix+Connect.Matrix')/max(Connect.Matrix(:));
    Names  = Connect.StructNames;
elseif isstruct(CMatrix)
    if isfield(CMatrix,'Matrix')
        C = CMatrix.Matrix;
        C = triu(C);
        C = (C+C')/max(C(:));
        if isfield(CMatrix,'StructNames')
            Names  = CMatrix.StructNames;
        else
            switch size(CMatrix,1)
                case 66
                    [Names,AcronSt,ColorsSt] = Load_Struct_Names(atype,66);
                case 68
                    [Names,AcronSt,ColorsSt] = Load_Struct_Names(atype,68);
                case 82
                    [Names,AcronSt,ColorsSt] = Load_Struct_Names(atype,82);
                otherwise
                    Names = num2str([1:size(C,1)]','%02d');
            end
        end
    else
        error('Connectivity Matrix does not exist');
        return;
    end
else
    C = CMatrix;
    C = triu(C);
    C = (C+C')/max(C(:));
    switch size(CMatrix,1)
        case 66
            [Names,AcronSt,ColorsSt] =  Load_Struct_Names(atype,66);
        case 68
            [Names,AcronSt,ColorsSt] = Load_Struct_Names(atype,68);
        case 82
            [Names,AcronSt,ColorsSt] = Load_Struct_Names(atype,82);
        otherwise
            Names = num2str([1:size(C,1)]','%02d');
    end
end
% =================== End of Detecting Connectivity type ================%

if nargin <2
    param.linethick = [2 5];      % Lines thickness
    param.connecthr = 0;        % Connectivity threshold
    param.blacklines = 0;         % Boolean variable to specified colored or black lines
    param.figcolor = 'white';     % Figure color
    param.textspace = 0.01;       % Text space for multiple colorbars
    param.newfigure = 1    ;      % Boolean variable to create new figure
    param.figurevisible = 'on';   % Figure visibility
    
    % Rings parameters
    param.N     = 10;             % Patches resolution
    param.rho   =[4 5];           % Interior and exterior rings radius
    param.tht   =[0 360];         % Rings theta limits
    param.tht0  =0;
    param.label =0;
    param.lines =[1];
    param.ringwidth = 0.15;       % Ring Width
    param.intringwidth = 0.025;   % Inter ring distance
else
    if ~isfield(param,'linethick')
        param.linethick = [2 5];      % Lines thickness
    end
    if ~isfield(param,'connecthr')
        param.connecthr = 0;        % Connectivity threshold
    end
    if ~isfield(param,'blacklines')
        param.blacklines = 0;         % Boolean variable to specified colored or black lines
    end
    if ~isfield(param,'figcolor')
        param.figcolor = 'black';     % Figure color
    end
    if ~isfield(param,'textspace')
        param.textspace = 0.01;       % Text space for multiple colorbars
    end
    if ~isfield(param,'newfigure')
        param.newfigure = 1   ;       % Boolean variable to create new figure
    end
    if ~isfield(param,'figurevisible')
        param.figurevisible = 'on';   % Figure visibility
    end
    
    % Rings parameters
    if ~isfield(param,'N ')
        param.N     = 10;             % Patches resolution
    end
    if ~isfield(param,'rho')
        param.rho   =[4 5];           % Interior and exterior rings radius
    end
    if ~isfield(param,'tht')
        param.tht   =[0 360];         % Rings theta limits
    end
    if ~isfield(param,'tht0')
        param.tht0  =0;
    end
    if ~isfield(param,'label')
        param.label =0;
    end
    if ~isfield(param,'lines')
        param.lines =[1];
    end
    if ~isfield(param,'ringwidth')
        param.ringwidth = 0.15;       % Ring Width
    end
    if ~isfield(param,'intringwidth')
        param.intringwidth = 0.025;   % Inter ring distance
    end
end

%% ================== End of Checking input parameters ================== %

%% ================= Colormaps for multiple colorbars =================== %
colorstrings = strvcat('hot(256)','bone(256)','hsv(256)','jet(256)','pink(256)','winter(256)','gray(256)');
%% ================ End of Colormaps for multiple colorbars ============= %

% ---------  XY coordinates of nodes in a circular layout  ---------------%
N = size(C,1);
ang = linspace(0,2*pi,N+1); ang = ang(2);
switch N
    case 66
        Ndel = [14 22 32];
        theta = linspace(pi/2+ang,3/2*pi-ang,(N/2+3))'; %theta(end) = [];
        theta(Ndel) = [];
        theta2 = linspace(pi/2-ang,-pi/2+ang,(N/2+3))'; %theta2(end) = [];
        theta2(Ndel) = [];
        theta = [theta;theta2];
    case 68
        Ndel = [14 22 32 37];
        theta = linspace(pi/2+ang,3/2*pi-ang,(N/2+4))'; %theta(end) = [];
        theta(Ndel) = [];
        theta2 = linspace(pi/2-ang,-pi/2+ang,(N/2+4))'; %theta2(end) = [];
        theta2(Ndel) = [];
        theta = [theta;theta2];
    case 82
        Nsubc = 7;
        Ndel = [8 22 30 40 45];
        theta = linspace(pi/2+ang,3/2*pi-ang,(N/2+5))'; %theta(end) = [];
        theta(Ndel) = [];
        theta2 = linspace(pi/2-ang,-pi/2+ang,(N/2+5))'; %theta2(end) = [];
        theta2(Ndel) = [];
        theta = [theta;theta2];
    otherwise
        theta = linspace(pi/2,5/2*pi,N+1)'; theta(end) = [];
        
end


%theta = linspace(pi/2,3/2*pi,N+1)'; theta(end) = [];
xy = param.rho(1).* [cos(theta) sin(theta)];

% -------------------- Network Labelling -------------------------------- %

% -------------------- Network Colors ----------------------------------- %
%LabNet = SubNetworks_Extraction(C,param.connecthr);

LabNet = Label_Graph_Components(C,param.connecthr);
%LabNet = SubNetworks_Labelling(C, param.connecthr);
% %
% % % Three Networks
% % LabNet(:,31) = LabNet(:,31)*2;
% % LabNet(31,:) = LabNet(31,:)*2;
% % % LabNet(:,19) = LabNet(:,19)*3;
% % % LabNet(19,:) = LabNet(19,:)*3;




disp(['Number of networks: ' num2str(max(LabNet(:))) '. Threshold = ' num2str(param.connecthr)]);
% LabNet = Label_Networks(C, param.connecthr);
if sum(LabNet(:))~=0
    LabNet = triu(LabNet);
    NNetworks = unique(nonzeros(LabNet(:)));
    NN = length(NNetworks);
    col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
    Ncolor = size(col,1);
    re = floor(NN/Ncolor);
    Col = repmat(col,[re+1 1]);
    c = accumarray(nonzeros(LabNet(:)),nonzeros(LabNet(:))*0+1); [a,b] = sort(c,'descend'); Col(b,:) = Col(1:length(c),:); % Selecting Colors
    Colors = Col(LabNet(find(LabNet)),:);
    if param.blacklines
        switch param.figcolor
            case 'black'
                Colors = Colors*0+1;
            case 'white'
                Colors = Colors*0;
        end
    end
    
    % -------------------- Points organization ------------------------------ %
    ind = find(LabNet ~=0);
    [xrow,xcol] = ind2sub(size(LabNet),ind);
    nxrow = xrow;
    nxcol = xcol;
    inddel = find((xcol-xrow)>N/2); % Detecting initial point
    nxrow(inddel) = xcol(inddel);
    nxcol(inddel) = xrow(inddel);
    xrow = nxrow;
    xcol = nxcol;
    
    % -------------------- Lines Thickness ---------------------------------- %
    [ut,i1t,it2] = unique(C(ind));
    LG = linspace(param.linethick(1),param.linethick(2),length(ut))';
    LG = LG(it2);
    
    % -----------------------  Bisector radius  ------------------------------%
    b = round(dist(xy')*10000)/10000;
    [ut,i1t,it2] = unique(b(:));
    bis_rad = linspace(0.8*param.rho(1),0.1*param.rho(1),length(ut));
    r_bis = reshape(bis_rad(it2),[size(b)]);
    r_bis = r_bis(ind);
    
    % -----------------------  Bisector angles  ------------------------------%
    theta_bis = theta(xrow) - (theta(xrow)-theta(xcol))/2;
    
    % -----------------------  Oposite angles  -------------------------------%
    
    if (N == 66)|(N == 68)|(N == 82)
        %indm = find((xrow-xcol)>0);
        indm = find(abs([theta(xrow) - theta(xcol) ])>pi);
        %     indr = find((xrow>=N/2)&(xcol<=(N/2+1)));
        %     indm = [indm;indr];
        if ~isempty(indm)
            theta_bis(indm) = theta_bis(indm)+pi;
        end
    end
    
    % ------------------- XY Bisector Coordinates  ---------------------------%
    xy_bis = [r_bis r_bis].* [cos(theta_bis) sin(theta_bis)];
end
% ------------------------ Labels of nodes -------------------------------%
if exist('Names','var')
    txt = cellstr(Names);
else
    txt = cellstr(num2str((1:N)','%02d'));
end
if param.newfigure
    switch param.figcolor
        case 'black'
            colordef black; hf = figure('numbertitle','off','name','Circular Graph','Color',[0 0 0],'Position',[0 0 1200 900],'Visible',param.figurevisible,'InvertHardcopy','off');
        case 'white'
            colordef white; hf = figure('numbertitle','off','name','Circular Graph','Color',[1 1 1],'Position',[0 0 1200 900],'Visible',param.figurevisible,'InvertHardcopy','off');
    end
end
% --------------- Show nodes and edges of the circular layout ----------- %
% if (N ~= 68)&(N ~= 82)
%     line(xy(:,1), xy(:,2), 'LineStyle','none', ...
%         'Marker','.', 'MarkerSize',30, 'Color','b')
% end

% --------------- Show nodes and edges of the bisectors ----------------- %
% hold on
% line(xy_bis(:,1), xy_bis(:,2), 'linestyle','none', ...
%     'marker','.', 'markersize',30, 'color','k')
%
% % ------------------ Show bisectors labels of the  ---------------------- %
% txt = cellstr(num2str((1:length(ind))','%02d'));
% h = text(xy_bis(:,1).*0.85, xy_bis(:,2).*0.85, txt, 'FontSize',8);
% set(h, {'Rotation'},num2cell(theta_bis*180/pi));

axis off;axis equal;
if sum(LabNet(:))~=0
    % ----------------- Points to create the parabolas -----------------------%
    Par2plot = [xy(xrow,:) xy_bis xy(xcol,:) ];
    hold on
    MaxLG = max(LG);
        for i = 1:size(Par2plot,1)
            x1 = Par2plot(i,1:2);
            x2 = Par2plot(i,3:4);
            x3 = Par2plot(i,5:6);
    
            % ------------------------  Rotation Matrix  -----------------------------%
            R = [cos(-theta_bis(i)) -sin(-theta_bis(i)); sin(-theta_bis(i)) cos(-theta_bis(i))];
    
            % ------------------------ Rotating the points  --------------------------%
            t = (R*[[x1(1); x2(1);x3(1)] [x1(2); x2(2);x3(2)] ]')';
    
            % -------------------- Parabola equation estimation ----------------------%
            Y=[t(1,1);t(2,1);t(3,1)];
            A=[t(1,2)^2 t(1,2) 1;t(2,2)^2 t(2,2) 1;t(3,2)^2 t(3,2) 1];
            X=inv(A)*Y;
    
            % ------------------------- Points evaluation ----------------------------%
            x=[0:t(3,2)/30:t(3,2)]';
            Y=X(1)*x.^2+ +X(2)*x+ X(3);
            Xr = [flipdim(Y,1);Y];
            Yr = [flipdim(x,1);-x];
    
            % ------------------- Inverse points rotation ----------------------------%
            xc = (inv(R)*[Xr(:) Yr(:)]')';
            hold on
            hline = line(xc(:,1),xc(:,2)','LineWidth',LG(i),'Color',(LG(i)/MaxLG)*Colors(i,:));
    
        end
    
    
    
%     % t = linspace(-pi,pi,length(adjacencyMatrix) + 1).'; % theta for each node
%     [row,col, v] = find(LabNet);
%     for i = 1:length(v)
%         if row(i) ~= col(i)
%             if abs(row(i) - col(i)) - length(LabNet)/2 == 0
%                 % points are diametric, so draw a straight line
%                 u = param.rho(1)*[cos(theta(row(i)));sin(theta(row(i)))];
%                 v = param.rho(1)*[cos(theta(col(i)));sin(theta(col(i)))];
%                 hline = line(...
%                     [u(1);v(1)],...
%                     [u(2);v(2)],...
%                     'LineWidth', LG(i),...
%                     'Color', (LG(i)/MaxLG)*Colors(i,:),...
%                     'PickableParts','none');
%             else % points are not diametric, so draw an arc
%                 u  = param.rho(1)*[cos(theta(row(i)));sin(theta(row(i)))];
%                 v  = param.rho(1)*[cos(theta(col(i)));sin(theta(col(i)))];
%                 x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
%                 y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
%                 r  = sqrt(x0^2 + y0^2 - 1);
%                 thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
%                 thetaLim(2) = atan2(v(2)-y0,v(1)-x0);
%                 
%                 if u(1) >= 0 && v(1) >= 0
%                     % ensure the arc is within the unit disk
%                     thetat = [linspace(max(thetaLim),pi,50),...
%                         linspace(-pi,min(thetaLim),50)].';
%                 else
%                     thetat = linspace(thetaLim(1),thetaLim(2)).';
%                 end
%                 
%                 hline = line(...
%                     r*cos(thetat)+x0,...
%                     r*sin(thetat)+y0,...
%                     'LineWidth', LG(i),...
%                     'Color', (LG(i)/MaxLG)*Colors(i,:),...
%                     'PickableParts','none');
%             end
%         end
%     end
    
    
    
end

rmaxLim = param.rho(1);
Nstruct = size(C,1);
%  Constructing Bullseye
% % %% =================== Reorganizing structures(Freesurfer) ===============%
% % if (Nstruct == 68)|(Nstruct == 66)|(Nstruct == 82)
% %     [GMcodes,NamesO, Colors] = Gray_Matter_codes('aparc+aseg');
% %     if Nstruct == 66
% %         Norg = [FroIds ParIds TempIds OccIds FroIdsR ParIdsR TempIdsR OccIdsR ];
% %         [a,it] = ismember(Norg,OldOrg);
% %     elseif Nstruct == 68
% %         Norg = [FroIds ParIds TempIds OccIds InsIds FroIdsR ParIdsR TempIdsR OccIdsR InsIdsR];
% %         [a,it] = ismember(Norg,OldOrg);
% %     elseif Nstruct == 82
% %         Norg = [SubCL  FroIds ParIds TempIds OccIds InsIds SubCR  FroIdsR ParIdsR TempIdsR OccIdsR  InsIdsR];
% %         [a,it] = ismember(Norg,OldOrg);
% %     end
% %
% % end
% % %% ============ End of Reorganizing structures(Freesurfer) ===============%



%% =========================== Drawing Data rings ======================= %
ColorTotal = [0 0 0];
if exist('CharFile','var')
    if exist(CharFile,'file')
        rorig = rmaxLim;
        % The matrices are defined from left to right...graph is CCW
        data = textread(CharFile,'%f','delimiter',';');
        fid = fopen(CharFile,'rt');
        line1 = fgetl(fid);
        ind = strfind(line1,';');
        Nc = length(ind) + 1;
        Nr = size(data,1)/Nc;
        data = reshape(data,[Nc Nr])';
        data = data./repmat(max(data),[size(data,1) 1]);
        Nrings = size(data,2);
        data = data(1:size(C,1),:);
        if Nstruct == 68|Nstruct == 82
            %data = data(it,:);
        end
        for i = 1: Nrings
            rminLim = (rorig + i*param.intringwidth + (i-1)*param.ringwidth);
            rmaxLim = (rorig + i*param.intringwidth + i*param.ringwidth);
            param.rho = [rminLim rmaxLim];
            minplaneimg = min(min(data(:,i))); % find minimum first.
            scaledimg = (floor(((data(:,i) - minplaneimg) ./ ...
                (max(max(data(:,i))) - minplaneimg)) * 255)); % perform scaling
            
            % convert the image to a true color image with the jet colormap.
            ColorTotal = [ColorTotal;eval(colorstrings(i,:))];
            colorimg = squeeze(ind2rgb(scaledimg,eval(colorstrings(i,:))));
            switch Nstruct
                case 66
                    %% ======================= Left Hemisphere =========================== %
                    % Frontal Left
                    Nfront = 13;
                    Ninterval = [1:Nfront];
                    param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h=patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Parietal Left
                    Nparie = 7;
                    Ninterval = [Nfront+1:Nfront+Nparie];
                    param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h=patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Temporal Left
                    Ntempo = 9;
                    Ninterval = [Nfront+Nparie+1:Nfront+Nparie+Ntempo];
                    param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h=patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Occipital Left
                    Nocc = 4;
                    Ninterval = [Nfront+Nparie+Ntempo+1:Nfront+Nparie+Ntempo+Nocc];
                    param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h=patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    
                    %% ====================== Right Hemisphere =========================== %
                    % Frontal Right
                    Nfront = 13;
                    Ninterval = [1:Nfront] + Nstruct/2;
                    param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h=patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Parietal Right
                    Nparie = 7;
                    Ninterval = [Nfront+1:Nfront+Nparie] + Nstruct/2;
                    param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h=patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Temporal Right
                    Ntempo = 9;
                    Ninterval = [Nfront+Nparie+1:Nfront+Nparie+Ntempo] + Nstruct/2;
                    param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h=patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Occipital Right
                    Nocc = 4;
                    Ninterval = [Nfront+Nparie+Ntempo+1:Nfront+Nparie+Ntempo+Nocc] + Nstruct/2;
                    param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h=patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                case 68
                    
                    %% ======================= Left Hemisphere =========================== %
                    % Frontal Left
                    Nfront = 13;
                    Ninterval = [1:Nfront];
                    param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h=patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Parietal Left
                    Nparie = 7;
                    Ninterval = [Nfront+1:Nfront+Nparie];
                    param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h=patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Temporal Left
                    Ntempo = 9;
                    Ninterval = [Nfront+Nparie+1:Nfront+Nparie+Ntempo];
                    param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h=patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Occipital Left
                    Nocc = 4;
                    Ninterval = [Nfront+Nparie+Ntempo+1:Nfront+Nparie+Ntempo+Nocc];
                    param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h=patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Insula Left
                    Nins = 1;
                    Ninterval = [Nfront+Nparie+Ntempo+Nocc+1:Nfront+Nparie+Ntempo+Nocc+Nins];
                    param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h=patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    %% ====================== Right Hemisphere =========================== %
                    % Frontal Right
                    Nfront = 13;
                    Ninterval = [1:Nfront] + Nstruct/2;
                    param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h=patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Parietal Right
                    Nparie = 7;
                    Ninterval = [Nfront+1:Nfront+Nparie] + Nstruct/2;
                    param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h=patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Temporal Right
                    Ntempo = 9;
                    Ninterval = [Nfront+Nparie+1:Nfront+Nparie+Ntempo] + Nstruct/2;
                    param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h=patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Occipital Right
                    Nocc = 4;
                    Ninterval = [Nfront+Nparie+Ntempo+1:Nfront+Nparie+Ntempo+Nocc] + Nstruct/2;
                    param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h=patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Insula Right
                    Nins = 1;
                    Ninterval = [Nfront+Nparie+Ntempo+Nocc+1:Nfront+Nparie+Ntempo+Nocc+Nins] + Nstruct/2;
                    param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h=patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                case 82
                    %% ======================= Left Hemisphere =========================== %
                    % Subcortical Left
                    Nsubc = 7;
                    Ninterval = [1:Nsubc];
                    param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h = patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Frontal Left
                    Nfront = 13;
                    Ninterval = [Nsubc+1:Nfront+Nsubc];
                    param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h = patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Parietal Left
                    Nparie = 7;
                    Ninterval = [Nsubc+Nfront+1:Nsubc+Nfront+Nparie];
                    param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h = patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Temporal Left
                    Ntempo = 9;
                    Ninterval = [Nsubc+Nfront+Nparie+1:Nsubc+Nfront+Nparie+Ntempo];
                    param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h = patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Occipital Left
                    Noccip = 4;
                    Ninterval = [Nsubc+Nfront+Nparie+Ntempo+1:Nsubc+Nfront+Nparie+Ntempo+Noccip];
                    param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h = patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Insula Left
                    Nins = 1;
                    Ninterval = [Nsubc+Nfront+Nparie+Ntempo+Noccip+1:Nsubc+Nfront+Nparie+Ntempo+Noccip+Nins];
                    param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h = patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    %% ====================== Right Hemisphere =========================== %
                    % Subcortical Right
                    Nsubc = 7;
                    Ninterval = [1:Nsubc]+ Nstruct/2;
                    param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h = patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Frontal Right
                    Nfront = 13;
                    Ninterval = [Nsubc+1:Nfront+Nsubc]+ Nstruct/2;
                    param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h = patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Parietal Right
                    Nparie = 7;
                    Ninterval = [Nsubc+Nfront+1:Nsubc+Nfront+Nparie]+ Nstruct/2;
                    param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h = patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Temporal Right
                    Ntempo = 9;
                    Ninterval = [Nsubc+Nfront+Nparie+1:Nsubc+Nfront+Nparie+Ntempo]+ Nstruct/2;
                    param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h = patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Occipital Right
                    Noccip = 4;
                    Ninterval = [Nsubc+Nfront+Nparie+Ntempo+1:Nsubc+Nfront+Nparie+Ntempo+Noccip]+ Nstruct/2;
                    param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h = patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                    
                    % Insula Right
                    Nins = 1;
                    Ninterval = [Nsubc+Nfront+Nparie+Ntempo+Noccip+1:Nsubc+Nfront+Nparie+Ntempo+Noccip+Nins] + Nstruct/2;
                    param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
                    [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
                    h = patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
                otherwise
                    param.tht   =[(theta(1)-ang/2)*180/pi (theta(end)+ang/2)*180/pi];
                    [X,Y]=gen_xy(data(:,i),param);
                    CDat = zeros(param.N*2,size(C,1),3);
                    CDat(:,:,1)=repmat(colorimg(:,1)',[param.N*2 1]);
                    CDat(:,:,2)=repmat(colorimg(:,2)',[param.N*2 1]);
                    CDat(:,:,3)=repmat(colorimg(:,3)',[param.N*2 1]);
                    h = patch(X,Y,repmat(lin(data(:,i)')',[param.N*2 1])*0,CDat);
            end
        end
        
    end
end
%% ==================== End of Drawing Data rings ======================= %

%% ========================== Creating Concatenated Colormaps =========== %
if size(ColorTotal,1)>1
    ColorTotal(1,:) = [];
    colormap(ColorTotal);
    h = colorbar;
    Nticks = 2*size(data,2) + 1; % Number of Ticks
    StartTicks = linspace(0,1,Nticks);
    Temp = [1:2:length(StartTicks)];Temp(1) = [];Temp(end) = [];
    Temp2 = [StartTicks(Temp) + param.textspace StartTicks(Temp) - param.textspace];StartTicks(Temp) = [];
    Ticks = sort([StartTicks Temp2]);
    set(h, 'YTick', Ticks);
    lab = [min(data);mean(data);max(data)];
    labels = num2str(lab(:));
    set(h, 'YTickLabel', labels);
    set(get(h,'title'),'string','Multiple Colorbars');
    
end
%% ==================== End of Creating Colormap ======================== %

%% ================= Drawing Rings for Structures ======================= %
param.rho = [rmaxLim + param.intringwidth*5 rmaxLim + param.intringwidth*5 + 0.4];
rmaxLim = rmaxLim+param.intringwidth+0.4;
switch Nstruct
    case 66
        
        %% ======================= Left Hemisphere =========================== %
        % Frontal Left
        Nfront = 13;
        Ninterval = [1:Nfront];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Parietal Left
        Nparie = 7;
        Ninterval = [Nfront+1:Nfront+Nparie];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Temporal Left
        Ntempo = 9;
        Ninterval = [Nfront+Nparie+1:Nfront+Nparie+Ntempo];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Occipital Left
        Nocc = 4;
        Ninterval = [Nfront+Nparie+Ntempo+1:Nfront+Nparie+Ntempo+Nocc];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        
        %% ====================== Right Hemisphere =========================== %
        % Frontal Right
        Nfront = 13;
        Ninterval = [1:Nfront] + Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Parietal Right
        Nparie = 7;
        Ninterval = [Nfront+1:Nfront+Nparie] + Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Temporal Right
        Ntempo = 9;
        Ninterval = [Nfront+Nparie+1:Nfront+Nparie+Ntempo] + Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Occipital Right
        Nocc = 4;
        Ninterval = [Nfront+Nparie+Ntempo+1:Nfront+Nparie+Ntempo+Nocc] + Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
    case 68
        
        %% ======================= Left Hemisphere =========================== %
        % Frontal Left
        Nfront = 13;
        Ninterval = [1:Nfront];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Parietal Left
        Nparie = 7;
        Ninterval = [Nfront+1:Nfront+Nparie];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Temporal Left
        Ntempo = 9;
        Ninterval = [Nfront+Nparie+1:Nfront+Nparie+Ntempo];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Occipital Left
        Nocc = 4;
        Ninterval = [Nfront+Nparie+Ntempo+1:Nfront+Nparie+Ntempo+Nocc];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Insula Left
        Nins = 1;
        Ninterval = [Nfront+Nparie+Ntempo+Nocc+1:Nfront+Nparie+Ntempo+Nocc+Nins];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        %% ====================== Right Hemisphere =========================== %
        % Frontal Right
        Nfront = 13;
        Ninterval = [1:Nfront] + Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Parietal Right
        Nparie = 7;
        Ninterval = [Nfront+1:Nfront+Nparie] + Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Temporal Right
        Ntempo = 9;
        Ninterval = [Nfront+Nparie+1:Nfront+Nparie+Ntempo] + Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Occipital Right
        Nocc = 4;
        Ninterval = [Nfront+Nparie+Ntempo+1:Nfront+Nparie+Ntempo+Nocc] + Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Insula Right
        Nins = 1;
        Ninterval = [Nfront+Nparie+Ntempo+Nocc+1:Nfront+Nparie+Ntempo+Nocc+Nins] + Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
    case 82
        
        %% ======================= Left Hemisphere =========================== %
        % Subcortical Left
        Nsubc = 7;
        Ninterval = [1:Nsubc];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Frontal Left
        Nfront = 13;
        Ninterval = [Nsubc+1:Nfront+Nsubc];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Parietal Left
        Nparie = 7;
        Ninterval = [Nsubc+Nfront+1:Nsubc+Nfront+Nparie];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Temporal Left
        Ntempo = 9;
        Ninterval = [Nsubc+Nfront+Nparie+1:Nsubc+Nfront+Nparie+Ntempo];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Occipital Left
        Noccip = 4;
        Ninterval = [Nsubc+Nfront+Nparie+Ntempo+1:Nsubc+Nfront+Nparie+Ntempo+Noccip];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Insula Left
        Nins = 1;
        Ninterval = [Nsubc+Nfront+Nparie+Ntempo+Noccip+1:Nsubc+Nfront+Nparie+Ntempo+Noccip+Nins];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        %% ====================== Right Hemisphere =========================== %
        % Subcortical Right
        Nsubc = 7;
        Ninterval = [1:Nsubc]+ Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Frontal Right
        Nfront = 13;
        Ninterval = [Nsubc+1:Nfront+Nsubc]+ Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Parietal Right
        Nparie = 7;
        Ninterval = [Nsubc+Nfront+1:Nsubc+Nfront+Nparie]+ Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Temporal Right
        Ntempo = 9;
        Ninterval = [Nsubc+Nfront+Nparie+1:Nsubc+Nfront+Nparie+Ntempo]+ Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Occipital Right
        Noccip = 4;
        Ninterval = [Nsubc+Nfront+Nparie+Ntempo+1:Nsubc+Nfront+Nparie+Ntempo+Noccip]+ Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Insula Right
        Nins = 1;
        Ninterval = [Nsubc+Nfront+Nparie+Ntempo+Noccip+1:Nsubc+Nfront+Nparie+Ntempo+Noccip+Nins]+ Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, ColorsSt, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
    otherwise
        % Unknown parcellation
        param.tht   =[(theta(1)-ang/2)*180/pi (theta(end)+ang/2)*180/pi];
        col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
        Ncolor = size(col,1);
        Nstruct = size(C,1);
        re = floor(size(C,1)/Ncolor);
        Col = repmat(col,[re+1 1]);
        Col = Col(1:Nstruct,:);
        [X,Y]=gen_xy([1:Nstruct]',param);
        CDat = zeros(param.N*2,Nstruct,3);
        CDat(:,:,1)=repmat(Col(:,1)',[param.N*2 1]);
        CDat(:,:,2)=repmat(Col(:,2)',[param.N*2 1]);
        CDat(:,:,3)=repmat(Col(:,3)',[param.N*2 1]);
        h = patch(X,Y,squeeze(CDat(:,:,1))*0,CDat);
end
%% ============== End of Drawing Rings for Structures ====================%
if Nstruct == 82|Nstruct == 68|Nstruct == 66
    txt = cellstr(AcronSt);
else
    txt = cellstr(num2str((1:length(xy))','%02d'));
end


% PLACE TEXT LABELS such that you always read 'left to right'
xy_text = rmaxLim.* [cos(theta) sin(theta)];
ipos       = xy_text(:,1) > 0;
s.t        = zeros(Nstruct,1);
s.t( ipos) = text(xy_text(ipos,1).*1.05, xy_text(ipos,2).*1.05, txt(ipos), 'FontSize',10,'FontName','Helvetica');
set(s.t( ipos), {'Rotation'},num2cell(theta(ipos)*180/pi));
s.t( ~ipos) = text(xy_text(~ipos,1).*1.05, xy_text(~ipos,2).*1.05, txt(~ipos), 'FontSize',10,'FontName','Helvetica');
set(s.t( ~ipos), {'Rotation'},num2cell(theta(~ipos)*180/pi-180),'Horiz','right');

% ADJUST FIGURE height width to fit text labels
xtn        = cell2mat(get(s.t,'extent'));
post       = cell2mat(get(s.t,'pos'));
sg         = sign(post(:,2));
posfa      = cell2mat(get([gcf gca],'pos'));
% Calculate xlim and ylim in data units as x (y) position + extension along x (y)
ylims      = post(:,2) + xtn(:,4).*sg;
ylims      = [min(ylims), max(ylims)];
xlims      = post(:,1) + xtn(:,3).*sg;
xlims      = [min(xlims), max(xlims)];
% Stretch figure
posfa(1,3) = (( diff(xlims)/2 - 1)*posfa(2,3) + 1) * posfa(1,3);
posfa(1,4) = (( diff(ylims)/2 - 1)*posfa(2,4) + 1) * posfa(1,4);
% Position it a bit lower (movegui slow)
posfa(1,2) = 100;

% Axis settings
set(gca, 'Xlim',xlims + [-1.2 1.2],'Ylim',ylims + [-1.6 1.2], 'color', 'k', 'layer','bottom', 'Xtick',[],'Ytick',[])
%set(gcf, 'pos' ,posfa(1,:),'Visible','on')
axis equal

%% ====================== Drawing Rings for Lobes ======================= %
param.rho = [rmaxLim + 2.1 rmaxLim + 2.6];
rmaxLim = rmaxLim+param.intringwidth+0.4;
switch Nstruct
    case 66
        
        %% ======================= Left Hemisphere =========================== %
        % Frontal Left
        Nfront = 13;
        Ninterval = [1:Nfront];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [12 73 153], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Parietal Left
        Nparie = 7;
        Ninterval = [Nfront+1:Nfront+Nparie];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [7 96 32], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Temporal Left
        Ntempo = 9;
        Ninterval = [Nfront+Nparie+1:Nfront+Nparie+Ntempo];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [210 113 26], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Occipital Left
        Nocc = 4;
        Ninterval = [Nfront+Nparie+Ntempo+1:Nfront+Nparie+Ntempo+Nocc];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [126 18 33], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);

        
        %% ====================== Right Hemisphere =========================== %
        % Frontal Right
        Nfront = 13;
        Ninterval = [1:Nfront] + Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [12 73 153], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Parietal Right
        Nparie = 7;
        Ninterval = [Nfront+1:Nfront+Nparie] + Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [7 96 32], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Temporal Right
        Ntempo = 9;
        Ninterval = [Nfront+Nparie+1:Nfront+Nparie+Ntempo] + Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [210 113 26], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Occipital Right
        Nocc = 4;
        Ninterval = [Nfront+Nparie+Ntempo+1:Nfront+Nparie+Ntempo+Nocc] + Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [126 18 33], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);

        
    case 68
        
        %% ======================= Left Hemisphere =========================== %
        % Frontal Left
        Nfront = 13;
        Ninterval = [1:Nfront];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [12 73 153], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Parietal Left
        Nparie = 7;
        Ninterval = [Nfront+1:Nfront+Nparie];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [7 96 32], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Temporal Left
        Ntempo = 9;
        Ninterval = [Nfront+Nparie+1:Nfront+Nparie+Ntempo];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [210 113 26], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Occipital Left
        Nocc = 4;
        Ninterval = [Nfront+Nparie+Ntempo+1:Nfront+Nparie+Ntempo+Nocc];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [126 18 33], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
%         % Insula Left
%         Nins = 1;
%         Ninterval = [Nfront+Nparie+Ntempo+Nocc+1:Nfront+Nparie+Ntempo+Nocc+Nins];
%         param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
%         [X,Y,CDat] = Create_Ring(Ninterval, [255 192  32],  param);
%         h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        %% ====================== Right Hemisphere =========================== %
        % Frontal Right
        Nfront = 13;
        Ninterval = [1:Nfront] + Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [12 73 153], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Parietal Right
        Nparie = 7;
        Ninterval = [Nfront+1:Nfront+Nparie] + Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [7 96 32], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Temporal Right
        Ntempo = 9;
        Ninterval = [Nfront+Nparie+1:Nfront+Nparie+Ntempo] + Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [210 113 26], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Occipital Right
        Nocc = 4;
        Ninterval = [Nfront+Nparie+Ntempo+1:Nfront+Nparie+Ntempo+Nocc] + Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [126 18 33], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
%         % Insula Right
%         Nins = 1;
%         Ninterval = [Nfront+Nparie+Ntempo+Nocc+1:Nfront+Nparie+Ntempo+Nocc+Nins] + Nstruct/2;
%         param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
%         [X,Y,CDat] = Create_Ring(Ninterval,[255 192  32], param);
%         h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
    case 82
        
        %% ======================= Left Hemisphere =========================== %
        % Subcortical Left
        Nsubc = 7;
        Ninterval = [1:Nsubc];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [187 205 31], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Frontal Left
        Nfront = 13;
        Ninterval = [Nsubc+1:Nfront+Nsubc];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [12 73 153], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Parietal Left
        Nparie = 7;
        Ninterval = [Nsubc+Nfront+1:Nsubc+Nfront+Nparie];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [7 96 32], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Temporal Left
        Ntempo = 9;
        Ninterval = [Nsubc+Nfront+Nparie+1:Nsubc+Nfront+Nparie+Ntempo];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [210 113 26], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Occipital Left
        Noccip = 4;
        Ninterval = [Nsubc+Nfront+Nparie+Ntempo+1:Nsubc+Nfront+Nparie+Ntempo+Noccip];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [126 18 33], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
%         % Insula Left
%         Nins = 1;
%         Ninterval = [Nsubc+Nfront+Nparie+Ntempo+Noccip+1:Nsubc+Nfront+Nparie+Ntempo+Noccip+Nins];
%         param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
%         [X,Y,CDat] = Create_Ring(Ninterval,[255 192  32], param);
%         h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        %% ====================== Right Hemisphere =========================== %
        % Subcortical Right
        Nsubc = 7;
        Ninterval = [1:Nsubc]+ Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [187 205 31], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Frontal Right
        Nfront = 13;
        Ninterval = [Nsubc+1:Nfront+Nsubc]+ Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [12 73 153], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Parietal Right
        Nparie = 7;
        Ninterval = [Nsubc+Nfront+1:Nsubc+Nfront+Nparie]+ Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [7 96 32], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Temporal Right
        Ntempo = 9;
        Ninterval = [Nsubc+Nfront+Nparie+1:Nsubc+Nfront+Nparie+Ntempo]+ Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [210 113 26], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
        % Occipital Right
        Noccip = 4;
        Ninterval = [Nsubc+Nfront+Nparie+Ntempo+1:Nsubc+Nfront+Nparie+Ntempo+Noccip]+ Nstruct/2;
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(1, [126 18 33], param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
        
%         % Insula Right
%         Nins = 1;
%         Ninterval = [Nsubc+Nfront+Nparie+Ntempo+Noccip+1:Nsubc+Nfront+Nparie+Ntempo+Noccip+Nins]+ Nstruct/2;
%         param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
%         [X,Y,CDat] = Create_Ring(Ninterval,[255 192  32], param);
%         h=patch(X,Y,CDat(:,:,1)*0,CDat/255);
end
%% ================ End ofDrawing Rings for Lobes ======================= %
if param.newfigure == 1
    camlight
end
return;

function [Names,Acronyms,Colors,Codes] = Load_Struct_Names(atype,nRegions);

% atype           : Atlas type
% nRegions        : Number of regions in the connectivity matrix. Three
%                   possible organizations are suported for aparc+aseg atlas.
%                   nRegions=66: Include just cortical regions
%                   nRegions=68: Include cortical regions + insula
%                   nRegions=82: Include cortical regions + insula +
%                   subcortical structures
switch atype
    case 'aparc+aseg'
        %% =============================== Left Hemisphere ===================================== %
        
        % ------ Subcortical Structures
        Struct.LHSubCort.Names = strvcat('Left-Thalamus-Proper',...
            'Left-Caudate',...
            'Left-Putamen',...
            'Left-Pallidum',...
            'Left-Hippocampus',...
            'Left-Amygdala',...
            'Left-Accumbens-area');
        Struct.LHSubCort.Colors = [0   118    14;...
            122   186   220;...
            236    13   176;...
            12    48   255;...
            220   216    20;...
            103   255   255;...
            255   165     0];
        Struct.LHSubCort.Acron = strvcat('LH-Thalam',...
            'LH-Caudt',...
            'LH-Putam',...
            'LH-Pallid',...
            'LH-Hippoc',...
            'LH-Amygd',...
            'LH-AccumbArea');
        Struct.LHSubCort.Codes = [10:13 17:18 26];
        
        % ------ Frontal Lobe
        Struct.LHFrontalCortex.Names = strvcat('ctx-lh-caudalanteriorcingulate',...
            'ctx-lh-caudalmiddlefrontal',...
            'ctx-lh-lateralorbitofrontal',...
            'ctx-lh-medialorbitofrontal',...
            'ctx-lh-paracentral',...
            'ctx-lh-parsopercularis',...
            'ctx-lh-parsorbitalis',...
            'ctx-lh-parstriangularis',...
            'ctx-lh-precentral',...
            'ctx-lh-rostralanteriorcingulate',...
            'ctx-lh-rostralmiddlefrontal',...
            'ctx-lh-superiorfrontal',...
            'ctx-lh-frontalpole');
        Struct.LHFrontalCortex.Colors = [125   100   160;...
            100    25     0;...
            35    75    50;...
            200    35    75;...
            60   220    60;...
            220   180   140;...
            20   100    50;...
            220    60    20;...
            60    20   220;...
            80    20   140;...
            75    50   125;...
            20   220   160;...
            100     0   100];
        Struct.LHFrontalCortex.Acron = strvcat('LH-CaudAntCing',...
            'LH-CaudMidFront',...
            'LH-LatOrbFront',...
            'LH-MedOrbFront',...
            'LH-Paracent',...
            'LH-ParsoPerc',...
            'LH-ParsOrb',...
            'LH-ParsTriang',...
            'LH-Precent',...
            'LH-RostAntCing',...
            'LH-RostMidFront',...
            'LH-SupFront',...
            'LH-FrontPole');
        Struct.LHFrontalCortex.Codes = [1002 1003 1012 1014 1017 1018 1019 1020 1024 1026 1027 1028 1032];
        
        % ------ Parietal Lobe
        Struct.LHParietalCortex.Names = strvcat('ctx-lh-inferiorparietal',...
            'ctx-lh-isthmuscingulate',...
            'ctx-lh-postcentral',...
            'ctx-lh-posteriorcingulate',...
            'ctx-lh-precuneus',...
            'ctx-lh-superiorparietal',...
            'ctx-lh-supramarginal');
        Struct.LHParietalCortex.Colors = [220    60   220;...
            140    20   140;...
            220    20    20;...
            220   180   220;...
            160   140   180;...
            20   180   140;...
            80   160    20];
        Struct.LHParietalCortex.Acron = strvcat('LH-InfPariet',...
            'LH-IsthCing',...
            'LH-PostCent',...
            'LH-PostCing',...
            'LH-Precuneus',...
            'LH-SupPariet',...
            'LH-Supramarg');
        Struct.LHParietalCortex.Codes = [1008 1010 1022 1023 1025 1029 1031];
        
        % ------ Temporal Lobe
        Struct.LHTemporalCortex.Names = strvcat('ctx-lh-bankssts',...
            'ctx-lh-entorhinal',...
            'ctx-lh-fusiform',...
            'ctx-lh-inferiortemporal',...
            'ctx-lh-middletemporal',...
            'ctx-lh-parahippocampal',...
            'ctx-lh-superiortemporal',...
            'ctx-lh-temporalpole',...
            'ctx-lh-transversetemporal');
        Struct.LHTemporalCortex.Colors = [25   100    40;...
            220    20    10;...
            180   220   140;...
            180    40   120;...
            160   100    50;...
            20   220    60;...
            140   220   220;...
            70    70    70;...
            150   150   200];
        Struct.LHTemporalCortex.Acron = strvcat('LH-BanksStS',...
            'LH-Entorh',...
            'LH-Fusif',...
            'LH-InfTemp',...
            'LH-MidTemp',...
            'LH-Parahip',...
            'LH-SupTemp',...
            'LH-TempPole',...
            'LH-TransvTemp');
        Struct.LHTemporalCortex.Codes = [1001 1006 1007 1009 1015 1016 1030 1033 1034];
        
        % ------ Occipital Lobe
        Struct.LHOccipitalCortex.Names = strvcat('ctx-lh-cuneus',...
            'ctx-lh-lateraloccipital',...
            'ctx-lh-lingual',...
            'ctx-lh-pericalcarine');
        Struct.LHOccipitalCortex.Colors = [220    20   100;...
            20    30   140;...
            225   140   140;...
            120   100    60];
        Struct.LHOccipitalCortex.Acron = strvcat('LH-Cuneus',...
            'LH-LatOccip',...
            'LH-Lingual',...
            'LH-PeriCalc');
        Struct.LHOccipitalCortex.Codes = [1001 1006 1007 1009 1015 1016 1030 1033 1034];
        
        % ------ Insula
        Struct.LHInsula.Names = 'ctx-lh-insula';
        Struct.LHInsula.Colors =  [255   192    32];
        Struct.LHInsula.Acron = 'LH-Insula';
        Struct.LHInsula.Codes = 1035;
        %% =========================== End of Left Hemisphere ================================== %
        
        
        %% =============================== Right Hemisphere ===================================== %
        
        % ------ Subcortical Structures
        Struct.RHSubCort.Names = strvcat('Right-Thalamus-Proper',...
            'Right-Caudate',...
            'Right-Putamen',...
            'Right-Pallidum',...
            'Right-Hippocampus',...
            'Right-Amygdala',...
            'Right-Accumbens-area');
        Struct.RHSubCort.Colors = [0   118    14;...
            122   186   220;...
            236    13   176;...
            12    48   255;...
            220   216    20;...
            103   255   255;...
            255   165     0];
        Struct.RHSubCort.Acron = strvcat('RH-Thalam',...
            'RH-Caudt',...
            'RH-Putam',...
            'RH-Pallid',...
            'RH-Hippoc',...
            'RH-Amygd',...
            'RH-AccumbArea');
        Struct.RHSubCort.Codes = [ 49:54 58];
        
        % ------ Frontal Lobe
        Struct.RHFrontalCortex.Names = strvcat('ctx-rh-caudalanteriorcingulate',...
            'ctx-rh-caudalmiddlefrontal',...
            'ctx-rh-lateralorbitofrontal',...
            'ctx-rh-medialorbitofrontal',...
            'ctx-rh-paracentral',...
            'ctx-rh-parsopercularis',...
            'ctx-rh-parsorbitalis',...
            'ctx-rh-parstriangularis',...
            'ctx-rh-precentral',...
            'ctx-rh-rostralanteriorcingulate',...
            'ctx-rh-rostralmiddlefrontal',...
            'ctx-rh-superiorfrontal',...
            'ctx-rh-frontalpole');
        Struct.RHFrontalCortex.Colors = [125   100   160;...
            100    25     0;...
            35    75    50;...
            200    35    75;...
            60   220    60;...
            220   180   140;...
            20   100    50;...
            220    60    20;...
            60    20   220;...
            80    20   140;...
            75    50   125;...
            20   220   160;...
            100     0   100];
        Struct.RHFrontalCortex.Acron = strvcat('RH-CaudAntCing',...
            'RH-CaudMidFront',...
            'RH-LatOrbFront',...
            'RH-MedOrbFront',...
            'RH-Paracent',...
            'RH-ParsoPerc',...
            'RH-ParsOrb',...
            'RH-ParsTriang',...
            'RH-Precent',...
            'RH-RostAntCing',...
            'RH-RostMidFront',...
            'RH-SupFront',...
            'RH-FrontPole');
        Struct.RHFrontalCortex.Codes = [2002 2003 2012 2014 2017 2018 2019 2020 2024 2026 2027 2028 2032];
        
        % ------ Parietal Lobe
        Struct.RHParietalCortex.Names = strvcat('ctx-rh-inferiorparietal',...
            'ctx-rh-isthmuscingulate',...
            'ctx-rh-postcentral',...
            'ctx-rh-posteriorcingulate',...
            'ctx-rh-precuneus',...
            'ctx-rh-superiorparietal',...
            'ctx-rh-supramarginal');
        Struct.RHParietalCortex.Colors = [220    60   220;...
            140    20   140;...
            220    20    20;...
            220   180   220;...
            160   140   180;...
            20   180   140;...
            80   160    20];
        Struct.RHParietalCortex.Acron = strvcat('RH-InfPariet',...
            'RH-IsthCing',...
            'RH-PostCent',...
            'RH-PostCing',...
            'RH-Precuneus',...
            'RH-SupPariet',...
            'RH-Supramarg');
        Struct.RHParietalCortex.Codes = [2008 2010 2022 2023 2025 2029 2031];
        
        % ------ Temporal Lobe
        Struct.RHTemporalCortex.Names = strvcat('ctx-rh-bankssts',...
            'ctx-rh-entorhinal',...
            'ctx-rh-fusiform',...
            'ctx-rh-inferiortemporal',...
            'ctx-rh-middletemporal',...
            'ctx-rh-parahippocampal',...
            'ctx-rh-superiortemporal',...
            'ctx-rh-temporalpole',...
            'ctx-rh-transversetemporal');
        Struct.RHTemporalCortex.Colors = [25   100    40;...
            220    20    10;...
            180   220   140;...
            180    40   120;...
            160   100    50;...
            20   220    60;...
            140   220   220;...
            70    70    70;...
            150   150   200];
        Struct.RHTemporalCortex.Acron = strvcat('RH-BanksStS',...
            'RH-Entorh',...
            'RH-Fusif',...
            'RH-InfTemp',...
            'RH-MidTemp',...
            'RH-Parahip',...
            'RH-SupTemp',...
            'RH-TempPole',...
            'RH-TransvTemp');
        Struct.RHTemporalCortex.Codes = [2001 2006 2007 2009 2015 2016 2030 2033 2034];
        
        % ------ Occipital Lobe
        Struct.RHOccipitalCortex.Names = strvcat('ctx-rh-cuneus',...
            'ctx-rh-lateraloccipital',...
            'ctx-rh-lingual',...
            'ctx-rh-pericalcarine');
        Struct.RHOccipitalCortex.Colors = [220    20   100;...
            20    30   140;...
            225   140   140;...
            120   100    60];
        Struct.RHOccipitalCortex.Acron = strvcat('RH-Cuneus',...
            'RH-LatOccip',...
            'RH-Lingual',...
            'RH-PeriCalc');
        Struct.RHOccipitalCortex.Codes = [2001 2006 2007 2009 2015 2016 2030 2033 2034];
        
        % ------ Insula
        Struct.RHInsula.Names = 'ctx-rh-insula';
        Struct.RHInsula.Colors =  [255   192    32];
        Struct.RHInsula.Acron = 'RH-Insula';
        Struct.RHInsula.Codes = 2035;
        
        %% =========================== End of Right Hemisphere ================================== %
        if nRegions == 66
            Names = strvcat(Struct.LHFrontalCortex.Names,Struct.LHParietalCortex.Names,Struct.LHTemporalCortex.Names,Struct.LHOccipitalCortex.Names,...
                Struct.RHFrontalCortex.Names,Struct.RHParietalCortex.Names,Struct.RHTemporalCortex.Names,Struct.RHOccipitalCortex.Names);
            
            Acronyms = strvcat(Struct.LHFrontalCortex.Acron,Struct.LHParietalCortex.Acron,Struct.LHTemporalCortex.Acron,Struct.LHOccipitalCortex.Acron,...
                Struct.RHFrontalCortex.Acron,Struct.RHParietalCortex.Acron,Struct.RHTemporalCortex.Acron,Struct.RHOccipitalCortex.Acron);
            
            Colors   =        [Struct.LHFrontalCortex.Colors;Struct.LHParietalCortex.Colors;Struct.LHTemporalCortex.Colors;Struct.LHOccipitalCortex.Colors;...
                Struct.RHFrontalCortex.Colors;Struct.RHParietalCortex.Colors;Struct.RHTemporalCortex.Colors;Struct.RHOccipitalCortex.Colors];
            
            Codes    =        [Struct.LHFrontalCortex.Codes Struct.LHParietalCortex.Codes Struct.LHTemporalCortex.Codes Struct.LHOccipitalCortex.Codes...
                Struct.RHFrontalCortex.Codes Struct.RHParietalCortex.Codes Struct.RHTemporalCortex.Codes Struct.RHOccipitalCortex.Codes]';
            
        elseif nRegions == 68
            Names = strvcat(Struct.LHFrontalCortex.Names,Struct.LHParietalCortex.Names,Struct.LHTemporalCortex.Names,Struct.LHOccipitalCortex.Names,Struct.LHInsula.Names,...
                Struct.RHFrontalCortex.Names,Struct.RHParietalCortex.Names,Struct.RHTemporalCortex.Names,Struct.RHOccipitalCortex.Names,Struct.RHInsula.Names);
            
            Acronyms = strvcat(Struct.LHFrontalCortex.Acron,Struct.LHParietalCortex.Acron,Struct.LHTemporalCortex.Acron,Struct.LHOccipitalCortex.Acron,Struct.LHInsula.Acron,...
                Struct.RHFrontalCortex.Acron,Struct.RHParietalCortex.Acron,Struct.RHTemporalCortex.Acron,Struct.RHOccipitalCortex.Acron,Struct.RHInsula.Acron);
            
            Colors   =        [Struct.LHFrontalCortex.Colors;Struct.LHParietalCortex.Colors;Struct.LHTemporalCortex.Colors;Struct.LHOccipitalCortex.Colors;Struct.LHInsula.Colors;...
                Struct.RHFrontalCortex.Colors;Struct.RHParietalCortex.Colors;Struct.RHTemporalCortex.Colors;Struct.RHOccipitalCortex.Colors;Struct.RHInsula.Colors];
            
            Codes    =        [Struct.LHFrontalCortex.Codes Struct.LHParietalCortex.Codes Struct.LHTemporalCortex.Codes Struct.LHOccipitalCortex.Codes Struct.LHInsula.Codes...
                Struct.RHFrontalCortex.Codes Struct.RHParietalCortex.Codes Struct.RHTemporalCortex.Codes Struct.RHOccipitalCortex.Codes Struct.RHInsula.Codes]';
            
        elseif nRegions == 82
            Names = strvcat(Struct.LHSubCort.Names,Struct.LHFrontalCortex.Names,Struct.LHParietalCortex.Names,Struct.LHTemporalCortex.Names,Struct.LHOccipitalCortex.Names,Struct.LHInsula.Names,...
                Struct.RHSubCort.Names,Struct.RHFrontalCortex.Names,Struct.RHParietalCortex.Names,Struct.RHTemporalCortex.Names,Struct.RHOccipitalCortex.Names,Struct.RHInsula.Names);
            
            Acronyms = strvcat(Struct.LHSubCort.Acron,Struct.LHFrontalCortex.Acron,Struct.LHParietalCortex.Acron,Struct.LHTemporalCortex.Acron,Struct.LHOccipitalCortex.Acron,Struct.LHInsula.Acron,...
                Struct.RHSubCort.Acron,Struct.RHFrontalCortex.Acron,Struct.RHParietalCortex.Acron,Struct.RHTemporalCortex.Acron,Struct.RHOccipitalCortex.Acron,Struct.RHInsula.Acron);
            
            Colors   =        [Struct.LHSubCort.Colors;Struct.LHFrontalCortex.Colors;Struct.LHParietalCortex.Colors;Struct.LHTemporalCortex.Colors;Struct.LHOccipitalCortex.Colors;Struct.LHInsula.Colors;...
                Struct.RHSubCort.Colors;Struct.RHFrontalCortex.Colors;Struct.RHParietalCortex.Colors;Struct.RHTemporalCortex.Colors;Struct.RHOccipitalCortex.Colors;Struct.RHInsula.Colors];
            
            Codes    =        [Struct.LHSubCort.Codes Struct.LHFrontalCortex.Codes Struct.LHParietalCortex.Codes Struct.LHTemporalCortex.Codes Struct.LHOccipitalCortex.Codes Struct.LHInsula.Codes...
                Struct.RHSubCort.Codes Struct.RHFrontalCortex.Codes Struct.RHParietalCortex.Codes Struct.RHTemporalCortex.Codes Struct.RHOccipitalCortex.Codes Struct.RHInsula.Codes]';
        end
        
end

return


%% Bullseye internal functions

function [X,Y]=gen_xy(data,param);
ind=1;
for j=1:size(data,1)
    for k=1:size(data,2)
        dtht=diff(param.tht)/size(data,1);
        tht1=dtht*(j-1)+param.tht(1)-param.tht0;
        tht2=dtht* j   +param.tht(1)-param.tht0;
        
        dr=diff(param.rho)/size(data,2);
        rho1=dr*(k-1)+param.rho(1);
        rho2=dr* k   +param.rho(1);
        
        ang = linspace(tht1/180*pi,tht2/180*pi,param.N);
        [arc1.x,arc1.y]=pol2cart(ang,rho1);
        [arc2.x,arc2.y]=pol2cart(ang,rho2);
        X(:,ind)=[arc1.x arc2.x(end:-1:1)];
        Y(:,ind)=[arc1.y arc2.y(end:-1:1)];
        ind=ind+1;
    end
end
return

function t=bull_labels(param);

% Labels are added CCW...this alone determines that the first ROW is INFEROLATERAL etc...
% label{1}='inferoseptal';
% label{2}='inferior';
% label{3}='inferolateral';
% label{4}='anterolateral';
% label{5}='anterior';
% label{6}='anteroseptal';
label{1}='anterior';
label{2}='anteroseptal';
label{3}='inferoseptal';
label{4}='inferior';
label{5}='inferolateral';
label{6}='anterolateral';

for k=1:length(label)
    %   ang=(60*(k-1)-param.tht0-30)*(pi/180);
    ang=(60*(k-1)-param.tht0-60)*(pi/180);
    %   [x,y]=pol2cart(ang,1.08*max(param.rho));
    [x,y]=pol2cart(ang,0.875*max(param.rho));
    t(k)=text(x,y,label{k});
    if rem((ang*180/pi),360)>=180
        set(t(k),'rotation',ang*(180/pi)+90,'HorizontalAlignment','Center','FontSize',14,'FontName','Courier','FontWeight','Bold');
    elseif (ang*180/pi)==0 | (ang*180/pi)==360
        set(t(k),'rotation',ang*(180/pi)+90,'HorizontalAlignment','Center','FontSize',14,'FontName','Courier','FontWeight','Bold');
    else
        set(t(k),'rotation',ang*(180/pi)-90,'HorizontalAlignment','Center','FontSize',14,'FontName','Courier','FontWeight','Bold');
    end
end

return

%  This function linearizes a matrix (m) of any dimension (eg M=m(:)).
%  If an index vector (ind) is given then the the ind entries of m(:) are
%  returned.
%
% SYNTAX: m=lin(m);
%         m=lin(m,ind);
%
% DBE 12/22/03

function m=lin(m,ind);

m=m(:);

if nargin==2
    m=m(ind);
end

function [X,Y,CDat] = Create_Ring(Ninterval, Colors, param);
%
% Syntax :
% [X,Y,CDat] = Create_Ring(Ninterval,theta, Colors, param);
%
% This function creates the coordinates and color data for bullseye
% plotting.
%
% Input Parameters:
%   Ninterval      : Structures interval
%   theta          : Theta angles (radians)
%   Colors         : Colors
%   param          : Parameters
%
%
% Output Parameters:
%   X              : X coordinates
%   Y              : Y coordinates
%   CDat           : Color Data
%
% Related references:
%
%
% See also: Plot_Circular_Graph
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% November 12th 2013
% Version $1.0

Nlength = length(Ninterval);
CDat = zeros(param.N*2,Nlength,3);
[X,Y]=gen_xy(Ninterval',param);
CDat(:,:,1)=repmat(Colors(Ninterval,1)',[param.N*2 1]);
CDat(:,:,2)=repmat(Colors(Ninterval,2)',[param.N*2 1]);
CDat(:,:,3)=repmat(Colors(Ninterval,3)',[param.N*2 1]);
return


function NETWORKS = Label_Networks(conn, threshold);
%
% Syntax :
%  Networks = Label_Networks(conn, threshold);
%
% This script labels the unconnected networks in a connectivity matrix
% provided they are higher than a certain threshold
%
% Input Parameters:
%   conn              :  Connectivity Matrix
%   threshold         :  Threshold of the minimum intensity in the
%                        connectivity matrix
%
% Output Parameters:
%   NETWORKS          :  Matrix that contains a label in every element
%                        in conn matrix according to the network they belong
%
% Related references:
%
%
% See also: graphshortestpath
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% Jul 27th 2012
% Version $1.0

% ======================= Main Program ================================== %
% conn = triu(conn);
%  conn_t = conn > threshold;
%  dims = size(conn);
%  NETWORKS = zeros(dims);
%  indices = find(conn_t>0);
%  num = 1;
%
%  [ix,iy] = ind2sub(size(conn_t),indices);
%  index = [ix iy];
%  Id = ix*0;
%  while sum(logical(Id)) ~= length(Id)
%      idf = index(1,:);
%      while ~isempty(idf)
%          indcon = find(sum(ismember(index,idf)') > 0);
%          idf = index(indcon,:);
%          Id(indcon) = num;
%          index(indcon,:) = [];
%          idf = unique(idf(:));
%      end
%      index(find(Id == num),:) = [];
%      num = num +1;
%  end
%      indices = find(conn_t>0);
% %
% %     elems = [];
%
%     while(length(elems)>0)
%         [X, Y] = ind2sub(dims, elems(1));
%         elems(1)=[];
%         NETWORKS(X,Y)=num;
%         NETWORKS(Y,X)=num;
%         conn_t(X,Y) = 0;
%         conn_t(Y,X) = 0;
%         a = find(conn_t(Y,:)>0);
%         for i = 1:length(a)
%             if (a(i)~=X)
%                 elems = [elems sub2ind(dims,Y,a(i))];
%             end
%         end
%         a = find(conn_t(X,:)>0);
%         for i = 1:length(a)
%             if (a(i)~=Y)
%                 elems = [elems sub2ind(dims,a(i),X)];
%             end
%         end
%
%     end
%     num = num+1;
%
%  end




conn_t = conn > threshold;
dims = size(conn);
NETWORKS = zeros(dims);
indices = find(conn_t>0);
num = 1;


while(sum(conn_t(:))>0)
    
    indices = find(conn_t>0);
    
    elems = [indices(1)];
    
    while(length(elems)>0)
        [X, Y] = ind2sub(dims, elems(1));
        elems(1)=[];
        NETWORKS(X,Y)=num;
        NETWORKS(Y,X)=num;
        conn_t(X,Y) = 0;
        conn_t(Y,X) = 0;
        a = find(conn_t(Y,:)>0);
        for i = 1:length(a)
            if (a(i)~=X)
                elems = [elems sub2ind(dims,Y,a(i))];
            end
        end
        a = find(conn_t(X,:)>0);
        for i = 1:length(a)
            if (a(i)~=Y)
                elems = [elems sub2ind(dims,a(i),X)];
            end
        end
        
    end
    num = num+1;
    
end

function [labid] = Recur_Corr(Surf,lab,labid,cont);
%
% Syntax :
% [labid] = Recur_Corr(Surf,lab,labid,cont);
%
% Recursive function for labelling correction.
%
% See also: Red_Surf Smooth_Surf Plot_Surf Surf_Comp Plot_oversurf
% Exp_Surf
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% December 1st 2007
% Version $1.0

ind = find(Surf.Is ==lab&(labid ==0));
A = unique(Surf.SurfData.faces(Surf.Tri(ind(1),3:2+Surf.Tri(ind(1),2)),:));A = A(:);
indt = find(Surf.Is(A)~= lab);
A(indt) = [];A(A == ind(1))= [];
T = unique([ind(1); A]);
labid(T) = cont;
An = rand(size(A));
while sum(A)~=sum(An)
    An = A;
    Neib = Surf.Tri(A,3:end); Neib = unique(nonzeros(Neib(:)));
    A = unique(Surf.SurfData.faces(Neib,:));A = A(:);
    indt = find(Surf.Is(A)~= lab);
    A(indt) = [];
    labid(A) = cont;
    T =unique([T;A]);
end
indn = find(Surf.Is ==lab&(labid ==0));
if ~isempty(indn)
    cont = cont+1;
    [labid] = Recur_Corr(Surf,lab,labid,cont);
else
    return;
end
if ~isempty(A)
    for i = 1:size(A,2)
        TA = unique(Surf.SurfData.faces(Surf.Tri(A(i),3:2+Surf.Tri(A(i),2)),:));
        indt = find(Surf.Is(TA(i))~= lab);
        TA(indt) = [] ;
        T = [T; A];
    end
else
    return;
end
return

function AcronSt  = Names_2_Acron(atype,Names);
%
% Syntax :
% AcronSt  = Names_2_Acron(atype,Names);
%
% This function creates anocryms for structures name according with a
% specified atlas type.
%
% Input Parameters:
%   atype          : Atlas type
%   Names          : Structure Names
%
% Output Parameters:
%     AcronSt      : Acronym names
%
% Related references:
%
%
% See also: Plot_Circular_Graph FreeSurf_2_Circos
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% November 12th 2013
% Version $1.0

switch atype
    case 'aparc+aseg'
        AcronSt = '';
        for i = 1:size(Names,1);
            
            % Left Hemisphere
            if strcmp(deblank(Names(i,:)),'Left-Thalamus')%|strcmp(deblank(Names(i,:)),'Right-Thalamus')
                AcronSt = strvcat(AcronSt,'LH-Thalam');
            elseif strcmp(deblank(Names(i,:)),'Left-Thalamus-Proper')%|strcmp(deblank(Names(i,:)),'Right-Thalamus-Proper')
                AcronSt = strvcat(AcronSt,'LH-ThalProp');
            elseif strcmp(deblank(Names(i,:)),'Left-Caudate')%|strcmp(deblank(Names(i,:)),'Right-Caudate')
                AcronSt = strvcat(AcronSt,'LH-Caudt');
            elseif strcmp(deblank(Names(i,:)),'Left-Putamen')%|strcmp(deblank(Names(i,:)),'Right-Putamen')
                AcronSt = strvcat(AcronSt,'LH-Putam');
            elseif strcmp(deblank(Names(i,:)),'Left-Pallidum')%|strcmp(deblank(Names(i,:)),'Right-Pallidum')
                AcronSt = strvcat(AcronSt,'LH-Pallid');
            elseif strcmp(deblank(Names(i,:)),'Left-Hippocampus')%|strcmp(deblank(Names(i,:)),'Right-Hippocampus')
                AcronSt = strvcat(AcronSt,'LH-Hippoc');
            elseif strcmp(deblank(Names(i,:)),'Left-Amygdala')%|strcmp(deblank(Names(i,:)),'Right-Amygdala')
                AcronSt = strvcat(AcronSt,'LH-Amygd');
            elseif strcmp(deblank(Names(i,:)),'Left-Accumbens-area')%|strcmp(deblank(Names(i,:)),'Right-Accumbens-area')
                AcronSt = strvcat(AcronSt,'LH-AccumbArea');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-bankssts')%|strcmp(deblank(Names(i,:)),'ctx-rh-bankssts')
                AcronSt = strvcat(AcronSt,'LH-BanksStS');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-caudalanteriorcingulate')%|strcmp(deblank(Names(i,:)),'ctx-rh-caudalanteriorcingulate')
                AcronSt = strvcat(AcronSt,'LH-CaudAntCing');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-caudalmiddlefrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-caudalmiddlefrontal')
                AcronSt = strvcat(AcronSt,'LH-CaudMidFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-cuneus')%|strcmp(deblank(Names(i,:)),'ctx-rh-cuneus')
                AcronSt = strvcat(AcronSt,'LH-Cuneus');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-entorhinal')%|strcmp(deblank(Names(i,:)),'ctx-rh-entorhinal')
                AcronSt = strvcat(AcronSt,'LH-Entorh');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-fusiform')%|strcmp(deblank(Names(i,:)),'ctx-rh-fusiform')
                AcronSt = strvcat(AcronSt,'LH-Fusif');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-inferiorparietal')%|strcmp(deblank(Names(i,:)),'ctx-rh-inferiorparietal')
                AcronSt = strvcat(AcronSt,'LH-InfPariet');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-inferiortemporal')%|strcmp(deblank(Names(i,:)),'ctx-rh-inferiortemporal')
                AcronSt = strvcat(AcronSt,'LH-InfTemp');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-isthmuscingulate')%|strcmp(deblank(Names(i,:)),'ctx-rh-isthmuscingulate')
                AcronSt = strvcat(AcronSt,'LH-IsthCing');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-lateraloccipital')%|strcmp(deblank(Names(i,:)),'ctx-rh-lateraloccipital')
                AcronSt = strvcat(AcronSt,'LH-LatOccip');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-lateralorbitofrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-lateralorbitofrontal')
                AcronSt = strvcat(AcronSt,'LH-LatOrbFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-lingual')%|strcmp(deblank(Names(i,:)),'ctx-rh-lingual')
                AcronSt = strvcat(AcronSt,'LH-Lingual');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-medialorbitofrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-medialorbitofrontal')
                AcronSt = strvcat(AcronSt,'LH-MedOrbFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-middletemporal')%|strcmp(deblank(Names(i,:)),'ctx-rh-middletemporal')
                AcronSt = strvcat(AcronSt,'LH-MidTemp');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-parahippocampal')%|strcmp(deblank(Names(i,:)),'ctx-rh-parahippocampal')
                AcronSt = strvcat(AcronSt,'LH-Parahip');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-paracentral')%|strcmp(deblank(Names(i,:)),'ctx-rh-paracentral')
                AcronSt = strvcat(AcronSt,'LH-Paracent');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-parsopercularis')%|strcmp(deblank(Names(i,:)),'ctx-rh-parsopercularis')
                AcronSt = strvcat(AcronSt,'LH-ParsoPerc');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-parsorbitalis')%|strcmp(deblank(Names(i,:)),'ctx-rh-parsorbitalis')
                AcronSt = strvcat(AcronSt,'LH-ParsOrb');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-parstriangularis')%|strcmp(deblank(Names(i,:)),'ctx-rh-parstriangularis')
                AcronSt = strvcat(AcronSt,'LH-ParsTriang');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-pericalcarine')%|strcmp(deblank(Names(i,:)),'ctx-rh-pericalcarine')
                AcronSt = strvcat(AcronSt,'LH-PeriCalc');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-postcentral')%|strcmp(deblank(Names(i,:)),'ctx-rh-postcentral')
                AcronSt = strvcat(AcronSt,'LH-PostCent');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-posteriorcingulate')%|strcmp(deblank(Names(i,:)),'ctx-rh-posteriorcingulate')
                AcronSt = strvcat(AcronSt,'LH-PostCing');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-precentral')%|strcmp(deblank(Names(i,:)),'ctx-rh-precentral')
                AcronSt = strvcat(AcronSt,'LH-Precent');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-precuneus')%|strcmp(deblank(Names(i,:)),'ctx-rh-precuneus')
                AcronSt = strvcat(AcronSt,'LH-Precuneus');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-rostralanteriorcingulate')%|strcmp(deblank(Names(i,:)),'ctx-rh-rostralanteriorcingulate')
                AcronSt = strvcat(AcronSt,'LH-RostAntCing');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-rostralmiddlefrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-rostralmiddlefrontal')
                AcronSt = strvcat(AcronSt,'LH-RostMidFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-superiorfrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-superiorfrontal')
                AcronSt = strvcat(AcronSt,'LH-SupFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-superiorparietal')%|strcmp(deblank(Names(i,:)),'ctx-rh-superiorparietal')
                AcronSt = strvcat(AcronSt,'LH-SupPariet');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-superiortemporal')%|strcmp(deblank(Names(i,:)),'ctx-rh-superiortemporal')
                AcronSt = strvcat(AcronSt,'LH-SupTemp');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-supramarginal')%|strcmp(deblank(Names(i,:)),'ctx-rh-supramarginal')
                AcronSt = strvcat(AcronSt,'LH-Supramarg');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-frontalpole')%|strcmp(deblank(Names(i,:)),'ctx-rh-frontalpole')
                AcronSt = strvcat(AcronSt,'LH-FrontPole');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-temporalpole')%|strcmp(deblank(Names(i,:)),'ctx-rh-temporalpole')
                AcronSt = strvcat(AcronSt,'LH-TempPole');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-transversetemporal')%|strcmp(deblank(Names(i,:)),'ctx-rh-transversetemporal')
                AcronSt = strvcat(AcronSt,'LH-TransvTemp');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-insula')%|strcmp(deblank(Names(i,:)),'ctx-rh-insula')
                AcronSt = strvcat(AcronSt,'LH-Insula');
                
                % Right Hemisphere
                
            elseif strcmp(deblank(Names(i,:)),'Right-Thalamus')%|strcmp(deblank(Names(i,:)),'Right-Thalamus')
                AcronSt = strvcat(AcronSt,'RH-Thalam');
            elseif strcmp(deblank(Names(i,:)),'Right-Thalamus-Proper')%|strcmp(deblank(Names(i,:)),'Right-Thalamus-Proper')
                AcronSt = strvcat(AcronSt,'RH-ThalProp');
            elseif strcmp(deblank(Names(i,:)),'Right-Caudate')%|strcmp(deblank(Names(i,:)),'Right-Caudate')
                AcronSt = strvcat(AcronSt,'RH-Caudt');
            elseif strcmp(deblank(Names(i,:)),'Right-Putamen')%|strcmp(deblank(Names(i,:)),'Right-Putamen')
                AcronSt = strvcat(AcronSt,'RH-Putam');
            elseif strcmp(deblank(Names(i,:)),'Right-Pallidum')%|strcmp(deblank(Names(i,:)),'Right-Pallidum')
                AcronSt = strvcat(AcronSt,'RH-Pallid');
            elseif strcmp(deblank(Names(i,:)),'Right-Hippocampus')%|strcmp(deblank(Names(i,:)),'Right-Hippocampus')
                AcronSt = strvcat(AcronSt,'RH-Hippoc');
            elseif strcmp(deblank(Names(i,:)),'Right-Amygdala')%|strcmp(deblank(Names(i,:)),'Right-Amygdala')
                AcronSt = strvcat(AcronSt,'RH-Amygd');
            elseif strcmp(deblank(Names(i,:)),'Right-Accumbens-area')%|strcmp(deblank(Names(i,:)),'Right-Accumbens-area')
                AcronSt = strvcat(AcronSt,'RH-AccumbArea');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-bankssts')%|strcmp(deblank(Names(i,:)),'ctx-rh-bankssts')
                AcronSt = strvcat(AcronSt,'RH-BanksStS');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-caudalanteriorcingulate')%|strcmp(deblank(Names(i,:)),'ctx-rh-caudalanteriorcingulate')
                AcronSt = strvcat(AcronSt,'RH-CaudAntCing');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-caudalmiddlefrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-caudalmiddlefrontal')
                AcronSt = strvcat(AcronSt,'RH-CaudMidFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-cuneus')%|strcmp(deblank(Names(i,:)),'ctx-rh-cuneus')
                AcronSt = strvcat(AcronSt,'RH-Cuneus');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-entorhinal')%|strcmp(deblank(Names(i,:)),'ctx-rh-entorhinal')
                AcronSt = strvcat(AcronSt,'RH-Entorh');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-fusiform')%|strcmp(deblank(Names(i,:)),'ctx-rh-fusiform')
                AcronSt = strvcat(AcronSt,'RH-Fusif');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-inferiorparietal')%|strcmp(deblank(Names(i,:)),'ctx-rh-inferiorparietal')
                AcronSt = strvcat(AcronSt,'RH-InfPariet');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-inferiortemporal')%|strcmp(deblank(Names(i,:)),'ctx-rh-inferiortemporal')
                AcronSt = strvcat(AcronSt,'RH-InfTemp');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-isthmuscingulate')%|strcmp(deblank(Names(i,:)),'ctx-rh-isthmuscingulate')
                AcronSt = strvcat(AcronSt,'RH-IsthCing');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-lateraloccipital')%|strcmp(deblank(Names(i,:)),'ctx-rh-lateraloccipital')
                AcronSt = strvcat(AcronSt,'RH-LatOccip');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-lateralorbitofrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-lateralorbitofrontal')
                AcronSt = strvcat(AcronSt,'RH-LatOrbFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-lingual')%|strcmp(deblank(Names(i,:)),'ctx-rh-lingual')
                AcronSt = strvcat(AcronSt,'RH-Lingual');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-medialorbitofrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-medialorbitofrontal')
                AcronSt = strvcat(AcronSt,'RH-MedOrbFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-middletemporal')%|strcmp(deblank(Names(i,:)),'ctx-rh-middletemporal')
                AcronSt = strvcat(AcronSt,'RH-MidTemp');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-parahippocampal')%|strcmp(deblank(Names(i,:)),'ctx-rh-parahippocampal')
                AcronSt = strvcat(AcronSt,'RH-Parahip');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-paracentral')%|strcmp(deblank(Names(i,:)),'ctx-rh-paracentral')
                AcronSt = strvcat(AcronSt,'RH-Paracent');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-parsopercularis')%|strcmp(deblank(Names(i,:)),'ctx-rh-parsopercularis')
                AcronSt = strvcat(AcronSt,'RH-ParsoPerc');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-parsorbitalis')%|strcmp(deblank(Names(i,:)),'ctx-rh-parsorbitalis')
                AcronSt = strvcat(AcronSt,'RH-ParsOrb');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-parstriangularis')%|strcmp(deblank(Names(i,:)),'ctx-rh-parstriangularis')
                AcronSt = strvcat(AcronSt,'RH-ParsTriang');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-pericalcarine')%|strcmp(deblank(Names(i,:)),'ctx-rh-pericalcarine')
                AcronSt = strvcat(AcronSt,'RH-PeriCalc');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-postcentral')%|strcmp(deblank(Names(i,:)),'ctx-rh-postcentral')
                AcronSt = strvcat(AcronSt,'RH-PostCent');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-posteriorcingulate')%|strcmp(deblank(Names(i,:)),'ctx-rh-posteriorcingulate')
                AcronSt = strvcat(AcronSt,'RH-PostCing');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-precentral')%|strcmp(deblank(Names(i,:)),'ctx-rh-precentral')
                AcronSt = strvcat(AcronSt,'RH-Precent');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-precuneus')%|strcmp(deblank(Names(i,:)),'ctx-rh-precuneus')
                AcronSt = strvcat(AcronSt,'RH-Precuneus');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-rostralanteriorcingulate')%|strcmp(deblank(Names(i,:)),'ctx-rh-rostralanteriorcingulate')
                AcronSt = strvcat(AcronSt,'RH-RostAntCing');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-rostralmiddlefrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-rostralmiddlefrontal')
                AcronSt = strvcat(AcronSt,'RH-RostMidFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-superiorfrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-superiorfrontal')
                AcronSt = strvcat(AcronSt,'RH-SupFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-superiorparietal')%|strcmp(deblank(Names(i,:)),'ctx-rh-superiorparietal')
                AcronSt = strvcat(AcronSt,'RH-SupPariet');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-superiortemporal')%|strcmp(deblank(Names(i,:)),'ctx-rh-superiortemporal')
                AcronSt = strvcat(AcronSt,'RH-SupTemp');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-supramarginal')%|strcmp(deblank(Names(i,:)),'ctx-rh-supramarginal')
                AcronSt = strvcat(AcronSt,'RH-Supramarg');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-frontalpole')%|strcmp(deblank(Names(i,:)),'ctx-rh-frontalpole')
                AcronSt = strvcat(AcronSt,'RH-FrontPole');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-temporalpole')%|strcmp(deblank(Names(i,:)),'ctx-rh-temporalpole')
                AcronSt = strvcat(AcronSt,'RH-TempPole');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-transversetemporal')%|strcmp(deblank(Names(i,:)),'ctx-rh-transversetemporal')
                AcronSt = strvcat(AcronSt,'RH-TransvTemp');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-insula')%|strcmp(deblank(Names(i,:)),'ctx-rh-insula')
                AcronSt = strvcat(AcronSt,'RH-Insula');
            end
        end
end
return;

function [GMcodes,Names, Colors] = Gray_Matter_codes(atlastype);
%
% Syntax :
% [GMcodes] = Gray_Matter_codes(atlastype);
%
% This function extract cortical gray matter codes from a specified atlas.
%
% Input Parameters:
%   atlastype     : Atlas type.
%
%
% Output Parameters:
%   GMcodes       : Gray matter codes.
%   Names         : Structure Names
%   Colors        : Structure Colors
%
% Related references:
%
%
% See also:
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% March 20th 2012
% Version $1.0
GMcodes = '';
Names = '';
cont = 0;
Colors = [0 0 0];
col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
Ncolor = size(col,1);
switch atlastype
    case 'aparc+aseg'
        [a,temp] = system('echo $FREESURFER_HOME');
        txtfile =  strcat(temp,filesep,'FreeSurferColorLUT.txt');
        fid = fopen(txtfile);
        cont = 0;lines = '';
        while 1
            cont = cont + 1;
            line = fgetl(fid);
            if ~ischar(line),   break,   end
            lines = strvcat(lines,line);
        end
        ind = find(sum(isspace(lines'))==size(lines,2));lines(ind,:) = [];
        Codes = '';
        Names = '';
        cont = 0;
        Colors = [0 0 0];
        for i = 1:size(lines,1)
            if ~strcmp(lower(deblank(lines(i,1))),'#');
                cont = cont + 1;
                temp = strread(deblank(lines(i,:)),'%s','delimiter',' ');
                code = temp{1};
                name = temp{2};
                Codes = strvcat(Codes,code);
                Names = strvcat(Names,name);
                Colors(cont,:) = [str2num(temp{3}) str2num(temp{4}) str2num(temp{5})];
            end
        end
        fclose(fid);
        GMcodes = [10:13 17:18 26 49:54 58 1001:1003 1005:1035 2001:2003 2005:2035];
        index = ismember(str2num(Codes),GMcodes);
        Names = Names(index,:);
        Colors = Colors(index,:);
    case 'a2009s+aseg'
        [a,temp] = system('echo $FREESURFER_HOME');
        txtfile =  strcat(temp,filesep,'FreeSurferColorLUT.txt');
        fid = fopen(txtfile);
        cont = 0;lines = '';
        while 1
            cont = cont + 1;
            line = fgetl(fid);
            if ~ischar(line),   break,   end
            lines = strvcat(lines,line);
        end
        ind = find(sum(isspace(lines'))==size(lines,2));lines(ind,:) = [];
        Codes = '';
        Names = '';
        cont = 0;
        Colors = [0 0 0];
        for i = 1:size(lines,1)
            if ~strcmp(lower(deblank(lines(i,1))),'#');
                cont = cont + 1;
                temp = strread(deblank(lines(i,:)),'%s','delimiter',' ');
                code = temp{1};
                name = temp{2};
                Codes = strvcat(Codes,code);
                Names = strvcat(Names,name);
                Colors(cont,:) = [str2num(temp{3}) str2num(temp{4}) str2num(temp{5})];
            end
        end
        fclose(fid);
        GMcodes = [10:13 17:18 26 49:54 58 11101:11175 12101:12175];
        index = ismember(str2num(Codes),GMcodes);
        Names = Names(index,:);
        Colors = Colors(index,:);
    case 'a2005s+aseg'
        [a,temp] = system('echo $FREESURFER_HOME');
        txtfile =  strcat(temp,filesep,'FreeSurferColorLUT.txt');
        fid = fopen(txtfile);
        cont = 0;lines = '';
        while 1
            cont = cont + 1;
            line = fgetl(fid);
            if ~ischar(line),   break,   end
            lines = strvcat(lines,line);
        end
        ind = find(sum(isspace(lines'))==size(lines,2));lines(ind,:) = [];
        Codes = '';
        Names = '';
        cont = 0;
        Colors = [0 0 0];
        for i = 1:size(lines,1)
            if ~strcmp(lower(deblank(lines(i,1))),'#');
                cont = cont + 1;
                temp = strread(deblank(lines(i,:)),'%s','delimiter',' ');
                code = temp{1};
                name = temp{2};
                Codes = strvcat(Codes,code);
                Names = strvcat(Names,name);
                Colors(cont,:) = [str2num(temp{3}) str2num(temp{4}) str2num(temp{5})];
            end
        end
        fclose(fid);
        GMcodes = [10:13 17:18 26 49:54 58 1102:1181 2102:2181];
        index = ismember(str2num(Codes),GMcodes);
        Names = Names(index,:);
        Colors = Colors(index,:);
    case 'yeoh7networks'
        
    case 'yeoh17networks'
        
    case 'parckmeans'
        GMcodes = [10:13 17:18 26 49:54 58 11101:11175 12101:12175];
        re = floor(90/Ncolor); col = repmat(col,[re+1 1]);
    case 'ibaspm116'
        txtfile = which('atlas116.cod');
        re = floor(90/Ncolor); col = repmat(col,[re+1 1]);
        [Codes,b] = textread(txtfile,'%u%s','delimiter','=');
        %         Codes = str2num(char(a));
        Names = char(b);
        GMcodes = [1:90];
        index = ismember(Codes,GMcodes);
        Names = Names(index,:);
        Colors = col(1:90,:);
    case 'ibaspm71'
        txtfile = which('atlas71.cod');
        col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
        Ncolor = size(col,1);
        [Codes,b] = textread(txtfile,'%u%s','delimiter','=');
        %         Codes = str2num(char(a));
        Names = char(b);
        GMcodes = [1     2     4     5     6     7     9    10    11    12    14    15    16    18    19    20 ...
            23    25    26    27    32    33    36    37    38    39    41    50    52    53    54    56    60 ...
            61    62    63    64    67    69    70    72    74    75    76    80    85    88    90    97    98 ...
            99   101   102   108   110   112   114   119   125   130   132   140   145   154   159   164   165 ...
            175   196   203   251];
        index = ismember(Codes,GMcodes);
        re = floor(length(index)/Ncolor); col = repmat(col,[re+1 1]);
        Names = Names(index,:);
        Colors = col(1:length(index),:);
end

return


function  Connect = Read_Connectivity_Matrix(ConnectFilename,vis,cl);
%
% Syntax :
%  Connect = Read_Connectivity_Matrix(ConnectFilename,vis,cl);
%
% This script reads the connectivity Matrix from a connectivity matrix filename.
%
% Input Parameters:
%   ConnectFilename   :  Text Connectivity Matrix Filename
%      vis            :  Boolean Variable to plot, or not, the Connectivity
%                        Matrix
%        cl           :  Colormap
%
% Output Parameters:
%   Connect           :  Connectivity Matrix
%
% Related references:
%
%
% See also:
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% May 18th 2012
% Version $1.0

%=====================Checking Input Parameters===========================%
if nargin ==0
    error('Please select a correct connectivity matrix');
    return;
end
if nargin <2
    vis = 1;
    cl = 'jet';
end
if nargin ==2
    cl = 'jet';
end
%=====================End of Checking Input Parameters====================%
%==========================  Main Program   ==============================%
if ~exist(ConnectFilename,'file')
    error('The specified Matrix file does not exist');
    return;
end
fid  = fopen(ConnectFilename,'rt');
Nstruct = str2num(char(textread(ConnectFilename,'%s',1,'headerlines',8)));
Snames = textread(ConnectFilename,'%s',Nstruct*2,'headerlines',9);
temp = char(Snames); names= temp(2:2:end,:);codes= temp(1:2:end,:);
temp = textread(ConnectFilename,'%f',Nstruct*Nstruct,'headerlines',Nstruct+11);
Connect.Matrix = reshape(temp,[Nstruct Nstruct]);
Connect.StructCodes = str2num(deblank(codes));
Connect.StructNames = deblank(names);
if vis
    [pth,nm,ext] = fileparts(ConnectFilename);
    colordef black;figure;set(gcf,'Color',[0 0 0]);
    imagesc(Connect.Matrix);hold on;plot([1:Nstruct],[1:Nstruct],'xw');colormap(cl);
    set(gca,'YTick',[1:Nstruct]);set(gca,'XTick',[1:Nstruct]); set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',Connect.StructNames);
    title([nm ' ===> ' num2str(Nstruct) ' Structures']);
    
    h = colorbar;
    if ~strcmp(cl,'spectral')
        colormap(cl);
    else
        colormap(spectral(64));
    end
    range = max(Connect.Matrix(:))-min(Connect.Matrix(:));
    values =  min(Connect.Matrix(:)):range/10:max(Connect.Matrix(:));
    set(h,'YTickLabel',num2str(values'));
    
end
%=============================End=========================================%
return;






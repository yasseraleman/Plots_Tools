function Plot_Connectivity_Graph(connectResults, connThresh)

%% ======================== Input Parameters =========================== %%
% Connectivity Results
connectResults = '/media/Data/PROCESSING_RESULTS/HCP/7-connectome/HCP100_Results_394Regions.mat';
outDir = '/media/UserBackup/Presentaciones_LIM/Presentación_BRAD_29092016/394-Regions/';
outName = 'HCP100_Results_394Regions';

% Connectivity Threshold employed to compute the BigClam results
connThresh = linspace(0.55,1,30); % Example

boolsave = 1; % Boolean Variable to save screenshots from the network configuration.

clMap = 'greenred'; % Colormap applied to the regions according to its strength in the BigClam Matrix
% none: No colormap is applied.
% Possible colormaps:  'jet','hot','bluered','greenred','hsv','cool','bone','pink','winter','autumn',
% 'summer','spring','copper','spectral','gray','cyanyellow','yellowcyan','green','red','blue','redblue','redgreen','belle','manitou',
% 'munising','saginaw','saugatuck','blackred','blackgreen','blackblue','blackyellow','blackcyan','blackmagenta','g-w-r','dkb-w-dkr',
% 'b-w-r','c-w-y'

%% ====================== End of Input Parameters ====================== %%

% Closing all figures
close all;

% Reading Connectivity Matrices, Coordinates and region names
tempVar = load(connectResults);
connecMatVar = eval(['tempVar.' char(fieldnames(tempVar))]);

regCoords = connecMatVar.RegionsMMCoord;          % Region Coordinates
regNames = connecMatVar.RegionsNames;             % Region Names
conMatrix = mean(connecMatVar.ConnectMatrices,3); % Mean Connectivity Matrix. THIS CAN BE CHANGED BY AN INDIVIDUAL CONNECTIVITY MATRIX
conMatrix = conMatrix./max(conMatrix(:)+eps);



%% ====================== Initial Parameters =========================== %%
param.figcolor = [1 1 1];     % Figure Color
param.figurevisible = 'on';   % Figure Visibility

set(0,'units','centimeters');
cm_screen = get(0,'screensize');  % Getting Screen Size

figPosition = [1 1 cm_screen(3)-3 cm_screen(4)-3*cm_screen(3)/cm_screen(4)]; % Figure Position
figUnits = 'centimeters'; % Figure Units

scalValue = 0.05;         % Transparency factor (0 - 1);
faceColors = [.8 .8 .8];  % Face color for the background Surfaces
lightIntensity = .25;     %  Light Intensity factor (0 - 1);
%% ================== End of Initial Parameters ======================== %%


%% ================= Reading Background Surface ======================== %%
filName = which('Plot_Cluster_Regions_over_Brain_BRADE');
[pth, nm, ext] = fileparts(filName);
load([pth filesep 'Background_Surfaces.mat']);

%% ================= End of Reading Background Surface================== %%

% Creating Figure;
FigID = figure('numbertitle','off','name','New Figure','Color',param.figcolor, 'units',figUnits,'Position',figPosition,'Visible',param.figurevisible,'InvertHardcopy','off');


%% ================= Plotting Sagittal View ============================ %%
hsagaxis = subplot(1,2,1); % Sagittal subplot

% Left Patch
hsag(1)=patch(Surfl.SurfData,'edgecolor','none','tag', 'patch','facelighting','phong');
set(hsag(1),'SpecularExponent',60); % Specular value

% Right Patch
hsag(2)=patch(Surfr.SurfData,'edgecolor','none','tag', 'patch','facelighting','phong');
set(hsag(2),'SpecularExponent',60); % Specular value
axis off;view([270 0]);set(gca,'Color',[1 1 1]);axis image % Some parameters for the current Axis

% Faces Colors
hsag(1).FaceColor = faceColors;
hsag(2).FaceColor = faceColors;

% Faces Transparency
hsag(1).FaceAlpha = scalValue;
hsag(2).FaceAlpha = scalValue;

%% =================== End of Plotting Sagittal View =================== %%

%% ================= Plotting Axial View ============================ %%
haxaxis = subplot(1,2,2); % Axial subplot

% Left Patch
hax(1)=patch(Surfl.SurfData,'edgecolor','none','tag', 'patch','facelighting','phong');
set(hax(1),'SpecularExponent',60); % Specular value

% Right Patch
hax(2)=patch(Surfr.SurfData,'edgecolor','none','tag', 'patch','facelighting','phong');
set(hax(2),'SpecularExponent',60); % Specular value
axis off;view([0 90]);set(gca,'Color',[1 1 1]);axis image % Some parameters for the current Axis

% Faces Colors
hax(1).FaceColor = faceColors;
hax(2).FaceColor = faceColors;

% Faces Transparency
hax(1).FaceAlpha = scalValue;
hax(2).FaceAlpha = scalValue;

%% =================== End of Plotting Sagittal View =================== %%

%% ===== Colors to plot different networks inside the same cluster ===== %%
cylcolors = [ 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9;0 1 0;];
cylcolors1 = [240,163,255;0,117,220;153,63,0;76,0,92;0,92,49;43,206,72;...
    255,204,153;148,255,181;143,124,0;157,204,0;194,0,136;...
    0,51,128;255,164,5;255,168,187;66,102,0;255,0,16;94,241,242;0,153,143;...
    224,255,102;116,10,255;153,0,0;255,255,128;255,255,0;255,80,5]/255;
cylcolors = [cylcolors;cylcolors1];

Ncolor = size(cylcolors,1);

%% == End of Colors to plot different networks inside the same cluster = %%


%% ==================== Plotting Axial View ============================ %%
mkdir([outDir filesep 'outFigs']);
vidObj = VideoWriter([outDir filesep 'outFigs' filesep outName '.avi']);
open(vidObj);

Nstruct = size(conMatrix,1); % Number of structures
Nthr = length(connThresh);  % Number of thresholds
for i = 1:Nthr
    
    % Detecting networks from structures that belong to the same cluster
    LabNet = Label_Graph_Components(conMatrix,connThresh(i));
    
    re = floor(max(LabNet(:))/Ncolor); cylcolors = repmat(cylcolors,[re+1 1]);

    
    % Strength Value According to degree of connectivity
    tempMat = conMatrix./max(conMatrix(:)+eps);
    tempMat(tempMat <connThresh(i)) = 0;
    strengthConnect = sum(tempMat,2);
    
    [X,Y] = find(triu(LabNet,1)); % Connected regions
    
    % Assigning colors
    sphereColors = repmat([0.7 0.7 0.7],[Nstruct 1]); % Same color for all the spheres
    ind = unique([X;Y]); % Structures with degree ~=0
    switch clMap
        case 'none'
            sphereColors(ind,:) = repmat([0.9 0.1 0.1],[length(ind) 1]);  % Different color for the spheres that belong to the same cluster
        otherwise
            % Colormap Gradient
            [ut,i1t,it2] = unique(strengthConnect(ind));
            
            if length(ut) >1
                [colMap] = colormaps_colors(clMap,size(ut,1));
                
                % Correcting for non linearities in the values scale
                m = size(colMap,1);
                X0 = linspace (1, size(colMap,1), size(ut,1));
                colMap = interp1(1:m,colMap,X0);
                
                Y0 = cumsum([0;diff(ut)/sum(diff(ut))]);
                X1 = linspace(0, 1, size(colMap,1));
                a = 1;b = -1;c = 0; % Line parameters
                Xo = X1(:);Yo = Y0(:);
                x = (b*(b*Xo-a*Yo) -a*c)./(a.^2+b.^2);
                y = (a*(-b*Xo+a*Yo) -b*c)./(a.^2+b.^2);
                
                indi = find(x >max(Xo));
                x(indi) = max(Xo);
                indi = find(x <min(Xo));
                x(indi) = min(Xo);
                
                indi = find(y >max(Yo));
                y(indi) = max(Yo);
                indi = find(y <min(Yo));
                y(indi) = min(Yo);
                colMap = abs(interp1(X1,colMap,x));
            else
                [colMap] = colormaps_colors(clMap,10);
                colMap = colMap(end,:);
            end
            sphereColors(ind,:) = colMap(it2,:);  % Color Gradient for the spheres that belong to the same cluster according to its strength value
    end
    
    
    % Assigning Transparency
    faceAlpha = 0.5*ones(Nstruct,1);  % Transparency = 0.5 for all the structures
    faceAlpha(ind,:) = 1*ones(length(ind),1);        % Transparency = 1 for the structures belonging to the same cluster
    
    % 1.---  Activating Sagittal subplot
    subplot(1,2,1);
    
    % Plotting Spheres located at the structures center of coordinates.
    hpsag = Plot_Points_as_Spheres(regCoords, 'Color', sphereColors,'Radius',2,'FigID',FigID,'FaceAlpha',faceAlpha);
    hsagaxis.Title.String = ['Sagittal View']; % Axis Title
    hsagaxis.FontSize = 15; % Title fontsize
    hsagaxis.Title.Color = [1 1 1] - param.figcolor; % Title color
    
    
    % 2.---  Activating Axial subplot
    subplot(1,2,2);
    
    % Plotting Spheres located at the structures center of coordinates.
    hpaxial = Plot_Points_as_Spheres(regCoords, 'Color', sphereColors,'Radius',2,'FigID',FigID,'FaceAlpha',faceAlpha);
    haxaxis.Title.String = [' Axial View']; % Axis Title
    haxaxis.FontSize = 15; % Title fontsize
    haxaxis.Title.Color = [1 1 1] - param.figcolor; % Title color
    
    
    % 3.--- Plotting Connected regions
    for j = 1:length(X)
        X1 = regCoords(X(j),:);
        X2 = regCoords(Y(j),:);
        
        subplot(1,2,1);hold on;
        color = cylcolors(LabNet(X(j),Y(j)),:);
        [hcyl(2*(j-1)+1) EndPlate1 EndPlate2] = Cylinder_test(X1,X2,conMatrix(X(j),Y(j)),20,color,0,0);
        subplot(1,2,2);hold on;
        [hcyl(2*j) EndPlate1 EndPlate2] = Cylinder_test(X1,X2,conMatrix(X(j),Y(j)),20,color,0,0);
        
    end
    
    % 4. --- Lights
    subplot(1,2,1);
    [x,y, z] = sph2cart(35,0,1);
    hl1 = light('Color',[1 1 1]*lightIntensity,'Position',[x y z]);
    
    [x,y, z] = sph2cart(125,0,1);
    hl2 = light('Color',[1 1 1]*lightIntensity,'Position',[x y z]);
    
    [x,y, z] = sph2cart(215,0,1);
    hl3 = light('Color',[1 1 1]*lightIntensity,'Position',[x y z]);
    
    [x,y, z] = sph2cart(305,0,1);
    hl4 = light('Color',[1 1 1]*lightIntensity,'Position',[x y z]);
    
    [x,y, z] = sph2cart(0,90,1);
    hl5 = light('Color',[1 1 1]*lightIntensity,'Position',[x y z]);
    
    [x,y, z] = sph2cart(0,270,1);
    hl6 = light('Color',[1 1 1]*lightIntensity,'Position',[x y z]);
    
    subplot(1,2,2);
    [x,y, z] = sph2cart(35,0,1);
    hl1 = light('Color',[1 1 1]*lightIntensity,'Position',[x y z]);
    
    [x,y, z] = sph2cart(125,0,1);
    hl2 = light('Color',[1 1 1]*lightIntensity,'Position',[x y z]);
    
    [x,y, z] = sph2cart(215,0,1);
    hl3 = light('Color',[1 1 1]*lightIntensity,'Position',[x y z]);
    
    [x,y, z] = sph2cart(305,0,1);
    hl4 = light('Color',[1 1 1]*lightIntensity,'Position',[x y z]);
    
    [x,y, z] = sph2cart(0,90,1);
    hl5 = light('Color',[1 1 1]*lightIntensity,'Position',[x y z]);
    
    [x,y, z] = sph2cart(0,270,1);
    hl6 = light('Color',[1 1 1]*lightIntensity,'Position',[x y z]);
    
    
    % 5.--- Creating Colorbar
    switch clMap
        case 'none'
        otherwise
            % ColorBar
            if length(ut) >1
                hc = colorbar('peer',gca);
                colormap(gca,colMap);
                [uord,ord] = unique(strengthConnect(ind));
                minVals = min(strengthConnect(ind));
                maxVals = max(strengthConnect(ind));
                if length(ord) < 11
                    nTicks = length(ord);
                    hc.Ticks = linspace(min(hc.Ticks),max(hc.Ticks),nTicks)
                else
                    nTicks = 11;
                end
                interValues = floor(linspace(1,length(ord),nTicks));
                values = uord(interValues);
                tempc = colMap(interValues,:);
                defvalues = linspace(values(1),values(end),nTicks); % Definitive Tick Values

                col = interp1(values,tempc,defvalues); % Definitive ColorValues
                range = max(defvalues)-min(defvalues);
                valuesInterp =  min(minVals):range/255:max(maxVals);
                colMap = interp1(defvalues,col,valuesInterp);
                colormap(gca,colMap);
                set(hc,'YTickLabel',num2str(defvalues(:)));
                %

            else
                hc = colorbar('peer',gca);
                hc.Ticks = max(hc.Ticks);
                colormap(gca,colMap);
                set(hc,'YTickLabel',num2str(unique(strengthConnect(ind))));

            end
                 figColor = get(FigID,'Color');
                hc.Title.String = ['Deg. of Conn.: Thresh = ' num2str(connThresh(i))];
                hc.FontName = 'Arial';
                hc.Color = ([1 1 1] - figColor)*.8;
                hc.FontSize = 15; % Title fontsize
    end
    pause(4);
    % 6.--- Saving Screenshot
    if boolsave
              mkdir([outDir filesep 'outFigs']);
              outFileName = [outDir filesep 'outFigs' filesep outName 'Thresh-' num2str(connThresh(i)) '.tiff'];
              figure(FigID);
              a = getframe(FigID);
              writeVideo(vidObj,a);
              imwrite(a.cdata,outFileName,'TIFF','Resolution',[100 100]);
    end
    
    % 6.--- Cleanning spheres and cylinders
    if i~= Nthr
        for j = 1:length(hpsag)
            delete(hpaxial(j));
            delete(hpsag(j));
        end
        for j = 1:2*length(X)
            delete(hcyl(j));
        end
    end
end
close(vidObj);
return;

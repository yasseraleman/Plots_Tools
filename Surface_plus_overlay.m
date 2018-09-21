function OutFigures = Surface_plus_overlay(SurfFile, ChFiles, OutDir, plotopts);
%
% Syntax :
% OutFigures = Surface_plus_overlay(SurfFile, ChFile, OutDir);
%
% The following program shows how overlay maps on top surfaces.
%
% Input Parameters:
%      SurfFile         : Surface File
%       ChFiles         : Characteristics Files
%         OutDir        : Output Directory
%         plotopts      : Plot options(ie. plotopts.cl : Plot Colormap)
%
%  Output Parameters:
%       OutFigures      : Figures Filenames
%
% Related references:
%
% See also:
%__________________________________________________
% Authors:  Yasser Aleman Gomez
% Unidad de Medicina y Cirugía Experimental, Hospital General Universitario Gregorio Marañón, Madrid
%
% November 1st 2012
% Version $1.0

%% ================== Cheking Input Parameters ========================= %%
% SurfFile = '/home/yaleman/Data/5-freesurfer_processing/ASPER_00001__101-20060510/surf/lh.inflated';
% ChFile = '/home/yaleman/Data/5-freesurfer_processing/ASPER_00001__101-20060510/surf/lh.curv';
% %ChFile = '/home/yaleman/Data/5-freesurfer_processing/ASPER_00001__101-20060510/label/lh.aparc.a2009s.annot';
[pth,name,ext] = fileparts(deblank(SurfFile));
if nargin < 3
    OutDir = pth;
end
if nargin < 4
    plotopts.cl = 'spectral';
    plotopts.figurevisible = 'off';
    plotopts.figcolor = 'black';
end
%% ================== End of Cheking Input Parameters ================== %%
%% ===================== Detecting Hemisphere ========================== %%
indl= strfind(SurfFile,'lh.');
indr= strfind(SurfFile,'rh.');
if ~isempty(indl)
    hemi = 'lh';clear indl;
elseif ~isempty(indr)
    hemi = 'rh';clear indr;
else
    warning('Hemispheres information could not be detected. ');
    hemi = '';
end
%% ===================== End of Detecting Hemisphere =================== %%
%Surf = read_surfreesurfer(SurfFile);
[OutFiles, SurfF] = Exp_Surf(SurfFile, '0', '','', 'imp','n'); Surf= SurfF{1}; % Reading Surface
if strfind(SurfFile,'.pial')
    name = 'pial';
elseif strfind(SurfFile, '.white');
    name = 'white';
elseif strfind(SurfFile, '.inflated');
    name = 'inflated';
elseif strfind(SurfFile, '.sphere');
    name = 'sphere';
elseif strfind(SurfFile, '.orig');
    name = 'orig';
end
OutFigures = '';
for k = 1:size(ChFiles,1)
    ChFile = deblank(ChFiles(k,:));
    [ptho,nameo,exto] = fileparts(deblank(ChFile));
    [txt,colortable] = read_cfiles(ChFile); % Reading Characteristics
    Surf.Is = txt;
    %% ===================== Painting Surface ============================== %%
    if colortable.table == 0
        [Colors] = Surf_Color(Surf,plotopts.cl);
        Surf.SurfData.FaceVertexCData = Colors;
    else
        sts = unique(txt);
        Ns = length(sts);
        Surf.SurfData.FaceVertexCData = zeros(size(Surf.SurfData.vertices,1),3);
        for i = 1:Ns
            ind = find(txt==sts(i));
            indc = find(colortable.table(:,5)==sts(i));
            if isempty(indc)
                Matname = 'unknown_structure';
                Color =   [1 1 1];      % Color
            else
                Matname = char(colortable.struct_names{indc});
                Color =  colortable.table(indc,1:3)/255;       % Color
            end
            Surf.SurfData.FaceVertexCData(ind,:) = repmat(Color,[length(ind) 1 ]);
        end
    end
    
    %% ===================== End of Painting Surface ======================= %%
    %% ===================== Creating Figures ============================== %%
    
    if strfind(ChFile,'.thickness')
        nameo = 'thickness';
    elseif strfind(ChFile, '.aparc.annot');
        nameo = 'aparc';
    elseif strfind(ChFile, '.a2009s.annot');
        nameo = 'a2009s';
    elseif strfind(ChFile, '.a2005s.annot');
        nameo = 'a2005s';
    elseif strfind(ChFile, '.curv');
        nameo = 'curv';
    elseif strfind(ChFile, '.sulc');
        nameo = 'sulc';
    end
    Basename = [name '_plus_' nameo];
    BasenameT = [name ' + ' nameo];
    if strcmp(hemi,'lh')
        % ----------- Plotting Surface  --------------------------------------%
        switch plotopts.figcolor
            case 'black'
                colordef black;
                colordef black; hf = figure('numbertitle','off','Color',[0 0 0],'Position',[0 0 1200 900],'Visible',plotopts.figurevisible,'InvertHardcopy','off');
            case 'white'
                colordef white;
                colordef white; hf = figure('numbertitle','off','Color',[1 1 1],'Position',[0 0 1200 900],'Visible',plotopts.figurevisible,'InvertHardcopy','off');
        end
        custom_plotsurf(Surf);
        view([270 0]);axis off;axis tight;axis equal; h=title(['Left Hemisphere. Lateral View: ' BasenameT]);set(h,'FontSize',15,'FontName','Arial');
        camlight;
        % ----------- Creating Legend ---------------------------------------%
        if sum(Surf.Is-floor(Surf.Is)) ~=0
            h = colorbar;
            if ~strcmp(plotopts.cl,'spectral')
                colormap(plotopts.cl);
            else
                colormap(spectral(64));
            end
            range = max(Surf(1).Is)-min(Surf(1).Is);values =  min(Surf(1).Is):range/10:max(Surf(1).Is);
            set(h,'YTickLabel',num2str(values'));
        end
        % ----------- Saving Figure ---------------------------------------%
        nametotal = ['LH_' Basename '-Lateral_View'];
        Figurename = [OutDir filesep nametotal '.jpg'];
        a = getframe(hf);
        imwrite(a.cdata,Figurename,'JPEG')
        OutFigures  = strvcat(OutFigures,Figurename);
        close(hf);
        
        % ----------- Plotting Surface  --------------------------------------%
        switch plotopts.figcolor
            case 'black'
                colordef black;
                colordef black; hf = figure('numbertitle','off','Color',[0 0 0],'Position',[0 0 1200 900],'Visible',plotopts.figurevisible,'InvertHardcopy','off');
            case 'white'
                colordef white;
                colordef white; hf = figure('numbertitle','off','Color',[1 1 1],'Position',[0 0 1200 900],'Visible',plotopts.figurevisible,'InvertHardcopy','off');
        end
        custom_plotsurf(Surf);
        view([90 0]);axis off;axis tight;axis equal; h=title(['Left Hemisphere. Medial View: ' BasenameT]);set(h,'FontSize',15,'FontName','Arial');
        camlight;
        
        % ----------- Creating Legend ---------------------------------------%
        if sum(Surf.Is-floor(Surf.Is)) ~=0
            h = colorbar;
            if ~strcmp(plotopts.cl,'spectral')
                colormap(plotopts.cl);
            else
                colormap(spectral(64));
            end
            range = max(Surf(1).Is)-min(Surf(1).Is);values =  min(Surf(1).Is):range/10:max(Surf(1).Is);
            set(h,'YTickLabel',num2str(values'));
        end
        
        % ----------- Saving Figure ---------------------------------------%
        nametotal = ['LH_' Basename '-Medial_View'];
        Figurename = [OutDir filesep nametotal '.jpg'];
        a = getframe(hf);
        imwrite(a.cdata,Figurename,'JPEG')
        OutFigures  = strvcat(OutFigures,Figurename);
        close(hf);
        
    elseif strcmp(hemi,'rh')
        
        % ----------- Plotting Surface  --------------------------------------%
        switch plotopts.figcolor
            case 'black'
                colordef black;
                colordef black; hf = figure('numbertitle','off','Color',[0 0 0],'Position',[0 0 1200 900],'Visible',plotopts.figurevisible,'InvertHardcopy','off');
            case 'white'
                colordef white;
                colordef white; hf = figure('numbertitle','off','Color',[1 1 1],'Position',[0 0 1200 900],'Visible',plotopts.figurevisible,'InvertHardcopy','off');
        end
        custom_plotsurf(Surf);
        view([90 0]);axis off;axis tight;axis equal; h=title(['Right Hemisphere. Lateral View: ' BasenameT]);set(h,'FontSize',15,'FontName','Arial');
        camlight;
        
        % ----------- Creating Legend ---------------------------------------%
        if sum(Surf.Is-floor(Surf.Is)) ~=0
            h = colorbar;
            if ~strcmp(plotopts.cl,'spectral')
                colormap(plotopts.cl);
            else
                colormap(spectral(64));
            end
            range = max(Surf(1).Is)-min(Surf(1).Is);values =  min(Surf(1).Is):range/10:max(Surf(1).Is);
            set(h,'YTickLabel',num2str(values'));
        end
        
        % ----------- Saving Figure ---------------------------------------%
        nametotal = ['RH_' Basename '-Lateral_View'];
        Figurename = [OutDir filesep nametotal '.jpg'];
        a = getframe(hf);
        imwrite(a.cdata,Figurename,'JPEG')
        OutFigures  = strvcat(OutFigures,Figurename);
        close(hf);
        
        % ----------- Plotting Surface  --------------------------------------%
        switch plotopts.figcolor
            case 'black'
                colordef black;
                colordef black; hf = figure('numbertitle','off','Color',[0 0 0],'Position',[0 0 1200 900],'Visible',plotopts.figurevisible,'InvertHardcopy','off');
            case 'white'
                colordef white;
                colordef white; hf = figure('numbertitle','off','Color',[1 1 1],'Position',[0 0 1200 900],'Visible',plotopts.figurevisible,'InvertHardcopy','off');
        end
        custom_plotsurf(Surf);
        view([270 0]);axis off;axis tight;axis equal; h=title(['Right Hemisphere. Medial View: ' BasenameT]);camlight;set(h,'FontSize',15,'FontName','Arial');
        
        % ----------- Creating Legend ---------------------------------------%
        if sum(Surf.Is-floor(Surf.Is)) ~=0
            h = colorbar;
            if ~strcmp(plotopts.cl,'spectral')
                colormap(plotopts.cl);
            else
                colormap(spectral(64));
            end
            range = max(Surf(1).Is)-min(Surf(1).Is);values =  min(Surf(1).Is):range/10:max(Surf(1).Is);
            set(h,'YTickLabel',num2str(values'));
        end
        
        % ----------- Saving Figure ---------------------------------------%
        nametotal = ['RH_' Basename '-Medial_View'];
        Figurename = [OutDir filesep nametotal '.jpg'];
        a = getframe(hf);
        imwrite(a.cdata,Figurename,'JPEG')
        OutFigures  = strvcat(OutFigures,Figurename);
        close(hf);
    else
        
        % ----------- Plotting Surface  --------------------------------------%
        switch plotopts.figcolor
            case 'black'
                colordef black;
                colordef black; hf = figure('numbertitle','off','Color',[0 0 0],'Position',[0 0 1200 900],'Visible',plotopts.figurevisible,'InvertHardcopy','off');
            case 'white'
                colordef white;
                colordef white; hf = figure('numbertitle','off','Color',[1 1 1],'Position',[0 0 1200 900],'Visible',plotopts.figurevisible,'InvertHardcopy','off');
        end
        custom_plotsurf(Surf);
        view([270 0]);axis off;axis tight;axis equal;
        camlight;
        
        % ----------- Creating Legend ---------------------------------------%
        if sum(Surf.Is-floor(Surf.Is)) ~=0
            h = colorbar;
            if ~strcmp(plotopts.cl,'spectral')
                colormap(plotopts.cl);
            else
                colormap(spectral(64));
            end
            range = max(Surf(1).Is)-min(Surf(1).Is);values =  min(Surf(1).Is):range/10:max(Surf(1).Is);
            set(h,'YTickLabel',num2str(values'));
        end
        
        % ----------- Saving Figure ---------------------------------------%
        nametotal = [Basename '-View_1'];
        Figurename = [OutDir filesep nametotal '.jpg'];
        a = getframe(hf);
        imwrite(a.cdata,Figurename,'JPEG')
        OutFigures  = strvcat(OutFigures,Figurename);
        close(hf);
        
        % ----------- Plotting Surface  --------------------------------------%
        switch plotopts.figcolor
            case 'black'
                colordef black;
                colordef black; hf = figure('numbertitle','off','Color',[0 0 0],'Position',[0 0 1200 900],'Visible',plotopts.figurevisible,'InvertHardcopy','off');
            case 'white'
                colordef white;
                colordef white; hf = figure('numbertitle','off','Color',[1 1 1],'Position',[0 0 1200 900],'Visible',plotopts.figurevisible,'InvertHardcopy','off');
        end
        custom_plotsurf(Surf);
        view([90 0]);axis off;axis tight;axis equal;
        camlight
        
        % ----------- Creating Legend ---------------------------------------%
        if sum(Surf.Is-floor(Surf.Is)) ~=0
            h = colorbar;
            if ~strcmp(plotopts.cl,'spectral')
                colormap(plotopts.cl);
            else
                colormap(spectral(64));
            end
            range = max(Surf(1).Is)-min(Surf(1).Is);values =  min(Surf(1).Is):range/10:max(Surf(1).Is);
            set(h,'YTickLabel',num2str(values'));
        end
        % ----------- Saving Figure ---------------------------------------%
        nametotal = [Basename '-View_1'];
        Figurename = [OutDir filesep nametotal '.jpg'];
        print( gcf,'-djpeg','-r300',Figurename);
%         a = getframe(hf);
%         imwrite(a.cdata,Figurename,'JPEG');
        OutFigures  = strvcat(OutFigures,Figurename);
        close(hf);
        
    end
end
return
%% ===================== End of Creating Figures ======================= %%


function custom_plotsurf(Surf);
Surf.SurfData.FaceColor = 'interp';
if isunix
    strsurf=patch(Surf.SurfData,'edgecolor','none', 'tag','patch','facelighting','gouraud');
else
    strsurf=patch(Surf.SurfData,'edgecolor','black','tag', 'patch','facelighting','gouraud');
end
return

function [Colors] = Surf_Color_custom(Surf,cl);
% Syntax :
% [Colors] = Surf_Color(Surf,cl);
%
% This function plots the surfaces contained in the cell file Surfa.
%
% Input Parameters:
%   Surf        : Surface variable.
%   cl          : Colormap used to see the results.
%
% Output Parameters:
%  Colors       : Output colormap matrix.
%
% Related references:
%
%
% See also: Smooth_Surf Surf_Comp Plot_Surf Plot_oversurf Exp_Surf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% September 29th 2007
% Version $1.0

%=========================Main program====================================%
txt=Surf.Is;
if sum(txt-floor(txt))
    [ut,i1t,it2] = unique(txt);
    switch cl
        case 'hot'
            col = hot(size(ut,1));
        case 'hsv'
            col = hsv(size(ut,1));
        case 'jet'
            col = jet(size(ut,1));
        case 'cool'
            col = cool(size(ut,1));
        case 'bone'
            col = bone(size(ut,1));
        case 'pink'
            col = pink(size(ut,1));
        case 'winter'
            col = winter(size(ut,1));
        case 'autumn'
            col = autumn(size(ut,1));
        case 'summer'
            col = summer(size(ut,1));
        case 'spring'
            col = spring(size(ut,1));
        case 'copper'
            col = copper(size(ut,1));
        case 'spectral'
            col = spectral(size(ut,1));
    end
    Colors = col(it2,:);
elseif sum(txt-floor(txt)) ==0
    %col = [213 221 227;175 206 227;149 196 228; 120 186 230;87 172 231;24 146 232;6 73 119;244 207 154;244 192 117;244 179 86;244 161 43; 212 133 20;158 101 19; 113 71 12]/255;
    col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
    Ncolor = size(col,1);
    ut = sort(unique(txt));
    re = floor(length(ut)/Ncolor); col = repmat(col,[re+1 1]);
    for j = 1:size(ut,1)
        indpos = find(txt == ut(j)); Colors(indpos,:) = repmat(col(j,:),[size(indpos,1) 1]);
    end
end
%========================End of main program==============================%
return

function s = spectral(m)
%SPECTRAL Black-purple-blue-green-yellow-red-white color map.
%
%         map = spectral(num_colors)
%
% SPECTRAL(M) returns an M-by-3 matrix containing a "spectral" colormap.
% SPECTRAL, by itself, is the same length as the current colormap.
%
% For example, to reset the colormap of the current figure:
%
%           colormap(spectral)
%
% See also HSV, GRAY, PINK, HOT, COOL, BONE, COPPER, FLAG,
%          COLORMAP, RGBPLOT.

if nargin < 1, m = size(get(gcf,'colormap'),1); end
base = [
    0.0000 0.0000 0.0000
    0.4667 0.0000 0.5333
    0.5333 0.0000 0.6000
    0.0000 0.0000 0.6667
    0.0000 0.0000 0.8667
    0.0000 0.4667 0.8667
    0.0000 0.6000 0.8667
    0.0000 0.6667 0.6667
    0.0000 0.6667 0.5333
    0.0000 0.6000 0.0000
    0.0000 0.7333 0.0000
    0.0000 0.8667 0.0000
    0.0000 1.0000 0.0000
    0.7333 1.0000 0.0000
    0.9333 0.9333 0.0000
    1.0000 0.8000 0.0000
    1.0000 0.6000 0.0000
    1.0000 0.0000 0.0000
    0.8667 0.0000 0.0000
    0.8000 0.0000 0.0000
    0.8000 0.8000 0.8000
    ];
n = length(base);
X0 = linspace (1, n, m);
s = interp1(1:n,base,X0);

return
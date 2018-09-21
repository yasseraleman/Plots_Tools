function OutFigures = Image_plus_surface(InputImage,Surfa, slices, OutDir, plotopts);
%
% Syntax :
% OutFigures = Image_plus_surface(InputImage,Surfa, slices, OutDir);
%
% The following program shows how overlay surface lines on top of anatomical 
% scans.
%
% Input Parameters:
%     InputImage       : Input Anatomical Image
%       Surfa          : Surfaces
%        slices        : Nx3 vector specifying wich slices will be ploted.
%                         ie. [Sagital_Slice Coronal_Slice  Axial_Slice];
%         OutDir        : Output Directory
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
warning off;
close all;
V = spm_vol(InputImage);
I = spm_read_vols(V);
[pth,name,ext] = fileparts(V(1).fname);
if nargin < 5
    plotopts.figurevisible = 'off';
    plotopts.figcolor = 'black';
end
if nargin < 4
    OutDir = pth;
end
if ischar(Surfa(1,:));
    Ns = size(Surfa,1)
    for i = 1:Ns
        [OutFiles, SurfF] = Exp_Surf(deblank(Surfa(i,:)), '0', '','', 'imp','n'); % Reading Surface;
        Surft(i) = SurfF{1}; 
    end
    Surfa = Surft;
end
Surfa = Surfa(:);
%% ============== End of Cheking Input Parameters ====================== %%
col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
Ncolor = size(col,1);
re = floor(size(Surfa,1)/Ncolor); 
col = repmat(col,[re+1 1]);
Ns = size(slices,1);
voxsize = sqrt(sum(V.mat(1:3,1:3).^2));
I = permute(I,[2 1 3]);
v = 90*eye(3);
[pth,nm,ext] = fileparts(deblank(InputImage));
OutFigures = '';
TitleF = ['Image + Surface: ' upper(name) ];
ATtile = [ upper(name) '_plus_Lines'];
for  t = 1:Ns
    slicet = slices(t,:);
    ind = find(slicet);
    for i = 1:length(ind)
        if ind(i) == 3
            nametotal = [ATtile '-Axial_' num2str(slicet(3))]; % Axial Slices
            DTtitle = [TitleF ' : Axial Slice Number ' num2str(slicet(3))];
        elseif ind(i) == 2
            nametotal = [ATtile '-Coronal_' num2str(slicet(2))]; % Coronal Slices
            DTtitle = [TitleF ' : Coronal Slice Number ' num2str(slicet(2))];
        elseif ind(i) == 1
            nametotal = [ATtile '-Sagital_' num2str(slicet(1))]; % Sagital Slices
            DTtitle = [TitleF ' : Sagital Slice Number ' num2str(slicet(1))];
        end
        slicei = [ 0 0 0];
        slicei(ind(i)) = slicet(ind(i));
        switch plotopts.figcolor
            case 'black'
                colordef black;
                colordef black; hf = figure('numbertitle','off','name',TitleF,'Color',[0 0 0],'Position',[0 0 1200 900],'Visible',plotopts.figurevisible,'InvertHardcopy','off');
            case 'white'
                colordef white;
                colordef white; hf = figure('numbertitle','off','name',TitleF,'Color',[1 1 1],'Position',[0 0 1200 900],'Visible',plotopts.figurevisible,'InvertHardcopy','off');
        end
        %% ============== Creating Slice =================================%
        fig = slice(I,slicei(1),slicei(2),slicei(3));
        slicei = logical(slicei);
        fign = fig(~slicei(1:3));
        figy = fig(slicei(1:3));
        if ~isempty(fign)&size(fign,1)==2
            set(fign(1),'Visible','off');
            set(fign(2),'Visible','off');
        elseif ~isempty(fign)&size(fign,1)==1
            set(fign(1),'Visible','off');
        end
        set(figy,'LineStyle','none');
        colormap gray; axis equal; hold on;
        if sum(slicei,2)==1
            view(v(slicei,[1 3]));
        else
            view(3)
        end
        axis off;
        %% ============== End of Creating Slice ==========================%
%         h1 =gca;
%         set(gca,'XTicklabel',num2str((str2num(get(gca,'XTicklabel')))-abs(V.mat(1,4))));set(h1,'XColor',[1 1 1]);xlabel(h1,' X axis(mm)');
%         set(gca,'YTicklabel',num2str((str2num(get(gca,'YTicklabel')))-abs(V.mat(2,4))));set(h1,'YColor',[1 1 1]);ylabel(h1,' Y axis(mm)');
%         set(gca,'ZTicklabel',num2str((str2num(get(gca,'ZTicklabel')))-abs(V.mat(3,4))));set(h1,'ZColor',[1 1 1]);zlabel(h1,' Z axis(mm)');
%         set(gca,'Color',[0 0 0]);
        title(DTtitle);
        plane = double(slicei);
        %% ============== Plotting Surface ===============================%
        plane = double(slicei);
        indt = find(plane == 1);
        for j = 1:size(Surfa,1)
            Surf = Surfa(j);
            if ~isfield(Surf,'color')
                Surf.color = col(j,:);
            end
            vertvox = (inv(V.mat)*[Surf.SurfData.vertices ones(size(Surf.SurfData.vertices),1)]')';
            Surf.SurfData.vertices = vertvox(:,1:3);
            for k = 1:size(indt,2)
                planet = zeros(1,4);
                planet(indt(k))=1;
                planet(4) = -slices(t,indt(k));
                [Xline,Yline,Zline] = cut_surf(planet,Surf.SurfData);
                line(Xline,Yline,Zline,'LineWidth',1.5,'Color',Surf.color);
            end
         
        end
        %% ============== Plotting Surface ===============================%
        Figurename = [OutDir filesep nametotal '.jpg'];
         print( gcf,'-djpeg','-r300',Figurename);
%         a = getframe(hf);
%         imwrite(a.cdata,Figurename,'JPEG')
        OutFigures  = strvcat(OutFigures,Figurename);
        close(hf);
    end
end
return;

function [Xline,Yline,Zline] = cut_surf(plane,fv)

% This function cuts a surface (defined by vertices and faces as represented
% by MATLAB patches) by a plane, leading to a curve in 3D space. The
% resulting curve is represented by a set of contigous lines in the space
%
% Syntax:
% [Xline,Yline,Zline] = cut_surf(plane,SurfData)
%
% INPUTS:
% plane : A 4-length vector with the parameters of the plane. If plane = [A
% B C D] then every 3D point P = (x,y,z) belonging to the plane satisfies
% plane*[P; 1]' = A*x + B*y + C*z + D = 0
% SurfData : surface structure as represented in MATLAB by patches:
%            SurfData.vertices
%            SurfData.faces
%
% OUTPUTS:
% Xline,Yline,Zline: Matrices with the line coordinates.
% The entire curve can be plotted by simply typing:
% line(Xline,Yline,Zline,'Properties1',Value1,...);
%
% Pedro Antonio Vald�s Hern�ndez
%
% October 29 2008

warning off %#ok
Xline = []; Yline = []; Zline = [];
oo = ones(size(fv.vertices,1),1);
maxdist = sqrt(max(sum((fv.vertices(fv.faces(:,1),:)-fv.vertices(fv.faces(:,2),:)).^2,2)));
maxdist = max(maxdist,sqrt(max(sum((fv.vertices(fv.faces(:,1),:)-fv.vertices(fv.faces(:,3),:)).^2,2))));
maxdist = max(maxdist,sqrt(max(sum((fv.vertices(fv.faces(:,2),:)-fv.vertices(fv.faces(:,3),:)).^2,2))));
vertx = find(abs(dot([fv.vertices oo],plane(oo,:),2)/norm(plane(1:3)))<maxdist);
indf = ismember(fv.faces,vertx);
[rindf,cindf] = find(indf); %#ok
rindf = unique(rindf);
Nf = length(rindf);
% h = waitbar(0,'cutting surface...');
for i = 1:Nf
    verts = fv.vertices(fv.faces(rindf(i),:),:);
    verts(:,4) = 1;
    diffv(1,:) = diff(verts([1 2],:));
    diffv(2,:) = diff(verts([2 3],:));
    diffv(3,:) = diff(verts([3 1],:));
    alpha = -verts*plane'./(diffv*plane');
    % NaN   : contains
    % < 0   : not contains down
    % -Inf  : parallel down
    % > 1   : not contains up
    % +Inf  : parallel up
    ind = find((alpha<1 & alpha >=0) | (alpha<=1 & alpha >0))  ;
    if ~isempty(ind) && length(ind)==2
        points = verts(ind,1:3) + alpha(ind,[1 1 1]).*diffv(ind,1:3);
        Xline = [Xline points(:,1)]; %#ok
        Yline = [Yline points(:,2)]; %#ok
        Zline = [Zline points(:,3)]; %#ok
    end
    %    waitbar(i/Nf,h);
end
% close(h);
warning on %#ok
return
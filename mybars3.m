function mybars3(Ms,Stds)
w = 0.3;
h = 8;
posx = 2;
posy = 3;
sh =1;
% Data = matres;
% options.linec = [0 0 0];
% options.line = 'on';
% options.facec = [1 0 0];
% Ms = squeeze(mean(Data,1));
% Stds = squeeze(std(Data,1));


fvt.faces = [0 0 0];
h = Ms(:)*100;
sh = Stds(:)*100;
[X,Y] = ind2sub(size(Ms),find(Ms));
for i = 1:length(X)
    if i ==1
        fv.vertices = ([1 -1 0;1 1 0;-1 1 0;-1 -1 0;1 -1 1;1 1 1;-1 1 1;-1 -1 1].*repmat([w/2 w/2 h(i)-sh(i)],[8 1]))+repmat([X(i) Y(i) 0],[8 1]);
        fv.faces= [1 2 3;1 3 4;1 2 6;1 6 5;2 3 6;3 7 6;1 4 8;1 8 5;3 4 7;4 7 8;5 6 7;5 7 8];
        fvss.vertices = ([1 -1 0;1 1 0;-1 1 0;-1 -1 0;1 -1 1;1 1 1;-1 1 1;-1 -1 1].*repmat([w/2 w/2 sh(i)],[8 1]))+repmat([X(i) Y(i) h(i)],[8 1]);
        fvss.faces= [1 2 3;1 3 4;1 2 6;1 6 5;2 3 6;3 7 6;1 4 8;1 8 5;3 4 7;4 7 8;5 6 7;5 7 8];
        fvsi.vertices = ([1 -1 0;1 1 0;-1 1 0;-1 -1 0;1 -1 1;1 1 1;-1 1 1;-1 -1 1].*repmat([w/2 w/2 sh(i)],[8 1]))+repmat([X(i) Y(i) h(i)-sh(i)],[8 1]);
        fvsi.faces= [1 2 3;1 3 4;1 2 6;1 6 5;2 3 6;3 7 6;1 4 8;1 8 5;3 4 7;4 7 8;5 6 7;5 7 8];
    else
        fv.faces=  [fv.faces;[1 2 3;1 3 4;1 2 6;1 6 5;2 3 6;3 7 6;1 4 8;1 8 5;3 4 7;4 7 8;5 6 7;5 7 8]+size(fv.vertices,1)];
        fv.vertices = [fv.vertices;(([1 -1 0;1 1 0;-1 1 0;-1 -1 0;1 -1 1;1 1 1;-1 1 1;-1 -1 1].*repmat([w/2 w/2 h(i)-sh(i)],[8 1]))+repmat([X(i) Y(i) 0],[8 1]))];
        fvss.faces=  [fvss.faces;[1 2 3;1 3 4;1 2 6;1 6 5;2 3 6;3 7 6;1 4 8;1 8 5;3 4 7;4 7 8;5 6 7;5 7 8]+size(fvss.vertices,1)];
        fvss.vertices = [fvss.vertices;(([1 -1 0;1 1 0;-1 1 0;-1 -1 0;1 -1 1;1 1 1;-1 1 1;-1 -1 1].*repmat([w/2 w/2 sh(i)],[8 1]))+repmat([X(i) Y(i) h(i)],[8 1]))];
        fvsi.faces=  [fvsi.faces;[1 2 3;1 3 4;1 2 6;1 6 5;2 3 6;3 7 6;1 4 8;1 8 5;3 4 7;4 7 8;5 6 7;5 7 8]+size(fvsi.vertices,1)];
        fvsi.vertices = [fvsi.vertices;(([1 -1 0;1 1 0;-1 1 0;-1 -1 0;1 -1 1;1 1 1;-1 1 1;-1 -1 1].*repmat([w/2 w/2 sh(i)],[8 1]))+repmat([X(i) Y(i) h(i)-sh(i)],[8 1]))];
    end
end
figure;
strsurf=patch(fv,'facecolor',options.facec,'edgecolor','none','tag','model0','facelighting','gouraud');
strsurfs=patch(fvss,'facecolor',[0 1 0],'edgecolor','none','tag','model0','facelighting','gouraud');
set(strsurfs,'FaceAlpha',[0.8]);
strsurfs=patch(fvsi,'facecolor',[0 0 1],'edgecolor','none','tag','model0','facelighting','gouraud');
set(strsurfs,'FaceAlpha',[0.8]);



hold on;
% for i = 1:6
%     line([fv.vertices(1:4,1); fv.vertices(1,1)],[fv.vertices(1:4,2); fv.vertices(1,2)],[fv.vertices(1:4,3); fv.vertices(1,3)],'Color',options.linec);
%     line([fv.vertices(5:8,1); fv.vertices(5,1)],[fv.vertices(5:8,2); fv.vertices(5,2)],[fv.vertices(5:8,3); fv.vertices(5,3)],'Color',options.linec);
%     line([fv.vertices(1,1); fv.vertices(5,1)],[fv.vertices(1,2); fv.vertices(1,2)],[fv.vertices(1,3); fv.vertices(1,3)],'Color',options.linec);
%     line([fv.vertices(4,1); fv.vertices(8,1)],[fv.vertices(4,2); fv.vertices(8,2)],[fv.vertices(4,3); fv.vertices(8,3)],'Color',options.linec);
%     line([fv.vertices(3,1); fv.vertices(7,1)],[fv.vertices(3,2); fv.vertices(7,2)],[fv.vertices(3,3); fv.vertices(7,3)],'Color',options.linec);
%     line([fv.vertices(2,1); fv.vertices(6,1)],[fv.vertices(2,2); fv.vertices(6,2)],[fv.vertices(2,3); fv.vertices(6,3)],'Color',options.linec);
% end
% fv.vertices = ([1 -1 0;1 1 0;-1 1 0;-1 -1 0;1 -1 1;1 1 1;-1 1 1;-1 -1 1].*repmat([w/2 w/2 sh],[8 1]))+repmat([posx posy h],[8 1]);
% fv.faces= [1 2 3;1 3 4;1 2 6;1 6 5;2 3 6;3 7 6;1 4 8;1 8 5;3 4 7;4 7 8;5 6 7;5 7 8];
% strsurf=patch(fv,'facecolor',[0 1 0],'edgecolor','none','tag','model0','facelighting','gouraud');
% 
% hold on;
% % for i = 1:6
% %     line([fv.vertices(1:4,1); fv.vertices(1,1)],[fv.vertices(1:4,2); fv.vertices(1,2)],[fv.vertices(1:4,3); fv.vertices(1,3)],'Color',options.linec);
% %     line([fv.vertices(5:8,1); fv.vertices(5,1)],[fv.vertices(5:8,2); fv.vertices(5,2)],[fv.vertices(5:8,3); fv.vertices(5,3)],'Color',options.linec);
% %     line([fv.vertices(1,1); fv.vertices(5,1)],[fv.vertices(1,2); fv.vertices(1,2)],[fv.vertices(1,3); fv.vertices(1,3)],'Color',options.linec);
% %     line([fv.vertices(4,1); fv.vertices(8,1)],[fv.vertices(4,2); fv.vertices(8,2)],[fv.vertices(4,3); fv.vertices(8,3)],'Color',options.linec);
% %     line([fv.vertices(3,1); fv.vertices(7,1)],[fv.vertices(3,2); fv.vertices(7,2)],[fv.vertices(3,3); fv.vertices(7,3)],'Color',options.linec);
% %     line([fv.vertices(2,1); fv.vertices(6,1)],[fv.vertices(2,2); fv.vertices(6,2)],[fv.vertices(2,3); fv.vertices(6,3)],'Color',options.linec);
% % end
% fv.vertices = ([1 -1 0;1 1 0;-1 1 0;-1 -1 0;1 -1 1;1 1 1;-1 1 1;-1 -1 1].*repmat([w/2 w/2 sh],[8 1]))+repmat([posx posy h-sh],[8 1]);
% fv.faces= [1 2 3;1 3 4;1 2 6;1 6 5;2 3 6;3 7 6;1 4 8;1 8 5;3 4 7;4 7 8;5 6 7;5 7 8];
% strsurf=patch(fv,'facecolor',[0 0 1],'edgecolor','none','tag','model0','facelighting','gouraud');
% set(strsurf,'FaceAlpha',[0.8]);
% hold on;
% for i = 1:6
%     line([fv.vertices(1:4,1); fv.vertices(1,1)],[fv.vertices(1:4,2); fv.vertices(1,2)],[fv.vertices(1:4,3); fv.vertices(1,3)],'Color',options.linec);
%     line([fv.vertices(5:8,1); fv.vertices(5,1)],[fv.vertices(5:8,2); fv.vertices(5,2)],[fv.vertices(5:8,3); fv.vertices(5,3)],'Color',options.linec);
%     line([fv.vertices(1,1); fv.vertices(5,1)],[fv.vertices(1,2); fv.vertices(1,2)],[fv.vertices(1,3); fv.vertices(1,3)],'Color',options.linec);
%     line([fv.vertices(4,1); fv.vertices(8,1)],[fv.vertices(4,2); fv.vertices(8,2)],[fv.vertices(4,3); fv.vertices(8,3)],'Color',options.linec);
%     line([fv.vertices(3,1); fv.vertices(7,1)],[fv.vertices(3,2); fv.vertices(7,2)],[fv.vertices(3,3); fv.vertices(7,3)],'Color',options.linec);
%     line([fv.vertices(2,1); fv.vertices(6,1)],[fv.vertices(2,2); fv.vertices(6,2)],[fv.vertices(2,3); fv.vertices(6,3)],'Color',options.linec);
% end
function Plot_Normals(Coord, Normals,headplot);
%
% Syntax :
% Plot_Normals(Coord, Normals,headplot);
%
% This function plots vectors over specified points
%
% Input Parameters:
%   Coord           : Points coordinates
%              
%   Normals         : Vectors to plot
%   headplot        : Boolean variable to plot or not the head of the
%                     vector
%
% Output Parameters:
%
%
% Related references:
%
%
% See also: Smooth_Surf Surf_Comp Red_Surf Plot_oversurf Atlas_Surf quiver3
% Exp_Surf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% December 30th 2008
% Version $1.0

%=====================Checking Input Parameters===========================%
if nargin<3
    headplot = 1;
end
%=====================End of Checking Input Parameters====================%
alpha = 0.33;
beta = 0.33;
scf = 0.1;
norma = sqrt(sum((Normals').^2))';
x = Coord(:,1); y = Coord(:,2); z = Coord(:,3);
u = Normals(:,1)./norma; v = Normals(:,2)./norma; w = Normals(:,3)./norma; u(isnan(u)) = 0;v(isnan(v)) = 0;w(isnan(w)) = 1;
if min(size(x))==1, n=sqrt(prod(size(x))); m=n; else [m,n]=size(x); end
delx = diff([min(x(:)) max(x(:))])/n;
dely = diff([min(y(:)) max(y(:))])/m;
delz = diff([min(z(:)) max(y(:))])/max(m,n);
del = sqrt(delx.^2 + dely.^2 + delz.^2);
if del>0
    len = sqrt((u/del).^2 + (v/del).^2 + (w/del).^2);
    maxlen = max(len(:));
else
    maxlen = 0;
end

if maxlen>0
    scf = scf*0.9 / maxlen;
else
    scf = scf*0.9;
end
col = abs([u v w]);
col = col./repmat(sqrt(sum((col').^2))',[1 3]);
u = u*scf; v = v*scf; w = w*scf;
beta = beta * sqrt(u.*u + v.*v + w.*w) ./ (sqrt(u.*u + v.*v) + eps);
hu = [x+u-alpha*(u+beta.*(v+eps)) x+u x+u-alpha*(u-beta.*(v+eps))];
hv = [y+v-alpha*(v-beta.*(u+eps)) y+v y+v-alpha*(v+beta.*(u+eps))];
hw = [z+w-alpha*w z+w z+w-alpha*w];
%figure;
for i = 1:size(u(:),1)
    line([x(i) x(i)+u(i)], [y(i) y(i)+v(i)],[z(i) z(i)+w(i)],'Color',col(i,:),'LineWidth',1);
    if headplot~=0
        line([hu(i,1) hu(i,2) hu(i,3)], [hv(i,1) hv(i,2) hv(i,3)],[hw(i,1) hw(i,2) hw(i,3)],'Color',col(i,:),'LineWidth',1);
    end
end
axis equal;view(3);

return
function OutCurve = my_3D_Fit(Coords,Np);

if nargin < 2
    Np = 100;
end



X = Coords(:,1);
Y = Coords(:,2);
Z = Coords(:,3);



% A = sqrt(cumsum([0;diff(X(:))].^2) + ...
%     cumsum([0;diff(Y(:))].^2) + ...
%     cumsum([0;diff(Z(:))].^2));
% 
% % Now, use interp1 on this, best with a
% % spline as the method.
% 
% a = linspace(A(1),A(end),1000);
% newxyz = interp1(A,[X(:),Y(:),Z(:)],a,'nearest');
% 
% plot3(X,Y,Z,'o')
% hold on
% plot3(newxyz(:,1),newxyz(:,2),newxyz(:,3),'-')
% 
% 
% 
% 
% 
% 
% 
% 
% % Xinterp = linspace(min(X),max(X),Np);
% % Yinterp = linspace(min(Y),max(Y),Np);
% % Zinterp = linspace(min(Z),max(Z),Np);
% 
 tInterp = linspace(1,length(X),Np);
% t = [1:length(X)]';
% 
% % Zinterp = foo([Xinterp(:) Yinterp(:)]);
% % 
% % Yinterp = foo(Xinterp);
% 
% 
% 
% Xinterp = linspace(min(X),max(X),100);
% Yinterp = linspace(min(Y),max(Y),100);
% Zinterp = linspace(min(Z),max(Z),100);
% 
% 
% 
% 
% 
% 
% 
% m=10; % number of data points
% n=3; % dimension
% P=rand(n,m);
% 
% F=spline((1:length(Coords)),Coords');
% 
% % Trajectory
% step=0.05;
% t=[1:step:length(Coords)];
% Ft=ppval(F,t);
% 
% plot3(Ft(1,:),Ft(2,:),Ft(3,:),'b', ...
%       Coords(:,1),Coords(:,1),Coords(:,3),'ro');
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
 sf = fit( [X, Y], Z, 'poly23', 'Robust', 'LAR' );
 Zp = sf(Xinterp,Yinterp);
% 
% sf = fit( [Y, Z], X, 'poly23', 'Robust', 'LAR' );
% Xp = sf(Xinterp,Yinterp);
% 
% sf = fit( [X, Z], Y, 'poly23', 'Robust', 'LAR' );
% Yp = sf(Xinterp,Yinterp);
% 
% 
% sf(Xinterp,Zinterp);


t = [1:length(X)]';
foo = fit(t, X, 'linearinterp');
Xinterp = foo(tInterp);

foo = fit(t, Y, 'linearinterp');
Yinterp = foo(tInterp);

foo = fit(t, Z, 'linearinterp');
Zinterp = foo(tInterp);


% 
% foo = fit(Y, Z, 'fourier2');
% Zy = foo(Yinterp);
% 
% 
% foo = fit(X, Y, 'fourier2');
% Yx = foo(Xinterp);
% 
% foo = fit(Z, Y, 'fourier2');
% Yz = foo(Zinterp);
% 
% 
% foo = fit(Y, X, 'fourier2');
% Xy = foo(Yinterp);
% 
% foo = fit(Z, X, 'fourier2');
% Xz = foo(Zinterp);

OutCurve = [Xinterp(:) Yinterp(:) Zinterp(:)];
% 
% OutCurve = [Xy(:)+Xz(:) Yx(:)+Yz(:) Zx(:)+Zy(:)];
return;


function Digit_image(Imfile);

Imfile  = '/media/UserBackup/Test.png';
idata = imread(Imfile);
figure;imshow(idata);

[x,y] = ginput(39);
x(1:end) =x(1:end)-x(1);
y(1:end) =y(1:end)-y(1);
%px =input('Entrar escala en X ');
px = .2;
%py =input('Entrar escala en Y ');
py = .1;

x(:) = x(:)*px/x(2);
y(:) = y(:)*py/y(2);
return

clear;
clc;

%% 球面グリッドの定義
theta = 0:pi/10:pi;                   % polar angle
phi = 0:pi/5:2*pi;                   % azimuth angle

[phi,theta] = meshgrid(phi,theta);    % define the grid

%% 球面調和関数の計算
%半径 5 の球の表面上で、度数 6、次数 1、振幅 0.5
degree = 1; %度数
order = 1; %次数
amplitude = 0.5; %振幅0.5
radius = 5; %半径

Ymn = legendre(degree,cos(theta(:,1)));
Ymn = Ymn(order+1,:)';
yy = Ymn;

for kk = 2: size(theta,1)
    yy = [yy Ymn];
end

yy = yy.*cos(order*phi);

order = max(max(abs(yy)));
rho = radius + amplitude*yy/order;

r = rho.*sin(theta);    % convert to Cartesian coordinates
x = r.*cos(phi);
y = r.*sin(phi);
z = rho.*cos(theta);

%% 球の表面での球面調和関数のプロット
figure
s = surf(x,y,z);

light               % add a light
lighting gouraud    % preferred lighting for a curved surface
axis equal off      % set axis equal and remove axis
view(40,30)         % set viewpoint
camzoom(1.5)        % zoom into scene

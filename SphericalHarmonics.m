
clear;
clc;

%% ���ʃO���b�h�̒�`
theta = 0:pi/10:pi;                   % polar angle
phi = 0:pi/5:2*pi;                   % azimuth angle

[phi,theta] = meshgrid(phi,theta);    % define the grid

%% ���ʒ��a�֐��̌v�Z
%���a 5 �̋��̕\�ʏ�ŁA�x�� 6�A���� 1�A�U�� 0.5
degree = 1; %�x��
order = 1; %����
amplitude = 0.5; %�U��0.5
radius = 5; %���a

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

%% ���̕\�ʂł̋��ʒ��a�֐��̃v���b�g
figure
s = surf(x,y,z);

light               % add a light
lighting gouraud    % preferred lighting for a curved surface
axis equal off      % set axis equal and remove axis
view(40,30)         % set viewpoint
camzoom(1.5)        % zoom into scene

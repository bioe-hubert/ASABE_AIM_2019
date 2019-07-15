% ------------------------------------------------------------------------------
%
% Dipole-array solution of steady saturated flow in heterogeneous medium
%
% configuration 1.
%
% ------------------------------------------------------------------------------

% build array of dipoles
rs 	  = 50;
[sx0,sy0] = meshgrid(-rs:4:rs,[-rs:4:rs]');
[sx1,sy1] = meshgrid(-(rs+2):4:(rs+2),[-(rs+2):4:(rs+2)]');
sz 	  = [sx0(:) ; sx1(:)] + [sy0(:) ; sy1(:)]*i;
sz 	  = sz(find(abs(sz) <= (rs*1.001)));

% plot dipole array with area of interest
figure
plot(real(sz),imag(sz),'k.')
hold on
tht = linspace(0,2*pi,101);
plot(rs*cos(tht),rs*sin(tht),'k')
plot([-2.25 2.25 2.25 -2.25 -2.25],[-0.25 -0.25 4.25 4.25 -0.25],'k')
hold off; axis equal; axis tight
title('Dipole Locations')

% compute dipole function around each dipole
[x,y] = meshgrid(linspace(-2.5,2.5,251),linspace(-0.5,4.5,251));
z = x + i*y;

% add contributions of all dipoles
omeg = -z;
for src = 1:length(sz)
    omeg = omeg - 1./(z - sz(src));
end

% set vector of values for contour plot
v = linspace(-4,4,49);

% set apparent porosity of medium
poros = 1 - 2*pi/16;

% compute flow velocities and adjust for zone inside circlar heterogeneities
[aqx,aqy] = gradient(real(omeg),5/250);
qx = -aqx(26:25:end-25,26:25:end-25);
qy = -aqy(26:25:end-25,26:25:end-25);
xq = x(26:25:end-25,26:25:end-25);
yq = y(26:25:end-25,26:25:end-25);
qx(find(abs(qx)>100)) = 0;
qy(find(abs(qy)>100)) = 0;
qx(find(xq.^2+yq.^2 < 1)) = 0;
qy(find(xq.^2+yq.^2 < 1)) = 0;
qx(find((xq-2).^2+(yq-2).^2 < 1)) = 0;
qy(find((xq-2).^2+(yq-2).^2 < 1)) = 0;
qx(find((xq+2).^2+(yq-2).^2 < 1)) = 0;
qy(find((xq+2).^2+(yq-2).^2 < 1)) = 0;
qx(find(xq.^2+(yq-4).^2 < 1)) = 0;
qy(find(xq.^2+(yq-4).^2 < 1)) = 0;

% plot velocities
figure
contour(x,y,real(omeg),v)
hold on
contour(x,y,imag(omeg),v,'r:')
quiver(xq,yq,qx,qy,0.75)
tht = linspace(0,2*pi,101)';
crc = [cos(tht) sin(tht)];
grey = [0.5 0.5 0.5];
fill(crc(:,1),crc(:,2),grey,crc(:,1)+2,crc(:,2)+2,grey, ...
    crc(:,1)-2,crc(:,2)+2,grey,crc(:,1),crc(:,2)+4,grey)
hold off; xlabel('x'); ylabel('y'); axis equal
axis([-2 2 0 4])
title('Flow Velocities')

% compute mean-removed flow velocities and adjust for zone inside heterogeneities
[aqx,aqy] = gradient(real(omeg+poros*z),5/250);
qx = -aqx(26:25:end-25,26:25:end-25);
qy = -aqy(26:25:end-25,26:25:end-25);
xq = x(26:25:end-25,26:25:end-25);
yq = y(26:25:end-25,26:25:end-25);
qx(find(xq.^2+yq.^2 <= 1)) = -poros;
qy(find(xq.^2+yq.^2 <= 1)) = 0;
qx(find((xq-2).^2+(yq-2).^2 <= 1)) = -poros;
qy(find((xq-2).^2+(yq-2).^2 <= 1)) = 0;
qx(find((xq+2).^2+(yq-2).^2 <= 1)) = -poros;
qy(find((xq+2).^2+(yq-2).^2 <= 1)) = 0;
qx(find(xq.^2+(yq-4).^2 <= 1)) = -poros;
qy(find(xq.^2+(yq-4).^2 <= 1)) = 0;

% plot mean-removed velocities
figure
yy   = linspace(0,2,101);
xx   = -1 - 0.27*cos(pi*yy/2);
xx2  = -1 - 0.285*cos(pi*yy/2);
xx3  = -1 - 0.29*cos(pi*yy/2);
xx4  = [xx2 -fliplr(xx2)];
yy4  = [yy fliplr(yy)];
gray = 0.9*[1 1 1];
fill(xx4,yy4,gray,2+xx4,2+yy4,gray,xx4-2,2+yy4,gray);
hold on
tht = linspace(0,2*pi,101)';
crc = [cos(tht) sin(tht)];
plot(crc(:,1),crc(:,2),'k',crc(:,1)+2,crc(:,2)+2,'k', ...
    crc(:,1)-2,crc(:,2)+2,'k',crc(:,1),crc(:,2)+4,'k')
quiver(xq,yq,qx,qy,0.75)
hold off; xlabel('x'); ylabel('y'); axis equal
axis([-2 2 0 4])
title('Mean-Removed Flow Velocities')





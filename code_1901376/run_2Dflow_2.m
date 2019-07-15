% ------------------------------------------------------------------------------
%
% Dipole-array solution of steady saturated flow in heterogeneous medium
%
% configuration 2.
%
% ------------------------------------------------------------------------------

% build array of dipoles
Lx = 3; Ly = 3;
rsx = 61*Lx; rsy = 61*Ly;
[sx0,sy0] = meshgrid(-rsx:Lx:rsx,[-rsy:Ly:rsy]');
sz = sx0(:) + sy0(:)*i;
rs = 61*max(Lx,Ly);
sz = sz(find(abs(sz - (Lx+Ly*i)/2) <= (rs*1.001)));

% plot dipole array
figure
plot(real(sz),imag(sz),'k.')
hold on
tht = linspace(0,2*pi,101);
plot(rs*cos(tht)+Lx/2,rs*sin(tht)+Ly/2,'k')
plot([-0.25 Lx+0.25 Lx+0.25 -0.25 -0.25],[-0.25 -0.25 Ly+0.25 Ly+0.25 -0.25],'k')
hold off; axis equal; axis tight
title('Dipole Locations')

% compute dipole function around each dipole
dx = Lx/8; dy = Ly/8; numsrc = 251;
[x,y] = meshgrid(linspace(-dx,Lx+dx,numsrc),linspace(-dy,Ly+dy,numsrc));
z = x + i*y;

% add contributions of all dipoles
omeg = -z;
for src = 1:length(sz)
    omeg = omeg - 1./(z - sz(src));
end

% set apparent porosity of medium
poros = 1 - pi/(Lx*Ly);

% compute flow velocities and adjust for zone inside circlar heterogeneities
[aqx,aqy] = gradient(real(omeg),(Lx+2*dx)/250);
xq = x(26:25:end-25,26:25:end-25);
yq = y(26:25:end-25,26:25:end-25);
qx = -aqx(26:25:end-25,26:25:end-25);
qy = -aqy(26:25:end-25,26:25:end-25);
qx(find(abs(qx)>100)) = 0;
qy(find(abs(qy)>100)) = 0;
qx(find(xq.^2+yq.^2 < 1)) = 0;
qy(find(xq.^2+yq.^2 < 1)) = 0;
qx(find((xq-Lx).^2+yq.^2 < 1)) = 0;
qy(find((xq-Lx).^2+yq.^2 < 1)) = 0;
qx(find(xq.^2+(yq-Ly).^2 < 1)) = 0;
qy(find(xq.^2+(yq-Ly).^2 < 1)) = 0;
qx(find((xq-Lx).^2+(yq-Ly).^2 < 1)) = 0;
qy(find((xq-Lx).^2+(yq-Ly).^2 < 1)) = 0;

% plot velocities
figure
dv  = (real(omeg(126,226)) - real(omeg(126,26))) / 20;
v   = linspace(real(omeg(126,26))-2*dv,real(omeg(126,226))+2*dv,25);
contour(x,y,real(omeg),v); hold on
vv  = imag(omeg(:,126));
dv  = (max(vv) - min(vv))/20;
nv1 = round(imag(omeg(51,26)) / dv);
nv2 = round(abs(min(vv)) / dv);
nv  = max(nv1,nv2);
v   = max(vv)+ [-nv*dv:dv:nv*dv];
contour(x,y,imag(omeg),v,'r:')
tht = linspace(0,2*pi,101)';
crc = [cos(tht) sin(tht)];
grey = 0.9*[1 1 1];
fill(crc(:,1),crc(:,2),grey,crc(:,1)+Lx,crc(:,2)+Ly,grey, ...
    crc(:,1)+Lx,crc(:,2),grey,crc(:,1),crc(:,2)+Ly,grey)
quiver(xq,yq,qx,qy,0.75)
hold off; xlabel('x'); ylabel('y'); axis equal
axis([0 Lx 0 Ly])
title('Flow Velocities')

% compute mean-removed flow velocities and adjust for zone inside heterogeneities
[aqx,aqy] = gradient(real(omeg+poros*z),(Lx+2*dx)/250);
xq = x(26:25:end-25,26:25:end-25);
yq = y(26:25:end-25,26:25:end-25);
qx = -aqx(26:25:end-25,26:25:end-25);
qy = -aqy(26:25:end-25,26:25:end-25);
qx(find(abs(qx)>100)) = 0;
qy(find(abs(qy)>100)) = 0;
qx(find(xq.^2+yq.^2 < 1)) = -poros;
qy(find(xq.^2+yq.^2 < 1)) = 0;
qx(find(xq.^2+(yq-Ly).^2 < 1)) = -poros;
qy(find(xq.^2+(yq-Ly).^2 < 1)) = 0;
qx(find((xq-Lx).^2+yq.^2 < 1)) = -poros;
qy(find((xq-Lx).^2+yq.^2 < 1)) = 0;
qx(find((xq-Lx).^2+(yq-Ly).^2 < 1)) = -poros;
qy(find((xq-Lx).^2+(yq-Ly).^2 < 1)) = 0;

% plot mean-removed velocities
figure
xb = linspace(-Lx/2,Lx/2,51);
yb = 1.336 - 0.1788*xb.^2 - 0.0448*xb.^4;
xc = linspace(0.971,Lx/2,21);
yc =  -0.1111 + 0.3732*xc.^2 + -0.0043*xc.^4;
xd = [-0.971 0.971];
yd = [0.237 0.237];
xx = [xd xc fliplr(xb)];
yy = [yd yc fliplr(yb)];
fill(xx,yy,grey,xx,3-yy,grey,3-xx,yy,grey,3-xx,3-yy,grey)
hold on
tht = linspace(0,2*pi,101)';
crc = [cos(tht) sin(tht)];
plot(crc(:,1),crc(:,2),'k',crc(:,1)+Lx,crc(:,2)+Ly,'k', ...
    crc(:,1)+Lx,crc(:,2),'k',crc(:,1),crc(:,2)+Ly,'k')
quiver(xq,yq,qx,qy,0.75,'b')
hold off; xlabel('x'); ylabel('y'); axis equal
axis([0 Lx 0 Ly])
title('Mean-Removed Flow Velocities')



function [t,XM] = ode45rk(F, tv, X0)
% Runge-Kutta Solution of a System of ODEs
%   F  = derivative function
%   tv = time vector
%   X0 = initial values
t       = tv(:);				% independent var output vector
dt      = t(2)-t(1);				% subinterval spacing
XM      = zeros(length(t),length(X0));		% initial solution matrix
XM(1,:) = X0;					% set initial values
% loop over independent variable values
for k = 1:length(t)-1
  tk   = t(k);					% current value of indep var
  Xk   = XM(k,:)';				% current dependent var vector
  dx1  = F(tk,      Xk);			% Runge-Kutta estimate 1
  dx2  = F(tk+dt/2, Xk+dx1*dt/2);		% Runge-Kutta estimate 2
  dx3  = F(tk+dt/2, Xk+dx2*dt/2);		% Runge-Kutta estimate 3
  dx4  = F(tk+dt,   Xk+dx3*dt);			% Runge-Kutta estimate 4
  XM(k+1,:) = XM(k,:) + (dx1+2*dx2+2*dx3+dx4)'*dt/6; % next solution value
end


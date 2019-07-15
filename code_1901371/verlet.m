function [t,XM] = verlet(F,tv,M,X,V)
% Verlet Symplectic Integrator (order 2)
%	F=force function, tv=time vector, M=mass vector
%	X=initial position matrix, V=initial velocity matrix 
% preliminary computations
p  = size(X,1);			% number of bodies
d  = size(X,2);			% number of spatial dimensions (2 or 3)
Md = repmat(M,1,d);		% replicated mass matrix
t  = tv(:);			% output time-vector
dt = t(2) - t(1);		% time step
A  = Ftot(F,M,X)./Md;		% initial acceleration
XM = zeros(length(t),d,p);	% initial solution matrix
XM(1,:,:) = X';			% store initial positions in matrix
% solution for successive times
for k = 2:length(t)
  X  = X + V*dt + 0.5*A*dt*dt;	% new positions
  XM(k,:,:) = X';		% store new positions in matrix
  An = Ftot(F,M,X)./Md;		% new accelerations
  V  = V + (An+A)*dt/2;		% new velocities
  A  = An;			% accelerations for next iter
end
% nested function to compute the force matrix
  function fout = Ftot(F,M,X)
    R    = repmat(reshape(X,p,1,d),[1,p,1])-repmat(reshape(X,1,p,d),[p,1,1]);
    dist = sqrt(sum(R.^2,3));	% distancess between all bodies
    Fod  = F(M*M',dist)./dist;	% forces between all bodies
    Fod(1:p+1:p*p) = 0;		% remove "self"-force
    fout = squeeze(sum(R.*repmat(Fod,[1,1,d]),2)); % total force on each body
  end
end


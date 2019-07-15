function [x,t,U] = FDM_sol(kvD,xv,tv,U0)
%	FDM_sol: Finite Difference solution of 1-D avection-difusion PDE
%	kvD    = vector of reaction rate, velocity and diffusion coefficient
%	xv, tv = vectors of spatial nodes and times for solution
%	U0     = initial conditions (IC) vector
x    = xv(:);				% spatial  node positions
t    = tv(:)';				% temporal node positions
nx   = length(x);			% number of spatial  nodes
nt   = length(t);			% number of temporal nodes
dx   = x(2)-x(1);			% spatial  spacing
dt   = t(2)-t(1);			% temporal spacing
lsd  =  kvD(2)/(2*dx) + kvD(3)/dx^2;	% left-side  spatial difference factor
csd  =  kvD(1) - 2*kvD(3)/dx^2;		% central    spatial difference factor
rsd  = -kvD(2)/(2*dx) + kvD(3)/dx^2;	% right-side spatial difference factor
SMAT = spdiags(repmat([lsd csd rsd],nx,1),[-1 0 1],nx,nx); % spatial matrix
SMAT(1,2)       = SMAT(1,2)+lsd;	% introduce BC at left-edge
SMAT(end,end-1) = SMAT(end,end-1)+rsd;	% introduce BC at right-edge
LHS  = speye(nx) - SMAT*dt/2;		% fully discrete left-hand-side  matrix
RHS  = speye(nx) + SMAT*dt/2;		% fully discrete right-hand-side matrix
U    = zeros(nx,nt);			% solution matrix (all nodes)
U(:,1) = U0;				% introduce IC into solution
% compute solution at successive times
for k = 2:nt
  U(:,k) = LHS\(RHS*U(:,k-1));		% solve for time-step k
end


% ------------------------------------------------------------------------------
%
% Demonstrate solution of 1-D advection-diffusion equation by Finite Differences
%
% ------------------------------------------------------------------------------

% define problem parameters ---
kvD = [-0.01 0.3 0.7];		% reaction rate, velocity and diffusion coefficient
nx  = 101;			% number of spatial nodes
nt  = 201;			% number of temporal nodes
xv  = linspace(0,10,nx);	% vector of spatial nodes
tv  = linspace(0,8,nt);		% vector of temporal nodes
U0  = zeros(nx,1);		% initial condition (IC) vector (all zeros)
dx  = xv(2)-xv(1);		% spacing between nodes
U0(round(1+4/dx)) = 3/dx;	% IC vector updated with Dirac pulse
[x,t,U] = FDM_sol(kvD,xv,tv,U0);% solve
% plot ------------------------
plot(x,U(:,1:40:end))
xlabel('x (units)'); ylabel('Concentration (units)')
axis([0 10 0 2]); legend([repmat('t = ',6,1) num2str(tv(1:40:end)',3)])




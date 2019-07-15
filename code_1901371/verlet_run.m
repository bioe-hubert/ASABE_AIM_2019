% ------------------------------------------------------------------------------
%
% Demonstrate solution of Hamiltonian system by Verlet method in 2-D
%
% ------------------------------------------------------------------------------

% masses and initial conditions for earth and one moon
tv = linspace(0,2.592e6,101);	% solution times (30 days, in seconds)
M  = [5.97e24 ; 7.35e22];	% masses of earth and moon (kg)
X0 = [ 0 0 ; 378000e3 0];	% initial positions (m)
V0 = [ 0 0 ; 0   1.02e3];	% initial velocities (m/s)

% force function (gravitational force)
F = @(MM,dist) -6.674e-11*MM./dist.^2;   % gravitational force (function)

% adjustment of earth velocity for zero total momentum (no space drift)
V0(1,:) = -sum(repmat(M,1,2).*V0)/M(1); % earth veloc adj for zero momentum

% solution (computation of positions)
[t,XM] = verlet(F,tv,M,X0,V0);

% plot results
XM = XM/1000;			% use Km as spatial units
plot(XM(:,1,1),XM(:,2,1),'b',XM(:,1,2),XM(:,2,2),'g')
axis equal
axis([-1 1 -1 1]*1.1*378000e3/1000)
xlabel('x (km)'); ylabel('y (km)'); legend('Earth','Moon')




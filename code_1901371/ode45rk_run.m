% ------------------------------------------------------------------------------
%
% Demonstrate use of Runge-Kutta method to solve Lotka-Volterra ODEs
%
% ------------------------------------------------------------------------------

% problem definition (ODEs, time vector, initial conditions)
F  = @(t,X) [0.5*X(1) - 0.01*prod(X) ; -0.4*X(2) + 0.001*prod(X)];
tv = linspace(0,100,501)';
X0 = [ 600 ; 10 ] ;

% solve the ODEs
[t,XM] = ode45rk(F,tv,X0);

% plot results
subplot(2,2,1)
plot(t,XM(:,1),'g',t,XM(:,2),'r')
legend('H','P'); xlabel('Time (years)'); ylabel('Population')
subplot(2,2,2)
plot(XM(:,1),XM(:,2))
xlabel('Herbivore Population'); ylabel('Predator Population')




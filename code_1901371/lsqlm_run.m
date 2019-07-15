% ------------------------------------------------------------------------------
%
% Demonstrate use of Levenberg-Marquardt to fit Langmuir isotherm to data
%
% ------------------------------------------------------------------------------

% adsorption data
c = [ 0.162 0.324 0.649 1.297 1.622 2.757 3.730 7.297 ]' ;
x = [ 158 200 265 366 417 500 583 692 ]' ;

% Langmuir model error function (returns residuals of fit)
E = @(U) U(2)*U(1)*c ./ (1+U(1)*c)  -  x;

% perform fit and show parameters
U     = lsqlm( E , [ 126 ; 1 ] );
K     = U(1)
alpha = U(2)

% plot result
cg = linspace(0,8,101);
P  = alpha*K*cg./(1+K*cg);
plot(c,x,'ro',cg,P,'b')
xlabel('c (\mu g/ml)'); ylabel('x (\mu g /g)')
text(3,300,'Isotherm Equation:')
text(3,250,['x = ' num2str(alpha,4) ' * ' num2str(K,4) ' c / (1 + ' num2str(K,4) ' c)'])



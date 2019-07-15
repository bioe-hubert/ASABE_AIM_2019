
% ------------------------------------------------------------------------------
%
% Demonstrate relationship between coefficient of determination
% and coefficient of correlation using streamflow hydrograph data.
%
% General case (eg. nonlinear regression):	R^2 < r^2
%
% Special case of linear regression:		R^2 = r^2
%
% ------------------------------------------------------------------------------

format compact

% streamflow hydrograph observations (mm/h) sampled every 4 hours
O   = [0.046523 0.122466 0.224406 0.296586 0.472074 0.677324 0.756003 ...
       0.735478 0.646536 0.543912 0.454970 0.379712 0.320189 0.279482 ...
       0.248010 0.220301 0.198408 0.181304 0.166936 0.157016 0.145727 ...
       0.137517 0.130676 0.120413]';
n   = length(O);		% number of observations
T   = 4*[0:n-1]';		% times at which observations were taken
% compute stats ----------------
mO  = mean(O)			% mean of observations
o   = O - mO;			% mean-removed observations
so2 = (o'*o)/(n-1)		% variance of observations
% plot ------------------------
figure(1); subplot(2,2,1); plot(T,O,'ro')
xlabel('Time (h)'); ylabel('O (mm/h)')
text(60,0.75,['M(O) = ' num2str(mO,4)])
text(60,0.65,['\sigma_O^2 = ' num2str(so2,4)])

% ------------------------------------------------------------------------------
%
% Use a line model: P = U1 + U2*T (linear least-square-error regression)
%
% ------------------------------------------------------------------------------
disp('------------- line model ------------------------------------')
mT  = mean(T)			% mean of independent variable (T)
t   = T - mT;			% mean-removed times
st2 = (t'*t)/(n-1)		% variance of times
sot = (o'*t)/(n-1)		% covariance between times and observations
U1  = mO - sot/st2*mT		% intercept of line model
U2  = sot/st2			% slope of least-square-error line model
P   = U1 + U2*T;		% model predictions
E   = P - O;			% model errors
% stats -----------------------
mP  = mean(P)			% mean prediction
p   = P - mP;			% mean-removed predictions
sp2 = (p'*p)/(n-1)		% variance of predictions
spo = (p'*o)/(n-1)		% covariance between predictions and observations
mE  = mean(E)			% mean error
e   = E - mE;			% mean-removed errors
se2 = (e'*e)/(n-1)		% variance of errors
R2  = 1 - sum((P-O).^2)/sum((O-mO).^2)	% coefficient of determination
r   = spo/sqrt(sp2*so2)			% coefficient of correlation
r2  = r^2				% coefficient of correlation squared
% plot ------------------------
TP  = linspace(0,100,201)';
figure(1); subplot(2,2,2); plot(T,O,'ro',TP,U1+U2*TP,'b')
xlabel('Time (h)'); ylabel('P and O (mm/h)'); % legend('O','P')
text(60,0.75,['M(P) = ' num2str(mP,4)],'fontsize',12)
text(60,0.65,['\sigma_p^2 = ' num2str(sp2,4)])
text(60,0.55,['\sigma_{po} = ' num2str(spo,4)])
text(60,0.45,['M(E) = ' num2str(mE,4)])
text(60,0.35,['\sigma_e^2 = ' num2str(se2,4)])
text(25,0.18,['R^2 = ' num2str(R2,4)])
text(25,0.08,['r^2 = ' num2str(r2,4)])
if U2 >= 0
  title(['P = ' num2str(U1,3) ' +  ' num2str(U2,3) ' T'])
else
  title(['P = ' num2str(U1,3) ' -  ' num2str(abs(U2),3) ' T'])
endif

% ------------------------------------------------------------------------------
%
% Use a dimensionless unit hydrograph model (nonlinear least squares)
%
%  parameters:   Qp = U1, tp = U2, m = U3
%  model:        P = Qp * exp(m ) * (T/tp)^m  * exp( -m*T/tp)
%                  = U1 * exp(U3) * (T/U2)^U3 * exp(-U3*T/U2)
%
% ------------------------------------------------------------------------------
disp('------------- original model ------------------------------------')
Qrs = @(U) U(1) * exp(U(3)) * (T/U(2)).^U(3) .* exp(-U(3)*T/U(2)) - O;
U   = lsqlm(Qrs, [0.8,30,4]')					% fit model to data
P   = U(1) * exp(U(3)) * (T/U(2)).^U(3) .* exp(-U(3)*T/U(2));	% model predictions
E   = P - O;							% prediction errors
% stats -----------------------
mP  = mean(P)			% mean prediction
p   = P - mP;			% mean-removed prediction
sp2 = (p'*p)/(n-1)		% variance of predictions
spo = (p'*o)/(n-1)		% covariance between predictions and observations
mE  = mean(E)			% mean error
e   = E - mE;			% mean-removed error
se2 = (e'*e)/(n-1)		% variance of errors
R2  = 1 - sum((P-O).^2)/sum((O-mO).^2)	% coefficient of determination
r   = spo/sqrt(sp2*so2)			% coefficient of correlation
r2  = r^2				% coefficient of correlation squared
% plot ------------------------
TP  = linspace(0,100,201)';
PTP = U(1) * (TP/U(2)).^U(3) .* exp(U(3)*(1-TP/U(2)));
figure(1); subplot(2,2,3); plot(T,O,'ro',TP,PTP,'b')
xlabel('Time (h)'); ylabel('P and O (mm/h)');
text(60,0.75,['M(P) = ' num2str(mP,4)],'fontsize',12)
text(60,0.65,['\sigma_p^2 = ' num2str(sp2,4)])
text(60,0.55,['\sigma_{po} = ' num2str(spo,4)])
text(60,0.45,['M(E) = ' num2str(mE,4)])
text(60,0.35,['\sigma_e^2 = ' num2str(se2,4)])
text(25,0.18,['R^2 = ' num2str(R2,4)])
text(25,0.08,['r^2 = ' num2str(r2,4)])
title(['P = ' num2str(U(1),3) ' ( T / ' num2str(U(2),3) ' )^{' num2str(U(3),3) '} ' ...
       'e^{ [ ' num2str(U(3),3) ' ( 1 - T / ' num2str(U(2),3) ' ) ] }' ])
% regression line -------------
aL  = mP - spo/so2*mO		% intercept
bL  = spo/so2			% slope
OL  = linspace(0,0.8,11)';
L   = aL + bL * OL;
subplot(2,2,4); plot(O,P,'bs',OL,L,'k:')
xlabel('O (mm/h)'); ylabel('P (mm/h)')
title(['Regression Line: L = ' num2str(aL,3) ' +  ' num2str(bL,3) ' O'])
text(0.5,0.18,['R^2 = ' num2str(R2,4)])
text(0.5,0.08,['r^2 = ' num2str(r2,4)])

% ------------------------------------------------------------------------------
%
% Use a linear adjustment of fitted predictions to account for baseflow
%
% ------------------------------------------------------------------------------
disp('------------- adjusted model ------------------------------------')
a   = mO - r*sqrt(so2/sp2)*mP	% intercept
b   = r*sqrt(so2/sp2)		% slope
M   = a + b*P;			% adjusted model output (predictions)
E   = M - O;			% adjusted model errors
% stats -----------------------
mM  = mean(M)			% mean of adjusted predictions
m   = M - mM;			% mean-removed adjusted predictions
sm2 = (m'*m)/(n-1)		% variance of adjusted predictions
smo = (m'*o)/(n-1)		% covariance between adjusted predictions and observations
mE  = mean(E)			% mean error
e   = E - mE;			% mean-removed error
se2 = (e'*e)/(n-1)		% variance of errors
R2  = 1 - sum((M-O).^2)/sum((O-mO).^2)	% coefficient of determination
r   = smo/sqrt(sm2*so2)			% coefficient of correlation
r2  = r^2				% coefficient of correlation squared
% plot ------------------------
TP  = linspace(0,100,201)';
MTP = a + b*PTP;
figure(2); subplot(2,2,1); plot(T,O,'ro',TP,MTP,'b')
xlabel('Time (h)'); ylabel('M and O (mm/h)');
text(60,0.75,['M(M) = ' num2str(mM,4)])
text(60,0.65,['\sigma_m^2 = ' num2str(sm2,4)])
text(60,0.55,['\sigma_{mo} = ' num2str(smo,4)])
text(60,0.45,['M(E) = ' num2str(mE,4)])
text(60,0.35,['\sigma_e^2 = ' num2str(se2,4)])
text(25,0.18,['R^2 = ' num2str(R2,4)])
text(25,0.08,['r^2 = ' num2str(r2,4)])
title(['M = ' num2str(a,3) ' + '  num2str(b,3) ' * ' ...
        num2str(U(1),3) ' ( T / ' num2str(U(2),3) ' )^{' num2str(U(3),3) '} ' ...
       'e^{ [ ' num2str(U(3),3) ' ( 1 - T / ' num2str(U(2),3) ' ) ] }' ])
% regression line -------------
aL  = (1 - r2)*mO		% intercept
bL  = r2			% slope
aL  = mM - smo/so2*mO		% intercept -- check
bL  = smo/so2			% slope -- check
OL  = linspace(0,0.8,11)';
L   = aL + bL * OL;
subplot(2,2,2); plot(O,M,'bs',OL,L,'k:')
xlabel('O (mm/h)'); ylabel('M (mm/h)')
title(['Regression Line: L = ' num2str(aL,3) ' +  ' num2str(bL,3) ' O'])
text(0.5,0.18,['R^2 = ' num2str(R2,4)])
text(0.5,0.08,['r^2 = ' num2str(r2,4)])

% ------------------------------------------------------------------------------
%
% dimensionless unit hydrograph model with baseflow
%  parameters:   Qp = U1, tp = U2, m = U3, a = U4, b = U5
%  model:        P = Qp * exp(m ) * (T/tp)^m  * exp( -m*T/tp) +  a +  b*T
%                  = U1 * exp(U3) * (T/U2)^U3 * exp(-U3*T/U2) + U4 + U5*T
%
% ------------------------------------------------------------------------------
disp('------------- modified model ------------------------------------')
Qrs = @(U) U(1)*exp(U(3))*(T/U(2)).^U(3) .* exp(-U(3)*T/U(2)) + U(4)+U(5)*T - O;
U   = lsqlm(Qrs, [0.8,30,4,0.05,0.01]')					% fit model to data
P   = U(1)*exp(U(3))*(T/U(2)).^U(3).*exp(-U(3)*T/U(2))+U(4)+U(5)*T;	% model predictions
E   = P - O;								% errors
% stats -----------------------
mP  = mean(P)			% mean of predictions
p   = P - mP;			% mean-removed predictions
sp2 = (p'*p)/(n-1)		% variance of predictions
spo = (p'*o)/(n-1)		% covariance between predictions and observations
mE  = mean(E)			% mean error
e   = E - mE;			% mean-removed errors
se2 = (e'*e)/(n-1)		% variance of errors
R2  = 1 - sum((P-O).^2)/sum((O-mO).^2)	% coefficient of determination
r   = spo/sqrt(sp2*so2)			% coefficient of correlation
r2  = r^2				% coefficient of correlation squared
% plot ------------------------
TP  = linspace(0,100,201)';
PTP = U(1) * exp(U(3)) * (TP/U(2)).^U(3) .* exp(-U(3)*TP/U(2)) + U(4) + U(5)*TP;
figure(2); subplot(2,2,3); plot(T,O,'ro',TP,PTP,'b')
xlabel('Time (h)'); ylabel('P and O (mm/h)');
text(60,0.75,['M(P) = ' num2str(mP,4)],'fontsize',12)
text(60,0.65,['\sigma_p^2 = ' num2str(sp2,4)])
text(60,0.55,['\sigma_{po} = ' num2str(spo,4)])
text(60,0.45,['M(E) = ' num2str(mE,4)])
text(60,0.35,['\sigma_e^2 = ' num2str(se2,4)])
text(25,0.18,['R^2 = ' num2str(R2,4)])
text(25,0.08,['r^2 = ' num2str(r2,4)])
title(['P = ' num2str(U(1),2) ' ( T / ' num2str(U(2),2) ' )^{' num2str(U(3),2) '} ' ...
       'e^{ [ ' num2str(U(3),2) ' ( 1 - T / ' num2str(U(2),2) ' ) ] }' ...
       ' + ' num2str(U(4),2) ' + ' num2str(U(5),2) ' T'])
% regression line -------------
aL  = mP - spo/so2*mO		% intercept
bL  = spo/so2			% slope
OL  = linspace(0,0.8,11)';
L   = aL + bL * OL;
subplot(2,2,4); plot(O,P,'bs',OL,L,'k:')
xlabel('O (mm/h)'); ylabel('P (mm/h)')
title(['Regression Line: L = ' num2str(aL,3) ' +  ' num2str(bL,3) ' O'])
text(0.5,0.18,['R^2 = ' num2str(R2,4)])
text(0.5,0.08,['r^2 = ' num2str(r2,4)])



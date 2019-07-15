% ------------------------------------------------------------------------------
%
% Demonstrate training and application of feed-forward ANN
%
% ------------------------------------------------------------------------------

% dataset to be interpolated
x = 	[-1.0  -0.8  -0.6  -0.4  -0.2  0.0  0.2  0.4  0.6  0.8  1.0 ];
y = 	[-0.95 -0.92 -0.85 -0.77 -0.51 0.01 0.09 0.19 0.26 0.35 0.42];

% training
alpha   = 0.05;
netgeom = [4 5 1];
X       = [x ; x.^2 ; x.^3 ; x.^4 ; x.^5];
Y       = y;
[W,rmse] = ANN_train(alpha,netgeom,X,Y);

% plot dataset and training result
subplot(2,2,1); plot(x,y,'ro')
axis([-1.5 1.5 -1 1]); xlabel('x'); ylabel('y')
subplot(2,2,2); semilogy(rmse)
axis([0 10000 1e-3 10]); xlabel('Iteration'); ylabel('RMSE')

% application
xa = linspace(-1.5,1.5,101);
Xa = [xa ; xa.^2 ; xa.^3 ; xa.^4 ; xa.^5];
Ya = ANN_apply(W,Xa);

% plot application result (interpolation of dataset)
subplot(2,2,3)
plot(x,y,'ro',xa,Ya,'b')
axis([-1.5 1.5 -1 1]); xlabel('x'); ylabel('y')
legend('data','5-4-5-1 tanh ANN')




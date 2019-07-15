function [W, rmse] = ANN_train(alph,nnl,X,Y)
%  Multi-layer acyclic feed-forward ANN training
%	alph = training speed, nnl = network geometry
%	X = trainng set inputs,  Y = training set outputs
%	W = cell array of trained layer weights
%	rmse = matrix of root mean square error training history
% initialize variables and parameters
fn   = @tanh;			% neuron tranfer function (tanh)
dfndx= @(fx) 1 - fx.^2;		% derivative of tranfer function
maxiter	= 10000;		% max number of training iterations
maxerr	= 1e-6;			% max allowable rmse for each output
nino = size(X,1);		% number of input  layer nodes
nono = size(Y,1);		% number of output layer nodes
ntrs = size(X,2);		% number of training samples
q    = length(nnl);		% number of layer (excluding input)
TX0  = [X ; ones(1,ntrs)];	% expanded input  samples (for bias)
TY   = [Y ; ones(1,ntrs)];	% expanded output samples (for bias)
W    = cell(q,1);		% layer weights, cell array (empty)
TX   = cell(q,1);		% layer outputs, cell array (empty)
rmse = maxerr*ones(maxiter,nono); % training error array (for output)
% initialize layer weights to random values (range = -1 to 1)
W{1} = [2*rand(nnl(1),nino+1)-1 ; zeros(1,nino) 100];
for k = 2:q
  W{k} = [2*rand(nnl(k),nnl(k-1)+1)-1 ; zeros(1,nnl(k-1)) 100];
end
W{k}(end,end) = 1;			% gamma corrected for output layer
% perform training iterations
iter  = 0;				% current iteration
irmse = 1+maxerr;			% initial error (greater than allowable)
while iter < maxiter & any(irmse > maxerr)
  iter = iter + 1;			% update current iteration
  % forward computation of layer outputs with new weights
  TX{1} = fn(W{1}*TX0);			% layer 1 output
  for k = 2:q-1
    TX{k} = fn(W{k}*TX{k-1});		% layer k output
  end
  TX{q} = W{q}*TX{q-1};			% layer q output (ANN output)
  % residuals and rmse
  R     = TX{q} - TY;			% residuals at ANN output
  irmse = sqrt(mean(R(1:end-1,:).^2,2))'; % rmse for this iteration
  rmse(iter,:) = irmse;			% store rmse (for output)
  % back-propagation of residuals and weight update
  if iter < maxiter & any(irmse > maxerr)
    fpR  = R;				% common update factor df/dx*R
    R    = W{q}'*fpR;			%' resids back-prop to layer q-1
    W{q} = W{q} - alph*fpR*TX{q-1}';	%' updated weights of layer q (output)
    for k = q-1:-1:2
      fpR  = dfndx(TX{k}).*R;		% common update factor df/dx*R
      R    = W{k}'*fpR;			%' resids back-prop to layer k-1
      W{k} = W{k} - alph*fpR*TX{k-1}';	%' updated weights of layer k
    end
    W{1} = W{1}-alph*(dfndx(TX{1}).*R)*TX0';% updated weights of layer 1
  end
end


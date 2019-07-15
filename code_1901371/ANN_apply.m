function Y = ANN_apply(W,X0)
%	ANN application (tanh activation function)
% 	W  = cell array of trained layer weights (matrices)
%	X0 = ANN input (vector or matrix)
% forward computation of ANN output
X = tanh(W{1}*[X0 ; ones(1,size(X0,2))]);% layer 1 output
for k = 2:length(W)-1
  X = tanh(W{k}*X);			% layer k output(s)
end
Y = W{length(W)}(1:end-1,:)*X;		% ANN output


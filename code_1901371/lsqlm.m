function  [U, diagnostic] = lsqlm(fres, U)
% Levenberg-Marquardt nonlinear Least Squares
% 	fres = handle to model residuals function
% 	U    = initial vector of model parameters
% (c) 2015, H.J. Montas (University of Maryland)

% constant problem parameters and termination criteria
n       = length(U);		% number of model parameters
rtol    = 1e-6;			% termination criterion on sqrt(SSE)
dUtol   = 1e-20;		% termination criterion on dU if U is zero
dUUtol  = 1e-8;			% termination criterion on dU/U
relstep = 1e-7;  		% relative step for finite differences
iter    = 250;			% max iterations to perform
% initial values for iterations
R       = fres(U);		% initial residuals
rsse    = norm(R);		% sqrt(SSE) for initial residuals
dU      = U;			% initial update vector (pseudo)
lambda  = 1;			% initial Levenberg-Marquardt parameter
% iterations ------------------
while rsse > rtol & norm(dU) > max(dUtol,dUUtol*norm(U)) & iter
  iter  = iter - 1;		% remaining number of iterations
  % Jacobian calculation (first-order finite differences)
  dUfd  = diag( relstep * max( relstep, abs(U) ) );
  for  k = 1:n
    J(:,k) = (fres(U+dUfd(:,k)) - R)/dUfd(k,k);
  end
  % update of model parameter vector and lambda
  JJ	= J'*J;			% Gauss-Newton matrix
  JR	= J'*R;			% descent vector
  orsse	= rsse;			% previous sqrt(SSE)
  while rsse >= orsse & norm(dU) > max(dUtol,dUUtol*norm(U))
    dU	   = (JJ+lambda*eye(n)) \ -JR;	% new parameter update vector
    R	   = fres(U+dU);	% updated residuals
    rsse   = norm(R);		% updated sqrt(SSE) of residuals
    lambda = 5*lambda;		% move away from Gauss-Newton
  end
  lambda = lambda/7.5;		% move back towards Gauss-Newton
  U 	 = U + dU;		% update parameter vector
end
% diagnostic: iterated out = 0, converged in SSE = 1, in dU = 2.
diagnostic = (rsse <= rtol) + 2*(norm(dU) <= max(dUtol,dUUtol*norm(U)));




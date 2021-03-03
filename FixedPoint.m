% Author: Jesse Sharp; Last Update: 03/03/2021
%  
function [root,Fevals] = FixedPoint(F,X0,tol,MaxFevals)
%Solves for a root, G(X)=0, for an equation in the form X = F(X), with
%initial guess X0. 
%Input tol specifies the absolute tolerance required for convergence
%Input MaxFevals specifies the maximum number of function evaluations to 
%perform before terminating the process if convergence is not achieved.

%% Set-up
Fevals = 0;
Err = inf;

%% Perform Fixed-point iteration
while Err > tol
    X1 = F(X0); 
    Fevals = Fevals+1;
    Err = norm(F(X1)-X1);
    X0 = X1;
    
    if Fevals > MaxFevals
    error('Maximum function evaluations of %i reached, current error of %e does not meet required tolerance of %e',MaxFevals,Err,tol)
    end
end
root = X1; 

% Author: Jesse Sharp; Last Update: 03/03/2021
%  
function [root,Fevals] = Aitken(F,X0,tol,MaxFevals,m)
%Solves for a root, G(X)=0, for an equation in the form X = F(X), with
%initial guess X0. 
%Input tol specifies the absolute tolerance required for convergence
%Input MaxFevals specifies the maximum number of function evaluations to 
%perform before terminating the process if convergence is not achieved.
%Input m specifies the desired dimension of the N x m difference matrices, 
%for "parital" Aitken method, m<N, where N is the system size; 
%m = N corresponds to the standard Aitken method. 

%% Set-up
Xhat=inf; %Arbitrary initialisation to ensure first error calculation continues loop
Err=inf;
Fevals = 0;

%% Perform Aitken method
%Generate initial function evaluation for Aitken method
X = zeros(length(X0),m+2);
X(:,1) = X0; 
X(:,2) = F(X0);
%We do not update Fevals here as it is updated immediately upon entering the while loop 

while Err > tol
    
    Fevals=Fevals+1;
    %Generate required function evaluations for one Aitken iteration
    for i = 2:m+1
        X(:,i+1) = F(X(:,i));
        Fevals = Fevals+1;
    end
    %Compute differences 
    for i = 1:m+1
        dX(:,i) = X(:,i+1)-X(:,i);
    end
    %Form difference matrices    
    DX0 = dX(:,1:end-1);
    DX1 = dX(:,2:end);
    %Form second difference matrix
    D2X0 = DX1-DX0;
    
    Xhatprev = Xhat;
    Xhat = X(:,1) - DX0*lsqminnorm(D2X0,dX(:,1));
    
    %Prepare for the next Aitken iteration.
    %The following block is located here for consistency with the Steffensen
    %method, however we move the Fevals=Fevals+1; line to the top
    %of the loop. This way it is appropriately not incurred on the iteration
    %where the method reaches convergence - since the Aitken method does not require 
    %the following function evaluation in the case where the Aitken method 
    %has already achieved an error sufficiently below tolerance. Unlike the
    %other methods considered, the error calculation for Aitken's method
    %does not involve a function evaluation. 
    X(:,1) = X(:,end);
    X(:,2) = F(X(:,1));
    
    Err = norm(Xhatprev-Xhat);
    
    if Fevals > MaxFevals
        error('Maximum function evaluations of %i reached, current error of %e does not meet required tolerance of %e',MaxFevals,Err,tol)
    end
end
root = Xhat; 


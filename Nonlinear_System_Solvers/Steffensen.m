% Author: Jesse Sharp; Last Update: 15/03/2021
%  
function [root,Fevals] = Steffensen(F,X0,tol,MaxFevals,m)
%Solves for a root, G(X)=0, for an equation in the form X = F(X), with
%initial guess X0. 
%Input tol specifies the absolute tolerance required for convergence
%Input MaxFevals specifies the maximum number of function evaluations to 
%perform before terminating the process if convergence is not achieved.
%Input m specifies the desired dimension of the N x m difference matrices, 
%for "parital" Steffensen method, m<N, where N is the system size; 
%m = N corresponds to the standard Steffensen method. 

%% Set-up
Err=inf;
Fevals = 0;

%% Perform Steffensen method
%Generate initial function evaluation for Steffensen method
X = zeros(length(X0),m+2);
X(:,1) = X0; 
X(:,2) = F(X0);
Fevals=Fevals+1;

while Err > tol
    
    %Generate required function evaluations for one Steffensen iteration
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
    
    Xhat = X(:,1) - DX0*lsqminnorm(D2X0,dX(:,1)); %Compute next Steffensen values
    %lsqminnorm is used rather than explicitly forming the Moore-Penrose
    %pseudoinverse and using Xhat = X(:,1) - DX0*pinv(D2X0)*dX(:,1);
    X(:,2) = F(Xhat);
    Fevals=Fevals+1;
    
    Err = norm(X(:,2)-Xhat);
    X(:,1) = Xhat;
    
    if Fevals > MaxFevals
        error('Maximum function evaluations of %i reached, current error of %e does not meet required tolerance of %e',MaxFevals,Err,tol)
    end    
end
root = Xhat;


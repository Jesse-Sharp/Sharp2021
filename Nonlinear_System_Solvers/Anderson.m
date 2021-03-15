% Author: Jesse Sharp; Last Update: 15/03/2021
%  
function [root,Fevals] = Anderson(F,X0,tol,MaxFevals,M,Droptol)
%Solves for a root, G(X)=0, for an equation in the form X = F(X), with
%initial guess X0.
%Input tol specifies the absolute tolerance required for convergence
%Input MaxFevals specifies the maximum number of function evaluations to 
%perform before terminating the process if convergence is not achieved. 
%Input M determines the maximum number of previous iterations to
%incorporate in each iteration.
%Input Droptol specifies the maximum accptable condition number of dG. 

%% Set-up
Fevals = 0;
Err = inf;

%% Perform Anderson method
%Generate initial function evaluations for Anderson method
X1 = F(X0);
X2 = F(X1);
Fevals = Fevals+2;
%Form matrix to store iterates
X = [X0,X1]; %Note: We only computed X2 so that we can form the 
%second difference matrix dG initially
%Compute initial resiuals
g0 = X1-X0;
g1 = X2-X1; 
%Form matrices
G = [g0,g1]; %Residuals
dX = X1-X0; %Differences (excl. X2)
%Form initial residual difference
dG = g1-g0;

m = 1; %initialise counter for matrix columns in Anderson method

while  Err > tol
    
    gamma = lsqminnorm(dG,G(:,end)); %Solve the least squares problem
    
    if m < M %Don't need to replace old matrix entries
        X = [X, X(:,end)+G(:,end)-(dX+dG)*gamma];  %Append updated solution
        dX = [dX,  X(:,end)-X(:,end-1)]; %Append updated differences
        G = [G, F(X(:,end))-X(:,end)]; %Append updated residual
        Fevals = Fevals+1;
        dG = [dG, G(:,end)-G(:,end-1)]; %Append updated residual differences 
        m = m+1;
    else %Need to replace old matrix entries
        Xn = X(:,end)+G(:,end)-(dX+dG)*gamma; %Compute updated iterate
        X = [X, Xn]; %Store updated iterate
        %Erase oldest stored values and include updated values.
        G(:,1:end-1) = G(:,2:end);
        G(:,end) = F(X(:,end))-X(:,end);
        Fevals = Fevals +1;
        dX(:,1:end-1) = dX(:,2:end);
        dX(:,end) = X(:,end)-X(:,end-1);
        dG(:,1:end-1) = dG(:,2:end);
        dG(:,end) = G(:,end)-G(:,end-1);
    end
    ConddG = cond(dG); %Check condition of the residual differences matrix
    %If matrix is poorly conditioned, remove oldest columns. 
    while m > 1 && ConddG>Droptol
        dG = dG(:,2:end);
        dX = dX(:,2:end);
        m = m - 1;
        fprintf('cond(D) = %e, reducing mAA to %d \n', ConddG, m);
        ConddG = cond(dG); %Check condition of updated residual differences matrix
    end
    Err = norm(G(:,end));
    
    if Fevals > MaxFevals
        error('Maximum function evaluations of %i reached, current error of %e does not meet required tolerance of %e',MaxFevals,Err,tol)
    end
end
root = X(:,end);





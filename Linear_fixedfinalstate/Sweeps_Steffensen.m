% Author: Jesse Sharp; Last Update: 03/03/2021
%  
%Performs the FBSM with the Steffensen method for a given guess theta of Lambda(end) 
%and model parameters. 
function [z,SweepFevals] = Sweeps_Steffensen(theta,parameters)
%% Set-up
parameters.theta = theta;
tol=parameters.tol;
MaxFevals = parameters.MaxFevals;
N = parameters.N;
omega = parameters.omega;
Tfinal = parameters.Tfinal;
t_y = linspace(0,Tfinal,N); %Time discretisation for plotting

%Steffensen parameter
m = parameters.m; 

y0=1; %Linear model initial conditon
X0 = zeros(N,1); %Initial guess for the control

SweepFevals = 0;
Err = inf;

%% Perform Steffensen method
%Generate initial control update for Steffensen method
X = zeros(length(X0),m+2);
X(:,1) = X0;
[~,~,X1] = FBSM(y0,X0,parameters);
X(:,2) = omega*X0+(1-omega)*X1; %Apply relaxation factor to aid convergence
SweepFevals = SweepFevals+1;

while Err > tol
    
    %Generate required function evaluations for one Steffensen iteration
    for i = 2:m+1
        [~,~,Xn] = FBSM(y0,X(:,i),parameters);
        X(:,i+1) = omega*X(:,i)+(1-omega)*Xn; %Apply relaxation factor to aid convergence
        SweepFevals = SweepFevals+1;
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
    [y,Lambda,Xn] = FBSM(y0,Xhat,parameters);
    X(:,2) = omega*Xhat+(1-omega)*Xn; %Apply relaxation factor to aid convergence
    SweepFevals=SweepFevals+1;
    
    Err = norm(X(:,2)-Xhat);
    X(:,1) = Xhat;
    
    if SweepFevals > MaxFevals
        error('Maximum function evaluations of %i reached, current error of %e does not meet required tolerance of %e',MaxFevals,Err,tol)
    end
end
z = [t_y',y,Xn,Lambda]; %Return time discretisation, state, control and co-state.

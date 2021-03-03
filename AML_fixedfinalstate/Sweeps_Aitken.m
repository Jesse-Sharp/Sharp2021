% Author: Jesse Sharp; Last Update: 03/03/2021
%   
%performs the FBSM for a given guess theta of Lambda(end,3) and model
%parameters
function [z,SweepFevals] = Sweeps_Aitken(theta,parameters) 
%% Set-up
parameters.theta = theta;
tol = parameters.tol;
MaxFevals = parameters.MaxFevals;
N = parameters.N;
omega = parameters.omega;
Tfinal = parameters.Tfinal;
t_y = linspace(0,Tfinal,N); %Time discretisation for plotting

%Aitken paremeter
m = parameters.m;

%AML model initial conditions
S(1) = 1-parameters.gs/parameters.ps;
A(1) = 0.3255;
L(1) = 0.3715;
y0 = [S(1),A(1),L(1)]; %State initial condition
X0 = zeros(N,1); %Initial guess for the control

Xhat = inf; %Arbitrary initialisation to ensure first error calculation continues loop

SweepFevals = 0; %Initialise Fevals for this sweep
Err = inf; %Initialise error

%% Perform Aitken method
%Generate initial control update for Aitken method
X = zeros(length(X0),m+2);
X(:,1) = X0;
[~,~,X1] = FBSM(y0,X0,parameters);
%We do not update SweepFevals here as it is updated immediately upon entering the while loop 
X(:,2) = omega*X0+(1-omega)*X1; %Apply relaxation factor to aid convergence

while Err > tol
    
    SweepFevals=SweepFevals+1;
    %Generate required function evaluations for one Aitken iteration
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
    
    Xhatprev = Xhat; %Store previous Aitken values
    Xhat = X(:,1) - DX0*lsqminnorm(D2X0,dX(:,1)); %Compute next Aitken values
    %lsqminnorm is used rather than explicitly forming the Moore-Penrose
    %pseudoinverse and using Xhat = X(:,1) - DX0*pinv(D2X0)*dX(:,1);
    
    %Prepare for the next Aitken iteration.
    %The following block is located here for consistency with the Steffensen
    %method, however we move the SweepFevals=SweepFevals+1; line to the top
    %of the loop. This way it is appropriately not incurred on the iteration
    %where the method reaches convergence - since the Aitken method does not require 
    %the following function evaluation in the case where the Aitken method 
    %has already achieved an error sufficiently below tolerance. Unlike the
    %other methods considered, the error calculation for Aitken's method
    %does not involve a function evaluation. 
    X(:,1) = X(:,end);
    [y,Lambda,Xn] = FBSM(y0,X(:,1),parameters); 
    X(:,2) = omega*X(:,1)+(1-omega)*Xn; %Apply relaxation factor to aid convergence
    
    Err = norm(Xhatprev-Xhat);
    
    if SweepFevals > MaxFevals
        error('Maximum function evaluations of %i reached, current error of %e does not meet required tolerance of %e',MaxFevals,Err,tol)
    end
end
z = [t_y',y,Xhat,Lambda]; %Return time discretisation, state, control and co-state.

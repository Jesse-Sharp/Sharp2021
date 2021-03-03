% Author: Jesse Sharp; Last Update: 03/03/2021
%  
%Performs the standard FBSM for a given guess theta of Lambda(end) 
%and model parameters. 
function [z,SweepFevals] = Sweeps(theta,parameters) 
%% Set-up
parameters.theta = theta;  
tol=parameters.tol;
MaxFevals = parameters.MaxFevals;
N = parameters.N;
omega = parameters.omega;
Tfinal = parameters.Tfinal;
t_y = linspace(0,Tfinal,N); %Time discretisation for plotting

%Linear model initial conditons
y0 = 1; %State initial condition
U = zeros(N,1); %Initial guess for the control

SweepFevals = 0; %Initialise Fevals for this sweep
Err = inf; %Initialise error

%% Perform FBSM
while Err > tol
    
    uprev = U; %Store control from previous iteration
    
    [y,Lambda,Uupdate] = FBSM(y0,U,parameters);
    SweepFevals = SweepFevals + 1;
    U = omega*U + (1-omega)*Uupdate; %Apply relaxation factor to aid convergence
    Err = norm(U-uprev);
    
    if SweepFevals > MaxFevals
        error('Maximum function evaluations of %i reached, current error of %e does not meet required tolerance of %e',MaxFevals,Err,tol)
    end
end
z = [t_y',y,U,Lambda]; %Return time discretisation, state, control and co-state

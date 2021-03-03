% Author: Jesse Sharp; Last Update: 03/03/2021
%  
%Performs the FBSM with Anderson acceleration for a given guess theta of Lambda(end,3) 
%and model parameters. 
function [z,SweepFevals] = Sweeps_Anderson(theta,parameters)
%% Set-up
parameters.theta = theta;  
tol=parameters.tol;
MaxFevals = parameters.MaxFevals;
N = parameters.N;
omega = parameters.omega;
Tfinal = parameters.Tfinal;
t_y = linspace(0,Tfinal,N); %Time discretisation for plotting

%Anderson parameters
M = parameters.M;
Droptol = parameters.Droptol;

%AML model initial conditions
S(1) = 1-parameters.gs/parameters.ps;
A(1) = 0.3255;
L(1) = 0.3715;
y0 = [S(1),A(1),L(1)]; %State initial condition
X0 = zeros(N,1); %Initial guess for the control

SweepFevals = 0;
Err = inf;

%% Perform Anderson method
%Generate initial control updates for Anderson method
[~,~,X1] = FBSM(y0,X0,parameters);
X1 = omega*X0 + (1-omega)*X1; %Apply relaxation factor to aid convergence
[~,~,X2] = FBSM(y0,X1,parameters);
X2 = omega*X1 + (1-omega)*X2; %Apply relaxation factor to aid convergence
SweepFevals = SweepFevals+2;
%Form matrix to store control guesses 
X = [X0,X1]; %Note: We only compute X2 so that we can form the 
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
        X = [X, X(:,end)+G(:,end)-(dX+dG)*gamma]; %Append updated control 
        dX = [dX,  X(:,end)-X(:,end-1)]; %Append updated differences
        [y,Lambda,Xn] = FBSM(y0,X(:,end),parameters); %Solve system with updated control to obtain next control value
        Xn = omega*X(:,end) + (1-omega)*Xn; %Apply relaxation factor to aid convergence
        G = [G, Xn-X(:,end)]; %Append updated residuals
        SweepFevals = SweepFevals+1;
        dG = [dG, G(:,end)-G(:,end-1)]; %Append updated residual differences
        m = m+1;
    else %Need to replace old matrix entries
        Xn = X(:,end)+G(:,end)-(dX+dG)*gamma; %Compute updated control
        X = [X, Xn]; %Store updated control
        [y,Lambda,Xn] = FBSM(y0,X(:,end),parameters); %Solve system with updated control
        SweepFevals = SweepFevals +1;
        Xn = omega*X(:,end) + (1-omega)*Xn; %Apply relaxation factor to aid convergence
        %Erase oldest stored values and include newest computed values.
        G(:,1:end-1) = G(:,2:end);
        G(:,end) = Xn-X(:,end);
        dX(:,1:end-1) = dX(:,2:end);
        dX(:,end) = X(:,end)-X(:,end-1);
        dG(:,1:end-1) = dG(:,2:end);
        dG(:,end) = G(:,end)-G(:,end-1);
    end
    ConddG = cond(dG); %Check condition of the residual differences matrix
    %If matrix is poorly conditioned, remove oldest columns
    while m > 1 && ConddG>Droptol 
        dG = dG(:,2:end);
        dX = dX(:,2:end);
        m = m - 1;
        fprintf('cond(D) = %e, reducing mAA to %d \n', ConddG, m);
        ConddG = cond(dG);
    end
    Err = norm(G(:,end));
    if SweepFevals > MaxFevals
        error('Maximum function evaluations of %i reached, current error of %e does not meet required tolerance of %e',MaxFevals,Err,tol)
    end
end
z = [t_y',y,Xn,Lambda]; %Return time discretisation, state, control and co-state.

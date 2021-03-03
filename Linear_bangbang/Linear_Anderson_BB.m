% Author: Jesse Sharp; Last Update: 03/03/2021
%  
function [Control,Fevals] = Linear_Anderson_BB(tol,MaxFevals,M,Droptol)
%Input tol specifies the absolute tolerance required for convergence
%Input MaxFevals specifies the maximum number of function evaluations to 
%perform before terminating the process if convergence is not achieved. 
%Input M determines the maximum number of previous iterations to
%incorporate in each iteration.
%Input Droptol specifies the maximum accptable condition number of dG.

%% Set-up
Tfinal = 1;  %Specified final time
parameters.dt = 2^(-8); %time-step
parameters.N = floor(Tfinal/parameters.dt+1); %Number of nodes in time discretisation
t_y = linspace(0,Tfinal,parameters.N); %Time discretisation for plotting
omega = 0; %Portion of previous iteration's control maintained when updating control

%Model parameters
parameters.Uupper = 2; %Control upper bound
parameters.Ulower = 0; %Control lower bound
parameters.Gamma = 0.5;
parameters.a = 1; %Weighting on state
parameters.b = 3; %Weighting on control

%Linear model initial conditions
y0=1; %State initial condition
X0 = zeros(parameters.N,1); %Initial guess for the control

Fevals = 0; %Initialise function evaluation count
Err = inf; %Initialise the error

%% Perform Anderson method
%Generate initial control updates for Anderson method
[~,~,X1] = FBSM_BB(y0,X0,parameters);
X1 = omega*X0 + (1-omega)*X1;
[~,~,X2] = FBSM_BB(y0,X1,parameters);
X2 = omega*X1 + (1-omega)*X2;
Fevals = Fevals+2;
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
    
    if m < M %don't need to replace old matrix entries
        X = [X, X(:,end)+G(:,end)-(dX+dG)*gamma]; %Append updated control 
        dX = [dX,  X(:,end)-X(:,end-1)]; %Append updated differences
        [~,~,Xn] = FBSM_BB(y0,X(:,end),parameters); %Solve system with updated control
        Xn = omega*X(:,end) + (1-omega)*Xn; %Apply relaxation factor to aid convergence
        G = [G, Xn-X(:,end)]; %Append updated residuals
        Fevals = Fevals+1;
        dG = [dG, G(:,end)-G(:,end-1)]; %Append updated residual differences
        m = m+1;
    else %Need to replace old matrix entries
        Xn = X(:,end)+G(:,end)-(dX+dG)*gamma; %Compute updated control
        X = [X Xn]; %Store updated control
        [~,~,Xn] = FBSM_BB(y0,X(:,end),parameters); %Solve system with updated control
        Fevals = Fevals +1;
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
        fprintf(' cond(D) = %e, reducing mAA to %d \n', ConddG, m);
        ConddG = cond(dG);
    end
    Err = norm(G(:,end));
    if Fevals > MaxFevals
        error('Maximum function evaluations of %i reached, current error of %e does not meet required tolerance of %e',MaxFevals,Err,tol)
    end    
    
    %%To view the solution as it converges, uncomment this block (line 96 to line 110)
    %%to produce interim figures
%     figure(1)
%     box on
%     [y,~,~] = FBSM_BB(y0,Xn,parameters); %Compute state ODE based on most recent control
%     line1 = plot(t_y,y,'b','LineWidth',2);
%     hold on
%     line2 = plot(t_y,Xn,'k--','LineWidth',2);
%     hL = legend([line1,line2],{'\it{x}','\it{u*}'},'Location','northeast');
%     ylabel('State','fontsize',18);
%     xlabel('\it{t}','fontsize',18);
%     axis([0,Tfinal,0,5])
%     text(0.2,4,sprintf('Fevals: %0.0f',Fevals));
%     xt = get(gca, 'XTick');
%     set(gca, 'FontSize', 18)
%     hold off
%     set(gca, 'FontName', 'Times New Roman')

end
Control = X(:,end);

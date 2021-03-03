% Author: Jesse Sharp; Last Update: 03/03/2021
%  
function [Control,Fevals] = AML_Anderson_BB(tol,MaxFevals,M,Droptol,omega)
%Input tol specifies the absolute tolerance required for convergence
%Input MaxFevals specifies the maximum number of function evaluations to 
%perform before terminating the process if convergence is not achieved. 
%Input M determines the maximum number of previous iterations to
%incorporate in each iteration.
%Input Droptol specifies the maximum accptable condition number of dG. 
%Input omega specifies the control update weighting to aid convergence. 

%% Set-up
Tfinal = 10;  %Specified final time
parameters.dt = 2^-11; %Time-step
parameters.N = Tfinal/parameters.dt+1; %Number of nodes in time discretisation
t_y = linspace(0,Tfinal,parameters.N); %Time discretisation for plotting

%Model parameters
parameters.ps = 0.5; %Proliferation of S
parameters.pa = 0.43; %Proliferation of A
parameters.pl = 0.27; %Proliferation of L
parameters.gs = 0.14; %Differentiation of S to A
parameters.ga = 0.44; %Differentiation of A to D
parameters.gl = 0.05; %Differentiation of L to T
parameters.ux = 0.275; %Migration of D into the blood stream
parameters.ut = 0.3; %Migration of T into the blood stream
parameters.K1 = 1; %Carrying capacity of the compartment with S
parameters.K2 = 1; %Carrying capacity of the compartment with SA + L
parameters.Alpha = 0.015; %Michaelis-Menten kinetic parameter alpha
parameters.Gamma = 0.1; %Michaelis-Menten kinetic parameter gamma
parameters.a1 = 1; %Weighting on negative impact of control
parameters.a2 = 2; %Weighting on negative impact of leukaemia
parameters.Uupper = 0.3; %Control upper bound
parameters.Ulower = 0; %Control lower bound

%AML model initial conditions
S(1) = 1-parameters.gs/parameters.ps; 
A(1) = 0.3255; 
L(1) = 0.3715; 
y0 = [S(1),A(1),L(1)]; %State initial condition
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
    
    %%To view the solution as it converges, uncomment this block (line 107 to line 122)
    %%to produce interim figures
%     figure(1)
%     box on
%     [y,~,~] = FBSM_BB(y0,X(:,end),parameters); %Compute state ODE based on most recent control
%     line1 = plot(t_y,y(:,1),'LineWidth',2);
%     hold on
%     line2 = plot(t_y,y(:,2),'LineWidth',2);
%     line3 = plot(t_y,y(:,3),'LineWidth',2);
%     line4 = plot(t_y,X(:,end),'k--','LineWidth',2);
%     hL = legend([line1,line2,line3,line4],{'\it{S}','\it{A}','\it{L}','\it{u*}'},'Location','northeast');
%     ylabel('State','fontsize',18);
%     xlabel('\it{t}','fontsize',18);
%     axis([0,Tfinal,0,1])
%     set(gca, 'FontSize', 18)
%     hold off
%     set(gca, 'FontName', 'Times New Roman')
%     drawnow
    
    if Fevals > MaxFevals
        error('Maximum function evaluations of %i reached, current error of %e does not meet required tolerance of %e',MaxFevals,Err,tol)
    end    

end
Control = X(:,end);

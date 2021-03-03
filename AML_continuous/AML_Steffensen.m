% Author: Jesse Sharp; Last Update: 03/03/2021
%   
function [Control,Fevals] = AML_Steffensen(tol,MaxFevals,m,omega)
%Input tol specifies the absolute tolerance required for convergence
%Input MaxFevals specifies the maximum number of function evaluations to 
%perform before terminating the process if convergence is not achieved.
%Input m specifies the desired dimension of the N x m difference matrices, 
%for "parital" Steffensen method, m<N, where N is the system size; 
%m = N corresponds to the standard Steffensen method. 
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

%AML model initial conditions
S(1) = 1-parameters.gs/parameters.ps; 
A(1) = 0.3255; 
L(1) = 0.3715; 
y0 = [S(1),A(1),L(1)]; %State initial condition
X0 = zeros(parameters.N,1); %Initial guess for the control

Fevals = 0; %Initialise function evaluation count
Err = inf; %Initialise the error

%% Perform Steffensen method
%Generate initial control update for Steffensen method
X = zeros(length(X0),m+2);
X(:,1) = X0;
[~,~,X1] = FBSM(y0,X0,parameters);
X(:,2) = omega*X0+(1-omega)*X1; %Apply relaxation factor to aid convergence
Fevals = Fevals+1;

while Err > tol
    
    %Generate required function evaluations for one Steffensen iteration
    for i = 2:m+1
        [~,~,Xn] = FBSM(y0,X(:,i),parameters);
        X(:,i+1) = omega*X(:,i)+(1-omega)*Xn; %Apply relaxation factor to aid convergence
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
    [~,~,Xn] = FBSM(y0,Xhat,parameters);
    X(:,2) = omega*Xhat+(1-omega)*Xn; %Apply relaxation factor to aid convergence
    Fevals=Fevals+1;
    
    Err = norm(X(:,2)-Xhat);
    X(:,1) = Xhat;
    
    %%To view the solution as it converges, uncomment this block (line 82 to line 97)
    %%to produce interim figures
%     figure(1)
%     box on
%     [y,~,~] = FBSM(y0,X(:,2),parameters); %Compute state ODE based on most recent control
%     line1 = plot(t_y,y(:,1),'LineWidth',2);
%     hold on
%     line2 = plot(t_y,y(:,2),'LineWidth',2);
%     line3 = plot(t_y,y(:,3),'LineWidth',2);
%     line4 = plot(t_y,X(:,2),'k--','LineWidth',2);
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
Control = X(:,2);
% Author: Jesse Sharp; Last Update: 03/03/2021
%  
function [Control,Fevals] = AML_Aitken(tol,MaxFevals,m,omega)
%Input tol specifies the absolute tolerance required for convergence
%Input MaxFevals specifies the maximum number of function evaluations to 
%perform before terminating the process if convergence is not achieved.
%Input m specifies the desired dimension of the N x m difference matrices, 
%for "parital" Aitken method, m<N, where N is the system size; 
%m = N corresponds to the standard Aitken method. 
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

Xhat = inf; %Arbitrary initialisation to ensure first error calculation continues loop
Fevals = 0; %Initialise function evaluation count
Err = inf; %Initialise the error

%% Perform Aitken method
%Generate initial control update for Aitken method
X = zeros(length(X0),m+2);
X(:,1) = X0;
[~,~,X1] = FBSM(y0,X0,parameters);
%We do not update SweepFevals here as it is updated immediately upon entering the while loop 
X(:,2) = omega*X0+(1-omega)*X1; %Apply relaxation factor to aid convergence

while Err > tol
    
    Fevals=Fevals+1;
    %Generate required function evaluations for one Aitken iteration
    for i = 2:m+1
        [~,~,Xn] = FBSM(y0,X(:,i),parameters);
        X(:,i+1) = omega*X(:,i)+(1-omega)*Xn;
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
    
    Xhatprev = Xhat; %Store previous Aitken values
    Xhat = X(:,1) - DX0*lsqminnorm(D2X0,dX(:,1)); %Compute next Aitken values
    %lsqminnorm is used rather than explicitly forming the Moore-Penrose
    %pseudoinverse and using Xhat = X(:,1) - DX0*pinv(D2X0)*dX(:,1);
    
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
    [~,~,Xn] = FBSM(y0,X(:,1),parameters);
    X(:,2) = omega*X(:,1)+(1-omega)*Xn; %Apply relaxation factor to aid convergence
    
    Err = norm(Xhatprev-Xhat);
    
    %%To view the solution as it converges, uncomment this block (line 94 to line 109)
    %%to produce interim figures
%     figure(1)
%     box on
%     [y,~,~] = FBSM(y0,Xhat,parameters); %Compute state ODE based on most recent control
%     line1 = plot(t_y,y(:,1),'LineWidth',2);
%     hold on
%     line2 = plot(t_y,y(:,2),'LineWidth',2);
%     line3 = plot(t_y,y(:,3),'LineWidth',2);
%     line4 = plot(t_y,Xhat,'k--','LineWidth',2);
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
Control = Xhat; %Return most recently computed Aitken value

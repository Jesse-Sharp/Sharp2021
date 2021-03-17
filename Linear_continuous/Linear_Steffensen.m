% Author: Jesse Sharp; Last Update: 03/03/2021
%    
function [Control,Fevals] = Linear_Steffensen(tol,MaxFevals,m)
%Input tol specifies the absolute tolerance required for convergence
%Input MaxFevals specifies the maximum number of function evaluations to 
%perform before terminating the process if convergence is not achieved.
%Input m specifies the desired dimension of the N x m difference matrices, 
%for "parital" Steffensen method, m<N, where N is the system size; 
%m = N corresponds to the standard Steffensen method. 

%% Set-up
Tfinal = 1;  %Specified final time
parameters.dt = 2^(-8); %time-step
parameters.N = floor(Tfinal/parameters.dt+1); %Number of nodes in time discretisation
t_y = linspace(0,Tfinal,parameters.N); %Time discretisation for plotting
omega = 0; %Portion of previous iteration's control maintained when updating control

%Model parameters
parameters.Gamma = 0.5;
parameters.a = 1; %Weighting on state
parameters.b = 1; %Weighting on control

%Linear model initial conditions
y0=1; %State initial condition
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
    
    if Fevals > MaxFevals
        error('Maximum function evaluations of %i reached, current error of %e does not meet required tolerance of %e',MaxFevals,Err,tol)
    end
    
    %%To view the solution as it converges, uncomment this block (line 72 to line 86)
    %%to produce interim figures
%     figure(1)
%     box on
%     [y,~,~] = FBSM(y0,X(:,2),parameters); %Compute state ODE based on most recent control
%     line1 = plot(t_y,y,'b','LineWidth',2);
%     hold on
%     line2 = plot(t_y,X(:,2),'k--','LineWidth',2);
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
Control = X(:,2);

%% Plot state with optimal control
figure
box on
[y,~,~] = FBSM(y0,Control,parameters); %Compute state ODE based on most recent control
line1 = plot(t_y,y,'LineWidth',2);
hold on
line2 = plot(t_y,Control,'k--','LineWidth',2);
hL = legend([line1,line2],{'\it{x}','\it{u*}'},'Location','northeast');
axis([0,Tfinal,0,5])
xt = get(gca, 'XTick');
set(gca, 'FontSize', 18)
hold off
set(gca, 'FontName', 'Times New Roman')
ylabel('State','fontsize',18);
xlabel('\it t','fontsize',18);

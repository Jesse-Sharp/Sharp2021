% Author: Jesse Sharp; Last Update: 03/03/2021
%  
%% Set-up
tol = 1e-10; %Absolute tolerance required for convergence
MaxFevals = 100; %Number of iterations to perform before giving up if convergence is not reached
Tfinal = 10;  %Specified final time
parameters.dt = 2^-11; %Time-step
parameters.N = Tfinal/parameters.dt+1; %Number of nodes in time discretisation
t_y = linspace(0,Tfinal,parameters.N); %Time discretisation for plots
omega = 0.55; %Portion of previous iteration's control maintained when updating control

%Model parameters
parameters.ps = 0.5; %Proliferation of S
parameters.pa = 0.43; %Proliferation of A
parameters.pl = 0.27; %Proliferation of L
parameters.ga = 0.44; %Differentiation of A to D
parameters.gl = 0.05; %Differentiation of L to T
parameters.gs = 0.14; %Differentiation of S to A
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
U = zeros(parameters.N,1); %Initial guess for the control

Fevals = 0; %Initialise function evaluation count
Err = inf; %Initialise the error

%% Forward-backward sweep
while Fevals<MaxFevals && Err > tol
    
    uprev=U; %Store control from previous iteration
    
    [y,Lambda,Uupdate] = FBSM(y0,U,parameters); %Forward Backward sweep (1 iteration)
    Fevals = Fevals+1;
    
    U = omega*U  + (1-omega)*Uupdate; %Update the value of the control and apply a relaxation factor.
    
    Err = norm(U-uprev);
    
    %%To view the solution as it converges, uncomment this block (line 53 to line 67)
    %%to produce interim figures
%     figure(1)
%     box on
%     line1 = plot(t_y,y(:,1),'LineWidth',2);
%     hold on
%     line2 = plot(t_y,y(:,2),'LineWidth',2);
%     line3 = plot(t_y,y(:,3),'LineWidth',2);
%     line4 = plot(t_y,U,'k--','LineWidth',2);
%     hL = legend([line1,line2,line3,line4],{'S','A','L','u*'},'Location','northeast');
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

%% Plot states with optimal control
figure
[y,~,~] = FBSM(y,U,parameters);  %Compute state ODE based on most recent control
box on
line1 = plot(t_y,y(:,1),'LineWidth',2);
hold on
line2 = plot(t_y,y(:,2),'LineWidth',2);
line3 = plot(t_y,y(:,3),'LineWidth',2);
line4 = plot(t_y,U,'k--','LineWidth',2);
hL = legend([line1,line2,line3,line4],{'\it{S}','\it{A}','\it{L}','\it{u*}'},'Location','northeast');
set(gca, 'FontName', 'Times New Roman')
ylabel('State','fontsize',18);
xlabel('\it t','fontsize',18);
axis([0,Tfinal,0,1])
set(gca, 'FontSize', 18)
hold off

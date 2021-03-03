% Author: Jesse Sharp; Last Update: 03/03/2021
%   
%% Set-up
tol = 1e-10; %Absolute tolerance required for convergence
MaxFevals = 100; %Number of function evaluations to perform before terminating the process if convergence is not reached
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
U = zeros(parameters.N,1); %Initial guess for the control

Fevals = 0; %Initialise function evaluation count
Err = inf; %Initialise the error

%% Forward-backward sweep
while Err > tol
    
    uprev = U;%Store control from previous iteration
    
    i = 0; %Initialise loop variable
    t = 0; %Initialise time for forward sweep
    
    %Forward sweep using fourth-order Runge-Kutta scheme
    [~,~,Uupdate] = FBSM(y0,U,parameters);
    
    Fevals = Fevals + 1;
    
    U = omega*U + (1-omega)*Uupdate; %Actual updated control after applying relaxation to aid convergence
    
    Err = norm(U-uprev);
    
    %%To view the solution as it converges, uncomment this block (line 43 to line 58)
    %%to produce interim figures
%     figure(1)
%     box on
%     [y,~,~] = FBSM(y0,U,parameters); %Compute state ODE based on most recent control
%     line1 = plot(t_y,y,'b','LineWidth',2);
%     hold on
%     line2 = plot(t_y,U,'k--','LineWidth',2);
%     hL = legend([line1,line2],{'\it{x}','\it{u*}'},'Location','northeast');
%     ylabel('State','fontsize',18);
%     xlabel('\it{t}','fontsize',18);
%     axis([0,Tfinal,0,5])
%     text(0.2,4,sprintf('Fevals: %0.0f',Fevals));
%     set(gca, 'FontSize', 18)
%     hold off
%     set(gca, 'FontName', 'Times New Roman')
%     drawnow
    
    if Fevals > MaxFevals
        error('Maximum function evaluations of %i reached, current error of %e does not meet required tolerance of %e',MaxFevals,Err,tol)
    end
    
end

%% Plot state with optimal control
figure
box on
[y,~,~] = FBSM(y0,U,parameters); %Compute state ODE based on most recent control
line1 = plot(t_y,y,'LineWidth',2);
hold on
line2 = plot(t_y,U,'k--','LineWidth',2);
hL = legend([line1,line2],{'\it{x}','\it{u*}'},'Location','northeast');
axis([0,Tfinal,0,5])
xt = get(gca, 'XTick');
set(gca, 'FontSize', 18)
hold off
set(gca, 'FontName', 'Times New Roman')
ylabel('State','fontsize',18);
xlabel('\it t','fontsize',18);
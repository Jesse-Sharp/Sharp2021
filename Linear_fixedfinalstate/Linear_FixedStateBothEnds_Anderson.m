% Author: Jesse Sharp; Last Update: 03/03/2021
%  
function [Sol,Fevals] = Linear_FixedStateBothEnds_Anderson(tol,MaxFevals,M,Droptol)
%Input tol specifies the absolute tolerance required for convergence
%Input MaxFevals specifies the maximum number of function evaluations to 
%perform for an inner problem before terminating the process if convergence 
%is not achieved. 
%Input M determines the maximum number of previous iterations to
%incorporate in each iteration.
%Input Droptol specifies the maximum acceptable condition number of dG. 

%% Set-up
parameters.Tfinal = 1; %Specified final time
Tfinal = parameters.Tfinal;  
parameters.dt = 2^(-8); %Time-step
parameters.N = floor(Tfinal/parameters.dt+1); %Number of nodes in time discretisation
parameters.omega = 0; %Portion of previous iteration's control maintained when updating control
parameters.tol = tol;
parameters.MaxFevals = MaxFevals;

%Anderson parameters
parameters.M = M;
parameters.Droptol = Droptol; 

%Model parameters
parameters.Gamma = 0.5; 
parameters.a = 1; %Weighting on state
parameters.b = 1; %Weighting on control

%Adapted FBSM inputs
SecantTol = 1e-10; %Convergence tolerance for secant method
s = 10; %Known fixed end point for state (i.e. y(end) = s)
parameters.theta1 = -10;
theta1 = parameters.theta1; %First guess for Lambda(end)
parameters.theta2  = 10;
theta2 = parameters.theta2; %Second guess for Lambda(end)

%Other set-up
CumulativeFevals = 0;

%Generate first secant method value using first initial guess for
%Lambda(end).
[z1,Fevals]  = Sweeps_Anderson(theta1,parameters);
CumulativeFevals = CumulativeFevals + Fevals;
Va = z1(end,2)-s; %Record difference between final value of y(end) and target, s.
%Note that z1,z2 and z have form z=[t_y',y,U,Lambda] so z(end,2) is y(end). 

%Plotting
figure(1)
line1=plot(z1(:,1),z1(:,2),'k');
axis([0,Tfinal,-20,30])
h = text(0.2,-1,sprintf('Cumulative Sweeps: %0.0f',CumulativeFevals));
hL = legend(line1,{'y1'},'Location','northwest');

%Generate second secant method value using second initial guess for
%Lambda(end).
[z2,Fevals] = Sweeps_Anderson(theta2,parameters);
CumulativeFevals = CumulativeFevals + Fevals;
Vb = z2(end,2)-s; %Record difference between final value of y(end) and target, s.

%Plotting
figure(1)
delete(h)
hold on
line1=plot(z2(:,1),z2(:,2),'k');
h = text(0.2,-1,sprintf('Cumulative Sweeps: %0.0f',CumulativeFevals));
hL = legend(line1,{'y1'},'Location','northwest');

%%Perform Secant method
SecantErr = inf; %Initialise error for secant method
while SecantErr > SecantTol
    
    %Choose best value to retain
    if (abs(Va) > abs(Vb))
        temp1 = theta1;
        theta1 = theta2;
        theta2 = temp1;
        temp2 = Va;
        Va = Vb;
        Vb = temp2;
    end
    
    d = Va*(theta2-theta1)/(Vb-Va); %Secant calculation
    
    %Updating stored values
    theta2 = theta1;
    Vb = Va;
    theta1 = theta1-d; %Update guess for Lambda(end)
    [z,Fevals] = Sweeps_Anderson(theta1,parameters); %Solve system with updated theta1
    Va = z(end,2) - s; %Compute how close y(end) is to target s
    SecantErr = abs(Va); %Error for secant method
    CumulativeFevals = CumulativeFevals + Fevals;
    
    %Plotting
    figure(1)
    delete(h)
    hold on
    line1=plot(z(:,1),z(:,2),'k');
    h = text(0.2,-1,sprintf('Cumulative Sweeps: %0.0f',CumulativeFevals));
    hL = legend(line1,{'y1'},'Location','northwest');
end
%Return function outputs
Sol=z;
Fevals = CumulativeFevals;


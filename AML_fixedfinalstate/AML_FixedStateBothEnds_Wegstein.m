% Author: Jesse Sharp; Last Update: 03/03/2021
%  
function [Sol,Fevals] = AML_FixedStateBothEnds_Wegstein(tol,MaxFevals,nth,bounding,lower,upper)
%Input tol specifies the absolute tolerance required for convergence
%Input MaxFevals specifies the maximum number of function evaluations to 
%perform for an inner problem before terminating the process if convergence 
%is not achieved. 
%Input nth specifies how frequently to update q; every nth iteration. 
%Input bounding acts as a switch, such that bounds are applied on q when
%bounding = 1.
%Inputs lower and upper specify the bounds to apply to q, if bounding = 1.

%% Set-up
parameters.Tfinal = 10; %Specified final time
Tfinal = parameters.Tfinal; 
parameters.dt = 2^-11; %Time-step
parameters.N = Tfinal/parameters.dt+1; %Number of nodes in time discretisation
parameters.omega = 0.55; %Portion of previous iteration's control maintained when updating control
parameters.tol = tol;
parameters.MaxFevals = MaxFevals;

%Wegstein parameters
parameters.nth = nth;
parameters.bounding = bounding;
parameters.lower = lower;
parameters.upper = upper;

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

%Adapted FBSM inputs
SecantTol = 1e-10; %Convergence tolerance for secant method
s = 0.05; %Known fixed end point for L (i.e. y(end,3) = s)
parameters.theta1 = 0;
theta1 = parameters.theta1; %First guess for Lambda(end,3)
parameters.theta2  = 10;
theta2 = parameters.theta2; %Second guess for Lambda(end,3)

%Other set-up
CumulativeFevals = 0;

%Generate first secant method value using first initial guess for
%Lambda(3,end).
[z1,Fevals]  = Sweeps_Wegstein(theta1,parameters); %Perform FBSM to convergence
CumulativeFevals = CumulativeFevals + Fevals; %Record number of times the bvp is solved.
Va = z1(end,4)-s; %Record difference between final value of y(end,3) and target, s.
%Note that z1,z2 and z have form z=[t_y',y,U,Lambda] so z(end,4) is y(end,3). 

%Plotting
figure(1)
line1=plot(z1(:,1),z1(:,4),'color',[237.15/255,175.95/255,33.15/255],'LineWidth',2);
axis([0,Tfinal,0,0.4])
h = text(4,0.3,sprintf('Cumulative Fevals: %0.0f',CumulativeFevals),'FontSize',18);
h2 = text(9.3,z1(end,4)+0.02,sprintf('%0.0f',CumulativeFevals),'FontSize',18);
hL = legend(line1,{'\it{L}'},'Location','northwest');
drawnow 

%Generate second secant method value using second initial guess for
%Lambda(3,end).
[z2,Fevals] = Sweeps_Wegstein(theta2,parameters);
CumulativeFevals = CumulativeFevals + Fevals;
Vb = z2(end,4)-s;

%Plotting
figure(1)
delete(h)
hold on
line1=plot(z2(:,1),z2(:,4),'color',[237.15/255,175.95/255,33.15/255],'LineWidth',2);
h = text(4,0.3,sprintf('Cumulative Fevals: %0.0f',CumulativeFevals),'FontSize',18);
h2 = text(9.3,z2(end,4)-0.01,sprintf('%0.0f',CumulativeFevals),'FontSize',18);
hL = legend(line1,{'\it{L}'},'Location','northwest');
drawnow 

%Perform Secant method
SecantErr = inf; %Initialise error for secant method
while SecantErr > SecantTol

if (abs(Va) > abs(Vb)) %Choose best value to retain
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
theta1 = theta1-d; %Update theta1 (guess for Lambda(3,end))
[z,Fevals] = Sweeps_Wegstein(theta1,parameters); %Solve system with updated theta1
Va = z(end,4) - s; %Compute how close y(end,3) is to target s
SecantErr = abs(Va); %Error for secant method
CumulativeFevals = CumulativeFevals + Fevals;

    %Plotting
    figure(1)
    delete(h)
    hold on
    line1=plot(z(:,1),z(:,4),'color',[237.15/255,175.95/255,33.15/255],'LineWidth',2);
    h = text(4,0.3,sprintf('Cumulative Fevals: %0.0f',CumulativeFevals),'FontSize',18);
    hL = legend(line1,{'\it{L}'},'Location','northwest');
    drawnow
end  
Sol=z;
Fevals = CumulativeFevals;


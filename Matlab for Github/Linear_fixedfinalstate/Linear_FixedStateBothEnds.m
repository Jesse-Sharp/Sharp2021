% Author: Jesse Sharp; Last Update: 03/03/2021
%  
%% Set-up  
parameters.Tfinal = 1;  %Specified final time
Tfinal = parameters.Tfinal; 
parameters.dt = 2^-8; %Time-step
parameters.N = floor(Tfinal/parameters.dt+1); %Number of nodes in time discretisation
parameters.omega = 0; %Portion of previous iteration's control maintained when updating control
parameters.tol = 1e-10; %Use for sweeps 
parameters.MaxFevals = 100; %Number of iterations to perform before giving up if convergence is not reached

%Model parameters
parameters.Gamma = 0.5;
parameters.a = 1; %Weighting on state
parameters.b = 1; %Weighting on control

%Adapted FBSM inputs
SecantTol = 1e-10; %Use for secant method
s = 10; %Known fixed end point for state (i.e. y(end) = s)
parameters.theta1 = -10;
theta1 = parameters.theta1; %First guess for Lambda(end)
parameters.theta2  = 10;
theta2 = parameters.theta2; %Second guess for Lambda(end)

%Other set-up
CumulativeFevals = 0;

%Generate first secant method value using first initial guess for
%Lambda(end).
[z1,Fevals]  = Sweeps(theta1,parameters);
CumulativeFevals = CumulativeFevals + Fevals;
Va = z1(end,2)-s; %Record difference between final value of y(end) and target, s.
%Note that z1,z2 and z have form z=[t_y',y,U,Lambda] so z(end,2) is y(end). 

%Plotting
figure(1)
line1=plot(z1(:,1),z1(:,2),'color',[0/255,114.75/255,188.7/255],'LineWidth',2);
axis([0,Tfinal,-20,30])
h = text(0.2,24,sprintf('Cumulative Fevals: %0.0f',CumulativeFevals),'FontSize',18);
h2 = text(0.91,z1(end,2)+4,sprintf('%0.0f',CumulativeFevals),'FontSize',18);
hL = legend(line1,{'\it{x}'},'Location','northwest');

%Generate second secant method value using second initial guess for
%Lambda(end).
[z2,Fevals] = Sweeps(theta2,parameters);
CumulativeFevals = CumulativeFevals + Fevals;
Vb = z2(end,2)-s; %Record difference between final value of y(end) and target, s.

%Plotting
figure(1)
delete(h)
hold on
line1=plot(z2(:,1),z2(:,2),'color',[0/255,114.75/255,188.7/255],'LineWidth',2);
h = text(0.2,24,sprintf('Cumulative Fevals: %0.0f',CumulativeFevals),'FontSize',18);
h2 = text(0.9,z2(end,2)+2,sprintf('%0.0f',CumulativeFevals),'FontSize',18);
hL = legend(line1,{'\it{x}'},'Location','northwest');
SecantErr = inf;
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
    [z,Fevals] = Sweeps(theta1,parameters); %Solve system with updated theta1
    Va = z(end,2) - s; %Compute how close y(end) is to target s
    SecantErr = abs(Va); %Error for secant method
    CumulativeFevals = CumulativeFevals + Fevals;
    
    %Plotting
    figure(1)
    delete(h)
    hold on
    line1=plot(z(:,1),z(:,2),'color',[0/255,114.75/255,188.7/255],'LineWidth',2);
    h = text(0.2,24,sprintf('Cumulative Fevals: %0.0f',CumulativeFevals),'FontSize',18);
    h2 = text(0.9,z(end,2)-3,sprintf('%0.0f',CumulativeFevals),'FontSize',18);
    hL = legend(line1,{'\it{x}'},'Location','northwest');
    
end

box on
set(gca, 'FontSize', 18)
set(gca, 'FontName', 'Times New Roman')
ylabel('State','fontsize',18);
xlabel('\it t','fontsize',18);

figure
box on
line1 = plot(z(:,1),z(:,2),'LineWidth',2);
hold on
line2 = plot(z(:,1),z(:,3),'k--','LineWidth',2);
h = text(0.2,12,sprintf('Cumulative Fevals: %0.0f',CumulativeFevals),'FontSize',18);
hL = legend([line1,line2],{'\it{x}','\it{u*}'},'Location','northeast');
axis([0,Tfinal,0,14])
set(gca, 'FontSize', 18)
set(gca, 'FontName', 'Times New Roman')
ylabel('State','fontsize',18);
xlabel('\it t','fontsize',18);
hold off

% Author: Jesse Sharp; Last Update: 03/03/2021
%  
function [Control,Fevals] = AML_Wegstein_BB(tol,MaxFevals,nth,bounding,lower,upper,omega)
%Input tol specifies the absolute tolerance required for convergence
%Input MaxFevals specifies the maximum number of function evaluations to 
%perform before terminating the process if convergence is not achieved.
%Input nth specifies how frequently to update q; every nth iteration. 
%Input bounding acts as a switch, such that bounds are applied on q when
%bounding = 1.
%Inputs lower and upper specify the bounds to apply to q, if bounding = 1.
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
iterations = 0; %Initialise iteration count

%% Perform Wegstein method
[~,~,X1] = FBSM_BB(y0,X0,parameters);
X1 = omega*X0+(1-omega)*X1;
[~,~,X2] = FBSM_BB(y0,X1,parameters);
X2 = omega*X1+(1-omega)*X2;
Fevals = Fevals + 2;

while Err > tol
    
    iterations = iterations+1;
    
    if iterations == 1 %Compute Wegstein's q on first iteration
        A = (X2-X1)./(X1-X0); %Compute A elementwise from the differences
        q = A./(A-1); %Compute q elementwise from A
        %Set NaN and Inf values to 0; choice of 0 is arbitrary, and could
        %be replaced with (for example) omega. 
        q(isnan(q))=0;
        q(isinf(q))=0;
        %Apply bounds to Wegstein q
        if bounding == 1
            if isreal(q)
                q = min(max(q,lower),upper); %Apply bounds on q
            else
                normq = norm(q);%norm of complex elements
                r = min(normq,abs(upper-lower)); %Set length between 0 and magnitude of difference of bounds. 
                %This magnitude is chosen arbitrarly, but appears to perform well
                phi = angle(q); %Compute angle of each complex number
                q = r.*exp(1i*phi); %Construct new (bounded) complex number maintaining original angle but with bounded length.
            end
        end
        %Update elementwise using Wegstein q as a weighting factor
        %to accelerate convergence
        Xn = q.*X0+(1-q).*X1;
        [~,~,Fn] = FBSM_BB(y0,Xn,parameters);
        Fevals = Fevals+1;
        Err = norm(Fn-Xn);
    else %If not first iteration, only update q if we are on the nth iteration
        if mod(iterations,nth) == 0 %Compute Wegstein's q
            A = (Fn - X2)./(Xn-X1); %Compute A elementwise from the differences
            q = A./(A-1); %Compute q elementwise from A
            q(isnan(q))=0;
            q(isinf(q))=0;
            %Apply bounds to Wegstein q
            if bounding == 1
                if isreal(q)
                    q = min(max(q,lower),upper); %Apply bounds on q
                else
                normq = norm(q); %Norm of complex elements
                r = min(normq,abs(upper-lower)); %Set length between 0 and magnitude of difference of bounds. 
                %This magnitude is chosen arbitrarly, but appears to perform well
                phi = angle(q); %Compute angle of each complex number
                q = r.*exp(1i*phi); %Construct new (bounded) complex number maintaining original angle but with bounded length. 
                end
            end
        end
        %Update stored values
        X1 = Xn;
        X2 = Fn;
        %Update elementwise using Wegstein q as a weighting factor
        %to accelerate convergence
        Xn = q.*Xn + (1-q).*Fn;
        [~,~,Fn] = FBSM_BB(y0,Xn,parameters); %Solve system with updated control to obtain next control value
        Fevals = Fevals+1;
        Err = norm(Fn-Xn);
    end
    if Fevals > MaxFevals
        error('Maximum function evaluations of %i reached, current error of %e does not meet required tolerance of %e',MaxFevals,Err,tol)
    end
    
    %%To view the solution as it converges, uncomment this block (line 119 to line 134)
    %%to produce interim figures
%     figure(1)
%     box on
%     [y,~,~] = FBSM_BB(y0,Fn,parameters); %Compute state ODE based on most recent control
%     line1 = plot(t_y,y(:,1),'LineWidth',2);
%     hold on
%     line2 = plot(t_y,y(:,2),'LineWidth',2);
%     line3 = plot(t_y,y(:,3),'LineWidth',2);
%     line4 = plot(t_y,Fn,'k--','LineWidth',2);
%     hL = legend([line1,line2,line3,line4],{'\it{S}','\it{A}','\it{L}','\it{u*}'},'Location','northeast');
%     ylabel('State','fontsize',18);
%     xlabel('\it{t}','fontsize',18);
%     axis([0,Tfinal,0,1])
%     set(gca, 'FontSize', 18)
%     hold off
%     set(gca, 'FontName', 'Times New Roman')
%     drawnow    
    
end
Control = Xn;
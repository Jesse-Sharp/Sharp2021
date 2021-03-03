% Author: Jesse Sharp; Last Update: 03/03/2021
%  
%Performs the FBSM with the Wegsetin method for a given guess theta of Lambda(end) 
%and model parameters. 
function [z,SweepFevals] = Sweeps_Wegstein(theta,parameters)
%% Set-up
parameters.theta = theta;
tol=parameters.tol;
MaxFevals = parameters.MaxFevals;
N = parameters.N;
omega = parameters.omega;
Tfinal = parameters.Tfinal;
t_y = linspace(0,Tfinal,N); %Time discretisation for plotting

%Wegstein parameters
nth = parameters.nth;
bounding = parameters.bounding;
lower = parameters.lower;
upper = parameters.upper;

y0 = 1; %Linear model initial conditon
X0 = zeros(N,1); %Initial guess for the control

SweepFevals = 0;
iterations = 0;
Err = inf;

%% Perform Wegstein method
%Generate initial control updates for Wegstein method
[~,~,X1] = FBSM(y0,X0,parameters);
X1 = omega*X0 + (1-omega)*X1; %Apply relaxation factor to aid convergence
[~,~,X2] = FBSM(y0,X1,parameters);
X2 = omega*X1 + (1-omega)*X2; %Apply relaxation factor to aid convergence
SweepFevals = SweepFevals+2;

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
            else %Bounding in the case that q is complex
                normq = norm(q);%norm of complex elements
                r = min(normq,abs(upper-lower)); %Set length between 0 and magnitude of difference of bounds. 
                %This magnitude is chosen arbitrarly, but appears to perform well
                phi = angle(q); %Compute angle of each complex number
                q = r.*exp(1i*phi); %Construct new (bounded) complex number maintaining original angle but with bounded length. 
            end
        end
        %Update control elementwise using Wegstein q as a weighting factor
        %to accelerate convergence
        Xn = q.*X0+(1-q).*X1;
        [y,Lambda,Fn] = FBSM(y0,Xn,parameters);
        SweepFevals = SweepFevals+1;
        Err = norm(Xn-Fn);
    else %If not first iteration, only update q if we are on the nth iteration
        if mod(iterations,nth) == 0 %Compute Wegstein's q
            A = (Fn - X2)./(Xn-X1); %Compute A elementwise from the differences
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
        %Update control elementwise using Wegstein q as a weighting factor
        %to accelerate convergence
        Xn = q.*Xn + (1-q).*Fn;
        [y,Lambda,Fn] = FBSM(y0,Xn,parameters); %Solve system with updated control to obtain next control value
        SweepFevals = SweepFevals+1;
        Err = norm(Xn-Fn);
    end
    
    if SweepFevals > MaxFevals
        error('Maximum function evaluations of %i reached, current error of %e does not meet required tolerance of %e',MaxFevals,Err,tol)
    end
end
z = [t_y',y,Fn,Lambda];
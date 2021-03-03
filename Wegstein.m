% Author: Jesse Sharp; Last Update: 03/03/2021
%  
function [root,Fevals] = Wegstein(F,X0,tol,MaxFevals,nth,bounding,lower,upper)
%Solves for a root, G(X)=0, for an equation in the form X = F(X), with
%initial guess X0. 
%Input tol specifies the absolute tolerance required for convergence
%Input MaxFevals specifies the maximum number of function evaluations to 
%perform before terminating the process if convergence is not achieved.
%Input nth specifies how frequently to update q; every nth iteration. 
%Input bounding acts as a switch, such that bounds are applied on q when
%bounding = 1.
%Inputs lower and upper specify the bounds to apply to q, if bounding = 1.

%% Set-up
Fevals = 0;
iterations = 0;
Err = inf;

%% Perform Wegstein method
%Generate initial function evaluations for Wegstein method
X1 = F(X0);
X2 = F(X1);
Fevals = Fevals + 2;

while Err > tol
    iterations = iterations+1;
        
    if iterations == 1
        A = (X2-X1)./(X1-X0); %Compute A elementwise from the differences
        q = A./(A-1); %Compute q elementwise from A
        %Set NaN and Inf values to 0; choice of 0 is arbitrary, and could
        %be replaced with (for example) omega. 
        q(isnan(q))=0; 
        q(isinf(q))=0;
        %Apply bounds to Wegstein q
        if bounding == 1
            if isreal(q)
                q = min(max(q,lower),upper); %Apply bounds to q
            else %Bounding in the case that q is complex
                normq = abs(q);%element-wise norm of complex elements
                r = min(normq,abs(upper-lower)); %Set length between 0 and magnitude of difference of bounds. 
                %This magnitude is chosen arbitrarly, but appears to perform well
                phi = angle(q); %Compute angle of each complex number
                q = r.*exp(1i*phi); %Construct new (bounded) complex number maintaining original angle but with bounded length. 
            end
        end
        %Update control elementwise using Wegstein q as a weighting factor
        %to accelerate convergence
        Xn = q.*X0+(1-q).*X1;
        Fn = F(Xn);
        Fevals = Fevals+1;
        Err = norm(Fn-Xn);
    else %If not first iteration, only update q if we are on the nth iteration, and we have not reached iteration Wegstop
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
                normq = abs(q);%element-wise norm of complex elements
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
        Fn = F(Xn);
        Fevals = Fevals+1;
        Err = norm(Fn-Xn);
    end
    
    if Fevals > MaxFevals
        error('Maximum function evaluations of %i reached, current error of %e does not meet required tolerance of %e',MaxFevals,Err,tol)
    end
end
root = Xn; 
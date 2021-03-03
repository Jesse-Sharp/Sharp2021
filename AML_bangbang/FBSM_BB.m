% Author: Jesse Sharp; Last Update: 03/03/2021
% 
function [y,Lambda,Uupdate] = FBSM_BB(y,U,parameters)
dt = parameters.dt;
N = parameters.N;

i = 0; %Initialise loop variable

%Forward sweep using fourth-order Runge-Kutta scheme
while i < N-1
    i=i+1;
    k1 = State(y(i,:),U(i),parameters);
    k2 = State(y(i,:)+dt*k1/2,0.5*(U(i)+U(i+1)),parameters);
    k3 = State(y(i,:)+dt*k2/2,0.5*(U(i)+U(i+1)),parameters);
    k4 = State(y(i,:)+dt*k3,U(i+1),parameters);
    y(i+1,:) = y(i,:) + (dt/6)*(k1+2*k2+2*k3+k4);
end

Lambda = zeros(length(y),3); %Initialise co-state
Lambda(end,:) = [0;0;0]; %Apply transversality conditions to obtain final time condition on Lambda (costate)
i = 0; %Initialise loop variable
j = N; %Initialise loop variable

%Backward sweep using fourth-order Runge-Kutta scheme 
while j > 1 
    i = i+1;
    j = N-i;
    k1 = Costate_BB(Lambda(j+1,:),y(j+1,:),U(j+1),parameters);
    k2 = Costate_BB(Lambda(j+1,:)-dt*k1/2,0.5*(y(j+1,:)+y(j+1,:)),0.5*(U(j)+U(j+1)),parameters);
    k3 = Costate_BB(Lambda(j+1,:)-dt*k2/2,0.5*(y(j+1,:)+y(j+1,:)),0.5*(U(j)+U(j+1)),parameters);
    k4 = Costate_BB(Lambda(j+1,:)-dt*k3,y(j+1,:),U(j),parameters);
    Lambda(j,:) = Lambda(j+1,:) - (dt/6)*(k1+2*k2+2*k3+k4);
end

Uupdate = Control_BB(y,Lambda,parameters); %Calculated updated control


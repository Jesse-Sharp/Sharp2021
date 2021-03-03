% Author: Jesse Sharp; Last Update: 03/03/2021
%  
function f = Control(y,Lambda,parameters)

a1 = parameters.a1;

f = y(:,3).*Lambda(:,3)/(2*a1);

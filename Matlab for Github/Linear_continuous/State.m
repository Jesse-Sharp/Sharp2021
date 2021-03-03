% Author: Jesse Sharp; Last Update: 03/03/2021
%  
function f = State(y,U,parameters)

Gamma = parameters.Gamma; 

f = Gamma*y+U;
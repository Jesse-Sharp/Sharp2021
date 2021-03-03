% Author: Jesse Sharp; Last Update: 03/03/2021
%  
function f = Costate(Lambda,y,parameters)

a = parameters.a;
Gamma = parameters.Gamma;

f = -2*a*y-Gamma*Lambda;
% Author: Jesse Sharp; Last Update: 03/03/2021
%  
function f = Control_BB(Lambda,parameters)

b = parameters.b;
Ulower = parameters.Ulower; 
Uupper = parameters.Uupper; 

f = max(Ulower*((sign(Lambda-b)-1)/(-2)),Uupper*((sign(Lambda-b)+1)/(2))); %Update control within bounds

% Author: Jesse Sharp; Last Update: 03/03/2021
%  
function f = Control_BB(y,Lambda,parameters)

a1 = parameters.a1;
Ulower = parameters.Ulower;
Uupper = parameters.Uupper;

f = max(Uupper*((sign(-y(:,3).*Lambda(:,3)+a1)-1)/(-2)),Ulower*((sign(-y(:,3).*Lambda(:,3)+a1)+1)/(2)));
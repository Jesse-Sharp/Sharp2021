% Author: Jesse Sharp; Last Update: 03/03/2021
%  
function f = State(y,U,parameters)

ps = parameters.ps; 
pa = parameters.pa; 
pl = parameters.pl; 
gs = parameters.gs; 
ga = parameters.ga; 
gl = parameters.gl; 
K1 = parameters.K1; 
K2 = parameters.K2; 
Alpha = parameters.Alpha; 
Gamma = parameters.Gamma;

f = [ps*y(1)*(K1-y(1))-gs*y(1),gs*y(1)+pa*y(2)*(K2-(y(2)+y(3)))-ga*y(2),pl*y(3)*(K2-(y(2)+y(3)))-gl*y(3)-Alpha*y(3)/(Gamma+y(3))-U*y(3)];